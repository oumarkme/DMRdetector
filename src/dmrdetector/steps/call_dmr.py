"""CALL_DMR (v2.2): DMR calling of one sample against the reference atlas.

The sites of the sample are joined to the smoothed atlas, restricted to panel
regions, filtered on |DIFF|, clustered, tested with a one-sample t-test
(mu = 0) and Benjamini-Hochberg adjusted per panel. The result is the v2.2 raw
table `dmr_raw_<id>_<type>_<scope>.tsv.gz`.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
import io
import os
import sys

import numpy as np
import pandas as pd

from dmrdetector import settings
from dmrdetector.fileio import bed
from dmrdetector.util import rstats


RAW_COLS = [
    "CHROM", "START", "END", "PANEL", "PANEL_SCOPE", "NUM_CPG", "NUM_READS",
    "MEAN_SAMPLE", "MEAN_REF", "MEAN_DIFF", "P_VALUE", "P_FDR", "DIRECTION"
]

# v2.2 defaults (nextflow.config), taken from the shared defaults module
DEFAULTS = dict(
    cpg_diff_threshold=settings.DMR_CPG_DIFF_THRESHOLD,
    cpg_min_sites=settings.DMR_CPG_MIN_SITES,
    icr_diff_threshold=settings.ICR_DIFF_THRESHOLD,
    icr_min_sites=settings.ICR_MIN_SITES,
    merge_distance=settings.DMR_MERGE_DISTANCE,
    fdr_threshold=settings.DMR_FDR_THRESHOLD,
    sex_f_diff_threshold=settings.SEX_F_DIFF_THRESHOLD,
    sex_f_min_sites=settings.SEX_F_MIN_SITES,
)


def _read_panel_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={0: str}, keep_default_na=True)


def _union_intervals(q: pd.DataFrame) -> dict:
    """Returns the merged inclusive START ranges per chromosome.

    The ranges are those covered by `tabix -R` on the BED regions [S, E): START
    in [S, E-1].
    """
    out = {}
    for chrom, g in q.groupby("CHROM", sort=False):
        lo = g["S"].to_numpy(np.int64)
        hi = g["E"].to_numpy(np.int64) - 1
        keep = hi >= lo
        lo, hi = lo[keep], hi[keep]
        if lo.size == 0:
            continue
        o = np.argsort(lo, kind="stable")
        lo, hi = lo[o], hi[o]
        cmax = np.maximum.accumulate(hi)
        new = np.ones(lo.size, bool)
        new[1:] = lo[1:] > cmax[:-1]
        starts = np.flatnonzero(new)
        ends = np.append(starts[1:], lo.size) - 1
        out[chrom] = (lo[starts], cmax[ends])
    return out


def _in_union(chrom: pd.Series, start: np.ndarray, union: dict) -> np.ndarray:
    ok = np.zeros(start.size, bool)
    chrom = chrom.to_numpy(dtype=object)
    for c in np.unique(chrom):
        if c not in union:
            continue
        lo, hi = union[c]
        m = chrom == c
        s = start[m]
        idx = np.searchsorted(lo, s, side="right") - 1
        good = idx >= 0
        good[good] = s[good] <= hi[idx[good]]
        ok[np.flatnonzero(m)] = good
    return ok


def _multiplicity(chrom: pd.Series, start: np.ndarray,
                  query: pd.DataFrame) -> np.ndarray:
    """Number of query lines [S, E) that contain each site START (0 when none).

    `tabix -R` (htslib 1.24, the version in the v2.2 baseline image) emits a
    record once per query line that overlaps it, including exact duplicate
    lines. v2.2 therefore sees duplicated sample and atlas rows wherever panel
    regions overlap (e.g. ICR rows with identical coordinates) and its merge()
    joins them m x m.
    """
    cnt = np.zeros(start.size, np.int64)
    chrom = chrom.to_numpy(dtype=object)
    for c, g in query.groupby("CHROM", sort=False):
        s_ = g["S"].to_numpy(np.int64)
        h_ = g["E"].to_numpy(np.int64) - 1
        keep = h_ >= s_
        s_, h_ = np.sort(s_[keep]), np.sort(h_[keep])
        m = np.flatnonzero(chrom == c)
        if m.size and s_.size:
            p = start[m]
            cnt[m] = np.searchsorted(s_, p, side="right") - np.searchsorted(
                h_, p, side="left")
    return cnt


def _read_fast(path: str, names: list, dtypes: dict, union: dict, allowed: set,
               threads: int):
    """Reads a BGZF file in memory (threaded inflate, numba parser).

    No index is needed. Only the columns used later are parsed (`_FAST_COLS`).

    Returns:
        The sites with START inside `union` as a DataFrame with CHROM and those
        columns, or None when the file is not BGZF or the fast reader cannot
        handle it.
    """
    if os.environ.get("DMRDETECTOR_NO_JIT") == "1":
        return None
    # numba is imported only when a file is read here (annotate does not need
    # it)
    from dmrdetector.fileio import bgzf
    from dmrdetector.fileio import numparse
    cols = _FAST_COLS[tuple(names)]
    spec = [(names.index(c), "i" if c in ("START", "COVERAGE",
                                          "N_METH") else "f") for c in cols]
    try:
        r = numparse.read_table(path, spec, threads)
    except bgzf.ParseError:
        return None
    if r is None:
        return None
    ints, floats, (rs, rnames) = r
    arrs = {}
    it, ft = iter(ints), iter(floats)
    for c, (_, k) in zip(cols, spec):
        arrs[c] = next(it) if k == "i" else next(ft)
    n = arrs["START"].size
    bounds = np.append(rs, n)
    norm = bed.normalize_primary_chrom(rnames).tolist() if len(rnames) else []
    keep_idx, chroms = [], []
    for i, c in enumerate(norm):
        if c not in allowed or c not in union:
            continue
        a, b = int(bounds[i]), int(bounds[i + 1])
        lo, hi = union[c]
        s_ = arrs["START"][a:b]
        j = np.searchsorted(lo, s_, side="right") - 1
        good = j >= 0
        good[good] = s_[good] <= hi[j[good]]
        ix = np.flatnonzero(good) + a
        keep_idx.append(ix)
        chroms.append(np.full(ix.size, c, dtype=object))
    if not keep_idx:
        return pd.DataFrame(
            {c: pd.Series(dtype=dtypes[c]) for c in ["CHROM"] + cols})
    ix = np.concatenate(keep_idx)
    df = pd.DataFrame({
        "CHROM": np.concatenate(chroms),
        **{
            c: arrs[c][ix] for c in cols
        }
    })
    return df.astype({c: dtypes[c] for c in df.columns})


# columns each table needs after reading (the rest of the 5 / 9 columns is never
# used by call_dmr)
_SAMPLE_COLS = ("CHROM", "START", "END", "COVERAGE", "N_METH")
_ATLAS_COLS = ("CHROM", "START", "END", "N", "COV", "PHI", "RAW_METH",
               "REF_METH", "REF_VAR")
_FAST_COLS = {
    _SAMPLE_COLS: ["START", "COVERAGE", "N_METH"],
    _ATLAS_COLS: ["START", "REF_METH", "REF_VAR"],
}


def _read_bed_restricted(path: str,
                         names: list,
                         dtypes: dict,
                         union: dict,
                         allowed: set,
                         query: pd.DataFrame | None = None,
                         threads: int | None = None) -> pd.DataFrame:
    """Reads a bgzip BED(-like) file and keeps the sites with START in `union`.

    This follows `tabix -R`. With `query`, each site is repeated once per query
    line containing it (tabix 1.24 behavior). BGZF files without an index, or
    with `threads` > 1, are read in memory by the numba reader (the index is not
    used). Indexed files with one thread, other files, a file the fast reader
    rejects and DMRDETECTOR_NO_JIT=1 use tabix region reads (indexed) or pandas
    (not indexed).
    """
    # Indexed file and one thread: the tabix region reads use less CPU and
    # memory than inflating the whole file. No index, or more than one thread:
    # read in memory.
    threads = threads or 1
    fast = tuple(names) in _FAST_COLS and (threads > 1 or
                                           not os.path.exists(path + ".tbi"))
    df = _read_fast(path, names, dtypes, union, allowed,
                    threads) if fast else None
    if df is not None:
        if query is not None and len(df):
            df = df.iloc[np.repeat(
                np.arange(len(df)),
                _multiplicity(df["CHROM"], df["START"].to_numpy(np.int64),
                              query))]
            df = df.reset_index(drop=True)
        return df
    tbi = path + ".tbi"
    frames = []
    use_tabix = os.path.exists(tbi)
    if use_tabix:
        import pysam
        tf = pysam.TabixFile(path)
        contigs = set(tf.contigs)
        for chrom, (lo, hi) in union.items():
            name = chrom if chrom in contigs else (
                chrom[3:] if chrom[3:] in contigs else None)
            if name is None:
                continue
            lines = []
            for a, b in zip(lo.tolist(), hi.tolist()):
                lines.extend(tf.fetch(name, a, b + 1))
            if lines:
                f = pd.read_csv(io.StringIO("\n".join(lines)),
                                sep="\t",
                                header=None,
                                names=names,
                                usecols=range(len(names)),
                                float_precision="round_trip")
                f["CHROM"] = chrom
                frames.append(f)
        tf.close()
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
            columns=names)
    else:
        df = pd.read_csv(path,
                         sep="\t",
                         header=None,
                         names=names,
                         usecols=range(len(names)),
                         comment="#",
                         dtype={"CHROM": str},
                         float_precision="round_trip")
    if len(df) == 0:
        return pd.DataFrame({n: pd.Series(dtype=dtypes[n]) for n in names})
    df["CHROM"] = bed.normalize_primary_chrom(df["CHROM"]).to_numpy()
    df = df[df["CHROM"].isin(allowed)]
    df = df[_in_union(df["CHROM"], df["START"].to_numpy(np.int64), union)]
    df = df.astype(dtypes).reset_index(drop=True)
    if query is not None and len(df):
        df = df.iloc[np.repeat(
            np.arange(len(df)),
            _multiplicity(df["CHROM"], df["START"].to_numpy(np.int64), query))]
        df = df.reset_index(drop=True)
    return df


def find_sites_in_regions(sites: pd.DataFrame,
                          regions: pd.DataFrame) -> pd.DataFrame:
    """Returns the overlaps of sites and regions, as data.table foverlaps().

    The site interval [START, START+1] and the region [REGION_START, REGION_END]
    are closed.

    Returns:
        One row per (site, region) overlap, ordered by site, with REGION_START
            and
        REGION_END.
    """
    has_core = "CORE_START" in regions.columns
    cols = list(sites.columns) + ["REGION_START", "REGION_END"
                                 ] + (["CORE_START"] if has_core else [])
    if len(sites) == 0 or len(regions) == 0:
        return pd.DataFrame({
            c: pd.Series(dtype=(sites[c].dtype if c in sites else "int64"))
            for c in cols
        })
    pieces = []
    spos = sites["START"].to_numpy(np.int64)
    schr = sites["CHROM"].to_numpy(dtype=object)
    for chrom, rg in regions.groupby("CHROM", sort=False):
        sidx_all = np.flatnonzero(schr == chrom)
        if sidx_all.size == 0:
            continue
        rs = rg["REGION_START"].to_numpy(np.int64)
        re_ = rg["REGION_END"].to_numpy(np.int64)
        core = rg["CORE_START"].to_numpy(np.int64) if has_core else rs
        lo = rs - 1  # START + 1 >= REGION_START
        # y is keyed by (CHROM, REGION_START, REGION_END)
        o = np.lexsort((re_, rs))
        lo, hi, rs, re_, core = lo[o], re_[o], rs[o], re_[o], core[o]
        rank = np.arange(lo.size)
        # split into layers of mutually non-overlapping regions
        layer_last: list[int] = []
        layer_of = np.empty(lo.size, np.int64)
        for i in np.argsort(lo, kind="stable"):
            for k, last in enumerate(layer_last):
                if last < lo[i]:
                    layer_last[k] = hi[i]
                    layer_of[i] = k
                    break
            else:
                layer_last.append(hi[i])
                layer_of[i] = len(layer_last) - 1
        sp = spos[sidx_all]
        for k in range(len(layer_last)):
            m = layer_of == k
            llo, lhi, lrs, lre, lcore, lrank = lo[m], hi[m], rs[m], re_[
                m], core[m], rank[m]
            ord_ = np.argsort(llo, kind="stable")
            llo, lhi, lrs, lre, lcore, lrank = llo[ord_], lhi[ord_], lrs[
                ord_], lre[ord_], lcore[ord_], lrank[ord_]
            idx = np.searchsorted(llo, sp, side="right") - 1
            good = idx >= 0
            good[good] = sp[good] <= lhi[idx[good]]
            if good.any():
                pieces.append((sidx_all[good], lrs[idx[good]], lre[idx[good]],
                               lcore[idx[good]], lrank[idx[good]]))
    if not pieces:
        return pd.DataFrame({
            c: pd.Series(dtype=(sites[c].dtype if c in sites else "int64"))
            for c in cols
        })
    si = np.concatenate([p[0] for p in pieces])
    rk = np.concatenate([p[4] for p in pieces])
    o = np.lexsort((rk, si))  # foverlaps: ordered by site, then by region key
    res = sites.iloc[si[o]].reset_index(drop=True)
    res["REGION_START"] = np.concatenate([p[1] for p in pieces])[o]
    res["REGION_END"] = np.concatenate([p[2] for p in pieces])[o]
    if has_core:
        res["CORE_START"] = np.concatenate([p[3] for p in pieces])[o]
    return res


def _cluster_ids(chrom: np.ndarray, start: np.ndarray, max_gap: float,
                 group: np.ndarray | None) -> np.ndarray:
    """Returns cumsum(gap > max_gap); the gap restarts per chrom or group."""
    n = chrom.size
    new = np.ones(n, bool)
    if n > 1:
        same = chrom[1:] == chrom[:-1]
        if group is not None:
            same &= group[1:] == group[:-1]
        gap = np.diff(start)
        new[1:] = ~same | (gap > max_gap)
    return np.cumsum(new) - 1


def _nanmean_by(vals: np.ndarray, starts: np.ndarray) -> np.ndarray:
    """Returns mean(x, na.rm = TRUE) per group.

    Accumulation is in long double with the second refinement pass of R.
    """
    v = vals.astype(np.longdouble)
    nn = ~np.isnan(vals)
    v0 = np.where(nn, v, np.longdouble(0))
    c = np.add.reduceat(nn.astype(np.int64), starts)
    cl = np.maximum(c, 1).astype(np.longdouble)
    s = np.add.reduceat(v0, starts) / cl
    d = np.where(nn, v - np.repeat(s, np.diff(np.r_[starts, len(v)])),
                 np.longdouble(0))
    s = s + np.add.reduceat(d, starts) / cl
    return np.where(c > 0, s.astype(float), np.nan)


def merge_sites(sites: pd.DataFrame,
                max_gap: float,
                by_region: bool = False,
                pos_col: str = "START") -> pd.DataFrame:
    """Clusters over-threshold sites and summarizes each cluster.

    A gap > max_gap starts a new cluster. With by_region=False the sites are
    sorted by CHROM, START and clustered per chromosome (CpG panel). With
    by_region=True they are grouped by (CHROM, REGION_START, REGION_END) in
    order of first appearance and clustered within the group (ICR panel).

    `pos_col` is the coordinate used for the gap, START (min) and END (max + 1).
    v2.2 quirk reproduced for the CpG panel: the foverlaps result has two
    columns named START (panel core START first, site position second). v2.2
    resolves `START` to the first one in gap, min and max, and its setorder() is
    a no-op because the table already carries the key (CHROM, START). CpG
    clusters are therefore built on the panel core START (pos_col="CORE_START")
    over rows in site order; the gap can be negative (no new cluster) when the
    core START drops back.
    """
    if len(sites) == 0:
        return pd.DataFrame({
            "CHROM": pd.Series(dtype=object),
            "START": pd.Series(dtype="int64"),
            "END": pd.Series(dtype="int64"),
            "NUM_CPG": pd.Series(dtype="int64"),
            "NUM_READS": pd.Series(dtype=float),
            "MEAN_SAMPLE": pd.Series(dtype=float),
            "MEAN_REF": pd.Series(dtype=float),
            "MEAN_DIFF": pd.Series(dtype=float),
            "DIFF_VALUES": pd.Series(dtype=object)
        })
    if by_region:
        gid = sites.groupby(["CHROM", "REGION_START", "REGION_END"],
                            sort=False).ngroup().to_numpy()
        o = np.lexsort((sites[pos_col].to_numpy(), gid))
        s = sites.iloc[o].reset_index(drop=True)
        group = gid[o]
    else:
        # no re-sorting for CORE_START (see above); START is already in site
        # order
        s = sites.reset_index(drop=True)
        group = None
    chrom = s["CHROM"].to_numpy(dtype=object)
    start = s[pos_col].to_numpy(np.int64)
    cl = _cluster_ids(chrom, start, max_gap, group)
    first = np.flatnonzero(np.r_[True, cl[1:] != cl[:-1]])
    last = np.r_[first[1:], len(s)] - 1
    diff = s["DIFF"].to_numpy(float)
    out = pd.DataFrame({
        "CHROM": chrom[first],
        "START": np.minimum.reduceat(start, first),
        "END": np.maximum.reduceat(start, first) + 1,
        "NUM_CPG": (last - first + 1).astype(np.int64),
        "NUM_READS": _nanmean_by(s["COVERAGE"].to_numpy(float), first),
        "MEAN_SAMPLE": _nanmean_by(s["SAMPLE_METH"].to_numpy(float), first),
        "MEAN_REF": _nanmean_by(s["REF_METH"].to_numpy(float), first),
        "MEAN_DIFF": _nanmean_by(diff, first),
    })
    out["DIFF_VALUES"] = [diff[a:b + 1] for a, b in zip(first, last)]
    return out


def finalise_dmrs(merged: pd.DataFrame, fdr_threshold: float,
                  min_sites: int) -> pd.DataFrame:
    """Tests each cluster and keeps the significant, large ones.

    It runs a t-test per cluster and BH over the non-NA p-values of the panel,
    and keeps FDR < `fdr_threshold` and NUM_CPG > `min_sites`.
    """
    if len(merged) == 0:
        out = merged.drop(columns="DIFF_VALUES").copy()
        out["P_VALUE"] = pd.Series(dtype=float)
        out["P_FDR"] = pd.Series(dtype=float)
        return out
    p = np.array([
        rstats.t_test_one_sample_p(d, 0.0) if len(d) >= 2 else np.nan
        for d in merged["DIFF_VALUES"]
    ])
    out = merged.drop(columns="DIFF_VALUES").copy()
    out["P_VALUE"] = p
    out["P_FDR"] = rstats.p_adjust_bh(p)
    keep = (out["P_FDR"] < fdr_threshold).to_numpy() & (out["NUM_CPG"]
                                                        > min_sites).to_numpy()
    return out[keep].reset_index(drop=True)


def resolve_thresholds(scope: str,
                       cpg_diff_threshold: float | None = None,
                       cpg_min_sites: int | None = None,
                       **kw: float) -> dict:
    """Returns the thresholds of a scope.

    As in the v2.2 workflow, female_X uses sex_f_diff_threshold and
    sex_f_min_sites; the other scopes use the dmr_cpg_* ones.
    """
    d = dict(DEFAULTS)
    d.update({k: v for k, v in kw.items() if v is not None})
    female_x = scope == "female_X"
    if cpg_diff_threshold is None:
        cpg_diff_threshold = d["sex_f_diff_threshold"] if female_x else d[
            "cpg_diff_threshold"]
    if cpg_min_sites is None:
        cpg_min_sites = d["sex_f_min_sites"] if female_x else d["cpg_min_sites"]
    d["cpg_diff_threshold"], d[
        "cpg_min_sites"] = cpg_diff_threshold, cpg_min_sites
    return d


def call_dmr(sample_bed: str,
             atlas: str,
             cpg_panel: str,
             icr_panel: str | None,
             scope: str = "autosomal",
             include_sex: bool = True,
             log: Callable[[str], None] = print,
             tabix_duplicates: bool = True,
             threads: int | None = None,
             stats: dict | None = None,
             **thr: float | None) -> pd.DataFrame:
    """Runs CALL_DMR for one sample and one scope.

    Args:
        sample_bed: Preprocessed sample BED.
        atlas: Smoothed reference atlas.
        cpg_panel: CpG island panel of the scope.
        icr_panel: ICR panel (autosomal scope only).
        scope: autosomal, male_X, female_X or Y.
        include_sex: Keep chrX and chrY.
        log: Receives progress messages.
        tabix_duplicates: True reproduces v2.2 (sites repeated once per
            overlapping panel query line, see `_multiplicity`); False reads
                every
            site once.
        threads: Reader threads.
        stats: If given, receives sample_sites and atlas_sites (CpGs inside the
            panel regions) and merged (sites present in both).
        **thr: Thresholds (see resolve_thresholds).

    Returns:
        The raw DMR table (RAW_COLS).
    """
    t = resolve_thresholds(scope, **thr)
    run_icr = scope == "autosomal" and icr_panel is not None and os.path.exists(
        icr_panel)
    allowed = {f"chr{i}" for i in range(1, 23)}
    if include_sex:
        allowed |= {"chrX", "chrY"}

    # ---- panels and tabix-equivalent query ----
    cpg = _read_panel_table(cpg_panel)
    cpg_cols = list(cpg.columns)
    if "START_EXTEND" in cpg_cols and "END_EXTEND" in cpg_cols:
        qcpg = pd.DataFrame({
            "CHROM": cpg["CHROM"],
            "S": cpg["START_EXTEND"],
            "E": cpg["END_EXTEND"]
        })
    else:
        qcpg = pd.DataFrame({
            "CHROM": cpg.iloc[:, 0],
            "S": cpg.iloc[:, 1],
            "E": cpg.iloc[:, 2]
        })
    qcpg["CHROM"] = bed.normalize_primary_chrom(qcpg["CHROM"]).to_numpy()
    qcpg = qcpg[qcpg["CHROM"].isin(allowed)]
    icr = None
    qs = [qcpg]
    if run_icr:
        icr = _read_panel_table(icr_panel).iloc[:, :3]
        icr.columns = ["CHROM", "START", "END"]
        icr["CHROM"] = bed.normalize_primary_chrom(icr["CHROM"]).to_numpy()
        icr = icr[icr["CHROM"].isin({f"chr{i}" for i in range(1, 23)})]
        qs.append(
            pd.DataFrame({
                "CHROM": icr["CHROM"],
                "S": icr["START"],
                "E": icr["END"]
            }))
    query = pd.concat(qs, ignore_index=True)
    query[["S", "E"]] = query[["S", "E"]].astype(np.int64)
    union = _union_intervals(query)

    log(f"Processing {os.path.basename(sample_bed)} | scope={scope} | "
        f"diff>{t['cpg_diff_threshold']:.2f}")
    sample = _read_bed_restricted(
        sample_bed, ["CHROM", "START", "END", "COVERAGE", "N_METH"], {
            "CHROM": object,
            "START": "int64",
            "END": "int64",
            "COVERAGE": "int64",
            "N_METH": "int64"
        }, union, allowed, query if tabix_duplicates else None, threads)
    sample["SAMPLE_METH"] = sample["N_METH"] / sample["COVERAGE"]
    atl = _read_bed_restricted(
        atlas, [
            "CHROM", "START", "END", "N", "COV", "PHI", "RAW_METH", "REF_METH",
            "REF_VAR"
        ], {
            "CHROM": object,
            "START": "int64",
            "END": "int64",
            "N": "int64",
            "COV": "int64",
            "PHI": float,
            "RAW_METH": float,
            "REF_METH": float,
            "REF_VAR": float
        }, union, allowed, query if tabix_duplicates else None, threads)
    log(f"  Sample sites in panel: {len(sample)}; atlas sites in panel: "
        f"{len(atl)}")

    merged = sample[["CHROM", "START", "COVERAGE", "N_METH", "SAMPLE_METH"
                    ]].merge(atl[["CHROM", "START", "REF_METH", "REF_VAR"]],
                             on=["CHROM", "START"],
                             how="inner")
    merged = merged.sort_values(["CHROM", "START"],
                                kind="stable").reset_index(drop=True)
    merged["DIFF"] = merged["SAMPLE_METH"] - merged["REF_METH"]
    log(f"  Merged: {len(merged)} sites")
    if stats is not None:
        stats.update(sample_sites=len(sample),
                     atlas_sites=len(atl),
                     merged=len(merged))

    # ---- CpG panel ----
    cpg = cpg[bed.normalize_primary_chrom(
        cpg["CHROM"]).isin(allowed).to_numpy()].copy()
    cpg["CHROM"] = bed.normalize_primary_chrom(cpg["CHROM"]).to_numpy()
    reg = pd.DataFrame({
        "CHROM": cpg["CHROM"],
        "REGION_START": cpg["START_EXTEND"].astype(np.int64),
        "REGION_END": cpg["END_EXTEND"].astype(np.int64),
        "CORE_START": cpg["START"].astype(np.int64)
    })
    cpg_sites = find_sites_in_regions(merged, reg)
    cpg_f = cpg_sites[np.abs(cpg_sites["DIFF"].to_numpy()) >
                      t["cpg_diff_threshold"]]
    cpg_dmrs = finalise_dmrs(
        merge_sites(cpg_f, t["merge_distance"], pos_col="CORE_START"),
        t["fdr_threshold"], t["cpg_min_sites"])
    cpg_dmrs["PANEL"] = "CpG"
    log(f"  CpG DMRs: {len(cpg_dmrs)} (from {len(cpg_sites)} sites)")

    # ---- ICR panel (autosomal scope only) ----
    icr_dmrs = cpg_dmrs.iloc[0:0].copy()
    if run_icr:
        ireg = icr.rename(columns={
            "START": "REGION_START",
            "END": "REGION_END"
        }).astype({
            "REGION_START": np.int64,
            "REGION_END": np.int64
        })
        icr_sites = find_sites_in_regions(merged, ireg)
        icr_f = icr_sites[np.abs(icr_sites["DIFF"].to_numpy()) >
                          t["icr_diff_threshold"]]
        icr_dmrs = finalise_dmrs(
            merge_sites(icr_f, t["merge_distance"], by_region=True),
            t["fdr_threshold"], t["icr_min_sites"])
        icr_dmrs["PANEL"] = "ICR"
        log(f"  ICR DMRs: {len(icr_dmrs)} (from {len(icr_sites)} sites)")

    allx = pd.concat([cpg_dmrs, icr_dmrs], ignore_index=True)
    allx["PANEL_SCOPE"] = scope
    allx["DIRECTION"] = np.where(allx["MEAN_DIFF"] > 0, "hyper", "hypo")
    allx = allx[RAW_COLS]
    allx = allx.sort_values(["CHROM", "START"],
                            kind="stable").reset_index(drop=True)
    log(f"Total DMRs: {len(allx)} (CpG: {len(cpg_dmrs)}, ICR: {len(icr_dmrs)})")
    return allx


def write_raw(df: pd.DataFrame, path: str) -> None:
    """Writes a raw DMR table as tab-separated text, as fwrite does."""
    bed.fwrite_table(df,
                     path,
                     int_cols=("START", "END", "NUM_CPG"),
                     float_cols=("NUM_READS", "MEAN_SAMPLE", "MEAN_REF",
                                 "MEAN_DIFF", "P_VALUE", "P_FDR"))


def add_parser(sub: argparse._SubParsersAction) -> None:
    """Adds the `call-dmr` command to `sub`."""
    p = sub.add_parser("call-dmr",
                       help="CALL_DMR: call DMRs for one sample and scope",
                       description="CALL_DMR: per-sample DMR calling against a "
                       "reference atlas (v2.2 algorithm).")
    p.add_argument(
        "--sample-bed",
        required=True,
        help="preprocessed sample BED.gz (chrom, start, end, coverage, "
        "n_meth); "
        ".tbi used if present")
    p.add_argument(
        "--atlas",
        required=True,
        help="smoothed reference atlas BED.gz (9 columns, narrow bandwidth)")
    p.add_argument(
        "--cpg-panel",
        required=True,
        help="CpG panel BED of the scope (DETECT_CPG_ISLANDS output)")
    p.add_argument("--icr-panel",
                   help="ICR panel BED (used for the autosomal scope only)")
    p.add_argument("--scope",
                   default="autosomal",
                   choices=["autosomal", "male_X", "female_X", "Y"])
    p.add_argument("--out",
                   required=True,
                   help="output dmr_raw_<id>_<type>_<scope>.tsv.gz")
    p.add_argument("--cpg-diff-threshold",
                   type=float,
                   help="default 0.40 (female_X: --sex-f-diff-threshold)")
    p.add_argument(
        "--cpg-min-sites",
        type=int,
        help="default 30 (female_X: --sex-f-min-sites); NUM_CPG must exceed it")
    p.add_argument("--sex-f-diff-threshold",
                   type=float,
                   default=DEFAULTS["sex_f_diff_threshold"])
    p.add_argument("--sex-f-min-sites",
                   type=int,
                   default=DEFAULTS["sex_f_min_sites"])
    p.add_argument("--icr-diff-threshold",
                   type=float,
                   default=DEFAULTS["icr_diff_threshold"])
    p.add_argument("--icr-min-sites",
                   type=int,
                   default=DEFAULTS["icr_min_sites"])
    p.add_argument("--merge-distance",
                   type=int,
                   default=DEFAULTS["merge_distance"])
    p.add_argument("--fdr-threshold",
                   type=float,
                   default=DEFAULTS["fdr_threshold"])
    p.add_argument(
        "--dedupe-sites",
        action="store_true",
        help="read each site once; default reproduces v2.2, where overlapping "
        "panel regions "
        "(e.g. identical ICR rows) make tabix return sites repeatedly and "
        "inflate NUM_CPG")
    p.add_argument("--no-sex-chromosomes",
                   action="store_true",
                   help="restrict to chr1-22 (include_sex_chromosomes=false)")
    p.add_argument(
        "--threads",
        type=int,
        default=bed.default_threads(4),
        help="threads for inflating and parsing the input BEDs (the result "
        "does "
        "not depend on it)")
    p.set_defaults(func=main)


def main(args: argparse.Namespace) -> int:
    """Runs the `call-dmr` command."""
    df = call_dmr(args.sample_bed,
                  args.atlas,
                  args.cpg_panel,
                  args.icr_panel,
                  scope=args.scope,
                  include_sex=not args.no_sex_chromosomes,
                  tabix_duplicates=not args.dedupe_sites,
                  threads=args.threads,
                  log=lambda m: print(m, file=sys.stderr),
                  cpg_diff_threshold=args.cpg_diff_threshold,
                  cpg_min_sites=args.cpg_min_sites,
                  sex_f_diff_threshold=args.sex_f_diff_threshold,
                  sex_f_min_sites=args.sex_f_min_sites,
                  icr_diff_threshold=args.icr_diff_threshold,
                  icr_min_sites=args.icr_min_sites,
                  merge_distance=args.merge_distance,
                  fdr_threshold=args.fdr_threshold)
    write_raw(df, args.out)
    return 0
