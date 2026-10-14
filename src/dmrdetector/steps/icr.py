"""DETECT_ICR_REGIONS (v2.2): the imprinted control region panel (autosomal).

The output has CHROM START END GENE ICR_ID NUM_CPG; END is the START of the last
hemi-methylated CpG plus 1. Only atlas sites that can matter are read: a tabix
index gives the merged candidate windows; otherwise `atlas_io` scans the whole
file with the same filters. The v2.2 multiplicity and clustering code then runs
on these few rows.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from io import StringIO
import os
import shutil
import tempfile
import zlib

import numpy as np
import pandas as pd

from dmrdetector import settings
from dmrdetector.fileio import atlas_io
from dmrdetector.fileio import bed
from dmrdetector.steps import cpg_islands


OUT_COLS = ["CHROM", "START", "END", "GENE", "ICR_ID", "NUM_CPG"]


def _empty() -> pd.DataFrame:
    return pd.DataFrame({
        c:
            pd.Series(dtype="object" if c in ("CHROM", "GENE",
                                              "ICR_ID") else "int64")
        for c in OUT_COLS
    })


def read_candidates(path: str) -> pd.DataFrame:
    """Reads the first 4 columns as CHROM START END GENE.

    The header is detected like fread does; only autosomal rows are kept.
    """
    with open(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")
    try:
        float(first[1])
        header = None
    except (ValueError, IndexError):
        header = 0
    df = pd.read_csv(path, sep="\t", header=header, dtype=str)
    df = df.iloc[:, :min(4, df.shape[1])]
    df.columns = ["CHROM", "START", "END", "GENE"][:df.shape[1]]
    df["CHROM"] = bed.normalize_primary_chrom(df["CHROM"])
    df = df[df["CHROM"].isin(bed.AUTOSOMES)].copy()
    df["START"] = df["START"].astype(float).astype(np.int64)
    df["END"] = df["END"].astype(float).astype(np.int64)
    return df.reset_index(drop=True)


def _core(a: pd.DataFrame, cand: pd.DataFrame, ref_lo, ref_hi, min_sites,
          merge_distance, v1_duplicates, log):
    """Clusters the hemi-methylated sites of the candidates into ICRs.

    Args:
        a: Autosomal sites (CHROM normalized, N filter applied) with CHROM START
            END REF_METH.
        cand: Candidate intervals with GENE.
        ref_lo: Smallest REF_METH of a hemi-methylated CpG.
        ref_hi: Largest REF_METH.
        min_sites: Minimum CpGs of an ICR.
        merge_distance: Largest gap (bp) inside one ICR.
        v1_duplicates: Count a site once per overlapping candidate (v2.2).
        log: Receives progress messages.
    """
    parts = []
    for chrom, c in cand.groupby("CHROM", sort=False):
        g = a[a["CHROM"] == chrom]
        if g.empty:
            continue
        s = g["START"].to_numpy(np.int64)
        e = g["END"].to_numpy(np.int64)
        # v2.2 pre-filters with `tabix -R` (BED, half-open), which returns a
        # record once per query row it overlaps (regions are not merged), then
        # joins with foverlaps (closed intervals). A site inside k overlapping
        # candidate intervals therefore enters the join k times and is counted
        # k times in NUM_CPG. Reproduced here; v1_duplicates=False counts once.
        qs = c["START"].to_numpy(np.int64)
        qe = c["END"].to_numpy(np.int64)
        mult = np.zeros(len(g), np.int64)
        for a0, b0 in zip(qs, qe):
            mult += (s < b0) & (e > a0)
        if not v1_duplicates:
            mult = np.minimum(mult, 1)
        meth = g["REF_METH"].to_numpy(float)
        ok = (meth >= ref_lo) & (meth <= ref_hi)
        for gene, a0, b0 in zip(c["GENE"], qs, qe):
            m = (mult > 0) & (s <= b0) & (e >= a0) & ok
            if m.any():
                parts.append(
                    pd.DataFrame({
                        "CHROM": chrom,
                        "START": np.repeat(s[m], mult[m]),
                        "GENE": gene
                    }))
    if not parts:
        log("  Hemi-methylated ICR CpGs: 0")
        return _empty()
    hemi = pd.concat(parts, ignore_index=True)
    log(f"  Hemi-methylated ICR CpGs: {len(hemi)}")
    hemi = hemi.sort_values(["GENE", "CHROM", "START"],
                            kind="stable").reset_index(drop=True)
    gap = hemi.groupby(["GENE", "CHROM"],
                       sort=False)["START"].diff().fillna(np.inf)
    hemi["cluster"] = (gap > merge_distance).astype(int).groupby(
        [hemi["GENE"], hemi["CHROM"]], sort=False).cumsum()
    icr = hemi.groupby(["GENE", "CHROM", "cluster"],
                       sort=False).agg(START=("START", "min"),
                                       END=("START", "max"),
                                       NUM_CPG=("START", "size")).reset_index()
    icr["END"] = icr["END"] + 1
    icr = icr[icr["NUM_CPG"] >= min_sites]
    if icr.empty:
        log("No ICR passed min_sites. Empty output.")
        return _empty()
    icr = icr.sort_values(["CHROM", "START", "GENE"],
                          kind="stable").reset_index(drop=True)
    idx = icr.groupby(["GENE", "CHROM"], sort=False).cumcount() + 1
    icr["ICR_ID"] = icr["GENE"] + "_" + icr["CHROM"] + "_icr" + idx.map(
        "{:02d}".format)
    out = icr[OUT_COLS]
    log(f"Total ICRs: {len(out)} across {out['GENE'].nunique()} genes")
    return out


def _overlaps_candidates(chrom: str, s: np.ndarray, e: np.ndarray,
                         cand: pd.DataFrame) -> np.ndarray:
    """Returns True where [s, e] (closed) overlaps a candidate of `chrom`.

    This is the same test as the kernel applies.
    """
    c = cand[cand["CHROM"] == chrom]
    if c.empty:
        return np.zeros(s.size, bool)
    a = c["START"].to_numpy(np.int64)
    b = c["END"].to_numpy(np.int64)
    o = np.argsort(a, kind="stable")
    a, pm = a[o], np.maximum.accumulate(b[o])
    k = np.searchsorted(a, e, side="right")
    return (k > 0) & (pm[np.maximum(k - 1, 0)] >= s)


def _windows_fetch(path: str, cand: pd.DataFrame, scratch_root: str | None):
    """Reads the atlas lines of the merged candidate windows through tabix.

    Returns:
        A DataFrame of the raw columns, or None when this route is not possible
        (no pysam, not indexable, records longer than 1 bp).
    """
    try:
        import pysam
    except ImportError:
        return None
    scratch = tempfile.mkdtemp(prefix="icr_tbi_", dir=scratch_root)
    try:
        try:
            ipath = bed.ensure_tabix(path, scratch)
            tbx = pysam.TabixFile(ipath)
        except (OSError, ValueError, RuntimeError, zlib.error):
            return None
        with tbx:
            by_norm: dict[str, list[str]] = {}
            names = list(tbx.contigs)
            for nm, norm in zip(
                    names,
                    bed.normalize_primary_chrom(pd.Series(names, dtype=str))):
                by_norm.setdefault(norm, []).append(nm)
            lines: list[str] = []
            for chrom, c in cand.groupby("CHROM", sort=False):
                # closed overlap needs START <= B and END >= A; query [A - 1, B
                # + 1) (half-open, 0-based)
                w = sorted(
                    zip((c["START"].to_numpy(np.int64) - 1).clip(0).tolist(),
                        (c["END"].to_numpy(np.int64) + 1).tolist()))
                merged: list[list[int]] = []
                for lo, hi in w:
                    if merged and lo <= merged[-1][1]:
                        merged[-1][1] = max(merged[-1][1], hi)
                    else:
                        merged.append([lo, hi])
                for contig in by_norm.get(chrom, []):
                    for lo, hi in merged:
                        lines.extend(tbx.fetch(contig, lo, hi))
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    if not lines:
        return pd.DataFrame({
            c: pd.Series(dtype="object" if c == "CHROM" else "float64")
            for c in atlas_io.ATLAS_COLS
        })
    df = pd.read_csv(StringIO("\n".join(lines) + "\n"),
                     sep="\t",
                     header=None,
                     names=atlas_io.ATLAS_COLS[:len(lines[0].split("\t"))],
                     dtype={"CHROM": str},
                     na_values=["NA"],
                     float_precision="round_trip")
    # a longer record could be returned by two windows
    if (df["END"] - df["START"] != 1).any():
        return None
    return df


def _windows_scan(path: str, cand: pd.DataFrame, ref_lo, ref_hi, min_samples):
    """Scans the whole file with the parallel reader.

    The rows are already filtered (chromosome, N, REF_METH, windows).
    """
    code = cand["CHROM"].map({
        n: i for i, n in enumerate(cpg_islands.NAMES[:22])
    }).to_numpy(np.int64)
    win = atlas_io.candidate_windows(code, cand["START"].to_numpy(np.int64),
                                     cand["END"].to_numpy(np.int64))
    code, start, end, _phi, meth, _cnt = atlas_io.read_columns(
        path,
        atlas_io.ICR,
        max_code=21,
        min_n=min_samples,
        strict=False,
        meth_lo=ref_lo,
        meth_hi=ref_hi,
        windows=win)
    return pd.DataFrame({
        "CHROM": np.array(cpg_islands.NAMES, dtype=object)[code],
        "START": start,
        "END": end,
        "REF_METH": meth
    })


def detect(atlas: atlas_io.LazyAtlas | pd.DataFrame | str,
           cand: pd.DataFrame,
           *,
           ref_lo: float = settings.ICR_REF_METH_MIN,
           ref_hi: float = settings.ICR_REF_METH_MAX,
           min_sites: int = settings.ICR_MIN_SITES,
           merge_distance: int = settings.DMR_MERGE_DISTANCE,
           min_samples: int = settings.CPG_MIN_SAMPLES_PER_SITE,
           v1_duplicates: bool = settings.ICR_V1_DUPLICATES,
           index: str = "auto",
           log: Callable[[str], None] = print) -> pd.DataFrame:
    """Detects the ICRs of the autosomal scope (v2.2 DETECT_ICR_REGIONS).

    Args:
        atlas: DataFrame with ATLAS_COLS, or a LazyAtlas or path (then only the
            candidate windows are read).
        cand: Candidates from `read_candidates`.
        ref_lo: Smallest REF_METH of a hemi-methylated CpG.
        ref_hi: Largest REF_METH.
        min_sites: Minimum CpGs of an ICR.
        merge_distance: Largest gap (bp) inside one ICR.
        min_samples: Sites need at least this many controls.
        v1_duplicates: Count a site once per overlapping candidate (v2.2).
        index: How a path is read. "auto" uses tabix windows when the file has a
            .tbi, else the parallel scan; "build" uses tabix windows and builds
                the
            index in scratch when missing; "scan" scans in parallel only.
        log: Receives progress messages.

    Returns:
        The ICRs (OUT_COLS).
    """
    if len(cand) == 0:
        log("No autosomal imprinting gene candidates. Empty output.")
        return _empty()
    path = atlas_io.unread_path(atlas)
    a = None
    if path is not None:
        no_jit = os.environ.get("DMRDETECTOR_NO_JIT") == "1"
        if index in ("build",
                     "auto") and (index == "build" or
                                  os.path.exists(path + ".tbi")) and not no_jit:
            a = _windows_fetch(path, cand, None)
        if a is None and not no_jit:
            try:
                a = _windows_scan(path, cand, ref_lo, ref_hi, min_samples)
                # the reader applied the N and REF_METH filters
                a["N"] = float(min_samples)
            except (atlas_io.Unsupported, ValueError, OSError, zlib.error):
                a = None
        if a is None:
            a = atlas_io.read_atlas_pandas(path)
    elif isinstance(atlas, atlas_io.LazyAtlas):
        a = atlas.frame()
    else:
        a = atlas
    a = a[["CHROM", "START", "END", "N", "REF_METH"]].copy()
    a["CHROM"] = bed.normalize_primary_chrom(a["CHROM"])
    a = a[a["CHROM"].isin(bed.AUTOSOMES)]
    a = a[a["N"] >= min_samples]  # >=, as v2.2 (CpG panel uses >)
    if len(a) == 0:
        log("No autosomal atlas sites after filtering. Empty output.")
        return _empty()
    # full table: keep only sites in a candidate window (as the readers do)
    if len(a) > 100_000:
        keep = np.zeros(len(a), bool)
        ch = a["CHROM"].to_numpy()
        for chrom in cand["CHROM"].unique():
            m = ch == chrom
            keep[m] = _overlaps_candidates(chrom,
                                           a["START"].to_numpy(np.int64)[m],
                                           a["END"].to_numpy(np.int64)[m], cand)
        a = a[keep]
    return _core(a, cand, ref_lo, ref_hi, min_sites, merge_distance,
                 v1_duplicates, log)


def main(args: argparse.Namespace) -> int:
    """Runs the `icr` command."""
    log = bed.warn
    log("=== Detect Imprinted Control Regions (ICR) ===")
    cand = read_candidates(args.icr_gene_list)
    out = detect(args.ref_atlas,
                 cand,
                 ref_lo=args.icr_ref_meth_min,
                 ref_hi=args.icr_ref_meth_max,
                 min_sites=args.icr_min_sites,
                 merge_distance=args.dmr_merge_distance,
                 min_samples=args.min_samples_per_site,
                 v1_duplicates=not args.no_v1_duplicates,
                 index=args.index,
                 log=log)
    out.to_csv(args.out, sep="\t", index=False)
    return 0


def add_parser(
        subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Adds the `icr` command to `subparsers`."""
    p = subparsers.add_parser(
        "icr", help="detect the ICR panel (v2.2 DETECT_ICR_REGIONS)")
    p.add_argument("--ref-atlas",
                   required=True,
                   help="autosomal narrow-sigma atlas (headerless BED, .gz ok)")
    p.add_argument("--icr-gene-list",
                   required=True,
                   help="imprinting gene BED (assets/icr_gene_list.bed)")
    p.add_argument("--out", default="icr_region.bed")
    p.add_argument("--icr-ref-meth-min",
                   type=float,
                   default=settings.ICR_REF_METH_MIN)
    p.add_argument("--icr-ref-meth-max",
                   type=float,
                   default=settings.ICR_REF_METH_MAX)
    p.add_argument("--icr-min-sites", type=int, default=settings.ICR_MIN_SITES)
    p.add_argument("--dmr-merge-distance",
                   type=int,
                   default=settings.DMR_MERGE_DISTANCE)
    p.add_argument("--min-samples-per-site",
                   type=int,
                   default=settings.CPG_MIN_SAMPLES_PER_SITE,
                   help="cpg_min_samples_per_site; atlas sites need N >= this")
    p.add_argument(
        "--no-v1-duplicates",
        action="store_true",
        default=not settings.ICR_V1_DUPLICATES,
        help=
        "count a site once even when it lies in several overlapping candidate "
        "intervals (v2.2 counts it once per interval; default reproduces v2.2)")
    p.add_argument(
        "--index",
        choices=["auto", "build", "scan"],
        default="auto",
        help=
        "atlas reading: auto = tabix windows if the atlas has a .tbi, else a "
        "parallel scan; "
        "build = tabix windows, indexing in a scratch directory when needed "
        "(the input is "
        "never modified); scan = parallel scan")
    p.set_defaults(func=main)
    return p
