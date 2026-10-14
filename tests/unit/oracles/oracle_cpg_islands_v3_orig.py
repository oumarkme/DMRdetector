"""DETECT_CPG_ISLANDS: CpG island panels from the narrow-sigma smoothed atlas.

Same algorithm, defaults, and output as v1 (main.nf). Output is a tab-separated
file with header CHROM START END START_EXTEND END_EXTEND LENGTH NUM_CPG.
Coordinates are the atlas START values as given (no 0/1-base conversion).
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

ATLAS_COLS = ["CHROM", "START", "END", "N", "COV", "PHI", "RAW_METH", "REF_METH", "REF_VAR"]
OUT_COLS = ["CHROM", "START", "END", "START_EXTEND", "END_EXTEND", "LENGTH", "NUM_CPG"]


def read_atlas(path: str) -> pd.DataFrame:
    """Headerless (optionally gzip/bgzip) 9-column atlas BED."""
    try:  # prefer the shared reader once io.py provides one
        from .io import read_atlas as _ra  # type: ignore
        return _ra(path)
    except (ImportError, AttributeError):
        pass
    df = pd.read_csv(path, sep="\t", header=None, names=ATLAS_COLS,
                     dtype={"CHROM": str}, na_values=["NA"])
    return df


def _empty() -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype="object" if c == "CHROM" else "int64")
                         for c in OUT_COLS})


def _cluster_ids(start: np.ndarray, max_gap: int) -> np.ndarray:
    """cumsum(d > max_gap) with d = START - shift(START, fill = START[1])."""
    d = np.diff(start, prepend=start[0])
    return np.cumsum(d > max_gap)


def _islands(core: pd.DataFrame, max_gap: int) -> pd.DataFrame:
    """Cluster CORE sites per chromosome; rows: CHROM, island_id, START, END, NUM_CPG."""
    rows = []
    for chrom, g in core.groupby("CHROM", sort=False):
        s = g["START"].to_numpy(np.int64)
        iid = _cluster_ids(s, max_gap)
        brk = np.flatnonzero(np.diff(iid, prepend=-1) != 0)
        ends = np.append(brk[1:], s.size)
        for b, e in zip(brk, ends):
            seg = s[b:e]
            rows.append((chrom, int(iid[b]), int(seg.min()), int(seg.max()), e - b))
    return pd.DataFrame(rows, columns=["CHROM", "island_id", "START", "END", "NUM_CPG"])


def detect(atlas: pd.DataFrame, scope: str = "autosomal", *, include_sex: bool = True,
           min_samples: int = 5, core_max_meth: float = 0.1, core_max_phi: float = 0.05,
           shore_max_meth: float = 0.3, female_min_meth: float = 0.40,
           female_max_meth: float = 0.60, max_gap: int = 150,
           min_length: int = 150, log=print) -> pd.DataFrame:
    a = atlas[["CHROM", "START", "N", "PHI", "REF_METH"]].copy()
    a["CHROM"] = a["CHROM"].astype(str)
    log(f"  Loaded {len(a)} sites")
    allowed = [f"chr{i}" for i in range(1, 23)] + (["chrX", "chrY"] if include_sex else [])
    a = a[a["CHROM"].isin(allowed)]
    log(f"  After chromosome filter (include_sex={include_sex}): {len(a)} sites")
    a = a[a["N"] > min_samples]  # strict, as v1
    log(f"  After N > {min_samples} filter: {len(a)} sites")
    # data.table key order: C-locale CHROM, then START
    a = a.sort_values(["CHROM", "START"], kind="stable").reset_index(drop=True)
    start = a["START"].to_numpy(np.int64)
    meth = a["REF_METH"].to_numpy(float)

    if scope == "female_X":
        log(f"  Female_X hemi-methylated model: REF_METH between {female_min_meth} and {female_max_meth}")
        core = (meth >= female_min_meth) & (meth <= female_max_meth)  # NaN -> False
        log(f"  Hemi-methylated sites: {int(core.sum())}")
        if core.sum() == 0:
            log("No hemi-methylated sites. Empty output.")
            return _empty()
        isl = _islands(a[core], max_gap)
        isl["LENGTH"] = isl["END"] - isl["START"] + 1
        isl = isl[isl["LENGTH"] >= min_length]
        log(f"  Islands after filtering: {len(isl)}")
        if len(isl) == 0:
            return _empty()
        isl["START_EXTEND"] = isl["START"]
        isl["END_EXTEND"] = isl["END"]
    else:
        phi = a["PHI"].to_numpy(float)
        core_na = np.isnan(meth) | np.isnan(phi)  # R: NA CORE
        core = (meth < core_max_meth) & (phi < core_max_phi)
        log(f"  CORE sites: {int(core.sum())}")
        if core.sum() == 0:
            log("No CORE sites. Empty output.")
            return _empty()
        # DIST to nearest CORE on the same chromosome (rolling join "nearest")
        dist = np.full(len(a), np.nan)
        chrom = a["CHROM"].to_numpy()
        bounds = np.flatnonzero(np.r_[True, chrom[1:] != chrom[:-1], True])
        for b0, b1 in zip(bounds[:-1], bounds[1:]):
            cs = start[b0:b1][core[b0:b1]]
            if cs.size == 0:
                continue
            x = start[b0:b1]
            i = np.searchsorted(cs, x)
            left = cs[np.clip(i - 1, 0, cs.size - 1)]
            right = cs[np.clip(i, 0, cs.size - 1)]
            dist[b0:b1] = np.minimum(np.abs(x - left), np.abs(x - right))
        with np.errstate(invalid="ignore"):
            shore = (meth < shore_max_meth) & (dist <= max_gap) & ~core & ~core_na
        log(f"  SHORE sites: {int(shore.sum())}")
        isl = _islands(a[core], max_gap)
        isl["LENGTH"] = isl["END"] - isl["START"] + 1
        isl = isl[(isl["LENGTH"] >= min_length) & (isl["NUM_CPG"] * 100 >= isl["LENGTH"])].copy()
        log(f"  Islands after filtering: {len(isl)}")
        if len(isl) == 0:
            return _empty()
        # extend with SHORE sites inside [START - max_gap, END + max_gap] (closed)
        s_ext = isl["START"].to_numpy(np.int64).copy()
        e_ext = isl["END"].to_numpy(np.int64).copy()
        sh_chrom = chrom[shore]
        sh_start = start[shore]
        for k, (c, s0, e0) in enumerate(zip(isl["CHROM"], isl["START"], isl["END"])):
            m = sh_chrom == c  # shores are sorted within chromosome
            ss = sh_start[m]
            lo = np.searchsorted(ss, s0 - max_gap, side="left")
            hi = np.searchsorted(ss, e0 + max_gap, side="right")
            if hi > lo:
                s_ext[k] = min(s0, ss[lo])
                e_ext[k] = max(e0, ss[hi - 1])
        isl["START_EXTEND"] = s_ext
        isl["END_EXTEND"] = e_ext

    res = isl[OUT_COLS].sort_values(["CHROM", "START"], kind="stable").reset_index(drop=True)
    log(f"Total islands: {len(res)}")
    return res


def main(args) -> int:
    log = lambda m: print(m, file=sys.stderr)
    log("=== Detect CpG Island Panel ===")
    atlas = read_atlas(args.ref_atlas)
    res = detect(atlas, args.scope, include_sex=args.include_sex_chromosomes,
                 min_samples=args.min_samples_per_site, core_max_meth=args.core_max_methylation,
                 core_max_phi=args.core_max_dispersion, shore_max_meth=args.shore_max_methylation,
                 female_min_meth=args.sex_f_ref_meth_min, female_max_meth=args.sex_f_ref_meth_max,
                 max_gap=args.max_gap_distance, min_length=args.min_island_length, log=log)
    res.to_csv(args.out, sep="\t", index=False)
    return 0


def add_parser(subparsers):
    p = subparsers.add_parser("cpg-islands", help="detect the CpG island panel (v1 DETECT_CPG_ISLANDS)")
    p.add_argument("--ref-atlas", required=True, help="narrow-sigma smoothed atlas (headerless BED, .gz ok)")
    p.add_argument("--scope", default="autosomal",
                   choices=["autosomal", "male_X", "female_X", "Y"])
    p.add_argument("--out", required=True, help="e.g. cpg_signal.bed or cpg_signal_<scope>.bed")
    p.add_argument("--include-sex-chromosomes", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--min-samples-per-site", type=int, default=5)
    p.add_argument("--core-max-methylation", type=float, default=0.1)
    p.add_argument("--core-max-dispersion", type=float, default=0.05)
    p.add_argument("--shore-max-methylation", type=float, default=0.3)
    p.add_argument("--sex-f-ref-meth-min", type=float, default=0.40)
    p.add_argument("--sex-f-ref-meth-max", type=float, default=0.60)
    p.add_argument("--max-gap-distance", type=int, default=150)
    p.add_argument("--min-island-length", type=int, default=150)
    p.set_defaults(func=main)
    return p

