"""DETECT_ICR_REGIONS: imprinted control region panel (autosomal only).

Same algorithm, defaults, and output as v1 (main.nf). Output header:
CHROM START END GENE ICR_ID NUM_CPG. END = last hemi-methylated CpG START + 1.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from oracle_cpg_islands_v3_orig import ATLAS_COLS, read_atlas

OUT_COLS = ["CHROM", "START", "END", "GENE", "ICR_ID", "NUM_CPG"]
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def _empty() -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype="object" if c in ("CHROM", "GENE", "ICR_ID") else "int64")
                         for c in OUT_COLS})


def normalize_primary_chrom(x: pd.Series) -> pd.Series:
    chrom = x.astype(str)
    base = chrom.str.replace(r"^chr", "", case=False, regex=True).str.upper()
    primary = base.isin([str(i) for i in range(1, 23)] + ["X", "Y"])
    return chrom.where(~primary, "chr" + base)


def read_candidates(path: str) -> pd.DataFrame:
    """First 4 columns as CHROM START END GENE; header auto-detected like fread."""
    with open(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")
    try:
        float(first[1])
        header = None
    except (ValueError, IndexError):
        header = 0
    df = pd.read_csv(path, sep="\t", header=header, dtype=str)
    df = df.iloc[:, : min(4, df.shape[1])]
    df.columns = ["CHROM", "START", "END", "GENE"][: df.shape[1]]
    df["CHROM"] = normalize_primary_chrom(df["CHROM"])
    df = df[df["CHROM"].isin(AUTOSOMES)].copy()
    df["START"] = df["START"].astype(float).astype(np.int64)
    df["END"] = df["END"].astype(float).astype(np.int64)
    return df.reset_index(drop=True)


def detect(atlas: pd.DataFrame, cand: pd.DataFrame, *, ref_lo: float = 0.40, ref_hi: float = 0.60,
           min_sites: int = 5, merge_distance: int = 150, min_samples: int = 5,
           v1_duplicates: bool = True, log=print) -> pd.DataFrame:
    if len(cand) == 0:
        log("No autosomal imprinting gene candidates. Empty output.")
        return _empty()
    a = atlas[["CHROM", "START", "END", "N", "REF_METH"]].copy()
    a["CHROM"] = normalize_primary_chrom(a["CHROM"])
    a = a[a["CHROM"].isin(AUTOSOMES)]
    a = a[a["N"] >= min_samples]  # >=, as v1 (CpG panel uses >)
    if len(a) == 0:
        log("No autosomal atlas sites after filtering. Empty output.")
        return _empty()

    parts = []
    for chrom, c in cand.groupby("CHROM", sort=False):
        g = a[a["CHROM"] == chrom]
        if g.empty:
            continue
        s = g["START"].to_numpy(np.int64)
        e = g["END"].to_numpy(np.int64)
        # v1 pre-filters with `tabix -R` (BED, half-open), which returns a
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
                parts.append(pd.DataFrame({"CHROM": chrom, "START": np.repeat(s[m], mult[m]),
                                           "GENE": gene}))
    if not parts:
        log("  Hemi-methylated ICR CpGs: 0")
        return _empty()
    hemi = pd.concat(parts, ignore_index=True)
    log(f"  Hemi-methylated ICR CpGs: {len(hemi)}")
    hemi = hemi.sort_values(["GENE", "CHROM", "START"], kind="stable").reset_index(drop=True)
    gap = hemi.groupby(["GENE", "CHROM"], sort=False)["START"].diff().fillna(np.inf)
    hemi["cluster"] = (gap > merge_distance).astype(int).groupby(
        [hemi["GENE"], hemi["CHROM"]], sort=False).cumsum()
    icr = hemi.groupby(["GENE", "CHROM", "cluster"], sort=False).agg(
        START=("START", "min"), END=("START", "max"), NUM_CPG=("START", "size")).reset_index()
    icr["END"] = icr["END"] + 1
    icr = icr[icr["NUM_CPG"] >= min_sites]
    if icr.empty:
        log("No ICR passed min_sites. Empty output.")
        return _empty()
    icr = icr.sort_values(["CHROM", "START", "GENE"], kind="stable").reset_index(drop=True)
    idx = icr.groupby(["GENE", "CHROM"], sort=False).cumcount() + 1
    icr["ICR_ID"] = icr["GENE"] + "_" + icr["CHROM"] + "_icr" + idx.map("{:02d}".format)
    out = icr[OUT_COLS]
    log(f"Total ICRs: {len(out)} across {out['GENE'].nunique()} genes")
    return out


def main(args) -> int:
    log = lambda m: print(m, file=sys.stderr)
    log("=== Detect Imprinted Control Regions (ICR) ===")
    cand = read_candidates(args.icr_gene_list)
    atlas = read_atlas(args.ref_atlas)
    out = detect(atlas, cand, ref_lo=args.icr_ref_meth_min, ref_hi=args.icr_ref_meth_max,
                 min_sites=args.icr_min_sites, merge_distance=args.dmr_merge_distance,
                 min_samples=args.min_samples_per_site,
                 v1_duplicates=not args.no_v1_duplicates, log=log)
    out.to_csv(args.out, sep="\t", index=False)
    return 0


def add_parser(subparsers):
    p = subparsers.add_parser("icr", help="detect the ICR panel (v1 DETECT_ICR_REGIONS)")
    p.add_argument("--ref-atlas", required=True, help="autosomal narrow-sigma atlas (headerless BED, .gz ok)")
    p.add_argument("--icr-gene-list", required=True, help="imprinting gene BED (assets/icr_gene_list.bed)")
    p.add_argument("--out", default="icr_region.bed")
    p.add_argument("--icr-ref-meth-min", type=float, default=0.40)
    p.add_argument("--icr-ref-meth-max", type=float, default=0.60)
    p.add_argument("--icr-min-sites", type=int, default=5)
    p.add_argument("--dmr-merge-distance", type=int, default=150)
    p.add_argument("--min-samples-per-site", type=int, default=5,
                   help="cpg_min_samples_per_site; atlas sites need N >= this")
    p.add_argument("--no-v1-duplicates", action="store_true",
                   help="count a site once even when it lies in several overlapping candidate "
                        "intervals (v1 counts it once per interval; default reproduces v1)")
    p.set_defaults(func=main)
    return p
