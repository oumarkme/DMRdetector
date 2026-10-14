"""DETECT_CPG_ISLANDS (v2.2): CpG island panel from the narrow-sigma atlas.

The output has CHROM START END START_EXTEND END_EXTEND LENGTH NUM_CPG, with the
atlas START values as given (no 0/1-base conversion). The atlas is read by
`atlas_io`; `_chrom_kernel` classifies CORE and SHORE sites, clusters by gap,
filters by length and density and extends the shores in one pass.
DMRDETECTOR_NO_JIT=1 runs the kernel as plain Python and reads with pandas.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
import os
import zlib

import numpy as np
import pandas as pd

from dmrdetector import settings
from dmrdetector.fileio import atlas_io
from dmrdetector.fileio import bed
from dmrdetector.util.jit import jit


OUT_COLS = [
    "CHROM", "START", "END", "START_EXTEND", "END_EXTEND", "LENGTH", "NUM_CPG"
]
# chromosome code (0..23) -> name; RANK = position in the C-locale (byte) sort
# of the names
NAMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
_BY_NAME = {n: i for i, n in enumerate(NAMES)}
# string order = C locale
RANK = np.argsort(np.argsort(np.array(NAMES, dtype=object))).astype(np.int8)
_BIG = np.int64(1) << 60


def _empty() -> pd.DataFrame:
    return pd.DataFrame({
        c: pd.Series(dtype="object" if c == "CHROM" else "int64")
        for c in OUT_COLS
    })


@jit
def _is_sorted(rank, start):
    for i in range(1, rank.size):
        if rank[i] < rank[i - 1] or (rank[i] == rank[i - 1] and
                                     start[i] < start[i - 1]):
            return False
    return True


@jit
def _chrom_kernel(start, meth, phi, female, core_max_meth, core_max_phi,
                  shore_max_meth, female_min, female_max, max_gap, min_length):
    """One chromosome (sites sorted by START).

    Returns (n_islands, START, END, START_EXTEND, END_EXTEND, NUM_CPG arrays,
    n_core, n_shore).

    CORE: REF_METH < core_max_meth and PHI < core_max_phi (female_X: female_min
    <= REF_METH <= female_max). SHORE: not CORE, REF_METH < shore_max_meth, no
    NA, and within max_gap of the nearest CORE (distance 0 for CORE). Islands:
    CORE sites clustered while the gap to the previous CORE is <= max_gap; kept
    when LENGTH = END - START + 1 >= min_length and NUM_CPG * 100 >= LENGTH
    (female_X: length only). Extension: first and last SHORE in [START -
    max_gap, END + max_gap].
    """
    n = start.size
    core = np.zeros(n, np.bool_)
    nc = 0
    for i in range(n):
        if female:
            c = meth[i] >= female_min and meth[i] <= female_max
        else:
            c = meth[i] < core_max_meth and phi[i] < core_max_phi
        if c:
            core[i] = True
            nc += 1
    cap = max(nc, 1)
    o_s = np.empty(cap, np.int64)
    o_e = np.empty(cap, np.int64)
    o_xs = np.empty(cap, np.int64)
    o_xe = np.empty(cap, np.int64)
    o_n = np.empty(cap, np.int64)
    ns = 0
    sh = np.empty(1, np.int64)
    if not female and nc > 0:
        dist = np.empty(n, np.int64)
        prev = -_BIG
        for i in range(n):
            if core[i]:
                prev = start[i]
            dist[i] = start[i] - prev if prev != -_BIG else _BIG
        nxt = _BIG
        for i in range(n - 1, -1, -1):
            if core[i]:
                nxt = start[i]
            if nxt != _BIG:
                d = nxt - start[i]
                if d < dist[i]:
                    dist[i] = d
        sh = np.empty(n, np.int64)
        for i in range(n):
            if ((not core[i]) and meth[i] < shore_max_meth and
                    dist[i] <= max_gap and
                    not (np.isnan(meth[i]) or np.isnan(phi[i]))):
                sh[ns] = start[i]
                ns += 1
    ni = 0
    cnt = 0
    first = np.int64(0)
    last = np.int64(0)
    for i in range(n + 1):
        close = False
        if i == n:
            close = cnt > 0
        elif core[i]:
            if cnt == 0:
                first = start[i]
                last = start[i]
                cnt = 1
            elif start[i] - last > max_gap:
                close = True
            else:
                last = start[i]
                cnt += 1
        if close:
            length = last - first + 1
            if female:
                ok = length >= min_length
            else:
                ok = length >= min_length and cnt * 100 >= length
            if ok:
                o_s[ni] = first
                o_e[ni] = last
                o_n[ni] = cnt
                xs = first
                xe = last
                if ns > 0:
                    lo = np.searchsorted(sh[:ns], first - max_gap, side="left")
                    hi = np.searchsorted(sh[:ns], last + max_gap, side="right")
                    if hi > lo:
                        xs = min(first, sh[lo])
                        xe = max(last, sh[hi - 1])
                o_xs[ni] = xs
                o_xe[ni] = xe
                ni += 1
            if i < n:  # the closing site starts the next cluster
                first = start[i]
                last = start[i]
                cnt = 1
            else:
                cnt = 0
    return ni, o_s, o_e, o_xs, o_xe, o_n, nc, ns


def _detect_arrays(code, start, meth, phi, scope, *, core_max_meth,
                   core_max_phi, shore_max_meth, female_min_meth,
                   female_max_meth, max_gap, min_length, log):
    """Returns the island table of the sites that passed the filters.

    `code` indexes NAMES; the sites passed the chromosome and N filters.
    """
    rank = RANK[code]
    # data.table key order: C-locale CHROM, then START; stable
    if not _is_sorted(rank, start):
        o = np.lexsort((start, rank))
        rank, start, meth, phi = rank[o], start[o], meth[o], phi[o]
    female = scope == "female_X"
    bounds = np.r_[0, np.flatnonzero(rank[1:] != rank[:-1]) + 1,
                   rank.size] if rank.size else np.zeros(1, np.int64)
    res, n_core, n_shore = [], 0, 0
    if female:
        log(f"  Female_X hemi-methylated model: REF_METH between "
            f"{female_min_meth} and {female_max_meth}")
    for b0, b1 in zip(bounds[:-1], bounds[1:]):
        r = _chrom_kernel(start[b0:b1], meth[b0:b1], phi[b0:b1], female,
                          core_max_meth, core_max_phi, shore_max_meth,
                          female_min_meth, female_max_meth, max_gap, min_length)
        ni = r[0]
        n_core += r[6]
        n_shore += r[7]
        res.append((int(rank[b0]), ni, r))
    if female:
        log(f"  Hemi-methylated sites: {n_core}")
        if n_core == 0:
            log("No hemi-methylated sites. Empty output.")
            return _empty()
    else:
        log(f"  CORE sites: {n_core}")
        if n_core == 0:
            log("No CORE sites. Empty output.")
            return _empty()
        log(f"  SHORE sites: {n_shore}")
    total = sum(x[1] for x in res)
    log(f"  Islands after filtering: {total}")
    if total == 0:
        return _empty()
    name_of_rank = np.empty(len(NAMES), object)
    name_of_rank[RANK] = NAMES
    chrom, cols = [], [[] for _ in range(5)]
    for rk, ni, r in res:
        if ni == 0:
            continue
        chrom.append(np.full(ni, name_of_rank[rk], dtype=object))
        for c, arr in zip(cols, r[1:6]):
            c.append(arr[:ni])
    s, e, xs, xe, nn = (np.concatenate(c) for c in cols)
    out = pd.DataFrame({
        "CHROM": np.concatenate(chrom),
        "START": s,
        "END": e,
        "START_EXTEND": xs,
        "END_EXTEND": xe,
        "LENGTH": e - s + 1,
        "NUM_CPG": nn
    })
    log(f"Total islands: {len(out)}")
    return out[OUT_COLS]


def detect(atlas: atlas_io.LazyAtlas | pd.DataFrame | str,
           scope: str = "autosomal",
           *,
           include_sex: bool = settings.INCLUDE_SEX_CHROMOSOMES,
           min_samples: int = settings.CPG_MIN_SAMPLES_PER_SITE,
           core_max_meth: float = settings.CPG_CORE_MAX_METHYLATION,
           core_max_phi: float = settings.CPG_CORE_MAX_DISPERSION,
           shore_max_meth: float = settings.CPG_SHORE_MAX_METHYLATION,
           female_min_meth: float = settings.SEX_F_REF_METH_MIN,
           female_max_meth: float = settings.SEX_F_REF_METH_MAX,
           max_gap: int = settings.CPG_MAX_GAP_DISTANCE,
           min_length: int = settings.CPG_MIN_ISLAND_LENGTH,
           log: Callable[[str], None] = print) -> pd.DataFrame:
    """Detects the CpG islands of one scope (v2.2 DETECT_CPG_ISLANDS).

    Args:
        atlas: DataFrame with ATLAS_COLS, or a LazyAtlas or path, which is read
            with the fast reader.
        scope: autosomal, male_X, female_X or Y.
        include_sex: Keep chrX and chrY.
        min_samples: Sites need more than this many controls.
        core_max_meth: CORE: REF_METH below this.
        core_max_phi: CORE: PHI below this.
        shore_max_meth: SHORE: REF_METH below this.
        female_min_meth: female_X: REF_METH at least this.
        female_max_meth: female_X: REF_METH at most this.
        max_gap: Largest gap (bp) inside one island.
        min_length: Shortest island (bp).
        log: Receives progress messages.

    Returns:
        The islands (OUT_COLS).
    """
    kw = dict(core_max_meth=core_max_meth,
              core_max_phi=core_max_phi,
              shore_max_meth=shore_max_meth,
              female_min_meth=female_min_meth,
              female_max_meth=female_max_meth,
              max_gap=max_gap,
              min_length=min_length,
              log=log)
    path = atlas_io.unread_path(atlas)
    max_code = 23 if include_sex else 21
    if path is not None and os.environ.get("DMRDETECTOR_NO_JIT") != "1":
        try:
            # sites that can be neither CORE nor SHORE never influence the
            # result: drop them while parsing
            mf = ((female_min_meth, female_max_meth,
                   False) if scope == "female_X" else
                  (-np.inf, max(core_max_meth, shore_max_meth), True))
            code, start, _end, phi, meth, cnt = atlas_io.read_columns(
                path,
                atlas_io.CPG,
                max_code=max_code,
                min_n=min_samples,
                strict=True,
                meth_lo=mf[0],
                meth_hi=mf[1],
                meth_hi_strict=mf[2])
        except (atlas_io.Unsupported, ValueError, OSError, zlib.error):
            code = None
        if code is not None:
            log(f"  Loaded {int(cnt[0])} sites")
            log(f"  After chromosome filter (include_sex={include_sex}): "
                f"{int(cnt[1])} sites")
            log(f"  After N > {min_samples} filter: {int(cnt[2])} sites")
            return _detect_arrays(code, start, meth, phi, scope, **kw)
    if path is not None:
        atlas = atlas_io.read_atlas_pandas(path)
    elif isinstance(atlas, atlas_io.LazyAtlas):
        atlas = atlas.frame()
    a = atlas[["CHROM", "START", "N", "PHI", "REF_METH"]]
    log(f"  Loaded {len(a)} sites")
    cat = pd.Categorical(a["CHROM"].astype(str))
    lut = np.array([_BY_NAME.get(c, -1) for c in cat.categories] + [-1],
                   np.int16)
    code = lut[np.where(cat.codes < 0, len(cat.categories), cat.codes)]
    keep = (code >= 0) & (code <= max_code)
    log(f"  After chromosome filter (include_sex={include_sex}): "
        f"{int(keep.sum())} sites")
    keep &= (a["N"].to_numpy(float) > min_samples)  # strict, as v2.2
    log(f"  After N > {min_samples} filter: {int(keep.sum())} sites")
    return _detect_arrays(code[keep].astype(np.int8),
                          a["START"].to_numpy(np.int64)[keep],
                          a["REF_METH"].to_numpy(float)[keep],
                          a["PHI"].to_numpy(float)[keep], scope, **kw)


def main(args: argparse.Namespace) -> int:
    """Runs the `cpg-islands` command."""
    log = bed.warn
    log("=== Detect CpG Island Panel ===")
    res = detect(args.ref_atlas,
                 args.scope,
                 include_sex=args.include_sex_chromosomes,
                 min_samples=args.min_samples_per_site,
                 core_max_meth=args.core_max_methylation,
                 core_max_phi=args.core_max_dispersion,
                 shore_max_meth=args.shore_max_methylation,
                 female_min_meth=args.sex_f_ref_meth_min,
                 female_max_meth=args.sex_f_ref_meth_max,
                 max_gap=args.max_gap_distance,
                 min_length=args.min_island_length,
                 log=log)
    res.to_csv(args.out, sep="\t", index=False)
    return 0


def add_parser(
        subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Adds the `cpg-islands` command to `subparsers`."""
    p = subparsers.add_parser(
        "cpg-islands",
        help="detect the CpG island panel (v2.2 DETECT_CPG_ISLANDS)")
    p.add_argument("--ref-atlas",
                   required=True,
                   help="narrow-sigma smoothed atlas (headerless BED, .gz ok)")
    p.add_argument("--scope",
                   default="autosomal",
                   choices=["autosomal", "male_X", "female_X", "Y"])
    p.add_argument("--out",
                   required=True,
                   help="e.g. cpg_signal.bed or cpg_signal_<scope>.bed")
    p.add_argument("--include-sex-chromosomes",
                   action=argparse.BooleanOptionalAction,
                   default=settings.INCLUDE_SEX_CHROMOSOMES)
    p.add_argument("--min-samples-per-site",
                   type=int,
                   default=settings.CPG_MIN_SAMPLES_PER_SITE)
    p.add_argument("--core-max-methylation",
                   type=float,
                   default=settings.CPG_CORE_MAX_METHYLATION)
    p.add_argument("--core-max-dispersion",
                   type=float,
                   default=settings.CPG_CORE_MAX_DISPERSION)
    p.add_argument("--shore-max-methylation",
                   type=float,
                   default=settings.CPG_SHORE_MAX_METHYLATION)
    p.add_argument("--sex-f-ref-meth-min",
                   type=float,
                   default=settings.SEX_F_REF_METH_MIN)
    p.add_argument("--sex-f-ref-meth-max",
                   type=float,
                   default=settings.SEX_F_REF_METH_MAX)
    p.add_argument("--max-gap-distance",
                   type=int,
                   default=settings.CPG_MAX_GAP_DISTANCE)
    p.add_argument("--min-island-length",
                   type=int,
                   default=settings.CPG_MIN_ISLAND_LENGTH)
    p.set_defaults(func=main)
    return p
