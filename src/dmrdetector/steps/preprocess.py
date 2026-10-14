"""PREPROCESS_SAMPLE and SPLIT_SEX_COMPONENTS (v2.2), input checks, the pool.

`filter` drops the CpGs on a heterozygous SNP of the sample VCF and keeps rows
with coverage > min_coverage; `split` writes the autosomal (with PAR1 and PAR2,
chrY PAR rows folded onto chrX), non-PAR chrX and non-PAR chrY files. Names
match as exact strings and intervals are 0-based half-open, as in v2.2.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from collections.abc import Iterable
import gzip
import os
import re
import sys
import time

import numpy as np
import pandas as pd

from dmrdetector import data
from dmrdetector import runtime
from dmrdetector import settings
from dmrdetector.fileio import bed
from dmrdetector.fileio import methbed
from dmrdetector.fileio import vcf as vcf_io


# chrY PAR2 (56,887,902-57,217,415) -> chrX PAR2 (155,701,382-156,030,895)
PAR2_SHIFT = 98_813_480
_AUTO_RE = re.compile(r"^(chr)?[0-9]+$")
_X_RE = re.compile(r"^(chr)?X$")
_Y_RE = re.compile(r"^(chr)?Y$")

# hg38 / GRCh38 primary assembly lengths (bp) of chr1-22, X and Y (UCSC
# hg38.chrom.sizes, Ensembl GRCh38)
HG38_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}
# bp by which an END may exceed the contig length before the file is rejected.
# Valid hg38 data never exceeds the length; the margin only keeps a stray
# off-by-one from stopping a run.
END_MARGIN = 1000

# Lab convention for calling the sex from the CpG records on non-PAR chrY: a
# male has at least this many (our ONT whole-blood controls: males 68,408 to
# 171,117, females 4,928 to 8,570 raw records). `run` counts them after the
# coverage filter, which lowers a male count and never raises a female one, so
# only a female sample above the cutoff is reported.
CHRY_MALE_MIN_CPG = 50_000

# Lab convention for calling the sex from the heterozygous SNPs on non-PAR chrX
# (atlas selection, doc/sex_chromosomes.md): a female has at least this many, a
# male fewer. Our ONT controls: females 60,362 to 83,844, males 6,532 to 19,199
# heterozygous SNVs by the lab's counting script. The count here is taken from
# the sites already in memory, which differs slightly (see
# `nonpar_chrx_het_count`).
FEMALE_CHRX_HET_MIN = 40_000

MODKIT_HEAD_LINES = 20000
PRIMARY_CONTIG = re.compile(r"^(chr)?([0-9]+|X|Y)$")


def check_hg38(max_end: dict[str, int], path: str) -> None:
    """Stops when a contig has a record ending beyond its hg38 length.

    The margin is END_MARGIN. Passing does not prove the file is hg38: hg19 and
    hg38 lengths differ by a few hundred kb to a few Mb, so a file is caught
    only when it has CpGs in the part of a contig that is longer in hg19 (for
    example chr1, chr9 and chr21).

    Args:
        max_end: Largest END seen per contig name (`chr1` or `1`); other contigs
            are ignored.
        path: File name for the message.

    Raises:
        runtime.RunError: A contig ends beyond its hg38 length.
    """
    for name, end in max_end.items():
        m = PRIMARY_CONTIG.match(name)
        length = HG38_LENGTHS.get("chr" + m.group(2)) if m else None
        if length is not None and end > length + END_MARGIN:
            raise runtime.RunError(
                f"{path}: {name} has a record ending at {end:,}, but {name} "
                f"is {length:,} bp in hg38/GRCh38. "
                "DMRdetector needs hg38/GRCh38 coordinates; hg19/GRCh37 and "
                "other assemblies are not supported. "
                "Align the reads to GRCh38 and call the methylation again.")


def _head_lines(path: str, limit: int):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    with (gzip.open(path, "rt") if gz else open(path)) as fh:
        n = 0
        for line in fh:
            if line.startswith("#") or line.startswith(
                    "chrom") or not line.strip():
                continue
            yield line.split()
            n += 1
            if n >= limit:
                return


def check_modkit(path: str, limit: int = MODKIT_HEAD_LINES) -> None:
    """Stops when a modkit bedMethyl file would count CpGs twice.

    The reader takes read depth (column 10) and modified reads (column 12) from
    every row and ignores the modification code (column 4) and the strand
    (column 6). A file with 5hmC (`h`) or 6mA (`a`) rows next to the 5mC (`m`)
    rows repeats each position, and a file written without `--combine-strands`
    has one row per strand (the minus-strand row sits on the G of the CpG).
    modkit writes the rows of one site next to each other, so the first `limit`
    data lines show either case; a file that changes only later is not detected.

    Raises:
        runtime.RunError: Such rows are in the head of the file.
    """
    # codes of rows that share a position with a row of another code, or are h /
    # a
    multi: set[str] = set()
    strands: set[str] = set()
    prev = None
    for f in _head_lines(path, limit):
        if len(f) < 18:
            continue
        strands.add(f[5])
        if f[3] in ("h", "a"):
            multi.add(f[3])
        if prev and prev[:2] == f[:2] and prev[3] != f[3]:
            multi.update((prev[3], f[3]))
        prev = f
    if multi or ("+" in strands and "-" in strands):
        why = (f"rows for more than one modification code "
               f"({', '.join(sorted(multi))}) in column 4"
               if multi else "one row per strand (+ and - in column 6)")
        raise runtime.RunError(
            f"{path}: this modkit bedMethyl file has {why}. DMRdetector "
            "counts every row as one CpG, so a CpG "
            "would count twice. Run `modkit pileup --cpg --combine-strands` "
            "and keep only the 5mC rows "
            "(column 4 = m; drop the h and a rows, for example with `awk '$4 "
            "== \"m\"'`).")


def contig_mismatch_warning(het: dict,
                            bed_contigs: Iterable[str]) -> str | None:
    """Returns a warning when no VCF contig name is a contig of the BED.

    Contig names compare as exact strings (as in v2.2), so `1` and `chr1` do
    not meet, and no CpG is removed for a heterozygous SNP.

    Args:
        het: {contig: heterozygous SNP positions} of the VCF.
        bed_contigs: Contig names of the methylation BED.
    """
    names = set(bed_contigs)
    vcf_names = [c for c, v in het.items() if len(v)]
    if not vcf_names or names.intersection(vcf_names):
        return None
    norm = set(bed.normalize_primary_chrom(sorted(names)))
    hint = (" (the same chromosomes with and without the 'chr' prefix)" if
            norm.intersection(bed.normalize_primary_chrom(vcf_names)) else "")
    return ("WARNING: the VCF has heterozygous SNPs on contigs "
            f"{', '.join(sorted(vcf_names)[:3])}"
            f"{' ...' if len(vcf_names) > 3 else ''}, but none of these names "
            "is a contig of the methylation BED "
            f"({', '.join(sorted(names)[:3])}"
            f"{' ...' if len(names) > 3 else ''}){hint}; no CpG is removed "
            "for a heterozygous SNP. Contig names must match exactly")


def wide_records_warning(start: np.ndarray, end: np.ndarray) -> str | None:
    """Returns a warning when a kept record is not one base long.

    A CpG record of the supported formats is [START, START + 1) (0-based,
    half-open). DMRdetector picks the records of a panel region by START, but
    `tabix -R` of v2.2 picks them by overlap, so a record with END - START != 1
    can be in one and not in the other. The check runs on the arrays that are
    already in memory.

    Args:
        start: START of the kept records (0-based).
        end: END of the kept records (exclusive).
    """
    n = int(np.count_nonzero(end - start != 1))
    if not n:
        return None
    return (f"WARNING: {n:,} of the kept methylation records have "
            "END - START != 1 (a "
            "CpG is one base); DMRdetector picks the records of a panel region "
            "by START, v2.2 (tabix) by overlap, so the calls near them can "
            "differ")


def vcf_warnings(n_samples: int, het: dict) -> list[str]:
    """Returns warnings for a VCF that cannot remove any CpG.

    Args:
        n_samples: Sample columns of the VCF (`vcf_io.sample_count`).
        het: {contig: heterozygous SNP positions} of the VCF.

    The two cases: no sample column (a sites-only VCF has no genotype), and a
    VCF with genotypes but no heterozygous SNP call. In both, the filter on
    heterozygous SNPs removes nothing. More than one sample column is reported
    by `preprocess_sample`.
    """
    if n_samples == 0:
        return [
            "WARNING: the VCF has no sample column (sites-only?), so it has "
            "no genotype; no CpG is removed for a heterozygous SNP"
        ]
    if not any(v.size for v in het.values()):
        return [
            "WARNING: the VCF has no heterozygous SNP call; no CpG is "
            "removed for a heterozygous SNP"
        ]
    return []


def nonpar_chrx_het_count(het: dict, par: dict) -> int | None:
    """Returns the number of heterozygous SNP sites on non-PAR chrX.

    The sites are those `vcf_io.het_snp_sites` found (a site counts when some
    ALT is a SNP and some sample is heterozygous; the lab's counting script
    takes the first sample and single-base ALT only), so multi-allelic sites
    count here. Contigs are matched by exact name (`chrX` or `X`).

    Args:
        het: {contig: sorted unique 0-based positions}.
        par: {contig: (starts, ends)} of the PARs, 0-based half-open
            (`read_par_bed`).

    Returns:
        The count, or None when a chrX contig of the VCF has no PAR entry (the
        PARs cannot be removed).
    """
    n = 0
    for c, pos in het.items():
        if not _X_RE.match(c):
            continue
        if c not in par:
            return None
        starts, ends = par[c]
        in_par = np.zeros(pos.size, dtype=bool)
        for a, b in zip(starts, ends):  # a union: PAR lines may overlap
            in_par |= (pos >= a) & (pos < b)
        n += int(np.count_nonzero(~in_par))
    return n


def sex_het_warning(sample_id: str, sex: str | None,
                    n_het: int | None) -> str | None:
    """Returns a warning when the sex label disagrees with the chrX SNPs.

    Args:
        sample_id: Sample ID.
        sex: M, F or empty.
        n_het: Heterozygous SNP sites on non-PAR chrX, or None when not
            counted.

    The lab convention is one cutoff, FEMALE_CHRX_HET_MIN: a female has at least
    that many, a male fewer. The check only warns.
    """
    if n_het is None:
        return None
    if sex == "M" and n_het >= FEMALE_CHRX_HET_MIN:
        return (f"{sample_id}: sex M, but the VCF has {n_het:,} heterozygous "
                "SNPs on chrX outside the PARs (a female has at least "
                f"{FEMALE_CHRX_HET_MIN:,} in our controls, a male fewer); "
                "check the sex of the sample")
    if sex == "F" and n_het < FEMALE_CHRX_HET_MIN:
        return (f"{sample_id}: sex F, but the VCF has {n_het:,} heterozygous "
                "SNPs on chrX outside the PARs (a female has at least "
                f"{FEMALE_CHRX_HET_MIN:,} in our controls); check the sex of "
                "the sample")
    return None


def preprocess_sample(meth_bed: str,
                      vcf: str,
                      min_coverage: float,
                      out_bed: str,
                      log: Callable[[str], None] = bed.warn,
                      fmt: str = "auto",
                      threads: int | None = None,
                      split_prefix: str | None = None,
                      par_bed: str | None = None) -> dict:
    """Runs v2.2 PREPROCESS_SAMPLE on one sample.

    Args:
        meth_bed: Methylation BED (5col, modkit or pbcpg).
        vcf: Sample VCF.
        min_coverage: Keep rows with coverage above this.
        out_bed: Output (bgzip plus tabix).
        log: Receives progress messages.
        fmt: "auto" (detect from the first data line), "5col", "modkit" or
            "pbcpg"; a forced format needs enough columns.
        threads: Threads of the byte-level reader and the BGZF writer (default
            DMRDETECTOR_THREADS, else 1).
        split_prefix: With `par_bed`, also write the three SPLIT_SEX_COMPONENTS
            files {split_prefix}_{autosomal,chrX,chrY}.bed.gz from the same
            arrays.
        par_bed: PAR coordinates for the split.

    Returns:
        Counts for logging and tests (with "split" counts when split, and
        "het_chrx", the heterozygous SNP sites on non-PAR chrX, when split and
        the VCF has one sample column, else None).

    Raises:
        ValueError: A forced format has too few columns.
        runtime.RunError: The input checks fail.
    """
    threads = threads or bed.default_threads()
    fmt, ncols, _ = methbed.scan_header(meth_bed, fmt)
    log(f"Detected {fmt} format")
    if fmt == "modkit":
        check_modkit(meth_bed)
    het = vcf_io.het_snp_sites(vcf, threads)
    n_vcf = vcf_io.sample_count(vcf)
    if n_vcf > 1:
        log(f"WARNING: the VCF has {n_vcf} sample columns; all samples are "
            "used, and a CpG is removed when any of them is heterozygous "
            "there")
    for w in vcf_warnings(n_vcf, het):
        log(w)
    # one sample column: with several, a site counts when any sample is
    # heterozygous, which is not the count of one person
    het_chrx = (nonpar_chrx_het_count(het, read_par_bed(par_bed))
                if split_prefix and par_bed and n_vcf == 1 else None)
    tab = None
    if ncols:
        tab = methbed.read_filtered(meth_bed, fmt, ncols, min_coverage, het,
                                    threads)
        if tab is not None and not tab.order_blocks():
            tab = None
    if tab is not None:
        _warn(log, contig_mismatch_warning(het, tab.names))
        check_hg38(
            {
                tab.names[g]: int(tab.end[a:b].max())
                for g, a, b in tab.blocks
                if b > a
            }, meth_bed)
        _warn(log, wide_records_warning(tab.start, tab.end))
    if tab is None:
        r = _preprocess_generic(meth_bed, het, fmt, min_coverage, out_bed, log)
        r["het_chrx"] = het_chrx
        if split_prefix:
            r["split"] = split_sex_components(out_bed, par_bed, split_prefix,
                                              log, threads)
        return r
    n = methbed.write_table(out_bed, tab, threads)
    if n == 0 and os.path.exists(out_bed + ".tbi"):
        os.remove(out_bed + ".tbi")
    r = {
        "format": fmt,
        "rows_in": tab.n_in,
        "het_sites": int(sum(v.size for v in het.values())),
        "removed_by_het": tab.n_drop,
        "rows_out": n,
        "het_chrx": het_chrx
    }
    if split_prefix:
        r["split"] = _split_table(tab, par_bed, split_prefix, log, threads)
    return r


def _warn(log: Callable[[str], None], message: str | None) -> None:
    if message:
        log(message)


def _preprocess_generic(meth_bed, het, fmt, min_coverage, out_bed, log) -> dict:
    """Runs the pandas implementation of preprocess_sample.

    It is the reference of the fast path and takes any file the reader accepts.
    """
    df, fmt = methbed.read_methylation_bed(meth_bed, fmt)
    n_in = len(df)
    _warn(log, contig_mismatch_warning(het, df["CHROM"].unique()))
    drop = np.zeros(n_in, dtype=bool)
    if het and n_in:
        chrom = df["CHROM"].to_numpy()
        for c, p0 in het.items():
            m = np.flatnonzero(chrom == c)
            if m.size:
                drop[m] = bed.sorted_point_overlap(df["START"].to_numpy()[m],
                                                   df["END"].to_numpy()[m], p0)
    keep = ~drop & (df["COVERAGE"].to_numpy() > min_coverage)
    df = df[keep].reset_index(drop=True)
    if len(df):
        check_hg38(df.groupby("CHROM")["END"].max().to_dict(), meth_bed)
        _warn(log, wide_records_warning(df["START"].to_numpy(),
                                        df["END"].to_numpy()))
    df["COV"] = bed.awk_number(df["COVERAGE"].to_numpy())
    df["NMOD"] = bed.awk_number(df["N_MOD"].to_numpy())
    if len(df):
        df = df.iloc[bed.bed_sort_order(df["CHROM"].to_numpy(),
                                        df["START"].to_numpy())]
        df = methbed.tie_sort(df)
    text = _fmt_df(df)
    bed.write_bed_bgzip_tabix(out_bed, text)
    return {
        "format": fmt,
        "rows_in": n_in,
        "het_sites": int(sum(v.size for v in het.values())),
        "removed_by_het": int(drop.sum()),
        "rows_out": len(df)
    }


def read_par_bed(path: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Returns the PAR interval starts and ends per exact chromosome name.

    BED, 0-based.
    """
    df = pd.read_csv(path,
                     sep="\t",
                     header=None,
                     usecols=[0, 1, 2],
                     names=["c", "s", "e"],
                     dtype={"c": str},
                     comment="#")
    return {
        c: (g["s"].to_numpy(np.int64), g["e"].to_numpy(np.int64))
        for c, g in df.groupby("c")
    }


def _in_par(chrom: np.ndarray, start: np.ndarray, end: np.ndarray,
            par) -> np.ndarray:
    """Returns the number of PAR intervals that overlap each row.

    The chromosome name must match exactly; this is the multiplicity of
    `bedtools intersect -wa`.
    """
    cnt = np.zeros(start.shape[0], dtype=np.int64)
    for c in np.unique(chrom):
        if c in par:
            m = chrom == c
            cnt[m] = bed.overlap_counts(start[m], end[m], *par[c])
    return cnt


def _fmt_df(df: pd.DataFrame) -> str:
    return bed.join_columns([
        df["CHROM"].tolist(), df["START"].astype(str).tolist(),
        df["END"].astype(str).tolist(), df["COV"].tolist(), df["NMOD"].tolist()
    ])


def split_sex_components(bed_path: str,
                         par_bed: str,
                         out_prefix: str,
                         log: Callable[[str], None] = bed.warn,
                         threads: int | None = None) -> dict:
    """Run v2.2 SPLIT_SEX_COMPONENTS on a preprocessed 5-column BED.

    Writes {out_prefix}_autosomal.bed.gz, _chrX.bed.gz, _chrY.bed.gz (+ .tbi).
    Coverage and N_MOD columns are carried as text, unchanged. A BED in the
    exact form the filter step writes (tab separated, canonical integers,
    sorted) takes the byte-level path; any other file the pandas one.
    """
    threads = threads or bed.default_threads()
    tab = methbed.read_filtered(bed_path,
                                "5col",
                                5,
                                -1.0, {},
                                threads,
                                strict=True)
    if tab is not None and tab.order_blocks():
        return _split_table(tab, par_bed, out_prefix, log, threads)
    return _split_generic(bed_path, par_bed, out_prefix, log)


def _split_table(tab, par_bed, out_prefix, log, threads) -> dict:
    counts = methbed.split_table(tab, read_par_bed(par_bed), out_prefix,
                                 _AUTO_RE, _X_RE, _Y_RE, PAR2_SHIFT, threads)
    for name in counts:
        f = f"{out_prefix}_{name}.bed.gz"
        if counts[name] == 0 and os.path.exists(f + ".tbi"):
            os.remove(f + ".tbi")
    log(f"split rows: {counts}")
    return counts


def _split_generic(bed_path: str,
                   par_bed: str,
                   out_prefix: str,
                   log=bed.warn) -> dict:
    """The pandas implementation of split_sex_components."""
    par = read_par_bed(par_bed)
    df = pd.read_csv(bed_path,
                     sep="\t",
                     header=None,
                     usecols=[0, 1, 2, 3, 4],
                     names=["CHROM", "START", "END", "COV", "NMOD"],
                     dtype={
                         "CHROM": str,
                         "START": np.int64,
                         "END": np.int64,
                         "COV": str,
                         "NMOD": str
                     },
                     keep_default_na=False)
    chrom = df["CHROM"].to_numpy(dtype=object)
    start = df["START"].to_numpy()
    end = df["END"].to_numpy()
    is_auto = np.fromiter((bool(_AUTO_RE.match(c)) for c in pd.unique(chrom)),
                          dtype=bool)
    auto_names = set(pd.unique(chrom)[is_auto])
    is_x = np.isin(chrom, [c for c in pd.unique(chrom) if _X_RE.match(c)])
    is_y = np.isin(chrom, [c for c in pd.unique(chrom) if _Y_RE.match(c)])

    core = df[np.isin(chrom, list(auto_names))]

    xi = np.flatnonzero(is_x)
    npar_x = _in_par(chrom[xi], start[xi], end[xi], par)
    par_x = df.iloc[np.repeat(xi, npar_x)].copy()
    par_x["CHROM"] = "chrX"

    yi = np.flatnonzero(is_y)
    npar_y = _in_par(chrom[yi], start[yi], end[yi], par)
    par_y = df.iloc[np.repeat(yi, npar_y)].copy()
    par_y["CHROM"] = "chrX"
    p2 = par_y["START"].to_numpy() >= 10_000_000  # awk: $2 < 10000000 -> PAR1
    par_y.loc[p2, "START"] += PAR2_SHIFT
    par_y.loc[p2, "END"] += PAR2_SHIFT

    auto = pd.concat([core, par_x, par_y], ignore_index=True)
    xnp = df.iloc[xi[npar_x == 0]]
    ynp = df.iloc[yi[npar_y == 0]]

    counts = {}
    for name, part in (("autosomal", auto), ("chrX", xnp), ("chrY", ynp)):
        part = part.reset_index(drop=True)
        if len(part):
            part = part.iloc[bed.bed_sort_order(
                part["CHROM"].to_numpy(dtype=object), part["START"].to_numpy())]
            part = methbed.tie_sort(part)
        bed.write_bed_bgzip_tabix(f"{out_prefix}_{name}.bed.gz", _fmt_df(part))
        counts[name] = len(part)
    log(f"split rows: {counts}")
    return counts


# steps of build and run (pool of samples, resume stamps)

# scope of each non-autosomal part, by the sex of the sample
_SEX_PARTS = {
    "M": (("chrX", "male_X"), ("chrY", "Y")),
    "F": (("chrX", "female_X"),)
}


def split_warnings(sample_id: str, sex: str | None, counts: dict | None,
                   rows_out: int) -> list[str]:
    """Returns warnings about the CpGs per part of a split sample.

    Args:
        sample_id: Sample ID.
        sex: M, F or empty.
        counts: Rows per part (autosomal, chrX, chrY) after filtering, or None
            when the sample was not split.
        rows_out: CpGs of the filtered sample, on any contig.

    The warnings cover an empty part that the sample is tested in, an empty
    autosomal part although CpGs remain on other contigs, and a female sample
    with as many chrY non-PAR CpGs as a male (CHRY_MALE_MIN_CPG).
    """
    if not counts:
        return []
    out = []
    if counts.get("autosomal", 0) == 0 and rows_out > 0:
        out.append(f"{sample_id}: no CpG on chr1-chr22 or the PARs, but "
                   f"{rows_out} on other contigs; are the contigs named chrN "
                   "or N (not RefSeq names such as NC_000001.11)?")
    for part, scope in _SEX_PARTS.get(sex or "", ()):
        if counts.get(part, 0) == 0 and rows_out > 0:
            out.append(f"{sample_id}: no CpG after filtering in the {scope} "
                       f"scope ({part} outside the PARs); it is empty for "
                       "this sample")
    if sex == "F" and counts.get("chrY", 0) >= CHRY_MALE_MIN_CPG:
        out.append(
            f"{sample_id}: sex F, but {counts['chrY']:,} CpGs remain on chrY "
            f"outside the PARs (a male has at least {CHRY_MALE_MIN_CPG:,} "
            "records before filtering in our controls); check the sex of the "
            "sample")
    return out


def filtered_path(ctx: runtime.Ctx, sid: str) -> str:
    """Returns the path of the filtered BED of a sample."""
    return str(ctx.workdir / "preprocessed" / f"{sid}_filtered.bed.gz")


def scope_bed(ctx: runtime.Ctx, sid: str, scope: str, split: bool) -> str:
    """Returns the BED a scope reads: the filtered file, or its split part."""
    if not split:
        return filtered_path(ctx, sid)
    return str(ctx.workdir / "split" /
               f"{sid}_{runtime.SCOPE_PART[scope]}.bed.gz")


def _t_preprocess(t: dict):
    out = t["out"]
    prefix = os.path.join(t["split_dir"], t["id"])
    parts = [f"{prefix}_{p}.bed.gz" for p in runtime.PARTS]
    state = {"split_done": False}

    def work():
        """Preprocesses the sample; returns its WARNING and NOTICE lines."""
        tp = os.path.join(t["split_dir"], f".tmp{os.getpid()}.{t['id']}")
        tp_parts = [f"{tp}_{p}.bed.gz" for p in runtime.PARTS]
        notes: list[str] = []

        def say(line: str) -> None:
            """Logs a line and keeps a warning, which a resume replays."""
            runtime.log(line)
            if "WARNING" in line or "NOTICE" in line:
                notes.append(line)

        with runtime.temp_output(out, *tp_parts) as tmp:
            try:
                # one pass writes the filtered file and, when a split is wanted,
                # the three split files
                r = preprocess_sample(t["meth_bed"],
                                      t["vcf"],
                                      t["min_cov"],
                                      tmp,
                                      log=lambda m: say(f"  {t['id']}: {m}"),
                                      fmt=t["fmt"],
                                      threads=t["threads"],
                                      split_prefix=tp if t["split"] else None,
                                      par_bed=t["par_bed"])
            except ValueError as e:
                if t["fmt"] == "auto":
                    raise
                if "needs at least" in str(e):
                    raise runtime.RunError(str(e)) from None
                raise runtime.RunError(
                    f"{t['meth_bed']}: cannot read as input_format {t['fmt']} "
                    f"({e}); "
                    "check input_format or use auto") from None
            runtime.commit(tmp, out)
            if t["split"]:
                for part, final in zip(tp_parts, parts):
                    runtime.commit(part, final)
                state["split_done"] = True
        if r["rows_in"] == 0:
            say(f"WARNING: {t['id']}: no CpG records in {t['meth_bed']}")
        elif r["rows_out"] == 0:
            say(f"WARNING: {t['id']}: no CpG is left after filtering "
                f"({r['rows_in']} in; {r['removed_by_het']} on het SNPs; "
                f"coverage must exceed {t['min_cov']}); every scope of this "
                "sample is empty")
        for w in split_warnings(t["id"], t.get("sex"), r.get("split"),
                                r["rows_out"]):
            say(f"WARNING: {w}")
        w = sex_het_warning(t["id"], t.get("sex"), r.get("het_chrx"))
        if w:
            say(f"WARNING: {w}")
        runtime.log(
            f"  {t['id']}: {r['format']}, {r['rows_in']} -> {r['rows_out']} "
            "CpGs "
            f"({r['removed_by_het']} on het SNPs, {r['het_sites']} het sites)")
        return notes

    # sex, split and the PARs decide which warnings the step stores, so a
    # resume with another sex label or scope setting runs the checks again
    key = dict(step="preprocess",
               min_cov=t["min_cov"],
               fmt=t["fmt"],
               meth=runtime.stat_of(t["meth_bed"]),
               vcf=runtime.stat_of(t["vcf"]),
               sex=t.get("sex"),
               split=t["split"],
               par=runtime.stat_of(t["par_bed"]) if t["par_bed"] else None)
    n1 = runtime.guarded(t["stamps"], f"preprocess_{t['id']}", key, [out], work)
    if not t["split"]:
        return [n1]

    def split():
        if state["split_done"]:  # written together with the filtered file
            return
        tp = os.path.join(t["split_dir"], f".tmp{os.getpid()}.{t['id']}")
        try:
            split_sex_components(out,
                                 t["par_bed"],
                                 tp,
                                 log=lambda m: runtime.log(f"  {t['id']}: {m}"),
                                 threads=t["threads"])
            for p, final in zip(runtime.PARTS, parts):
                runtime.commit(f"{tp}_{p}.bed.gz", final)
        finally:
            runtime.cleanup(*[f"{tp}_{p}.bed.gz" for p in runtime.PARTS])

    key2 = dict(step="split",
                bed=runtime.stat_of(out),
                par=runtime.stat_of(t["par_bed"]) if t["par_bed"] else None)
    return [
        n1,
        runtime.guarded(t["stamps"], f"split_{t['id']}", key2, parts, split)
    ]


def preprocess_samples(ctx: runtime.Ctx, samples: list[runtime.Sample],
                       min_cov: float, split: bool,
                       par_bed: str | None) -> None:
    """Preprocesses samples in the process pool, with resume stamps."""
    tasks = [
        dict(id=s.sample_id,
             meth_bed=s.meth_bed,
             vcf=s.vcf,
             min_cov=min_cov,
             fmt=ctx.cfg["input_format"],
             out=filtered_path(ctx, s.sample_id),
             stamps=ctx.stamps,
             split=split,
             split_dir=str(ctx.workdir / "split"),
             sex=s.sex,
             par_bed=par_bed) for s in samples
    ]
    runtime.log(f"preprocess: {len(tasks)} samples (coverage > {min_cov}, "
                f"input_format {ctx.cfg['input_format']})")
    t0 = time.time()
    runtime.pool_map(ctx, _t_preprocess, tasks, runtime.TASK_GB["preprocess"],
                     "preprocess")
    runtime.log(f"preprocess finished ({time.time() - t0:.0f} s)")


def _main_filter(args) -> int:
    out = args.out or f"{args.sample_id}_filtered.bed.gz"
    r = preprocess_sample(args.meth_bed,
                          args.vcf,
                          args.min_coverage,
                          out,
                          fmt=args.input_format,
                          threads=args.threads)
    print(f"{args.sample_id}: {r}", file=sys.stderr)
    return 0


def _main_split(args) -> int:
    par = args.par_bed or data.asset("par_hg38.bed")
    split_sex_components(args.bed,
                         par,
                         os.path.join(args.outdir, args.sample_id),
                         threads=args.threads)
    return 0


def main(args: argparse.Namespace) -> int:
    """Runs the `preprocess filter` or `preprocess split` subcommand."""
    return args.subfunc(args)


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    """Adds `preprocess` with `filter` and `split` to `subparsers`."""
    p = subparsers.add_parser("preprocess",
                              help="filter one sample and split it by scope")
    p.set_defaults(func=main)
    sub = p.add_subparsers(dest="preprocess_cmd", required=True)

    f = sub.add_parser(
        "filter",
        help="drop het-SNP CpGs and low coverage; write {id}_filtered.bed.gz")
    f.add_argument("--sample-id", required=True)
    f.add_argument("--meth-bed",
                   required=True,
                   help="5-col, modkit bedMethyl, or pb-CpG-tools BED (.gz ok)")
    f.add_argument("--vcf", required=True, help="sample VCF (.vcf or .vcf.gz)")
    f.add_argument("--min-coverage",
                   type=float,
                   default=min(settings.ATLAS_MIN_READ_COVERAGE,
                               settings.DMR_MIN_COVERAGE),
                   help="keep rows with coverage > this (v2.2: "
                   "min(atlas_min_read_coverage, "
                   "dmr_min_coverage)) "
                   "[%(default)s]")
    f.add_argument("--input-format",
                   choices=settings.INPUT_FORMATS,
                   default=settings.INPUT_FORMAT,
                   help="methylation BED format [%(default)s]")
    f.add_argument("--out", help="default {sample-id}_filtered.bed.gz")
    _add_threads(f)
    f.set_defaults(subfunc=_main_filter)

    s = sub.add_parser(
        "split",
        help="split into autosomal(+PAR) / non-PAR chrX / non-PAR chrY")
    s.add_argument("--sample-id", required=True)
    s.add_argument("--bed", required=True, help="{id}_filtered.bed.gz")
    s.add_argument("--par-bed", help="default assets/par_hg38.bed")
    s.add_argument("--outdir", default=".")
    _add_threads(s)
    s.set_defaults(subfunc=_main_split)


def _add_threads(p):
    p.add_argument("--threads",
                   type=int,
                   default=None,
                   help="threads of the reader and the BGZF writer (default "
                   "$DMRDETECTOR_THREADS, else 1)")
