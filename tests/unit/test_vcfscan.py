"""Tests of the VCF byte scan (numba) against the pysam reference loop.

The VCFs are synthetic.
"""
import os
import random
import sys

import numpy as np
import pytest


sys.path.insert(0, os.path.dirname(__file__))
from dmrdetector.fileio import vcf as vcfscan  # noqa: E402


N = int(os.environ.get("VCFSCAN_N", "60000"))
HEAD = (
    "##fileformat=VCFv4.2\n" +
    "".join(f"##contig=<ID={c}>\n"
            for c in ("chr1", "chr2", "chrX", "1", "chrUn_x")) +
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">\n'
)
PAIRS = [("A", "G"), ("AC", "A"), ("A", "AT"), ("CGTGTGT", "C,TGTGTGT"),
         ("AT", "GT"), ("AT", "GC"), ("a", "g"), ("A", "<DEL>"), ("A", "*"),
         ("A", "G,*"), ("ACGT", "ACCT"), ("A", "A"), ("A", "X"),
         ("AAA", "AAAA,AGA"), ("C", "T,G"), ("AGC", "A[chr1:5[")]
GTS = [
    "0/1", "1|0", "1/1", "0/0", "./.", "./1", "0/.", "1", ".", "0|2", "1/2",
    "0/1/1", "0/0/1", "2|2"
]


def make(path, n, seed=1, samples=1, fmt="GT:DP", sample_fmt=None):
    rnd = random.Random(seed)
    lines = [
        HEAD + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" +
        "".join(f"\tS{i}" for i in range(samples)) + "\n"
    ]
    for i in range(n):
        ref, alt = rnd.choice(PAIRS)
        nalt = alt.count(",") + 1
        row = [
            rnd.choice(("chr1", "chr2", "chrX", "1", "chrUn_x")),
            str(rnd.randint(1, 5000)), ".", ref, alt, "30", "PASS", "DP=3", fmt
        ]
        for _ in range(samples):
            gt = rnd.choice(GTS)
            if any(
                    int(a) > nalt
                    for a in gt.replace("|", "/").split("/")
                    if a.isdigit()):
                gt = "0/1" if nalt >= 1 else gt
            row.append(f"{gt}:7" if fmt == "GT:DP" else (
                f"7:{gt}" if fmt == "DP:GT" else gt))
        lines.append("\t".join(row) + "\n")
    with open(path, "w") as fh:
        fh.write("".join(lines))


def same(a, b):
    assert sorted(a) == sorted(b)
    for c in a:
        assert np.array_equal(a[c], b[c]), c


@pytest.mark.parametrize("samples,fmt", [(1, "GT:DP"), (3, "GT:DP"),
                                         (2, "DP:GT"), (1, "GT")])
def test_matches_pysam(tmp_path, samples, fmt):
    p = str(tmp_path / "t.vcf")
    make(p, N, samples=samples, fmt=fmt)
    got = vcfscan.scan_het_sites(p, threads=3)
    ref = vcfscan._het_snp_sites_pysam(p)
    # GT not first: htslib does not expose it, the scan defers to pysam
    if fmt == "DP:GT":
        assert got is None and ref == {}
        return
    assert got is not None
    same(got, ref)
    assert sum(v.size for v in got.values()) > 0


def test_gzip_and_bgzf_and_stitching(tmp_path):
    import gzip

    import pysam
    p = str(tmp_path / "t.vcf")
    make(p,
         N * 3 if os.environ.get("DMRDETECTOR_NO_JIT") != "1" else N,
         seed=5,
         samples=2)
    ref = vcfscan._het_snp_sites_pysam(p)
    g = str(tmp_path / "t.vcf.gz")
    with open(p, "rb") as fi, gzip.open(g, "wb") as fo:
        fo.write(fi.read())
    same(vcfscan.scan_het_sites(g, threads=2), ref)
    b = str(tmp_path / "t.bg.vcf.gz")
    with open(p, "rb") as fi, pysam.BGZFile(b, "wb") as fo:
        fo.write(fi.read())
    same(vcfscan.scan_het_sites(b, threads=4), ref)
    same(vcfscan.het_snp_sites(b), ref)


@pytest.mark.parametrize(
    "line",
    [
        "chr1\t10\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\r",  # CR
        "chr1\t10\t.\tA\tG\t30\tPASS\t.\tGT\t0|a",  # non-numeric allele
        "chr1\t10\t.\tA\tG\t30\tPASS\t.\tGT\t0/5",  # allele beyond ALT
        "chr1\t10\t.\tA\tG,,T\t30\tPASS\t.\tGT\t0/1",  # empty allele
        "chr1\tx\t.\tA\tG\t30\tPASS\t.\tGT\t0/1",  # bad POS
    ])
def test_anomaly_falls_back(tmp_path, line):
    p = str(tmp_path / "t.vcf")
    with open(p, "w", newline="") as fh:
        fh.write(HEAD +
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" +
                 line + "\n")
    assert vcfscan.scan_het_sites(p) is None


def test_no_alt_no_format_no_samples(tmp_path):
    p = str(tmp_path / "t.vcf")
    with open(p, "w") as fh:
        fh.write(HEAD + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 "chr1\t10\t.\tA\tG\t30\tPASS\t.\n"
                 "chr1\t11\t.\tA\t.\t30\tPASS\t.\n"
                 "chr1\t12\t.\tA\tG\t30\tPASS\t.\tGT\n")
    same(vcfscan.scan_het_sites(p), vcfscan._het_snp_sites_pysam(p))
    assert vcfscan.scan_het_sites(p) == {}


def test_empty_and_header_only(tmp_path):
    p = str(tmp_path / "t.vcf")
    with open(p, "w") as fh:
        fh.write(HEAD +
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    assert vcfscan.scan_het_sites(p) == {}
