"""Unit tests for io.py, preprocess.py, atlas.py (synthetic data).

PYTHONPATH=src python -m pytest tests/unit/test_preprocess_atlas.py Comparisons
with the real v1 R/bash code are in tests/cmp_preprocess.sh and cmp_atlas.sh.
"""
import gzip
import os

import numpy as np
import pandas as pd
import pytest

from dmrdetector import data
from dmrdetector.fileio import bed as dio  # noqa: E402
from dmrdetector.fileio import methbed  # noqa: E402
from dmrdetector.fileio import vcf as vcf_io  # noqa: E402
from dmrdetector.steps import atlas
from dmrdetector.steps import preprocess

from helpers import generate  # noqa: E402
from helpers import read_bed5  # noqa: E402


PAR = data.asset("par_hg38.bed")


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    d = tmp_path_factory.mktemp("synth")
    generate(str(d), n=4, seed=3)
    return str(d)


def zlines(p):
    with gzip.open(p, "rt") as fh:
        return fh.read().splitlines()


# ---- formatting ----------------------------------------------------------


def test_fmt_double_matches_fwrite_notation():
    # values and expected strings taken from data.table::fwrite (R 4.4)
    x = np.array([
        0, 1, 0.1, 1e-4, 2e-4, 1e-5, 1.5e-5, 123456, 1e5, 1e6, 1.2e6, 1e15,
        0.0001234, 100000.5, 1 / 3, 2 / 3, 0.001, 0.0012
    ])
    exp = [
        "0", "1", "0.1", "1e-04", "2e-04", "1e-05", "1.5e-05", "123456",
        "1e+05", "1e+06", "1200000", "1e+15", "0.0001234", "100000.5",
        "0.333333333333333", "0.666666666666667", "0.001", "0.0012"
    ]
    assert dio.fmt_doubles(x) == exp


def test_awk_number():
    assert dio.awk_number(np.array([5.0, 12.5, 1234567.0, 0.1234567
                                   ])) == ["5", "12.5", "1234567", "0.123457"]


def test_normalize_primary_chrom():
    out = dio.normalize_primary_chrom(
        ["1", "chr2", "CHRX", "y", "chrUn_x", "chrM"])
    assert list(out) == ["chr1", "chr2", "chrX", "chrY", "chrUn_x", "chrM"]


def test_c_order():
    assert dio.c_order(["chr2", "chr10", "chrX",
                        "chr1"]) == ["chr1", "chr10", "chr2", "chrX"]


# ---- readers -------------------------------------------------------------


def test_three_formats_agree(synth):
    ref, f0 = methbed.read_methylation_bed(
        os.path.join(synth, "S01.5col.bed.gz"))
    assert f0 == "5col"
    for fmt in ("modkit", "pbcpg"):
        df, f = methbed.read_methylation_bed(
            os.path.join(synth, f"S01.{fmt}.bed.gz"))
        assert f == fmt
        pd.testing.assert_frame_equal(df, ref, check_dtype=False)


# ---- preprocess ----------------------------------------------------------


def test_filter_removes_het_snps_only_and_applies_coverage(synth, tmp_path):
    out = str(tmp_path / "S01_filtered.bed.gz")
    r = preprocess.preprocess_sample(os.path.join(synth, "S01.5col.bed.gz"),
                                     os.path.join(synth, "S01.vcf"),
                                     15,
                                     out,
                                     log=lambda m: None)
    assert os.path.exists(out + ".tbi")
    raw, _ = methbed.read_methylation_bed(os.path.join(synth,
                                                       "S01.5col.bed.gz"))
    het = vcf_io.het_snp_sites(os.path.join(synth, "S01.vcf"))
    # contig "21" is not in the BED: removes nothing;
    # hom/missing/FILTERed-but-het follow bcftools
    assert "21" in het and "chr21" in het
    got = pd.read_csv(out, sep="\t", header=None)
    assert (got[3] > 15).all()
    assert r["rows_out"] == len(got) < len(raw)
    # no remaining row overlaps a het site
    for c, p0 in het.items():
        g = got[got[0] == c]
        assert not dio.sorted_point_overlap(g[1].to_numpy(), g[2].to_numpy(),
                                            p0).any()


def test_split_par_folding(synth, tmp_path):
    f = str(tmp_path / "S02_filtered.bed.gz")
    preprocess.preprocess_sample(os.path.join(synth, "S02.5col.bed.gz"),
                                 os.path.join(synth, "S02.vcf"),
                                 5,
                                 f,
                                 log=lambda m: None)
    preprocess.split_sex_components(f,
                                    PAR,
                                    str(tmp_path / "S02"),
                                    log=lambda m: None)
    src = pd.read_csv(f, sep="\t", header=None, dtype={0: str})
    auto = pd.read_csv(str(tmp_path / "S02_autosomal.bed.gz"),
                       sep="\t",
                       header=None,
                       dtype={0: str})
    x = pd.read_csv(str(tmp_path / "S02_chrX.bed.gz"), sep="\t", header=None)
    y = pd.read_csv(str(tmp_path / "S02_chrY.bed.gz"), sep="\t", header=None)
    sx = src[src[0] == "chrX"]
    sy = src[src[0] == "chrY"]

    def in_par(d):
        return (((d[1] < 2781479) & (d[2] > 10000)) | ((d[1] < 156030895) &
                                                       (d[2] > 155701382)))

    def in_par_y(d):
        return (((d[1] < 2781479) & (d[2] > 10000)) | ((d[1] < 57217415) &
                                                       (d[2] > 56887902)))

    assert len(x) == (~in_par(sx)).sum()
    assert len(y) == (~in_par_y(sy)).sum()
    n_auto_core = src[~src[0].isin(["chrX", "chrY"])].shape[0]
    assert len(auto) == n_auto_core + in_par(sx).sum() + in_par_y(sy).sum()
    # chrY PAR2 rows are shifted onto chrX PAR2
    y2 = sy[(sy[1] >= 56887902 - 1) & in_par_y(sy) & (sy[1] > 10_000_000)]
    shifted = y2[1] + preprocess.PAR2_SHIFT
    assert set(shifted).issubset(set(auto[auto[0] == "chrX"][1]))
    assert (auto[0] != "chrY").all()
    # sorted by (chrom text, start)
    key = list(zip(auto[0], auto[1]))
    assert key == sorted(key)


# ---- atlas build ---------------------------------------------------------


def _brute_atlas(beds, min_samples, w):
    frames = [read_bed5(b) for b in sorted(beds)]
    d = pd.concat(frames, ignore_index=True)
    d["ratio"] = d["N_METH"] / d["COVERAGE"]
    rows = []
    for (c, s, e), g in d.groupby(["CHROM", "START", "END"], sort=False):
        n = len(g)
        if n < min_samples:
            continue
        cov = g["COVERAGE"].sum()
        raw = g["N_METH"].sum() / cov
        s2 = np.var(g["ratio"].to_numpy(), ddof=1)
        nbar = cov / n
        mu = raw * (1 - raw)
        num, den = s2 * nbar - mu, mu * (nbar - 1)
        phi = num / den if (den > 1e-9 and num > 0) else 0.0
        rows.append((c, s, e, n, cov, min(phi, 1.0), raw))
    t = pd.DataFrame(
        rows, columns=["CHROM", "START", "END", "N", "COV", "phi_raw", "RAW"])
    trend = np.median(t.phi_raw[t.phi_raw > 0]) if (t.phi_raw
                                                    > 0).any() else 0.01
    t["PHI"] = (1 - w) * t.phi_raw + w * trend
    return t.sort_values(["CHROM", "START"],
                         kind="stable").reset_index(drop=True), trend


def test_build_matches_brute_force(synth, tmp_path):
    beds = []
    for k in (1, 2, 3, 4):
        f = str(tmp_path / f"S0{k}_filtered.bed.gz")
        preprocess.preprocess_sample(os.path.join(synth, f"S0{k}.5col.bed.gz"),
                                     os.path.join(synth, f"S0{k}.vcf"),
                                     5,
                                     f,
                                     log=lambda m: None)
        preprocess.split_sex_components(f,
                                        PAR,
                                        str(tmp_path / f"S0{k}"),
                                        log=lambda m: None)
        beds.append(str(tmp_path / f"S0{k}_autosomal.bed.gz"))
    out = str(tmp_path / "raw.bed.gz")
    trend = atlas.build_atlas(beds,
                              out,
                              3,
                              0.1,
                              threads=2,
                              scratch=str(tmp_path / "tbi"),
                              window=500_000,
                              log=lambda m: None)
    exp, etrend = _brute_atlas(beds, 3, 0.1)
    assert trend == pytest.approx(etrend, rel=1e-12)
    got = pd.read_csv(out,
                      sep="\t",
                      header=None,
                      names=["CHROM", "START", "END", "N", "COV", "PHI", "RAW"],
                      dtype={"CHROM": str})
    assert len(got) == len(exp)
    for c in ("CHROM", "START", "END", "N", "COV"):
        assert (got[c].to_numpy() == exp[c].to_numpy()).all()
    np.testing.assert_allclose(got.RAW, exp.RAW, rtol=1e-12)
    np.testing.assert_allclose(got.PHI, exp.PHI, rtol=1e-10, atol=1e-14)


def test_build_skips_tiny_files_and_errors_without_input(tmp_path):
    empty = str(tmp_path / "e.bed.gz")
    dio.write_bed_bgzip_tabix(empty, "")
    assert os.path.getsize(empty) <= 50
    with pytest.raises(SystemExit):
        atlas.build_atlas([empty],
                          str(tmp_path / "o.bed.gz"),
                          scratch=str(tmp_path),
                          log=lambda m: None)


# ---- smoothing -----------------------------------------------------------


def test_smooth_matches_direct_convolution():
    rng = np.random.default_rng(0)
    start = np.unique(rng.integers(1000, 30000, 600))
    end = start + 1
    cov = rng.integers(5, 60, start.size).astype(float)
    raw = rng.random(start.size)
    phi = rng.random(start.size) * 0.2
    sigma = 40
    k = atlas.gaussian_kernel(sigma)
    ref, var = atlas.smooth_chromosome(start, end, cov, raw, phi, k)
    span = int(end.max() - start.min() + 1)
    vm = np.zeros(span)
    vc = np.zeros(span)
    idx = start - start.min()
    vm[idx], vc[idx] = raw * cov, cov
    sm = np.convolve(vm, k, mode="same")
    sc = np.convolve(vc, k, mode="same")
    exp = np.where(sc[idx] > 1e-9, sm[idx] / sc[idx], 0)
    np.testing.assert_allclose(ref, exp, rtol=1e-10)
    np.testing.assert_allclose(var,
                               exp * (1 - exp) * phi,
                               rtol=1e-10,
                               atol=1e-14)


def test_smooth_chromosome_filter_and_order(tmp_path):
    rows = []
    for c in ("chrUn_x", "22", "chrY", "chr2", "chrX", "chr10"):
        for s in (100, 200, 400):
            rows.append((c, s, s + 1, 5, 100, 0.01, 0.5))
    raw = str(tmp_path / "raw.bed.gz")
    dio.write_bgzip_text(raw, [
        dio.join_columns(
            [[r[i] if i == 0 else str(r[i]) for r in rows] for i in range(7)])
    ])
    out = str(tmp_path / "ref.bed.gz")
    atlas.smooth_atlas(raw,
                       10,
                       out,
                       include_sex=True,
                       threads=2,
                       log=lambda m: None)
    got = pd.read_csv(out, sep="\t", header=None, dtype={0: str})
    assert got.shape[1] == 9
    assert list(dict.fromkeys(
        got[0])) == ["chr10", "chr2", "chr22", "chrX", "chrY"]
    atlas.smooth_atlas(raw,
                       10,
                       out,
                       include_sex=False,
                       threads=1,
                       log=lambda m: None)
    got = pd.read_csv(out, sep="\t", header=None, dtype={0: str})
    assert list(dict.fromkeys(got[0])) == ["chr10", "chr2", "chr22"]
    assert np.allclose(got[7], 0.5)
