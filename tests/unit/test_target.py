import gzip
import os

import numpy as np
import pandas as pd
import pysam
import pytest

from dmrdetector.fileio import bed as dio
from dmrdetector.steps import target


def _bgzip_bed(path, rows):
    rows = sorted(rows, key=lambda r: (r[0], r[1]))
    dio.write_bed_bgzip_tabix(
        str(path), "".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


def _build(tmp_path, n=40, sample_shift=0.0, seed=1):
    rng = np.random.default_rng(seed)
    pos = 1000 + 10 * np.arange(n)
    ref = np.full(n, 0.5)
    atlas = [("chr1", int(p), int(p) + 1, 8, 100, 0.01, 0.5, 0.5,
              0.0004 + 1e-5 * i) for i, p in enumerate(pos)]
    # second region on chr2: few CpGs (< 30) ; third on chrX (skipped)
    pos2 = 5000 + 10 * np.arange(10)
    atlas += [("chr2", int(p), int(p) + 1, 8, 100, 0.01, 0.5, 0.5, 0.0004)
              for p in pos2]
    cov = rng.integers(10, 30, n)
    meth = np.clip(ref + sample_shift + rng.normal(0, 0.03, n), 0, 1)
    smp = [("chr1", int(p), int(p) + 1, int(c), int(round(m * c)))
           for p, c, m in zip(pos, cov, meth)]
    smp += [("chr1", 990, 991, 20, 10),
            ("chr1", 1000 + 10 * n + 5, 1000 + 10 * n + 6, 20, 10)
           ]  # outside region
    smp += [("chr2", int(p), int(p) + 1, 20, 18) for p in pos2]
    a = _bgzip_bed(tmp_path / "atlas.bed.gz", atlas)
    s = _bgzip_bed(tmp_path / "S1.combined.bed.gz", smp)
    tb = tmp_path / "target.bed"
    tb.write_text("CHROM\tSTART\tEND\tREGION_NAME\tGENE_SYMBOL\tDISEASE_NAME\n"
                  f"chr1\t1000\t{1000 + 10 * n}\tR1\tG1\tDisease A (hypo)\n"
                  "chr2\t5000\t5100\tR2\tG2;G3\tDisease B (hyper);Disease C\n"
                  "chrX\t100\t200\tRX\tGX\tDisease X\n")
    return a, s, str(tb), pos


def test_hypo_region_and_not_evaluable(tmp_path):
    a, s, tb, _ = _build(tmp_path, sample_shift=-0.45)
    msgs = []
    df = target.target_report([{
        "sample_id": "S1",
        "type": "CASE",
        "merged": s
    }],
                              a,
                              tb,
                              str(tmp_path / "o.tsv"),
                              log=msgs.append)
    assert any("RX" in m for m in msgs)  # chrX skipped with a message
    assert list(df.columns) == target.OUT_COLS and "P_VALUE" not in df.columns
    assert set(df.REGION_NAME) == {"R1", "R2"} and set(
        df.HAPLOTYPE) == {"merged"}
    r1 = df[df.REGION_NAME == "R1"].iloc[0]
    assert r1.NUM_CPG == 40  # outside-region CpGs excluded
    assert r1.DIFF < -0.4 and r1.DIRECTION == "hypo" and r1.P_FDR < 1e-6
    r2 = df[df.REGION_NAME == "R2"].iloc[0]
    assert r2.NUM_CPG == 10 and pd.isna(r2.DIRECTION) and np.isnan(
        r2.P_FDR)  # constant diffs: t-test undefined -> not evaluable
    txt = (tmp_path / "o.tsv").read_text().splitlines()
    assert txt[0].split("\t") == target.OUT_COLS
    assert txt[2].split("\t")[-1] == "NA" and txt[2].split("\t")[-2] == "NA"


def test_normal_and_icr_threshold(tmp_path):
    a, s, tb, pos = _build(tmp_path, sample_shift=-0.30)
    icr = tmp_path / "icr.bed"
    icr.write_text("CHROM\tSTART\tEND\tGENE\n" + "chr1\t900\t1100\tG1\n")
    # diff -0.30: below the non-ICR threshold 0.40 -> normal; above ICR 0.25 ->
    # hypo
    r = target.target_report([{
        "sample_id": "S1",
        "type": "CASE",
        "merged": s
    }], a, tb)
    assert r[r.REGION_NAME == "R1"].iloc[0].DIRECTION == "normal"
    r = target.target_report([{
        "sample_id": "S1",
        "type": "CASE",
        "merged": s
    }],
                             a,
                             tb,
                             icr_region=str(icr))
    assert r[r.REGION_NAME == "R1"].iloc[0].DIRECTION == "hypo"


def test_values_match_manual(tmp_path):
    a, s, tb, pos = _build(tmp_path, sample_shift=0.1)
    df = target.target_report([{
        "sample_id": "S1",
        "type": "CASE",
        "merged": s
    }],
                              a,
                              tb,
                              min_coverage=12)
    sm = pd.read_csv(s,
                     sep="\t",
                     header=None,
                     names=["c", "s", "e", "cov", "n"])
    sm = sm[(sm.s >= 1000) & (sm.s < 1400) & (sm["cov"] > 12)]
    r = df[df.REGION_NAME == "R1"].iloc[0]
    assert r.NUM_CPG == len(sm)
    assert r["SAMPLE_METH (%)"] == pytest.approx(
        round((sm.n / sm["cov"]).mean() * 100, 1))
    assert r["COVERAGE"] == pytest.approx(round(sm["cov"].mean(), 1))
    assert r.DIFF == pytest.approx((sm.n / sm["cov"]).mean() - 0.5)
    idx = ((sm.s - 1000) // 10).to_numpy()
    assert r["REF_SD (%)"] == pytest.approx(
        round(np.sqrt(np.mean(0.0004 + 1e-5 * idx)) * 100, 1))


def test_haplotypes_het_snp_and_order(tmp_path):
    a, s, tb, pos = _build(tmp_path, sample_shift=-0.45)
    # hap1: same CpGs, coverage 5 (> 4 passes, <= 8 would fail merged); hap2
    # missing
    hap = [("chr1", int(p), int(p) + 1, 5, 1) for p in pos]
    h1 = _bgzip_bed(tmp_path / "S1.hap1.bed.gz", hap)
    vcf = tmp_path / "s.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                   f"chr1\t{int(pos[3]) + 1}\t.\tA\tG\t50\tPASS\t.\tGT\t0|1\n")
    smp = [{
        "sample_id": "S2",
        "type": "CTRL",
        "merged": s,
        "hap1": h1,
        "vcf": str(vcf)
    }, {
        "sample_id": "S1",
        "type": "CASE",
        "merged": s,
        "hap1": h1
    }]
    df = target.target_report(smp, a, tb)
    assert list(dict.fromkeys(df.SAMPLE_ID)) == ["S1", "S2"]  # CASE first
    s1 = df[(df.SAMPLE_ID == "S1") & (df.REGION_NAME == "R1")]
    s2 = df[(df.SAMPLE_ID == "S2") & (df.REGION_NAME == "R1")]
    assert list(s1.HAPLOTYPE) == ["merged", "hap1"]
    h = s1[s1.HAPLOTYPE == "hap1"].iloc[0]
    assert h.NUM_CPG == 40 and np.isnan(h.DIFF) and np.isnan(
        h.P_FDR) and pd.isna(h.DIRECTION)
    assert s2[s2.HAPLOTYPE ==
              "merged"].iloc[0].NUM_CPG == 39  # het-SNP CpG removed
    # hap1 R2 has no CpG: NUM_CPG 0, NA values
    r2 = df[(df.SAMPLE_ID == "S1") & (df.REGION_NAME == "R2") &
            (df.HAPLOTYPE == "hap1")].iloc[0]
    assert r2.NUM_CPG == 0 and np.isnan(r2["SAMPLE_METH (%)"])


def test_bh_over_evaluable_regions_only(tmp_path):
    a, s, tb, _ = _build(tmp_path, sample_shift=-0.2)
    df = target.target_report([{
        "sample_id": "S1",
        "type": "CASE",
        "merged": s
    }], a, tb)
    r1 = df[df.REGION_NAME == "R1"].iloc[0]
    # only one evaluable region: BH adjusted p equals the raw p
    from dmrdetector.util.rstats import t_test_one_sample_p
    sm = pd.read_csv(s,
                     sep="\t",
                     header=None,
                     names=["c", "s", "e", "cov", "n"])
    sm = sm[(sm.s >= 1000) & (sm.s < 1400) & (sm["cov"] > 8)]
    assert r1.P_FDR == pytest.approx(
        t_test_one_sample_p((sm.n / sm["cov"]).to_numpy() - 0.5))


def test_samplesheet_haplotypes_come_from_the_columns_only(tmp_path):
    for t in ("combined", "hap1"):
        (tmp_path / f"X.{t}.bed.gz").write_text("")
    f = tmp_path / "s.csv"
    f.write_text(
        f"sample_id,type,meth_bed,vcf,meth_bed_hap1,meth_bed_hap2\nA,CASE,{tmp_path}/X.combined.bed.gz,,,\n"
        f"B,CASE,{tmp_path}/X.combined.bed.gz,,{tmp_path}/X.hap1.bed.gz,\n")
    a, b = target.read_samplesheet(str(f))
    assert "hap1" not in a and b["hap1"].endswith(
        "X.hap1.bed.gz") and "hap2" not in b  # no file-name convention


def test_cli_registered():
    import argparse
    p = argparse.ArgumentParser()
    sub = p.add_subparsers()
    target.add_parser(sub)
    a = p.parse_args([
        "target", "--samplesheet", "s.csv", "--atlas", "a", "--target-bed", "t",
        "--dmr-fdr-threshold", "0.05"
    ])
    assert a.fdr_threshold == 0.05 and a.func is target.main


def test_few_cpg_regions_callable(tmp_path):
    # no minimum-CpG requirement: 28 and 3 CpGs are called whenever the t-test
    # is defined
    a, s, _, _ = _build(tmp_path, sample_shift=-0.45)
    tb = tmp_path / "t28.bed"
    tb.write_text("chr1\t1000\t1280\tR28\tG\tDisease A (hypo)\n"
                  "chr1\t1000\t1030\tR3\tG\tDisease A (hypo)\n"
                  "chr1\t1000\t1010\tR1\tG\tDisease A (hypo)\n")
    df = target.target_report([{
        "sample_id": "S1",
        "type": "CASE",
        "merged": s
    }],
                              a,
                              str(tb),
                              icr_gene_list=os.devnull).set_index("REGION_NAME")
    assert df.loc["R28", "NUM_CPG"] == 28 and df.loc["R28",
                                                     "DIRECTION"] == "hypo"
    assert df.loc["R3", "NUM_CPG"] == 3 and df.loc[
        "R3", "DIRECTION"] == "hypo" and df.loc["R3", "P_FDR"] < 0.01
    # one CpG: t-test undefined -> NA
    assert df.loc["R1", "NUM_CPG"] == 1 and pd.isna(
        df.loc["R1", "DIRECTION"]) and np.isnan(df.loc["R1", "P_FDR"])


def test_region_without_disease_direction_still_gets_direction(tmp_path):
    a, s, _, _ = _build(tmp_path, sample_shift=-0.45)
    tb = tmp_path / "td.bed"
    tb.write_text(
        "chr1\t1000\t1400\tRD\tG\tAngelman syndrome (locus, methylation not "
        "diagnostic)\n"
        "chr1\t1000\t1200\tRU\tG\tDisease B (hyper);Disease C (hypo)\n")
    df = target.target_report([{
        "sample_id": "S1",
        "type": "CASE",
        "merged": s
    }],
                              a,
                              str(tb),
                              icr_gene_list=os.devnull).set_index("REGION_NAME")
    assert df.loc["RD", "DIRECTION"] == "hypo" and df.loc["RD", "P_FDR"] < 0.01
    assert df.loc["RU", "DIRECTION"] == "hypo"


def test_no_usable_cpg_writes_report_and_warns_once(tmp_path):
    # a BED on a contig that the atlas and the sample lack used to stop with
    # "You are trying to merge on object and float64 columns for key 'CHROM'"
    a, s, _, _ = _build(tmp_path)
    tb = tmp_path / "absent.bed"
    tb.write_text("chr5\t1000\t2000\tR5\tG5\tDisease A (hypo)\n"
                  "chr6\t1000\t2000\tR6\tG6\tDisease B (hyper)\n")
    out = tmp_path / "t.tsv"
    msgs = []
    df = target.target_report([{
        "sample_id": "S1",
        "type": "CASE",
        "merged": s
    }],
                              a,
                              str(tb),
                              str(out),
                              icr_gene_list=os.devnull,
                              log=msgs.append)
    assert list(df.columns) == target.OUT_COLS
    assert set(df.REGION_NAME) == {"R5", "R6"} and (df.NUM_CPG == 0).all()
    assert df.DIRECTION.isna().all() and df["P_FDR"].isna().all()
    lines = out.read_text().splitlines()
    assert lines[0] == "\t".join(target.OUT_COLS) and len(lines) == 3
    warns = [m for m in msgs if m.startswith("WARNING")]
    assert len(warns) == 1 and "no target region" in warns[0]
    assert "S1" in warns[0]


def test_atlas_without_cpg_in_regions_has_typed_columns():
    atl = target._unique_atlas([])
    assert atl["CHROM"].dtype == object and atl["START"].dtype == np.int64


_PBCPG_HEAD = ("chrom\tbegin\tend\tmod_score\ttype\tcov\test_mod_count\t"
               "est_unmod_count\tdiscretized_mod_score\n")


def _fmt_line(fmt, r):
    c, s, e, cov, n = r
    if fmt == "pbcpg":
        return (f"{c}\t{s}\t{e}\t{100 * n / cov:.2f}\tTotal\t{cov}\t{n}\t"
                f"{cov - n}\t{n / cov:.3f}\n")
    if fmt in ("modkit", "modkit_sp"):
        # modkit writes the columns from the 10th on after spaces
        sep = " " if fmt == "modkit_sp" else "\t"
        tail = sep.join(
            map(str, [cov, f"{100 * n / cov:.2f}", n, cov - n, 0, 0, 0, 0, 0]))
        return (f"{c}\t{s}\t{e}\tm\t{cov}\t+\t{s}\t{e}\t255,0,0\t"
                f"{tail}\n")
    return f"{c}\t{s}\t{e}\t{cov}\t{n}\n"


def _write_track(path, src, fmt, packing="bgzf"):
    """Writes the rows of the 5-column BED `src` in `fmt` without an index."""
    rows = [(c, int(s), int(e), int(v), int(n)) for c, s, e, v, n in (
        ln.split("\t") for ln in dio.open_text(src).read().splitlines())]
    text = (_PBCPG_HEAD if fmt == "pbcpg" else "") + "".join(
        _fmt_line(fmt, r) for r in rows)
    path = str(path)
    if packing == "bgzf":
        dio.write_bgzip_text(path, [text])
    elif packing == "gzip":
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)
    return path


def _report(tmp_path, a, tb, merged, hap1=None, hap2=None, vcf=None):
    rec = {"sample_id": "S1", "type": "CASE", "merged": merged}
    for k, v in (("hap1", hap1), ("hap2", hap2), ("vcf", vcf)):
        if v:
            rec[k] = v
    return target.target_report([rec],
                                a,
                                tb,
                                scratch=str(tmp_path / "scratch"),
                                log=lambda m: None)


@pytest.mark.parametrize("route", ["default", "tabix"])
@pytest.mark.parametrize("fmt,packing",
                         [("pbcpg", "bgzf"), ("pbcpg", "gzip"),
                          ("pbcpg", "plain"), ("modkit", "bgzf"),
                          ("modkit", "gzip"), ("modkit", "plain"),
                          ("modkit_sp", "bgzf"), ("modkit_sp", "plain"),
                          ("5col", "gzip"), ("5col", "plain")])
def test_formats_give_the_5col_table(tmp_path, monkeypatch, route, fmt,
                                     packing):
    a, s, tb, _ = _build(tmp_path, sample_shift=-0.45)
    want = _report(tmp_path, a, tb, s)
    if route == "tabix":
        monkeypatch.setenv("DMRDETECTOR_NO_JIT", "1")
    m = _write_track(tmp_path / f"m.{fmt}.{packing}.bed.gz", s, fmt, packing)
    got = _report(tmp_path, a, tb, m)
    assert want.NUM_CPG.max() == 40
    pd.testing.assert_frame_equal(got, want)
    if packing == "bgzf":  # a user-indexed file: the reads go through its .tbi
        if fmt == "pbcpg":
            pysam.tabix_index(m,
                              seq_col=0,
                              start_col=1,
                              end_col=2,
                              zerobased=True,
                              line_skip=1,
                              force=True,
                              keep_original=True)
        else:
            dio.tabix_bed(m)
        pd.testing.assert_frame_equal(_report(tmp_path, a, tb, m), want)


@pytest.mark.parametrize("route", ["default", "tabix"])
@pytest.mark.parametrize("fmt", ["pbcpg", "modkit"])
def test_haplotype_tracks_in_any_format(tmp_path, monkeypatch, route, fmt):
    a, s, tb, pos = _build(tmp_path, sample_shift=-0.45)
    # same file names in different folders: the reads must not share an index
    h = [("chr1", int(p), int(p) + 1, 5 + i % 3, 1) for i, p in enumerate(pos)]
    for d in ("hap1", "hap2", "m5", "mf", "c1", "c2"):
        (tmp_path / d).mkdir()
    h1 = _bgzip_bed(tmp_path / "hap1" / "S1.bed.gz", h)
    h2 = _bgzip_bed(tmp_path / "hap2" / "S1.bed.gz", h[::-1][:20][::-1])
    m5 = _bgzip_bed(tmp_path / "m5" / "S1.bed.gz", [
        (c, int(a), int(b), int(v), int(n)) for c, a, b, v, n in (
            ln.split("\t") for ln in dio.open_text(s).read().splitlines())
    ])
    want = _report(tmp_path, a, tb, m5, h1, h2)
    assert set(want.HAPLOTYPE) == {"merged", "hap1", "hap2"}
    if route == "tabix":
        monkeypatch.setenv("DMRDETECTOR_NO_JIT", "1")
    conv = {
        k: _write_track(tmp_path / d / "S1.bed.gz", v, fmt)
        for k, d, v in (("hap1", "c1", h1), ("hap2", "c2", h2))
    }
    mf = _write_track(tmp_path / "mf" / "S1.bed.gz", s, fmt)
    pd.testing.assert_frame_equal(_report(tmp_path, a, tb, mf, **conv), want)
    # mixed: a 5-column combined BED with pb-CpG-tools haplotype BEDs
    pd.testing.assert_frame_equal(_report(tmp_path, a, tb, m5, **conv), want)


@pytest.mark.parametrize("route", ["default", "tabix"])
def test_track_without_data_lines_gives_empty_regions(tmp_path, monkeypatch,
                                                      route):
    a, _, tb, _ = _build(tmp_path)
    if route == "tabix":
        monkeypatch.setenv("DMRDETECTOR_NO_JIT", "1")
    m = tmp_path / "empty.bed.gz"
    dio.write_bgzip_text(str(m), [_PBCPG_HEAD])
    df = _report(tmp_path, a, tb, str(m))
    assert set(df.REGION_NAME) == {"R1", "R2"} and (df.NUM_CPG == 0).all()
