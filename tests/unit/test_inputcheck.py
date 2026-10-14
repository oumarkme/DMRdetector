"""Input checks: hg38 coordinates and modkit rows that count a CpG twice.

Synthetic data only. The modkit tests also document what the reader does with
such rows: it takes columns 10 and 12 of every row and ignores the modification
code (column 4) and the strand (column 6).
"""

import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from dmrdetector import cli
from dmrdetector import run
from dmrdetector.fileio import bed as dio
from dmrdetector.fileio import methbed
from dmrdetector.runtime import RunError
from dmrdetector.steps import preprocess
from dmrdetector.steps import target

from helpers import read_bed5


ROOT = Path(__file__).resolve().parents[2]
SRC = os.environ.get("DMRDETECTOR_SRC", str(ROOT / "src"))
POS = [10468 + 7 * i for i in range(40)]  # C positions (0-based) of 40 CpGs
COV, MOD = 20, 15


def modkit_row(chrom, start, code="m", strand="+", cov=COV, mod=MOD):
    """Returns one bedMethyl row.

    It has 9 BED columns, then Nvalid_cov, the fraction modified, Nmod,
    Ncanonical, Nother_mod, ...
    """
    return "\t".join([
        chrom,
        str(start),
        str(start + 1), code,
        str(cov), strand,
        str(start),
        str(start + 1), "255,0,0",
        str(cov), f"{100 * mod / cov:.2f}",
        str(mod),
        str(cov - mod), "0", "0", "0", "0", "0"
    ])


def write(path, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("".join(r + "\n" for r in rows))
    return str(path)


def vcf(path):
    Path(path).write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )
    return str(path)


def pre(bed, tmp_path, name="o"):
    out = str(tmp_path / f"{name}.bed.gz")
    r = preprocess.preprocess_sample(bed,
                                     vcf(tmp_path / "s.vcf"),
                                     5,
                                     out,
                                     log=lambda m: None,
                                     threads=2)
    return r, read_bed5(out)


@pytest.fixture
def combined(tmp_path):
    return write(tmp_path / "combined.bed.gz",
                 [modkit_row("chr1", p, strand=".") for p in POS])


def test_reader_counts_5hmc_rows_as_extra_sites(tmp_path):
    """Documents that an `h` row at the same position becomes a second row."""
    rows = []
    for p in POS:
        rows += [modkit_row("chr1", p), modkit_row("chr1", p, code="h", mod=2)]
    df, fmt = methbed.read_methylation_bed(write(tmp_path / "mh.bed.gz", rows))
    assert fmt == "modkit" and len(df) == 2 * len(POS)
    assert df["START"].duplicated().sum() == len(POS)
    with pytest.raises(RunError, match=r"code.*\(h, m\).*--combine-strands"):
        preprocess.check_modkit(str(tmp_path / "mh.bed.gz"))


def test_reader_counts_strand_split_rows_as_extra_sites(tmp_path):
    """Documents the rows of a file written without --combine-strands.

    The minus-strand row sits on the G (START + 1) with its own, partial counts.
    """
    rows = []
    for p in POS:
        rows += [
            modkit_row("chr1", p, strand="+", cov=10, mod=7),
            modkit_row("chr1", p + 1, strand="-", cov=10, mod=8)
        ]
    f = write(tmp_path / "split.bed.gz", rows)
    df, _ = methbed.read_methylation_bed(f)
    assert len(df) == 2 * len(POS) and set(df["COVERAGE"]) == {10.0}
    with pytest.raises(RunError, match="one row per strand.*--combine-strands"):
        preprocess.check_modkit(f)


def test_preprocess_stops_on_both_cases_and_writes_nothing(tmp_path):
    f = write(tmp_path / "mh.bed.gz", [
        r for p in POS
        for r in (modkit_row("chr1", p), modkit_row("chr1", p, code="h"))
    ])
    with pytest.raises(RunError, match="modkit pileup --cpg --combine-strands"):
        pre(f, tmp_path)
    assert not (tmp_path / "o.bed.gz").exists()


def test_correct_modkit_files_pass_and_match_the_5col_table(tmp_path, combined):
    """Checks that correct modkit files pass and give the 5-column table.

    Combined strands, a single `+` strand (as the test files) and a repeated `m`
    row all pass.
    """
    plus = write(tmp_path / "plus.bed.gz", [modkit_row("chr1", p) for p in POS])
    five = write(tmp_path / "five.bed.gz",
                 [f"chr1\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS])
    ref = pre(five, tmp_path, "five_out")[1]
    for f in (combined, plus):
        preprocess.check_modkit(f)
        r, got = pre(f, tmp_path)
        assert r["format"] == "modkit" and got.equals(ref)
    dup = write(tmp_path / "dup.bed.gz",
                [modkit_row("chr1", p, strand=".") for p in POS] +
                [modkit_row("chr1", POS[0], strand=".")])
    preprocess.check_modkit(
        dup)  # same code twice is not evidence of 5hmC or strand rows


def test_modkit_check_reads_plain_text_and_skips_headers(tmp_path):
    p = tmp_path / "h.bed"
    p.write_text("#chrom\tstart\n" +
                 "".join(modkit_row("chr1", q, strand=".") + "\n" for q in POS))
    preprocess.check_modkit(str(p))


def test_hg38_check():
    preprocess.check_hg38(
        {
            "chr1": 248956422,
            "1": 248956422,
            "chrX": 156040895,
            "chrM": 99999999,
            "chr1_KI270706v1_random": 10**9
        }, "f")
    preprocess.check_hg38({"chr9": 138394717 + preprocess.END_MARGIN}, "f")
    for name, end in (("chr1", 249_240_000), ("9", 141_000_000),
                      ("chrY", 59_000_000)):  # hg19-sized ends
        with pytest.raises(RunError,
                           match=r"hg38/GRCh38.*hg19/GRCh37.*not supported"):
            preprocess.check_hg38({name: end}, "f.bed.gz")


def test_hg38_table_matches_the_pipeline_contigs():
    assert set(preprocess.HG38_LENGTHS) == {
        f"chr{c}" for c in [*range(1, 23), "X", "Y"]
    }


@pytest.mark.parametrize("fast", [True, False])
def test_preprocess_stops_on_hg19_coordinates(tmp_path, monkeypatch, fast):
    from dmrdetector.fileio import methbed
    if not fast:
        monkeypatch.setattr(methbed, "read_filtered", lambda *a, **k: None)
    rows = [f"chr9\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS
           ] + [f"chr9\t141100000\t141100001\t{COV}\t{MOD}"]
    with pytest.raises(
            RunError,
            match=r"chr9 has a record ending at 141,100,001.*138,394,717"):
        pre(write(tmp_path / "hg19.bed.gz", rows), tmp_path)
    assert not (tmp_path / "o.bed.gz").exists()


def test_build_prints_the_error_from_a_worker_pool(tmp_path):
    """Checks the error of a worker in the preprocessing pool of `build`.

    The message must reach the user as `ERROR: ...` with exit status 1.
    """
    v = vcf(tmp_path / "s.vcf")
    good = write(tmp_path / "good.bed.gz",
                 [f"chr1\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS])
    bad = write(tmp_path / "bad.bed.gz", [
        r for p in POS
        for r in (modkit_row("chr1", p), modkit_row("chr1", p, code="h"))
    ])
    sheet = tmp_path / "s.csv"
    sheet.write_text("sample_id,sex,meth_bed,vcf\n" +
                     "".join(f"{s},{x},{b},{v}\n"
                             for s, x, b in (("A", "F", good), ("B", "M", good),
                                             ("C", "F", good),
                                             ("D", "M", bad))))
    r = subprocess.run([
        sys.executable, "-m", "dmrdetector.cli", "build", "--input",
        str(sheet), "--outdir",
        str(tmp_path / "ref"), "--threads", "8", "--set",
        "include_sex_chromosomes=false"
    ],
                       capture_output=True,
                       text=True,
                       env=dict(os.environ, PYTHONPATH=SRC))
    assert r.returncode == 1, r.stderr[-1500:]
    assert 'ERROR:' in r.stderr, r.stderr[-1500:]
    assert '--combine-strands' in r.stderr, r.stderr[-1500:]
    assert 'Traceback' not in r.stderr, r.stderr[-1500:]


# ---- checks of the haplotype tracks and the stand-alone target command -----


def _target_inputs(tmp_path):
    from test_target import _build
    return _build(tmp_path, sample_shift=-0.2)


def _bgzip_rows(path, rows, index=True):
    """Writes sorted text rows as BGZF, with a .tbi unless `index` is False."""
    text = "".join(r + "\n" for r in sorted(
        rows, key=lambda r: (r.split("\t")[0], int(r.split("\t")[1]))))
    if index:
        dio.write_bed_bgzip_tabix(str(path), text)
    else:
        dio.write_bgzip_text(str(path), [text])
    return str(path)


def _target_report(atlas, merged, tb, hap1, **kw):
    return target.target_report([{
        "sample_id": "S1",
        "type": "CASE",
        "merged": merged,
        "hap1": hap1
    }],
                                atlas,
                                tb,
                                log=lambda m: None,
                                threads=2,
                                **kw)


@pytest.mark.parametrize("index", [True, False])
def test_a_haplotype_track_with_5hmc_rows_is_refused(tmp_path, index):
    atlas, merged, tb, pos = _target_inputs(tmp_path)
    rows = [
        modkit_row("chr1", int(p), code=c, strand=".")
        for p in pos
        for c in ("m", "h")
    ]
    hap1 = _bgzip_rows(tmp_path / "S1.hap1.bed.gz", rows, index)
    with pytest.raises(RunError, match=r"\(h, m\).*--combine-strands"):
        _target_report(atlas, merged, tb, hap1)


@pytest.mark.parametrize("index", [True, False])
def test_a_haplotype_track_beyond_the_hg38_length_is_refused(tmp_path, index):
    atlas, merged, tb, pos = _target_inputs(tmp_path)
    rows = [f"chr9\t{int(p)}\t{int(p) + 1}\t20\t10" for p in pos
           ] + ["chr9\t141100000\t141100001\t20\t10"]
    hap1 = _bgzip_rows(tmp_path / "S1.hap1.bed.gz", rows, index)
    with pytest.raises(RunError,
                       match=r"chr9 has a record ending at 141,100,001"):
        _target_report(atlas, merged, tb, hap1)


def test_a_valid_modkit_haplotype_track_passes_the_checks(tmp_path):
    atlas, merged, tb, pos = _target_inputs(tmp_path)
    five = _bgzip_rows(tmp_path / "five.bed.gz",
                       [f"chr1\t{int(p)}\t{int(p) + 1}\t20\t12" for p in pos])
    modkit = _bgzip_rows(
        tmp_path / "S1.hap1.bed.gz",
        [modkit_row("chr1", int(p), strand=".", cov=20, mod=12) for p in pos])
    a = _target_report(atlas, merged, tb, five)
    b = _target_report(atlas, merged, tb, modkit)
    assert a.equals(b)


def test_the_target_command_prints_the_check_error_and_exits_with_1(
        tmp_path, capsys):
    atlas, merged, tb, pos = _target_inputs(tmp_path)
    rows = [
        modkit_row("chr1", int(p), code=c, strand=".")
        for p in pos
        for c in ("m", "h")
    ]
    hap1 = _bgzip_rows(tmp_path / "S1.hap1.bed.gz", rows)
    sheet = tmp_path / "s.csv"
    sheet.write_text("sample_id,type,meth_bed,vcf,meth_bed_hap1\n"
                     f"S1,CASE,{merged},,{hap1}\n")
    rc = cli.main([
        "target", "--samplesheet",
        str(sheet), "--atlas", atlas, "--target-bed", tb, "--out",
        str(tmp_path / "t.tsv")
    ])
    assert rc == 1 and "ERROR:" in capsys.readouterr().err
    assert not (tmp_path / "t.tsv").exists()


def test_the_target_report_honors_input_format(tmp_path):
    atlas, merged, tb, pos = _target_inputs(tmp_path)
    five = [f"chr1\t{int(p)}\t{int(p) + 1}\t20\t12" for p in pos]
    # eight columns: auto takes the layout of pb-CpG-tools (columns 6 and 7)
    eight = [f"{r}\t100\t0\t1" for r in five]
    t5 = _bgzip_rows(tmp_path / "t5.bed.gz", five, index=False)
    t8 = _bgzip_rows(tmp_path / "t8.bed.gz", eight, index=False)
    want = _target_report(atlas, merged, tb, t5)
    auto = _target_report(atlas, merged, tb, t8)
    forced = _target_report(atlas, merged, tb, t8, input_format="5col")
    h = want["HAPLOTYPE"] == "hap1"
    assert forced[h].equals(want[h])
    assert not auto[auto["HAPLOTYPE"] == "hap1"].equals(want[h])
    assert target._cov_cols(8,
                            "5col") == (3, 4) and target._cov_cols(8) == (5, 6)
    with pytest.raises(RunError, match="cannot read as input_format modkit"):
        _target_report(atlas, merged, tb, t5, input_format="modkit")


def test_a_stale_index_of_a_track_is_rebuilt_in_scratch(tmp_path):
    atlas, merged, tb, pos = _target_inputs(tmp_path)
    path = tmp_path / "S1.hap1.bed.gz"
    _bgzip_rows(path, [f"chr1\t{int(p)}\t{int(p) + 1}\t20\t12" for p in pos])
    scratch = str(tmp_path / "scratch")
    assert target._indexable(str(path), scratch) == str(path)  # fresh
    old = os.path.getmtime(str(path)) - 100
    os.utime(str(path) + ".tbi", (old, old))  # index older than the data
    got = target._indexable(str(path), scratch)
    assert got != str(path) and got.startswith(scratch)
    assert os.path.exists(got + ".tbi")


def test_scratch_links_of_equal_names_do_not_collide(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir()
    b.mkdir()
    pa = _bgzip_rows(a / "x.bed.gz", ["chr1\t1\t2\t5\t1"], index=False)
    pb = _bgzip_rows(b / "x.bed.gz", ["chr2\t1\t2\t5\t1"], index=False)
    scratch = str(tmp_path / "s")
    la, lb = dio.ensure_tabix(pa, scratch), dio.ensure_tabix(pb, scratch)
    assert la != lb and os.path.basename(la) == os.path.basename(
        lb) == "x.bed.gz"
    import pysam
    assert pysam.TabixFile(la).contigs == ["chr1"]
    assert pysam.TabixFile(lb).contigs == ["chr2"]


# ---- warnings about input that leaves nothing to compare ---------------------


def _task(tmp_path, bed_path, sex, split=True, min_cov=5):
    (tmp_path / "split").mkdir(exist_ok=True)
    (tmp_path / "stamps").mkdir(exist_ok=True)
    return dict(id="R",
                meth_bed=bed_path,
                vcf=vcf(tmp_path / "s.vcf"),
                min_cov=min_cov,
                fmt="auto",
                out=str(tmp_path / "R.bed.gz"),
                stamps=str(tmp_path / "stamps"),
                split=split,
                split_dir=str(tmp_path / "split"),
                sex=sex,
                par_bed=str(preprocess.data.asset("par_hg38.bed")),
                threads=1)


def test_warning_when_no_cpg_is_left_after_filtering(tmp_path, capsys):
    rows = [f"chr1\t{p}\t{p + 1}\t3\t1" for p in POS]  # coverage 3 <= 5
    preprocess._t_preprocess(
        _task(tmp_path, write(tmp_path / "low.bed.gz", rows), "F"))
    err = capsys.readouterr().err
    assert ("WARNING: R: no CpG is left after filtering (40 in; 0 on het "
            "SNPs; coverage must exceed 5)") in err


def test_warning_for_refseq_contig_names(tmp_path, capsys):
    rows = [f"NC_000001.11\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS]
    preprocess._t_preprocess(
        _task(tmp_path, write(tmp_path / "refseq.bed.gz", rows), "F"))
    err = capsys.readouterr().err
    assert "WARNING: R: no CpG on chr1-chr22 or the PARs, but 40 on other " \
        "contigs" in err and "NC_000001.11" in err


def test_warning_for_an_empty_sex_scope(tmp_path, capsys):
    rows = [f"chr1\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS]
    f = write(tmp_path / "auto.bed.gz", rows)
    preprocess._t_preprocess(_task(tmp_path, f, "M"))
    err = capsys.readouterr().err
    assert "no CpG after filtering in the male_X scope" in err
    assert "no CpG after filtering in the Y scope" in err
    assert "no CpG on chr1-chr22" not in err


def test_split_warnings_and_the_sex_check():
    n = preprocess.CHRY_MALE_MIN_CPG
    counts = {"autosomal": 100, "chrX": 10, "chrY": n}
    f = preprocess.split_warnings("S", "F", counts, 100 + 10 + n)
    assert len(f) == 1 and f"{n:,} CpGs remain on chrY" in f[0]
    assert "check the sex" in f[0]
    counts["chrY"] = n - 1
    assert preprocess.split_warnings("S", "F", counts, 100 + 10 + n - 1) == []
    counts["chrY"] = n
    assert preprocess.split_warnings("S", "M", counts, 100 + 10 + n) == []
    assert preprocess.split_warnings("S", "", counts, 100 + 10 + n) == []
    assert preprocess.split_warnings("S", "F", None, 5) == []


def test_no_match_with_the_atlas_is_an_error_in_the_autosomal_scope():
    stats = dict(sample_sites=0, atlas_sites=57, merged=0)
    msg = run.check_matched("S1", "autosomal", stats)
    assert "S1 autosomal: 0 sample CpGs match the reference atlas" in msg
    assert "sample CpGs in the panel regions: 0" in msg
    assert "atlas sites in the panel regions: 57" in msg
    assert "RefSeq" in run.NO_MATCH_CAUSES and "hg38" in run.NO_MATCH_CAUSES
    assert run.check_matched("S1", "Y", dict(stats, merged=3)) is None


def _no_match_task(tmp_path, scope, contig):
    from helpers import write_call_dmr_inputs
    d = tmp_path / "syn"
    write_call_dmr_inputs(str(d))
    sample = f"{d}/{scope}.sample.bed.gz"
    rows = read_bed5(sample)
    rows["CHROM"] = contig
    rows = rows.drop_duplicates("START")
    new = str(tmp_path / f"{scope}.other.bed.gz")
    _bgzip_rows(new,
                ["\t".join(map(str, r)) for r in rows.itertuples(index=False)])
    (tmp_path / "stamps").mkdir(exist_ok=True)
    return dict(id="S1",
                scope=scope,
                bed=new,
                atlas=f"{d}/{scope}.atlas.bed.gz",
                panel=f"{d}/{scope}.cpg.bed",
                icr_panel=f"{d}/icr.bed" if scope == "autosomal" else None,
                out=str(tmp_path / f"{scope}.raw.tsv.gz"),
                include_sex=True,
                tabix_duplicates=True,
                stamps=str(tmp_path / "stamps"),
                threads=1,
                thr={})


def test_run_stops_when_no_cpg_matches_the_autosomal_atlas(tmp_path):
    t = _no_match_task(tmp_path, "autosomal", "NC_000001.11")
    with pytest.raises(RunError) as e:
        run._t_call_dmr(t)
    assert "S1 autosomal: 0 sample CpGs match the reference atlas" in str(
        e.value)
    assert "RefSeq" in str(e.value)


def test_run_warns_when_no_cpg_matches_a_sex_atlas(tmp_path, capsys):
    t = _no_match_task(tmp_path, "male_X", "NC_000023.11")
    steps = run._t_call_dmr(t)
    assert steps[0][1] == "done"
    err = capsys.readouterr().err
    assert "WARNING: S1 male_X: 0 sample CpGs match the reference atlas" in err


# ---- the VCF and the BED must name contigs alike -----------------------------


def _vcf_text(contig, samples=("S",), pos=POS[3] + 1):
    gts = "\t".join("0|1" for _ in samples)
    return ("##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
            "\t".join(samples) + "\n"
            f"{contig}\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t{gts}\n")


@pytest.mark.parametrize("fast", [True, False])
def test_warning_when_the_vcf_names_contigs_the_bed_does_not(
        tmp_path, monkeypatch, fast):
    if not fast:
        monkeypatch.setattr(methbed, "read_filtered", lambda *a, **k: None)
    f = write(tmp_path / "chr.bed.gz",
              [f"chr1\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS])
    msgs = []

    def go(contig):
        v = tmp_path / "s.vcf"
        v.write_text(_vcf_text(contig))
        return preprocess.preprocess_sample(f,
                                            str(v),
                                            5,
                                            str(tmp_path / "o.bed.gz"),
                                            log=msgs.append,
                                            threads=1)

    r = go("1")
    assert r["rows_out"] == len(POS) and r["removed_by_het"] == 0  # as v2.2
    assert any("none of these names is a contig of the methylation BED" in m and
               "with and without the 'chr' prefix" in m for m in msgs)
    msgs.clear()
    r = go("chr1")
    assert r["removed_by_het"] == 1  # the SNP sits on the fourth CpG
    assert not any("none of these names" in m for m in msgs)


def test_warning_for_a_vcf_with_several_samples(tmp_path):
    f = write(tmp_path / "chr.bed.gz",
              [f"chr1\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS])
    msgs = []
    for samples, warned in ((("A",), False), (("A", "B"), True)):
        v = tmp_path / "s.vcf"
        v.write_text(_vcf_text("chr1", samples))
        msgs.clear()
        preprocess.preprocess_sample(f,
                                     str(v),
                                     5,
                                     str(tmp_path / "o.bed.gz"),
                                     log=msgs.append,
                                     threads=1)
        hits = [m for m in msgs if "sample columns" in m]
        assert bool(hits) == warned
        if warned:
            assert "WARNING: the VCF has 2 sample columns; all samples are " \
                "used" in hits[0]


HEAD = ("##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n')


@pytest.mark.parametrize("fast", [True, False])
def test_warning_for_a_vcf_that_cannot_remove_any_cpg(tmp_path, monkeypatch,
                                                      fast):
    if not fast:
        monkeypatch.setattr(methbed, "read_filtered", lambda *a, **k: None)
    f = write(tmp_path / "chr.bed.gz",
              [f"chr1\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS])
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    rec = f"chr1\t{POS[3] + 1}\t.\tA\tG\t50\tPASS\t."
    cases = {  # name: (VCF text, warning that must appear, rows removed)
        "sites-only": (HEAD + cols + "\n" + rec + "\n", "no sample column", 0),
        "no heterozygous call":
            (HEAD + cols + "\tFORMAT\tS\n" + rec + "\tGT\t1/1\n" +
             rec.replace("A\tG", "AT\tA") + "\tGT\t0/1\n",
             "no heterozygous SNP call", 0),
        "a heterozygous SNP": (HEAD + cols + "\tFORMAT\tS\n" + rec +
                               "\tGT\t0|1\n", None, 1),
    }
    for name, (text, want, removed) in cases.items():
        msgs = []
        v = tmp_path / "s.vcf"
        v.write_text(text)
        r = preprocess.preprocess_sample(f,
                                         str(v),
                                         5,
                                         str(tmp_path / "o.bed.gz"),
                                         log=msgs.append,
                                         threads=1)
        assert r["removed_by_het"] == removed, name
        hits = [m for m in msgs if "WARNING: the VCF has no" in m]
        if want is None:
            assert not hits, name
        else:
            assert len(hits) == 1 and want in hits[0], (name, msgs)


def _chrx_vcf(path, n_nonpar, n_par):
    """Writes a one-sample VCF with heterozygous SNPs on chrX."""
    lines = [HEAD + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
    # 1-based positions: non-PAR after PAR1 (ends at 2,781,479), PAR1 itself
    pos = [3_000_000 + 7 * i for i in range(n_nonpar)]
    pos += [100_000 + 10 * i for i in range(n_par)]
    lines += [f"chrX\t{p}\t.\tC\tT\t50\tPASS\t.\tGT\t0|1" for p in pos]
    Path(path).write_text("\n".join(lines) + "\n")
    return str(path)


def test_nonpar_chrx_het_count_follows_the_par_bed():
    par = preprocess.read_par_bed(preprocess.data.asset("par_hg38.bed"))
    pos = np.array([
        1000, 10_000, 2_781_478, 2_781_479, 5_000_000, 155_701_381, 155_701_382,
        156_030_894, 156_030_895
    ])  # 0-based; PAR1 is [10000, 2781479), PAR2 [155701382, 156030895)
    assert preprocess.nonpar_chrx_het_count({
        "chrX": pos,
        "chr1": pos
    }, par) == 5
    assert preprocess.nonpar_chrx_het_count({"X": pos}, par) == 5
    assert preprocess.nonpar_chrx_het_count({"chr1": pos}, par) == 0
    # overlapping or repeated PAR lines remove a site once
    twice = {
        "chrX": (np.repeat(par["chrX"][0], 2), np.repeat(par["chrX"][1], 2))
    }
    assert preprocess.nonpar_chrx_het_count({"chrX": pos}, twice) == 5
    assert preprocess.nonpar_chrx_het_count({"chrx": pos}, par) == 0  # not X
    # a chrX name that the PAR file does not know: the count is not made
    assert preprocess.nonpar_chrx_het_count({"chrX": pos},
                                            {"X": par["X"]}) is None


def test_sex_label_check_uses_one_cutoff_in_both_directions():
    w = preprocess.sex_het_warning
    cut = preprocess.FEMALE_CHRX_HET_MIN
    assert cut == 40_000
    assert w("S", "M", cut - 1) is None and "sex M" in w("S", "M", cut)
    assert w("S", "F", cut) is None and "sex F" in w("S", "F", cut - 1)
    assert w("S", "", 0) is None and w("S", "F", None) is None
    assert "41,000 heterozygous SNPs on chrX" in w("S", "M", 41_000)


def _pre_task(tmp_path, bed, vcf, sex):
    for d in ("stamps", "split"):
        (tmp_path / d).mkdir(exist_ok=True)
    return dict(id="S1",
                meth_bed=bed,
                vcf=vcf,
                min_cov=5,
                fmt="auto",
                out=str(tmp_path / "S1_filtered.bed.gz"),
                stamps=str(tmp_path / "stamps"),
                split=True,
                split_dir=str(tmp_path / "split"),
                sex=sex,
                par_bed=preprocess.data.asset("par_hg38.bed"),
                threads=1)


@pytest.mark.parametrize("fast", [True, False])
def test_sex_label_is_checked_against_the_chrx_snps(tmp_path, monkeypatch,
                                                    capsys, fast):
    if not fast:
        monkeypatch.setattr(methbed, "read_filtered", lambda *a, **k: None)
    bed = write(tmp_path / "x.bed.gz", [
        f"chrX\t{p}\t{p + 1}\t{COV}\t{MOD}"
        for p in range(3_000_000, 3_000_300, 7)
    ])
    female = _chrx_vcf(tmp_path / "f.vcf", 41_000, 300)  # PAR sites not counted
    male = _chrx_vcf(tmp_path / "m.vcf", 19_000, 500)
    r = preprocess.preprocess_sample(
        bed,
        female,
        5,
        str(tmp_path / "o.bed.gz"),
        log=lambda m: None,
        threads=1,
        split_prefix=str(tmp_path / "sp"),
        par_bed=preprocess.data.asset("par_hg38.bed"))
    assert r["het_chrx"] == 41_000
    # without a split the count is not made (no sex is needed then)
    r = preprocess.preprocess_sample(bed,
                                     male,
                                     5,
                                     str(tmp_path / "o2.bed.gz"),
                                     log=lambda m: None,
                                     threads=1)
    assert r["het_chrx"] is None
    for sex, vcf, warned in (("F", female, False), ("M", female, True),
                             ("M", male, False), ("F", male, True)):
        preprocess._t_preprocess(_pre_task(tmp_path, bed, vcf, sex))
        err = capsys.readouterr().err
        assert (f"S1: sex {sex}, but the VCF has" in err) == warned, (sex, err)
    # a VCF with two samples is not counted: a site counts for any sample
    two = tmp_path / "two.vcf"
    two.write_text(
        Path(female).read_text().replace("\tS\n", "\tS\tT\n").replace(
            "\t0|1", "\t0|1\t0|0"))
    r = preprocess.preprocess_sample(
        bed,
        str(two),
        5,
        str(tmp_path / "o3.bed.gz"),
        log=lambda m: None,
        threads=1,
        split_prefix=str(tmp_path / "sp3"),
        par_bed=preprocess.data.asset("par_hg38.bed"))
    assert r["het_chrx"] is None


@pytest.mark.parametrize("order", ["descending", "interleaved"])
def test_an_unsorted_track_is_a_clear_error_of_the_target_report(
        tmp_path, order):
    atlas, merged, tb, pos = _target_inputs(tmp_path)
    rows = [f"chr1\t{int(p)}\t{int(p) + 1}\t20\t12\n" for p in pos]
    if order == "descending":
        rows.reverse()
    else:  # chr1, chr2, chr1: the contig is not in one block
        rows[1:1] = ["chr2\t5000\t5001\t20\t12\n"]
    track = str(tmp_path / "S1.hap1.bed.gz")
    dio.write_bgzip_text(track, ["".join(rows)])
    with pytest.raises(RunError,
                       match="cannot index the track.*sorted by contig and "
                       "START"):
        _target_report(atlas, merged, tb, track)


def test_a_stale_index_gives_the_same_report(tmp_path):
    atlas, merged, tb, pos = _target_inputs(tmp_path)
    hap1 = _bgzip_rows(tmp_path / "S1.hap1.bed.gz",
                       [f"chr1\t{int(p)}\t{int(p) + 1}\t20\t12" for p in pos])
    fresh = _target_report(atlas, merged, tb, hap1)
    old = os.path.getmtime(hap1) - 100
    for f in (hap1, merged, atlas):
        os.utime(f + ".tbi", (old, old))  # every index older than its data
    stale = _target_report(atlas, merged, tb, hap1)
    assert fresh.equals(stale)


def test_build_reads_the_control_rows_of_a_sheet_with_a_type_column(tmp_path):
    """Checks that 12 CTRL and 2 CASE rows give an atlas of 12 controls (v2.2).

    The CASE rows hold the same CpGs, so a build that read them would give a
    maximum N of 14.
    """
    v = vcf(tmp_path / "s.vcf")
    bed = write(tmp_path / "s.bed.gz",
                [f"chr1\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS])
    rows = [(f"C{i:02d}", "CTRL", "FM"[i % 2]) for i in range(12)]
    rows += [("P01", "CASE", "F"), ("P02", "case", "M")]
    sheet = tmp_path / "s.csv"
    sheet.write_text("sample_id,type,sex,meth_bed,vcf\n" +
                     "".join(f"{s},{t},{x},{bed},{v}\n" for s, t, x in rows))
    out = tmp_path / "ref"
    r = subprocess.run([
        sys.executable, "-m", "dmrdetector.cli", "build", "--input",
        str(sheet), "--outdir",
        str(out), "--threads", "8", "--set", "include_sex_chromosomes=false"
    ],
                       capture_output=True,
                       text=True,
                       env=dict(os.environ, PYTHONPATH=SRC))
    assert r.returncode == 0, r.stderr[-1500:]
    assert "NOTICE" in r.stderr and "2 of 14 rows are dropped" in r.stderr
    with gzip.open(out / "ref_atlas_150.bed.gz", "rt") as fh:
        n = {int(ln.split("\t")[3]) for ln in fh}
    assert n == {12}
    info = json.loads((out / "reference.json").read_text())
    assert [s["sample_id"] for s in info["samples"]
           ] == [f"C{i:02d}" for i in range(12)]
    assert info["scopes"]["autosomal"]["n_ctrl"] == 12


def test_records_longer_than_one_base_give_one_warning(tmp_path, monkeypatch):
    """Checks the count of END - START != 1 on the fast and the pandas path."""
    rows = [f"chr1\t{p}\t{p + (2 if i in (3, 8, 9) else 1)}\t{COV}\t{MOD}"
            for i, p in enumerate(POS)]
    wide = write(tmp_path / "wide.bed.gz", rows)
    one = write(tmp_path / "one.bed.gz",
                [f"chr1\t{p}\t{p + 1}\t{COV}\t{MOD}" for p in POS])
    for generic in (False, True):
        with monkeypatch.context() as m:
            if generic:
                m.setattr(methbed, "read_filtered", lambda *a, **k: None)
            for path, expected in ((wide, 3), (one, 0)):
                msgs = []
                out = str(tmp_path / "o.bed.gz")
                preprocess.preprocess_sample(path,
                                             vcf(tmp_path / "s.vcf"),
                                             5,
                                             out,
                                             log=msgs.append,
                                             threads=2)
                got = [x for x in msgs if "END - START" in x]
                assert len(got) == (expected > 0)
                if expected:
                    assert got[0].startswith(f"WARNING: {expected} of the kept")


def test_only_the_kept_records_are_counted_as_longer_than_one_base():
    start = np.array([1, 10, 20])
    assert preprocess.wide_records_warning(start, start + 1) is None
    assert preprocess.wide_records_warning(start[:0], start[:0]) is None
    assert preprocess.wide_records_warning(
        start, start + [1, 3, 1]).startswith("WARNING: 1 of the kept")


# Functions whose pandas reads hold only text and integer columns (panels, gene
# lists, PAR intervals, BED coordinates): no decimal is parsed there.
_READS_WITHOUT_DECIMALS = {
    ("steps/annotate.py", "override_icr_genes"),
    ("steps/call_dmr.py", "_read_panel_table"),
    ("steps/icr.py", "read_candidates"),
    ("steps/preprocess.py", "read_par_bed"),
    ("steps/preprocess.py", "_split_generic"),
    ("steps/target.py", "_read_intervals"),
}


def test_every_pandas_read_of_decimals_is_a_round_trip_parse():
    """Checks that no read_csv that parses decimals uses the default parser.

    The default parser of pandas can be one ulp off the correctly rounded value,
    and the fast readers are correctly rounded; the fallbacks must match them
    (the compiled and the plain-Python runs must write the same text).
    """
    import ast
    src = Path(SRC) / "dmrdetector"
    bad = []
    for path in sorted(src.rglob("*.py")):
        tree = ast.parse(path.read_text())
        rel = path.relative_to(src).as_posix()
        for fn in ast.walk(tree):
            if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for call in ast.walk(fn):
                if not (isinstance(call, ast.Call) and
                        isinstance(call.func, ast.Attribute) and
                        isinstance(call.func.value, ast.Name) and
                        call.func.value.id == "pd" and
                        call.func.attr in ("read_csv", "read_table")):
                    continue
                exact = any(
                    k.arg == "float_precision" and
                    getattr(k.value, "value", None) == "round_trip"
                    for k in call.keywords)
                if not exact and (rel, fn.name) not in _READS_WITHOUT_DECIMALS:
                    bad.append(f"{rel}:{call.lineno} in {fn.name}")
    assert not bad, bad


def test_resume_with_another_sex_label_checks_the_label_again(tmp_path, capsys):
    bed = write(tmp_path / "x.bed.gz", [
        f"chrX\t{p}\t{p + 1}\t{COV}\t{MOD}"
        for p in range(3_000_000, 3_000_300, 7)
    ])
    female = _chrx_vcf(tmp_path / "f.vcf", 41_000, 300)
    t = _pre_task(tmp_path, bed, female, "M")
    assert preprocess._t_preprocess(t)[0][1] == "done"
    assert "S1: sex M, but the VCF has 41,000" in capsys.readouterr().err
    # the same command: skipped, and the warning of the first run is repeated
    assert preprocess._t_preprocess(t)[0][1] == "skipped"
    assert "S1: sex M, but the VCF has 41,000" in capsys.readouterr().err
    # the label was corrected: the step runs again and the warning is gone
    t["sex"] = "F"
    assert preprocess._t_preprocess(t)[0][1] == "done"
    assert "sex M" not in capsys.readouterr().err
