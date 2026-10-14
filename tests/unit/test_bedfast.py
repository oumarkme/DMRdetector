"""Byte-level preprocess (methbed) against the pandas code and pysam's tabix.

Run `pytest tests/unit/test_bedfast.py`, also with DMRDETECTOR_NO_JIT=1
(interpreted kernels; slow: set BEDFAST_N=1).
"""
import gzip
import os

import pytest

from dmrdetector import data
from dmrdetector.fileio import bgzf as bz
from dmrdetector.fileio import methbed
from dmrdetector.steps import preprocess  # noqa: E402

from helpers import generate  # noqa: E402
from helpers import parse_tbi  # noqa: E402


PAR = data.asset("par_hg38.bed")
N = int(os.environ.get("BEDFAST_N", "3"))


def gz_bytes(path):
    with gzip.open(path, "rb") as fh:
        return fh.read()


def run(sample_bed,
        vcf,
        mc,
        out,
        fast=True,
        fmt="auto",
        threads=2,
        split=True,
        monkeypatch=None):
    if not fast:
        monkeypatch.setattr(methbed, "read_filtered", lambda *a, **k: None)
    r = preprocess.preprocess_sample(sample_bed,
                                     vcf,
                                     mc,
                                     out + "_filtered.bed.gz",
                                     log=lambda m: None,
                                     fmt=fmt,
                                     threads=threads)
    if split:
        preprocess.split_sex_components(out + "_filtered.bed.gz",
                                        PAR,
                                        out,
                                        log=lambda m: None,
                                        threads=threads)
    return r


def same_outputs(a, b, split=True):
    for part in ("filtered", "autosomal", "chrX",
                 "chrY") if split else ("filtered",):
        pa, pb = f"{a}_{part}.bed.gz", f"{b}_{part}.bed.gz"
        assert gz_bytes(pa) == gz_bytes(pb), part
        assert os.path.exists(pa + ".tbi") == os.path.exists(pb + ".tbi"), part


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    d = tmp_path_factory.mktemp("synth")
    generate(str(d), n=N, seed=3)
    return str(d)


def clean_5col(synth, path, shuffle_contigs=True):
    """Writes a clean S01 5-column table and returns its path.

    It has no duplicate starts and contigs in natural (not C) order, as MethBat
    writes it.
    """
    rows = [
        line.split("\t")
        for line in gz_bytes(f"{synth}/S01.5col.bed.gz").decode().splitlines()
    ]
    seen, out = set(), []
    for r in rows:
        k = (r[0], r[1])
        if k not in seen:
            seen.add(k)
            out.append(r)
    order = {"chr21": 0, "chr22": 1, "chrX": 2, "chrY": 3}
    out.sort(key=lambda r: (order[r[0]], int(r[1])))
    with gzip.open(path, "wt") as fh:
        fh.write("".join("\t".join(r) + "\n" for r in out))
    return path


@pytest.mark.parametrize("fmt", ["5col", "modkit", "pbcpg"])
@pytest.mark.parametrize("mc", [5, 15])
def test_formats_with_ties_match_generic(synth, tmp_path, monkeypatch, fmt, mc):
    f = f"{synth}/S01.{fmt}.bed.gz"
    v = f"{synth}/S01.vcf"
    r1 = run(f, v, mc, str(tmp_path / "fast"))
    with monkeypatch.context() as m:
        r2 = run(f, v, mc, str(tmp_path / "gen"), fast=False, monkeypatch=m)
    assert r1 == r2
    same_outputs(str(tmp_path / "fast"), str(tmp_path / "gen"))


def test_clean_file_takes_fast_path_and_matches_generic(synth, tmp_path,
                                                        monkeypatch):
    bed = clean_5col(synth, str(tmp_path / "clean.bed.gz"))
    v = f"{synth}/S01.vcf"
    called = []
    orig = methbed.write_table
    monkeypatch.setattr(methbed, "write_table",
                        lambda *a, **k: called.append(1) or orig(*a, **k))
    r1 = run(bed, v, 10, str(tmp_path / "fast"), threads=3)
    assert called, "fast path not taken"
    with monkeypatch.context() as m:
        m.setattr(methbed, "read_filtered", lambda *a, **k: None)
        r2 = run(bed, v, 10, str(tmp_path / "gen"))
    assert r1 == r2
    same_outputs(str(tmp_path / "fast"), str(tmp_path / "gen"))
    assert r1["removed_by_het"] > 0 and r1["rows_out"] > 100


def test_tbi_equals_pysam(synth, tmp_path):
    import pysam
    bed = clean_5col(synth, str(tmp_path / "clean.bed.gz"))
    out = str(tmp_path / "o")
    run(bed, f"{synth}/S01.vcf", 10, out, threads=3)
    for part in ("filtered", "autosomal", "chrX", "chrY"):
        mine = f"{out}_{part}.bed.gz"
        ref = str(tmp_path / f"ref_{part}.bed.gz")
        with open(mine, "rb") as fi, open(ref, "wb") as fo:
            fo.write(fi.read())
        pysam.tabix_index(ref, preset="bed", force=True, keep_original=True)
        assert parse_tbi(mine + ".tbi") == parse_tbi(ref + ".tbi"), part
        with pysam.TabixFile(mine) as t:
            assert sum(1 for _ in t.fetch(t.contigs[0])) > 0


def test_bgzf_blocks_valid_and_threads_independent(synth, tmp_path):
    bed = clean_5col(synth, str(tmp_path / "clean.bed.gz"))
    outs = []
    for th in (1, 4):
        o = str(tmp_path / f"t{th}.bed.gz")
        preprocess.preprocess_sample(bed,
                                     f"{synth}/S01.vcf",
                                     5,
                                     o,
                                     log=lambda m: None,
                                     threads=th)
        outs.append(o)
    assert open(outs[0], "rb").read() == open(
        outs[1], "rb").read()  # identical bytes for any thread count
    raw = open(outs[0], "rb").read()
    assert raw.endswith(bz.EOF_BLOCK)
    import pysam
    with pysam.BGZFile(outs[0], "rb") as fh:  # htslib reads it
        assert fh.read() == gz_bytes(outs[0])


def test_header_spaces_and_odd_files_fall_back(synth, tmp_path, monkeypatch):
    # '#' header, space separated, CRLF, decimals, a mid-file '#': generic
    # result either way
    lines = [
        "#chrom\tstart\tend\tcov\tmod", "chr21 5010000 5010001 30 10",
        "chr21\t5010100\t5010101\t30\t5", "chr22\t5\t6\t40\t2"
    ]
    for name, text in {
            "hdr": "\n".join(lines) + "\n",
            "crlf": "\r\n".join(lines[1:]) + "\r\n",
            "dec": "\n".join(lines[1:3] + ["chr22\t5\t6\t40.5\t2"]) + "\n",
            "mid": "\n".join(lines[1:3] + ["#c", lines[3]]) + "\n"
    }.items():
        p = str(tmp_path / f"{name}.bed")
        open(p, "w", newline="").write(text)
        try:
            r1 = run(p, f"{synth}/S01.vcf", 0, str(tmp_path / f"f{name}"))
        except Exception as e:  # same failure as the generic code
            with monkeypatch.context() as m:
                with pytest.raises(type(e)):
                    run(p,
                        f"{synth}/S01.vcf",
                        0,
                        str(tmp_path / f"g{name}"),
                        fast=False,
                        monkeypatch=m)
            continue
        with monkeypatch.context() as m:
            r2 = run(p,
                     f"{synth}/S01.vcf",
                     0,
                     str(tmp_path / f"g{name}"),
                     fast=False,
                     monkeypatch=m)
        assert r1 == r2
        same_outputs(str(tmp_path / f"f{name}"), str(tmp_path / f"g{name}"))


def test_empty_inputs(synth, tmp_path):
    p = str(tmp_path / "e.bed.gz")
    with gzip.open(p, "wt") as fh:
        fh.write("")
    r = run(p, f"{synth}/S01.vcf", 5, str(tmp_path / "e"))
    assert r["rows_out"] == 0 and not os.path.exists(
        str(tmp_path / "e_filtered.bed.gz.tbi"))
    assert gz_bytes(str(tmp_path / "e_chrY.bed.gz")) == b""


# ---- pipeline step: one pass for filter + split, resume stamps (`build`)
# -----------------------------------


def _ctx(tmp_path, fmt="auto"):
    from dmrdetector import runtime
    cfg = {"input_format": fmt}
    ctx = runtime.Ctx("run", cfg, tmp_path / "o", tmp_path / "o" / "work", 2,
                      64.0)
    runtime.prepare_dirs(ctx, "preprocessed", "split")
    return runtime, ctx


def _stamp_status(runtime, ctx):
    return {k: v["status"] for k, v in ctx.timings.items()}


def test_pipeline_one_pass_equals_two_steps_and_resumes_from_stamps(
        synth, tmp_path):
    """One pass writes the files of the two steps; stamps decide what reruns.

    Covers the filter-and-split pass of `preprocess_samples`, the skip of a
    finished sample, the rerun of the split alone when only its stamp is
    missing, and the full pass again when both stamps are missing.
    """
    runtime, ctx = _ctx(tmp_path)
    clean = clean_5col(synth, str(tmp_path / "clean.bed.gz"))
    s = runtime.Sample("S01", "M", clean, f"{synth}/S01.vcf")
    preprocess.preprocess_samples(ctx, [s], 10, True, PAR)
    assert _stamp_status(runtime, ctx) == {
        "preprocess_S01": "done",
        "split_S01": "done"
    }
    ref = {
        p: gz_bytes(preprocess.scope_bed(ctx, "S01", sc, True))
        for p, sc in (("a", "autosomal"), ("x", "male_X"), ("y", "Y"))
    }
    ref["f"] = gz_bytes(preprocess.filtered_path(ctx, "S01"))
    assert all(
        os.path.exists(preprocess.scope_bed(ctx, "S01", sc, True) + ".tbi")
        for sc in ("autosomal", "male_X", "Y"))
    # the one pass equals the two separate steps
    out = str(tmp_path / "sep")
    preprocess.preprocess_sample(clean,
                                 f"{synth}/S01.vcf",
                                 10,
                                 out + "_filtered.bed.gz",
                                 log=lambda m: None)
    preprocess.split_sex_components(out + "_filtered.bed.gz",
                                    PAR,
                                    out,
                                    log=lambda m: None)
    assert gz_bytes(out + "_autosomal.bed.gz") == ref["a"] and gz_bytes(
        out + "_chrX.bed.gz") == ref["x"]
    assert gz_bytes(out + "_chrY.bed.gz") == ref["y"] and gz_bytes(
        out + "_filtered.bed.gz") == ref["f"]
    # second run: everything skipped
    ctx.timings.clear()
    preprocess.preprocess_samples(ctx, [s], 10, True, PAR)
    assert _stamp_status(runtime, ctx) == {
        "preprocess_S01": "skipped",
        "split_S01": "skipped"
    }
    # split stamp removed: only the standalone split runs, outputs unchanged
    os.remove(os.path.join(ctx.stamps, "split_S01.json"))
    ctx.timings.clear()
    preprocess.preprocess_samples(ctx, [s], 10, True, PAR)
    assert _stamp_status(runtime, ctx) == {
        "preprocess_S01": "skipped",
        "split_S01": "done"
    }
    assert gz_bytes(preprocess.scope_bed(ctx, "S01", "autosomal",
                                         True)) == ref["a"]
    # both stamps removed: one pass again
    for n in ("preprocess_S01", "split_S01"):
        os.remove(os.path.join(ctx.stamps, n + ".json"))
    ctx.timings.clear()
    preprocess.preprocess_samples(ctx, [s], 10, True, PAR)
    assert _stamp_status(runtime, ctx) == {
        "preprocess_S01": "done",
        "split_S01": "done"
    }
    assert gz_bytes(preprocess.scope_bed(ctx, "S01", "Y", True)) == ref["y"]


def test_pipeline_forced_format_errors(synth, tmp_path):
    runtime, ctx = _ctx(tmp_path, "modkit")
    s = runtime.Sample("S01", "M", f"{synth}/S01.5col.bed.gz",
                       f"{synth}/S01.vcf")
    with pytest.raises(runtime.RunError, match="needs at least 12"):
        preprocess.preprocess_samples(ctx, [s], 10, False, None)
