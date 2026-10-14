"""Real `build` runs: the log shows names, and a resume repeats the warnings.

The data are the controls of the smoke test, copied to folders with names that
appear nowhere else. One control has no CpG record, which gives a warning that
carries the path of its BED. Each check reads the files that the build wrote.
"""

import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pysam
import pytest

from dmrdetector import build
from dmrdetector import data


ROOT = Path(__file__).resolve().parents[2]
SRC = os.environ.get("DMRDETECTOR_SRC", str(ROOT / "src"))
EMPTY_WARNING = "WARNING: CTRL01: no CpG records in CTRL01.bed.gz"
NAME_DIR = "cohort_zq7"  # appears only in paths of the samples


def _smoke():
    spec = importlib.util.spec_from_file_location(
        "run_smoke", ROOT / "tests" / "smoke" / "run_smoke.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    """Returns (samplesheet, sample folder, assets folder, parent folder)."""
    sm = _smoke()
    base = tmp_path_factory.mktemp("logs")
    sm.write_data(base)
    folder = base / NAME_DIR / "sub_k3"
    folder.mkdir(parents=True)
    for f in (base / "data" / "5col").iterdir():
        shutil.copy(f, folder / f.name)
    (folder / "none.txt").write_text("")
    pysam.tabix_compress(str(folder / "none.txt"),
                         str(folder / "CTRL01.bed.gz"),
                         force=True)  # a control without a CpG record
    sheet = sm.sheet(base / "s.csv", base / "data")
    sheet.write_text(sheet.read_text().replace(f"{base}/data/5col",
                                               str(folder)))
    assets = base / "assets_zq7"
    assets.mkdir()
    shutil.copy(data.asset("par_hg38.bed"), assets / "par_hg38.bed")
    shutil.copy(data.asset("icr_gene_list.bed"), assets / "icr_genes.bed")
    return sheet, folder, assets, base


def _build(sheet, out, *extra):
    cmd = [
        sys.executable, "-m", "dmrdetector.cli", "build", "--input",
        str(sheet), "--outdir",
        str(out), "--threads", "6", "--keep-tmp", "--set",
        "cpg_min_samples_per_site=2", "--set", "cpg_min_island_length=100",
        "--set", "icr_min_sites=5", *extra
    ]
    return subprocess.run(cmd,
                          capture_output=True,
                          text=True,
                          env=dict(os.environ, PYTHONPATH=SRC))


def _section(log, name):
    """Returns the lines of one [section] of a log."""
    lines = log.splitlines()
    start = next(
        i for i, ln in enumerate(lines) if ln.startswith(f"[{name}]")) + 1
    end = next((i for i in range(start, len(lines))
                if re.match(r"\[[a-z]+\]", lines[i])), len(lines))
    return lines[start:end]


@pytest.fixture(scope="module")
def first_build(cohort, tmp_path_factory):
    sheet, _, assets, _ = cohort
    out = tmp_path_factory.mktemp("ref")
    r = _build(sheet, out, "--set", f"par_bed={assets / 'par_hg38.bed'}",
               "--set", f"icr_gene_list={assets / 'icr_genes.bed'}")
    assert r.returncode == 0, r.stderr[-2000:]
    return out, r


def test_build_log_and_reference_json_hold_no_folder_of_the_samples(
        cohort, first_build):
    _, folder, assets, base = cohort
    out, r = first_build
    log = (out / "build.log").read_text()
    ref = (out / "reference.json").read_text()
    # the warning with the path of a sample is there, by name
    assert any(EMPTY_WARNING in ln for ln in _section(log, "warnings"))
    # the progress log of the work folder keeps the path (it is not shared)
    progress = (out / "work" / "progress.log").read_text()
    assert str(folder / "CTRL01.bed.gz") in progress
    for text in (log, ref):
        assert str(base.parent) not in text
        assert str(folder) not in text and str(assets) not in text
    params = json.loads(ref)["parameters"]
    assert params["par_bed"] == "par_hg38.bed"
    assert params["icr_gene_list"] == "icr_genes.bed"


def test_resumed_build_repeats_the_warnings_of_the_skipped_steps(
        cohort, first_build, tmp_path):
    sheet, _, assets, _ = cohort
    out, _ = first_build
    first = (out / "build.log").read_text()
    warned = [ln for ln in _section(first, "warnings") if ln]
    assert any(EMPTY_WARNING in ln for ln in warned)
    assert any("CTRL05: sex F, but the VCF has 0 heterozygous SNPs" in ln
               for ln in warned)
    # the same command again: the work folder was kept, so every step is
    # skipped, and the log of this run must still show the warnings
    r = _build(sheet, out, "--set", f"par_bed={assets / 'par_hg38.bed'}",
               "--set", f"icr_gene_list={assets / 'icr_genes.bed'}")
    assert r.returncode == 0, r.stderr[-2000:]
    second = (out / "build.log").read_text()
    steps = _section(second, "steps")
    assert "preprocess_CTRL01: skipped, 0.0" in steps
    assert not any(": done" in ln for ln in steps)
    # the second build also warns that it overwrites the first reference
    again = [
        ln for ln in _section(second, "warnings")
        if ln and "already holds a reference" not in ln
    ]
    assert _untimed(again) == _untimed(warned)


def _untimed(lines):
    """Returns the sorted lines without their time stamp."""
    return sorted(ln.split("] ", 1)[1] for ln in lines)


def _error_build(cohort, out, *extra):
    """Runs a build in which one BED has a record beyond the end of chr22."""
    sheet, folder, _, base = cohort
    bad = base / NAME_DIR / "bad"
    bad.mkdir(parents=True, exist_ok=True)
    for f in folder.iterdir():
        shutil.copy(f, bad / f.name)
    with gzip.open(folder / "CTRL02.bed.gz", "rt") as fh:
        text = fh.read() + "chr22\t60000000\t60000001\t20\t5\n"
    (bad / "CTRL02.txt").write_text(text)
    pysam.tabix_compress(str(bad / "CTRL02.txt"),
                         str(bad / "CTRL02.bed.gz"),
                         force=True)
    bad_sheet = base / "bad.csv"
    bad_sheet.write_text(sheet.read_text().replace(str(folder), str(bad)))
    return bad, _build(bad_sheet, out, *extra)


def test_build_error_block_names_the_file_not_its_folder(cohort, tmp_path):
    bad, r = _error_build(cohort, tmp_path / "ref")
    assert r.returncode == 1 and "CTRL02.bed.gz: chr22 has a record" in r.stderr
    log = (tmp_path / "ref" / "build.log").read_text()
    err = [ln for ln in log.splitlines() if ln.startswith("error:")]
    assert len(err) == 1
    assert "CTRL02.bed.gz: chr22 has a record ending at 60,000,001" in err[0]
    assert str(bad) not in log and str(bad.parent) not in log
    # the same build with full paths shows the folder: the check is not empty
    _, r = _error_build(cohort, tmp_path / "ref_full", "--set",
                        "log_input_paths=full")
    assert r.returncode == 1
    full = (tmp_path / "ref_full" / "build.log").read_text()
    assert f"{bad}/CTRL02.bed.gz: chr22 has a record" in full


def test_reference_json_records_path_options_by_name_unless_full():
    cfg = dict(par_bed="/a/b/par.bed",
               icr_gene_list=None,
               smooth_sigma=150,
               threads=4,
               tmp_dir="/t",
               keep_tmp=False,
               log_input_paths="name")
    got = build._reference_parameters(cfg)
    assert got == dict(par_bed="par.bed", icr_gene_list=None, smooth_sigma=150)
    full = build._reference_parameters(dict(cfg, log_input_paths="full"))
    assert full["par_bed"] == "/a/b/par.bed" and "log_input_paths" not in full
