"""Checks that the results of `build` and `run` do not depend on --threads.

The synthetic smoke data gives byte-identical atlases (with their .tbi), panels
and annotated DMR tables, identical text of the raw DMR tables (the gzip header
holds the time of writing) and of the target report at 2, 4 and 36 threads.
"""
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from dmrdetector import runtime


ROOT = Path(__file__).resolve().parents[2]
THREADS = (2, 4, 36)


def _smoke():
    spec = importlib.util.spec_from_file_location(
        "run_smoke_inv", ROOT / "tests" / "smoke" / "run_smoke.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _cli(*args, env):
    return subprocess.run(
        [sys.executable, "-m", "dmrdetector.cli", *map(str, args)],
        capture_output=True,
        text=True,
        env=env)


@pytest.fixture(scope="module")
def outputs(tmp_path_factory):
    sm = _smoke()
    base = tmp_path_factory.mktemp("inv")
    sm.write_data(base)
    ctrl = sm.sheet(base / "ctrl.csv", base / "data")
    gtf = sm.small_gtf(base / "genes.gtf.gz")
    tb = base / "target.bed"
    tb.write_text(
        f"chr22\t{sm.ICR_ISLAND}\t{sm.ICR_ISLAND + sm.ISLAND_LEN}\tT1\t"
        "DGCR6L\tTest syndrome (hypo)\n")
    env = dict(os.environ,
               PYTHONPATH=os.environ.get("DMRDETECTOR_SRC", str(ROOT / "src")))
    out = {}
    for th in THREADS:
        ref, res = base / f"ref{th}", base / f"run{th}"
        r = _cli("build",
                 "--input",
                 ctrl,
                 "--outdir",
                 ref,
                 "--threads",
                 th,
                 "--keep-tmp",
                 *sm.BUILD_OPTS,
                 env=env)
        assert r.returncode == 0, r.stderr[-1500:]
        n, mem, _ = runtime.resolve_threads(th)  # lowered on a small host
        assert f"threads: {n}  memory budget: {mem:.0f} GB" in r.stderr
        for sid, sex in sm.CASE:
            r = _cli("run",
                     *sm.sample_args(base / "data", sid, sex),
                     "--outdir",
                     res,
                     "--reference",
                     ref,
                     "--gtf",
                     gtf,
                     "--threads",
                     th,
                     *sm.RUN_OPTS,
                     "--target-region",
                     tb,
                     env=env)
            assert r.returncode == 0, r.stderr[-1500:]
        out[th] = (ref, res)
    return out


def _files(ref: Path, res: Path) -> dict:
    found = {}
    pats = [(ref, [
        "ref_atlas_150*.bed.gz", "ref_atlas_150*.bed.gz.tbi", "cpg_signal*.bed",
        "icr_region.bed", "work/raw_atlas/*.bed.gz"
    ]), (res, ["*.dmr.tsv", "*.target_region.tsv"])]
    for d, globs in pats:
        for g in globs:
            for p in sorted(d.glob(g)):
                data = p.read_bytes()
                if g.startswith(
                        "work/raw_atlas"
                        # plain gzip: its header holds the time of writing, so
                        # compare the text
                ):
                    data = gzip.decompress(data)
                found[f"{d.name[:3]}/{p.relative_to(d)}"] = data
    return found


def test_outputs_are_byte_identical_at_2_4_and_36_threads(outputs):
    base = _files(*outputs[2])
    assert len(base) >= 20 and any(k.endswith(".tbi") for k in base) and any(
        "CASE01.dmr.tsv" in k for k in base)
    for th in (4, 36):
        got = _files(*outputs[th])
        assert set(got) == set(base)
        assert [k for k in base if got[k] != base[k]] == [], th


def test_reference_json_is_identical_apart_from_the_creation_time(outputs):

    def load(th):
        d = json.loads((outputs[th][0] / "reference.json").read_text())
        d.pop("created")
        return d

    assert "threads" not in load(2)["parameters"]
    assert load(4) == load(2) and load(36) == load(2)


def test_log_records_the_effective_threads(outputs):
    for th in THREADS:
        for sid in ("CASE01", "CASE02"):
            log = (outputs[th][1] / f"{sid}.log").read_text()
            n, mem, _ = runtime.resolve_threads(th)
            assert f'threads: {n}, memory budget {mem:.0f} GB' in log
            assert f'threads = {n}  # flag' in log
        assert f"threads = {n}  # flag" in (outputs[th][0] /
                                            "build.log").read_text()
