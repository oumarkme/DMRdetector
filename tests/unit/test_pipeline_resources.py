"""Tests of the thread budget, the temporary folder and the atlas index.

They cover `build` and `run`, the gene table hand-over of `run` and the tabix
index of the atlas that `build` writes.
"""
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from dmrdetector import data
from dmrdetector import run
from dmrdetector import runtime
from dmrdetector.steps import annotate


ROOT = Path(__file__).resolve().parents[2]


def _ctx(tmp_path, threads, mem=None):
    cfg = dict(threads=threads, input_format="auto", tmp_dir=None)
    return runtime.Ctx("build", cfg, tmp_path / "out",
                       tmp_path / "out" / "work", threads,
                       8.0 * threads if mem is None else mem)


def test_thread_plan_never_exceeds_the_thread_or_memory_budget(tmp_path):
    rng = np.random.default_rng(0)
    for _ in range(300):
        threads = int(rng.integers(2, 51))
        n = int(rng.integers(1, 15))
        gbs = [float(x) for x in rng.choice([0.5, 2, 5, 14], n)]
        w = [float(x) for x in rng.choice([1, 1, 6, 10], n)]
        ctx = _ctx(tmp_path, threads)
        workers, th = runtime.thread_plan(ctx, gbs, w)
        assert 1 <= workers <= min(max(1, threads // 2), n) and all(
            t >= 2 for t in th)
        # the `workers` tasks that run together (worst case: the ones with most
        # threads) fit in the thread budget
        assert sum(sorted(th, reverse=True)[:workers]) <= max(threads, workers)
        assert max(th) <= threads
        # and the largest estimates of that many tasks fit in 8 GB x threads
        assert sum(sorted(gbs,
                          reverse=True)[:workers]) <= ctx.mem_gb or workers == 1


def test_five_panel_tasks_are_not_five_times_eight_threads(tmp_path):
    workers, th = runtime.thread_plan(_ctx(tmp_path, 30), [2.0] * 5, [1.0] * 5)
    assert workers == 5 and th == [6] * 5  # 30 threads, not 40
    workers, th = runtime.thread_plan(_ctx(tmp_path, 12), [2.0] * 12,
                                      [1.0] * 12)
    assert workers == 6 and th == [
        2
    ] * 12  # a worker counts its main thread: never 12 workers of 1 + main
    workers, th = runtime.thread_plan(_ctx(tmp_path, 30), [20.0, 8.0, 8.0, 3.0],
                                      [10.0, 1.0, 1.0, 1.0])
    assert workers == 4 and th[0] > 4 * th[1] and sum(
        th) <= 30  # autosomal gets the bulk


def test_memory_limits_the_workers_below_the_thread_count(tmp_path):
    workers, th = runtime.thread_plan(_ctx(tmp_path, 8), [20.0] * 6,
                                      [1.0] * 6)  # 64 GB: three 20 GB tasks
    assert workers == 3 and th == [2 + int(2 / 3)] * 6
    with pytest.raises(runtime.RunError, match="exceeds the budget"):
        runtime.thread_plan(_ctx(tmp_path, 2), [20.0], [1.0])


def test_scope_plan_small_budgets_give_the_autosomal_scope_every_thread(
        tmp_path):
    """Checks the scope plan below and above 16 threads.

    Below 16 threads the scopes run in turn with the whole budget (a sex scope
    beside the autosomal one would fix its thread count at a fraction of the
    threads for the whole smoothing); from 16 threads they run together.
    """
    from dmrdetector import build
    scopes = ["autosomal", "male_X", "female_X", "Y"]
    gbs = [20.0, 9.0, 9.0, 3.5]
    for threads in (2, 4, 8, 15):
        w, mw = build.scope_plan(threads, scopes)
        workers, th = runtime.thread_plan(_ctx(tmp_path, threads, 1000.0), gbs,
                                          w, mw)
        assert workers == 1 and th == [threads] * 4, (threads, th)
    for threads in (16, 24, 30, 50):
        w, mw = build.scope_plan(threads, scopes)
        workers, th = runtime.thread_plan(_ctx(tmp_path, threads), gbs, w, mw)
        assert workers == 4 and th[1:] == [
            2, 2, 2
        ] and th[0] >= threads - 9 and sum(th) <= threads, (threads, th)


def test_scope_memory_shares_fit_in_the_budget():
    from dmrdetector import build
    needs = [20.0, 9.0, 9.0, 3.5]
    for workers in (1, 2, 3, 4):
        for mem in (40.0, 96.0, 288.0):
            share = build.scope_memory(needs, workers, mem)
            assert share[1:] == needs[1:] and share[0] >= needs[0]
            assert share[0] == max(
                20.0, mem - sum(sorted(needs[1:], reverse=True)[:workers - 1]))


def test_memgate_admits_in_order_while_the_estimates_fit():
    """Checks that MemGate starts FFT tasks in order while their sum fits."""
    import threading
    import time

    from dmrdetector.steps import atlas
    gate = atlas.MemGate(10.0)
    live, peak, lock = [0.0], [0.0], threading.Lock()

    def job(i, gb):
        gate.acquire(i, gb)
        with lock:
            live[0] += gb
            peak[0] = max(peak[0], live[0])
        time.sleep(0.02)
        with lock:
            live[0] -= gb
        gate.release(gb)

    gbs = [6.0, 4.0, 4.0, 3.0, 3.0, 1.0]
    ts = [threading.Thread(target=job, args=(i, g)) for i, g in enumerate(gbs)]
    [t.start() for t in ts]
    [t.join() for t in ts]
    assert 0 < peak[0] <= 10.0


def _probe(t):
    return [(
        f"{t['i']}|{tempfile.gettempdir()}|{os.environ['TMPDIR']}|{os.environ['DMRDETECTOR_THREADS']}|{t['threads']}",
        "done", 0.0)]


def test_tmp_dir_and_threads_reach_every_worker(tmp_path, monkeypatch):
    monkeypatch.delenv("DMRDETECTOR_THREADS", raising=False)
    old = (os.environ.get("TMPDIR"), tempfile.tempdir)
    try:
        ctx = _ctx(tmp_path, 8)
        runtime.prepare_dirs(ctx)
        assert tempfile.gettempdir() == str(ctx.scratch) == os.environ["TMPDIR"]
        for threads, n in ((8, 4), (2,
                                    1)):  # a fork pool and the in-process path
            ctx = _ctx(tmp_path, threads)
            ctx.timings.clear()
            runtime.pool_map(ctx, _probe, [dict(i=i) for i in range(n)], 4.0,
                             "probe")
            names = list(ctx.timings)
            assert len(names) == n
            for nm in names:
                _, tmpd, env, thr, t = nm.split("|")
                assert tmpd == env == str(
                    ctx.scratch) and thr == t and int(t) == max(
                        1, threads // min(threads, n) - 1)
        assert "DMRDETECTOR_THREADS" not in os.environ  # restored in the parent
    finally:
        tempfile.tempdir = old[1]
        if old[0] is None:
            os.environ.pop("TMPDIR", None)
        else:
            os.environ["TMPDIR"] = old[0]


def test_annotate_gets_the_gene_table_as_an_argument(tmp_path, monkeypatch):
    seen = {}

    def fake_annotate(raw,
                      gtf,
                      icr_panel=None,
                      icr_gene_list=None,
                      *,
                      genes=None):
        seen["genes"] = genes
        import pandas as pd
        return pd.DataFrame(columns=annotate.OUTPUT_COLS)

    monkeypatch.setattr(annotate, "annotate", fake_annotate)
    monkeypatch.setattr(annotate, "write_annotated",
                        lambda df, path: Path(path).write_text("x"))
    orig = annotate.load_genes
    monkeypatch.setattr(run, "_GENES", "TABLE")
    for n in ("a", "g", "l"):
        (tmp_path / n).write_text("x")
    t = dict(id="S1",
             raw=[str(tmp_path / "a")],
             gtf=str(tmp_path / "g"),
             icr_panel=None,
             icr_gene_list=str(tmp_path / "l"),
             stamps=str(tmp_path),
             out=str(tmp_path / "o.tsv"))
    run._t_annotate(t)
    assert seen[
        # nothing is patched
        "genes"] == "TABLE" and annotate.load_genes is orig


# ---- the atlas written by `build` has its index
# -------------------------------------------------------------


def _smoke():
    spec = importlib.util.spec_from_file_location(
        "run_smoke", ROOT / "tests" / "smoke" / "run_smoke.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    sm = _smoke()
    base = tmp_path_factory.mktemp("build")
    sm.write_data(base)
    sheet = sm.sheet(base / "full.csv", base / "data")
    out = base / "ref"
    r = subprocess.run([
        sys.executable, "-m", "dmrdetector.cli", "build", "--input",
        str(sheet), "--outdir",
        str(out), "--threads", "6", "--keep-tmp", *sm.BUILD_OPTS
    ],
                       capture_output=True,
                       text=True,
                       env=dict(os.environ,
                                PYTHONPATH=os.environ.get(
                                    "DMRDETECTOR_SRC", str(ROOT / "src"))))
    assert r.returncode == 0, r.stderr[-2000:]
    return sm, out, base


def test_build_writes_the_atlas_index_with_the_atlas(built):
    sm, out, _ = built
    atlases = sorted(out.glob("ref_atlas_150*.bed.gz"))
    assert len(atlases) == 4
    missing = [a.name for a in atlases if not Path(str(a) + ".tbi").is_file()]
    # fails when the smooth step leaves the index behind
    assert not missing, f"atlas without .tbi: {missing}"
    assert sm.atlas_index_problems(out) == []
    # the index is not older than its data, so htslib trusts it
    for a in atlases:
        assert os.stat(str(a) + ".tbi").st_mtime_ns >= os.stat(a).st_mtime_ns


def test_build_leaves_no_temporary_files_and_no_stray_index(built):
    _, out, base = built
    stray = [p.name for p in out.rglob("*") if p.name.startswith(".tmp")]
    assert not stray, stray
    # nothing re-indexed later: the .tbi is exactly the one the smooth step
    # wrote (an index of the final file)
    a = out / "ref_atlas_150.bed.gz"
    import pysam
    with pysam.TabixFile(str(a)) as tbx:
        assert tbx.contigs == ["chr22"]
        assert sum(1 for _ in tbx.fetch("chr22")) == len(
            gzip.open(a).read().splitlines())


def test_build_records_the_md5_of_the_icr_gene_list(built):
    _, out, _ = built
    info = json.loads((out / "reference.json").read_text())
    assert info["icr_gene_list_md5"] == runtime.md5_file(
        data.asset("icr_gene_list.bed"))


def test_a_smooth_step_without_index_fails_the_build(tmp_path, monkeypatch):
    """Checks that a smooth step without an index is an error, not a gap.

    The step declares `<atlas>.tbi` as an output.
    """
    from dmrdetector import build
    from dmrdetector.steps import atlas
    monkeypatch.setattr(atlas, "build_atlas",
                        lambda beds, out, *a, **k: Path(out).write_bytes(b"x"))

    def smooth_without_index(raw, sigma, out, *a, **k):
        Path(out).write_bytes(b"x")  # no out + ".tbi"

    monkeypatch.setattr(atlas, "smooth_atlas", smooth_without_index)
    t = dict(scope="autosomal",
             beds=[],
             raw=str(tmp_path / "raw.bed.gz"),
             final=str(tmp_path / "final.bed.gz"),
             scratch=str(tmp_path),
             sigma=150,
             include_sex=True,
             stamps=str(tmp_path),
             threads=1,
             min_samples_atlas=3,
             shrinkage=0.1,
             panels=[])
    with pytest.raises(runtime.RunError, match="did not produce"):
        build._t_scope(t)


def test_cli_and_annotate_import_without_numba():
    """Checks that annotate and the option handling load without numba.

    numba is imported only when a step that needs it runs.
    """
    code = ("import sys; sys.modules['numba'] = None; "
            "import dmrdetector.fileio.bed, dmrdetector.steps.annotate, "
            "dmrdetector.cli, dmrdetector.build, "
            "dmrdetector.run, dmrdetector.runtime, dmrdetector.settings, "
            "dmrdetector.steps.call_dmr")
    r = subprocess.run([sys.executable, "-c", code],
                       capture_output=True,
                       text=True,
                       env=dict(os.environ,
                                PYTHONPATH=os.environ.get(
                                    "DMRDETECTOR_SRC", str(ROOT / "src"))))
    assert r.returncode == 0, r.stderr[-1000:]


def test_tmp_dir_option_keeps_every_temporary_file_out_of_tmpdir(
        built, tmp_path):
    """Checks that every temporary file goes below --tmp-dir, not TMPDIR.

    `build` and `run` (with the target report) run with TMPDIR at an empty
    sentinel directory. The sentinel stays empty (a file made before set_tmp
    would land in it), the work files are below --tmp-dir, not in --outdir, and
    the run leaves --tmp-dir empty.
    """
    sm, ref, base = built
    sentinel, tmpd = tmp_path / "sentinel", tmp_path / "scratch"
    sentinel.mkdir()
    env = dict(os.environ,
               PYTHONPATH=os.environ.get("DMRDETECTOR_SRC", str(ROOT / "src")),
               TMPDIR=str(sentinel))
    res = ["--threads", "6", "--tmp-dir", str(tmpd)]
    out_b = tmp_path / "ref"
    sheet = sm.sheet(tmp_path / "ctrl.csv", base / "data")
    r = subprocess.run([
        sys.executable, "-m", "dmrdetector.cli", "build", "--input",
        str(sheet), "--outdir",
        str(out_b), *res, "--keep-tmp", *sm.BUILD_OPTS
    ],
                       capture_output=True,
                       text=True,
                       env=env)
    assert r.returncode == 0, r.stderr[-1500:]
    tb = tmp_path / "target.bed"
    tb.write_text(
        f"chr22\t{sm.ICR_ISLAND}\t{sm.ICR_ISLAND + sm.ISLAND_LEN}\tT1\t"
        "DGCR6L\tTest syndrome (hypo)\n")
    gtf = sm.small_gtf(tmp_path / "genes.gtf.gz")
    out_r = tmp_path / "run"
    r = subprocess.run([
        sys.executable, "-m", "dmrdetector.cli", "run",
        *sm.sample_args(base / "data", "CASE01", "F"), "--outdir",
        str(out_r), "--reference",
        str(out_b), "--gtf",
        str(gtf), *res, *sm.RUN_OPTS, "--target-region",
        str(tb)
    ],
                       capture_output=True,
                       text=True,
                       env=env)
    assert r.returncode == 0, r.stderr[-1500:]
    assert (out_r / "CASE01.target_region.tsv").is_file()
    assert list(sentinel.iterdir()) == [], list(sentinel.iterdir())
    assert not (out_b / "work").exists() and not (out_r / "work").exists()
    kept = runtime.build_workdir(str(tmpd), out_b.resolve())
    assert (kept / "stamps").is_dir() and (
        out_b / "ref_atlas_150.bed.gz.tbi").is_file()
    assert [p.name for p in tmpd.iterdir()
           ] == [kept.name]  # the run removed its own folder only


def test_resume_with_a_missing_index_reruns_only_that_scope(built, tmp_path):
    """Checks that a missing index reruns only that scope.

    The smooth step declares the index as an output: after one .tbi is deleted,
    a rerun redoes that scope's smoothing and CpG panel (the panel reads the
    atlas) and nothing else (the work folder is kept with --keep-tmp).
    """
    sm, _, base = built
    out = tmp_path / "ref"
    sheet = sm.sheet(tmp_path / "ctrl.csv", base / "data")
    env = dict(os.environ,
               PYTHONPATH=os.environ.get("DMRDETECTOR_SRC", str(ROOT / "src")))
    cmd = [
        sys.executable, "-m", "dmrdetector.cli", "build", "--input",
        str(sheet), "--outdir",
        str(out), "--threads", "6", "--keep-tmp", *sm.BUILD_OPTS
    ]
    assert subprocess.run(cmd, capture_output=True, text=True,
                          env=env).returncode == 0
    (out / "ref_atlas_150_Y.bed.gz.tbi").unlink()
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr[-1500:]
    done = sorted(
        ln.split("] ", 1)[1].split(":")[0]
        for ln in r.stderr.splitlines()
        if ": done" in ln)
    assert done == ["smooth_atlas_Y"
                   ] or done == ["cpg_islands_Y", "smooth_atlas_Y"], done
    assert (out / "ref_atlas_150_Y.bed.gz.tbi").is_file()


def test_smooth_atlas_refuses_a_budget_below_its_estimate(built, tmp_path):
    from dmrdetector.steps import atlas
    sm, out, _ = built
    raw = out / "work" / "raw_atlas" / "raw_atlas.bed.gz"
    with pytest.raises(RuntimeError, match="smoothing needs about"):
        atlas.smooth_atlas(str(raw),
                           150,
                           str(tmp_path / "s.bed.gz"),
                           True,
                           2,
                           log=lambda m: None,
                           mem_gb=0.5)
    a, b = tmp_path / "a.bed.gz", tmp_path / "b.bed.gz"
    atlas.smooth_atlas(str(raw),
                       150,
                       str(a),
                       True,
                       2,
                       log=lambda m: None,
                       mem_gb=16.0)
    atlas.smooth_atlas(str(raw), 150, str(b), True, 4, log=lambda m: None)
    assert a.read_bytes() == b.read_bytes()


def test_scope_need_fits_the_minimum_budget_for_any_cohort_size(tmp_path):
    """Checks that a scope worker needs no more memory for a larger cohort.

    The raw atlas reads windows and the smoothing holds arrays that follow the
    number of sites (12 samples of 0.21 GB each; 40 and 130 samples
    are sparse files of the same size). All fit the 16 GB of `--threads 2`.
    """
    from dmrdetector import build
    from dmrdetector import runtime
    for n in (12, 40, 130):
        beds = []
        for i in range(n):
            p = tmp_path / f"s{n}_{i}_autosomal.bed.gz"
            with open(p, "wb") as fh:
                fh.truncate(int(2.51e9 / 12))
            beds.append(str(p))
        need = build.scope_need_gb(beds, "autosomal")
        assert need <= 16.0, (n, need)
        assert need > runtime.smooth_gb(
            runtime.SCOPE_LARGEST_SPAN["autosomal"]
        )  # the smoothing of chr1 is the largest term
        assert runtime.raw_atlas_min_gb(sum(os.path.getsize(b) for b in beds),
                                        runtime.SCOPE_EXTENT["autosomal"],
                                        3) < 6.5
    ctx = _ctx(tmp_path, 2)
    assert runtime.thread_plan(ctx,
                               [build.scope_need_gb(beds, "autosomal"), 2.0],
                               [1.0, 1.0])[0] == 1
