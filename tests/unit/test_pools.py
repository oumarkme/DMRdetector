"""Worker pools: a killed worker is an error, never a hang.

`multiprocessing.Pool` waits forever for the result of a worker that the
out-of-memory killer ended; a ProcessPoolExecutor raises BrokenProcessPool,
which `runtime.fork_pool` turns into a RunError. Each check runs in a
subprocess with a timeout, so a regression fails the test instead of hanging
the suite.
"""

import os
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = os.environ.get("DMRDETECTOR_SRC", str(ROOT / "src"))
TIMEOUT = 120  # seconds; a healthy script needs a few
MESSAGE = "a worker process was killed, probably out of memory; use fewer"


def _script(code: str, **env: str) -> subprocess.CompletedProcess:
    base = {
        k: v
        for k, v in os.environ.items()
        if k not in ("DMRDETECTOR_NO_JIT", "DMRDETECTOR_EXECUTOR")
    }
    return subprocess.run([sys.executable, "-c", code],
                          capture_output=True,
                          text=True,
                          timeout=TIMEOUT,
                          env=dict(base, PYTHONPATH=SRC, **env))


KILL = """
import os, signal
def die(*a, **k):
    os.kill(os.getpid(), signal.SIGKILL)
"""

ATLAS = """
from dmrdetector.fileio import bed
def atlas_file(path, ncol, n=600):
    # CHROM START END N COV PHI RAW_METH [REF_METH REF_VAR]
    rows = [("chr1", 1000 + 50 * i, 1001 + 50 * i, 5, 20, 0.1, 0.5, 0.5, 0.0)
            for i in range(n)]
    bed.write_bgzip_text(path, [bed.join_columns(
        [[str(r[i]) for r in rows] for i in range(ncol)])])
"""


def test_no_module_uses_multiprocessing_pool():
    found = []
    for path in (ROOT / "src").rglob("*.py"):
        text = path.read_text()
        if re.search(r"\.Pool\(|import Pool\b|multiprocessing\.pool", text):
            found.append(path.name)
    assert not found, found


def test_fork_pool_turns_a_killed_worker_into_a_run_error():
    r = _script(KILL + """
from dmrdetector import runtime
def task(i):
    if i == 2:
        die()
    return i
try:
    with runtime.fork_pool(2, "demo") as ex:
        list(ex.map(task, range(6)))
except runtime.RunError as e:
    print("RUNERROR", e)
""")
    assert r.returncode == 0, r.stderr
    assert f"RUNERROR demo: {MESSAGE}" in r.stdout


def test_fork_pool_keeps_the_order_and_the_other_errors():
    r = _script("""
from dmrdetector import runtime
def task(i):
    if i == 3:
        raise ValueError("task three")
    return i * i
with runtime.fork_pool(3, "demo") as ex:
    print(list(ex.map(task, range(3))))
try:
    with runtime.fork_pool(3, "demo") as ex:
        list(ex.map(task, range(6)))
except ValueError as e:
    print("VALUEERROR", e)
""")
    assert r.returncode == 0, r.stderr
    assert "[0, 1, 4]" in r.stdout and "VALUEERROR task three" in r.stdout


def test_atlas_reader_survives_a_killed_worker_as_an_error(tmp_path):
    r = _script(KILL + ATLAS + f"""
from dmrdetector import runtime
from dmrdetector.fileio import atlas_io, bgzf
atlas_file({str(tmp_path / 'a.bed.gz')!r}, 9)
bgzf.blocks_starting_in = die
try:
    atlas_io.read_columns({str(tmp_path / 'a.bed.gz')!r}, atlas_io.CPG,
                          max_code=23, min_n=0, strict=True, threads=2)
except runtime.RunError as e:
    print("RUNERROR", e)
""")
    assert r.returncode == 0, r.stderr
    assert f"RUNERROR atlas reader: {MESSAGE}" in r.stdout


def test_atlas_reader_forks_and_gives_the_rows_of_the_thread_route(tmp_path):
    code = ATLAS + f"""
import hashlib
from dmrdetector import runtime
from dmrdetector.fileio import atlas_io
atlas_file({str(tmp_path / 'a.bed.gz')!r}, 9)
forks = []
real = runtime.fork_pool
def spy(*a, **k):
    forks.append(a)
    return real(*a, **k)
runtime.fork_pool = spy
cols = atlas_io.read_columns({str(tmp_path / 'a.bed.gz')!r}, atlas_io.CPG,
                             max_code=23, min_n=0, strict=True, threads=2)
h = hashlib.md5(b"".join(c.tobytes() for c in cols[:5])).hexdigest()
print("FORKS", len(forks), h)
"""
    forked = _script(code).stdout.split()
    threaded = _script(code, DMRDETECTOR_EXECUTOR="thread").stdout.split()
    assert forked[:2] == ["FORKS", "1"] and threaded[:2] == ["FORKS", "0"]
    assert forked[2] == threaded[2]


def test_tabix_atlas_build_surfaces_a_killed_worker(tmp_path):
    r = _script(KILL + f"""
import numpy as np
import pysam
from dmrdetector import runtime
from dmrdetector.steps import atlas
files = []
for k in range(3):
    plain = {str(tmp_path)!r} + f"/s{{k}}.bed"
    with open(plain, "w") as fh:
        for i in range(200):
            p = 1000 + 10 * i
            fh.write(f"chr1\\t{{p}}\\t{{p + 1}}\\t20\\t5\\n")
    pysam.tabix_compress(plain, plain + ".gz", force=True)
    pysam.tabix_index(plain + ".gz", preset="bed", force=True)
    files.append(plain + ".gz")
atlas.aggregate_sites = die
try:
    atlas.build_atlas(files, {str(tmp_path / 'raw.bed.gz')!r}, 2, 0.5,
                      threads=2, scratch={str(tmp_path)!r}, window=1000,
                      log=lambda m: None)
except runtime.RunError as e:
    print("RUNERROR", e)
""",
                DMRDETECTOR_NO_JIT="1")
    assert r.returncode == 0, r.stderr
    assert f"RUNERROR atlas: {MESSAGE}" in r.stdout


def test_tabix_atlas_formatting_pool_surfaces_a_killed_worker(tmp_path):
    r = _script(KILL + f"""
import pysam
from dmrdetector import runtime
from dmrdetector.steps import atlas
files = []
for k in range(3):
    plain = {str(tmp_path)!r} + f"/s{{k}}.bed"
    with open(plain, "w") as fh:
        for i in range(200):
            p = 1000 + 10 * i
            fh.write(f"chr1\\t{{p}}\\t{{p + 1}}\\t20\\t5\\n")
    pysam.tabix_compress(plain, plain + ".gz", force=True)
    pysam.tabix_index(plain + ".gz", preset="bed", force=True)
    files.append(plain + ".gz")
atlas._format_raw_window = die  # the second pool, forked after the first
try:
    atlas.build_atlas(files, {str(tmp_path / 'raw.bed.gz')!r}, 2, 0.5,
                      threads=2, scratch={str(tmp_path)!r}, window=1000,
                      log=lambda m: None)
except runtime.RunError as e:
    print("RUNERROR", e)
""",
                DMRDETECTOR_NO_JIT="1")
    assert r.returncode == 0, r.stderr
    assert f"RUNERROR atlas: {MESSAGE}" in r.stdout


def test_pandas_smoothing_surfaces_a_killed_worker(tmp_path):
    r = _script(KILL + ATLAS + f"""
from dmrdetector import runtime
from dmrdetector.steps import atlas
atlas_file({str(tmp_path / 'raw.bed.gz')!r}, 7)
atlas.smooth_chromosome = die
try:
    atlas.smooth_atlas({str(tmp_path / 'raw.bed.gz')!r}, 40,
                       {str(tmp_path / 'ref.bed.gz')!r}, True, 2,
                       log=lambda m: None)
except runtime.RunError as e:
    print("RUNERROR", e)
""",
                DMRDETECTOR_NO_JIT="1")
    assert r.returncode == 0, r.stderr
    assert f"RUNERROR smoothing: {MESSAGE}" in r.stdout
