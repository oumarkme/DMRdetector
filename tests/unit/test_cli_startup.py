"""Checks that command registration stays cheap.

Building the main parser imports no numba, scipy or pandas. The stage commands
(call-dmr, atlas, ...) load their modules only when they are used. The checks
run in a subprocess so that modules imported by other tests do not count.
"""

import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import tempfile
import threading

import pytest


ROOT = Path(__file__).resolve().parents[2]
SRC = os.environ.get("DMRDETECTOR_SRC", str(ROOT / "src"))
HEAVY = ("numba", "llvmlite", "scipy", "pandas")


def _py(code: str, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-c", code, *args],
                          capture_output=True,
                          text=True,
                          env=dict(os.environ, PYTHONPATH=SRC))


def test_build_parser_imports_no_heavy_module():
    r = _py("import sys, dmrdetector.cli as c; c.build_parser(); "
            f"print(sorted(m for m in {HEAVY!r} if m in sys.modules))")
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "[]"


@pytest.mark.parametrize("argv",
                         [["--version"], ["--help"], ["run", "--help"],
                          ["build", "--help"], ["config", "show", "run"]])
def test_light_commands_import_no_heavy_module(argv):
    r = _py(
        "import sys, dmrdetector.cli as c\n"
        "try:\n    c.main(sys.argv[1:])\nexcept SystemExit:\n    pass\n"
        f"print('HEAVY', sorted(m for m in {HEAVY!r} if m in sys.modules))",
        *argv)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().splitlines()[-1] == "HEAVY []"


def test_hidden_stage_commands_stay_hidden_but_parse():
    r = _py(
        "import dmrdetector.cli as c; p = c.build_parser(); h = "
        "p.format_help(); "
        "print('call-dmr' in h, 'get-genes' in h); "
        "a = p.parse_args(['call-dmr', '--sample-bed', 'a', '--atlas', 'b', "
        "'--cpg-panel', 'c', '--out', 'd']); "
        "print(a.command, a.scope, a.func.__name__)")
    assert r.returncode == 0, r.stderr
    assert r.stdout.split("\n")[:2] == ["False True", "call-dmr autosomal main"]


def test_stage_help_loads_the_real_parser():
    r = subprocess.run(
        [sys.executable, "-m", "dmrdetector.cli", "atlas", "smooth", "--help"],
        capture_output=True,
        text=True,
        env=dict(os.environ, PYTHONPATH=SRC))
    assert r.returncode == 0, r.stderr
    assert 'usage: dmrdetector atlas smooth' in r.stdout
    assert '--raw-atlas' in r.stdout


def _user_folder(parent):
    """Returns the per-user cache folder for the seed of the environment."""
    from dmrdetector.util import jit
    return Path(parent) / jit.cache_folder_name(os.geteuid())


def _cache_env(monkeypatch, tmp_path):
    from dmrdetector.util import jit
    monkeypatch.setenv("NUMBA_CACHE_DIR", "x")
    monkeypatch.delenv("NUMBA_CACHE_DIR")  # restored after the test
    monkeypatch.delenv("DMRDETECTOR_NO_JIT", raising=False)
    monkeypatch.setenv(jit.SEED_ENV, "")  # not the cache of the image
    jit._seed_once.cache_clear()
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    return jit


def test_numba_cache_defaults_to_a_private_folder_of_the_user(
        tmp_path, monkeypatch):
    jit = _cache_env(monkeypatch, tmp_path)
    want = _user_folder(tmp_path)
    assert want.name == f"dmrdetector-numba-{os.geteuid()}-noseed"
    assert jit.set_default_cache_dir() == str(want)
    assert os.environ["NUMBA_CACHE_DIR"] == str(want)
    assert stat.S_IMODE(want.stat().st_mode) == 0o700
    # a second call finds the variable set and changes nothing
    assert jit.set_default_cache_dir() is None


def test_numba_cache_variable_set_by_the_caller_is_kept(tmp_path, monkeypatch):
    jit = _cache_env(monkeypatch, tmp_path)
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "mine"))
    assert jit.set_default_cache_dir() is None
    assert os.environ["NUMBA_CACHE_DIR"] == str(tmp_path / "mine")
    assert not _user_folder(tmp_path).exists()


def test_numba_cache_folder_readable_by_others_is_closed(tmp_path, monkeypatch):
    jit = _cache_env(monkeypatch, tmp_path)
    want = _user_folder(tmp_path)
    want.mkdir(mode=0o755)
    os.chmod(want, 0o755)
    assert jit.set_default_cache_dir() == str(want)
    assert stat.S_IMODE(want.stat().st_mode) == 0o700


def test_numba_cache_folder_writable_by_group_or_world_is_not_trusted(
        tmp_path, monkeypatch):
    # numba unpickles its cache files, so a folder that another user may have
    # written to is neither repaired nor used
    jit = _cache_env(monkeypatch, tmp_path)
    want = _user_folder(tmp_path)
    for mode in (0o775, 0o757):
        want.mkdir(mode=mode, exist_ok=True)
        os.chmod(want, mode)
        monkeypatch.delenv("NUMBA_CACHE_DIR", raising=False)
        got = jit.set_default_cache_dir()
        assert got != str(want) and Path(got).parent == tmp_path
        assert stat.S_IMODE(Path(got).stat().st_mode) == 0o700
        assert stat.S_IMODE(want.stat().st_mode) == mode  # not repaired


def test_numba_cache_folder_of_another_user_is_refused(tmp_path):
    from dmrdetector.util import jit
    folder = tmp_path / "mine"
    assert jit._is_private_folder(str(folder), os.geteuid()) is True
    assert jit._is_private_folder(str(folder), os.geteuid() + 1) is False


def test_numba_cache_parent_writable_by_others_needs_the_sticky_bit(
        tmp_path, monkeypatch):
    jit = _cache_env(monkeypatch, tmp_path)
    home = tmp_path / "home"
    (home / ".cache").mkdir(parents=True, mode=0o700)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
    want = _user_folder(tmp_path)
    os.chmod(tmp_path, 0o777)  # anyone could rename a folder away
    got = jit.set_default_cache_dir()
    # neither the per-user folder nor a mkdtemp folder may sit there: the
    # replacement attack works on any name; the user cache folder is used
    assert not want.exists() and Path(got).parent == home / ".cache"
    monkeypatch.delenv("NUMBA_CACHE_DIR")
    os.chmod(home / ".cache", 0o777)  # no safe parent left: numba decides
    assert jit.set_default_cache_dir() is None
    assert "NUMBA_CACHE_DIR" not in os.environ
    os.chmod(tmp_path, 0o1777)  # like /tmp
    assert jit.set_default_cache_dir() == str(want)


def test_numba_cache_never_follows_a_link_or_a_file(tmp_path, monkeypatch):
    jit = _cache_env(monkeypatch, tmp_path)
    want = _user_folder(tmp_path)
    target = tmp_path / "other"
    target.mkdir()
    want.symlink_to(target)  # someone else's folder behind a link
    got = jit.set_default_cache_dir()
    assert got != str(want) and Path(got).parent == tmp_path
    assert stat.S_IMODE(Path(got).stat().st_mode) == 0o700
    assert list(target.iterdir()) == []


def test_command_line_sets_the_numba_cache_before_any_kernel_is_built(tmp_path):
    env = {
        k: v
        for k, v in os.environ.items()
        if k not in ("NUMBA_CACHE_DIR", "DMRDETECTOR_NO_JIT")
    }
    r = subprocess.run([
        sys.executable, "-c", "import os, dmrdetector.cli as c; "
        "c.main(['config', 'show', 'run']); "
        "print(os.environ['NUMBA_CACHE_DIR'])"
    ],
                       capture_output=True,
                       text=True,
                       env=dict(env, PYTHONPATH=SRC, TMPDIR=str(tmp_path)))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().splitlines()[-1] == str(_user_folder(tmp_path))


def test_version_and_help_make_no_cache_folder(tmp_path):
    for flag in ("--version", "--help"):
        r = subprocess.run([sys.executable, "-m", "dmrdetector", flag],
                           capture_output=True,
                           text=True,
                           env=dict(os.environ,
                                    PYTHONPATH=SRC,
                                    TMPDIR=str(tmp_path)))
        assert r.returncode == 0, r.stderr
    assert list(tmp_path.iterdir()) == []


def test_jit_sets_the_numba_cache_for_callers_that_skip_the_command_line(
        tmp_path):
    env = {
        k: v
        for k, v in os.environ.items()
        if k not in ("NUMBA_CACHE_DIR", "DMRDETECTOR_NO_JIT")
    }
    r = subprocess.run([
        sys.executable, "-c", "import os, numba; "
        "from dmrdetector.util import jit; jit.jit(cache=False); "
        "assert numba.config.CACHE_DIR == os.environ['NUMBA_CACHE_DIR']; "
        "print(os.environ['NUMBA_CACHE_DIR'])"
    ],
                       capture_output=True,
                       text=True,
                       env=dict(env, PYTHONPATH=SRC, TMPDIR=str(tmp_path)))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == str(_user_folder(tmp_path))


def _seed(root, name="seed 1"):
    """Writes a fake seed: two kernel files, a link and a stray file."""
    sub = root / "fileio_abc"
    sub.mkdir(parents=True)
    (root / "SEED").write_text(name)
    (sub / "f.f-1.py312.nbi").write_bytes(b"index " + name.encode())
    (sub / "f.f-1.py312.0.nbc").write_bytes(b"data " + name.encode())
    (sub / "stray.txt").write_text("x")
    (sub / "link.nbc").symlink_to(root / "SEED")
    return root


def test_seed_copies_the_kernels_once_with_private_modes(
        tmp_path, monkeypatch):
    from dmrdetector.util import jit
    seed, dest = _seed(tmp_path / "seed"), tmp_path / "cache"
    dest.mkdir(mode=0o700)
    monkeypatch.setenv(jit.SEED_ENV, str(seed))
    assert jit.seed_cache(str(dest)) is True
    got = sorted(p.name for p in (dest / "fileio_abc").iterdir())
    assert got == ["f.f-1.py312.0.nbc", "f.f-1.py312.nbi"]  # no link, no stray
    index = dest / "fileio_abc" / "f.f-1.py312.nbi"
    assert index.read_bytes() == b"index seed 1"
    assert stat.S_IMODE(index.stat().st_mode) == 0o600
    assert (dest / ".seeded").read_text() == "seed 1"
    assert jit.seed_cache(str(dest)) is False  # same seed: nothing to do
    (seed / "SEED").write_text("seed 2")
    (seed / "fileio_abc" / "f.f-1.py312.nbi").write_bytes(b"index seed 2")
    assert jit.seed_cache(str(dest)) is True  # another seed replaces the files
    assert index.read_bytes() == b"index seed 2"


def test_seed_failure_leaves_numba_to_compile(tmp_path, monkeypatch):
    from dmrdetector.util import jit
    seed = _seed(tmp_path / "seed")
    blocked = tmp_path / "blocked"
    blocked.write_text("a file, not a folder")
    monkeypatch.setenv(jit.SEED_ENV, str(seed))
    assert jit.seed_cache(str(blocked)) is False
    assert not list(tmp_path.glob("**/*.tmp*"))
    monkeypatch.setenv(jit.SEED_ENV, str(tmp_path / "no_such_seed"))
    assert jit.seed_cache(str(tmp_path / "cache")) is False
    monkeypatch.setenv(jit.SEED_ENV, "")
    assert jit.seed_cache(str(tmp_path / "cache")) is False
    assert not (tmp_path / "cache").exists()


def test_default_and_caller_folders_are_seeded_once_per_process(
        tmp_path, monkeypatch):
    jit = _cache_env(monkeypatch, tmp_path)
    monkeypatch.setenv(jit.SEED_ENV, str(_seed(tmp_path / "seed")))
    want = _user_folder(tmp_path)
    assert jit.set_default_cache_dir() == str(want)
    assert (want / ".seeded").is_file()
    mine = tmp_path / "mine"
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(mine))
    jit._seed_once.cache_clear()
    assert jit.set_default_cache_dir() is None  # the variable stays
    assert os.environ["NUMBA_CACHE_DIR"] == str(mine)
    assert (mine / "fileio_abc" / "f.f-1.py312.nbi").is_file()
    (mine / ".seeded").unlink()
    jit.set_default_cache_dir()  # second call of the process: no new copy
    assert not (mine / ".seeded").exists()


def test_seeded_cache_is_used_by_a_fresh_process(tmp_path):
    """A kernel compiled in one cache folder loads from a copy of it."""
    kernel = tmp_path / "kern.py"
    kernel.write_text(
        "import numpy as np\n"
        "from dmrdetector.util.jit import jit\n"
        "@jit\n"
        "def total(x):\n"
        "    s = 0.0\n"
        "    for v in x:\n"
        "        s += v * v\n"
        "    return s\n"
        "print('TOTAL', total(np.arange(5.0)))\n")
    base = {
        k: v
        for k, v in os.environ.items()
        if k not in ("NUMBA_CACHE_DIR", "DMRDETECTOR_NO_JIT",
                     "NUMBA_DEBUG_CACHE")
    }

    def run(cache, seed, debug):
        env = dict(base, PYTHONPATH=SRC, NUMBA_CACHE_DIR=str(cache),
                   DMRDETECTOR_NUMBA_SEED=str(seed))
        if debug:
            env["NUMBA_DEBUG_CACHE"] = "1"
        r = subprocess.run([sys.executable, str(kernel)], capture_output=True,
                           text=True, env=env)
        assert r.returncode == 0, r.stderr
        assert "TOTAL 30.0" in r.stdout
        return r.stdout

    first = tmp_path / "first"
    run(first, "", False)  # compiles and saves
    (first / "SEED").write_text("test seed")
    second = tmp_path / "second"
    out = run(second, first, True)
    assert "data loaded" in out and "data saved" not in out


def test_seed_identity_names_the_cache_folder_of_each_build(
        tmp_path, monkeypatch):
    jit = _cache_env(monkeypatch, tmp_path)
    a = _seed(tmp_path / "a", "dmrdetector 3.0.10 numba cpu generic\nm 1")
    b = _seed(tmp_path / "b", "dmrdetector 3.0.10 numba cpu generic\nm 2")
    ia, ib = jit.seed_identity(str(a)), jit.seed_identity(str(b))
    assert re.fullmatch("[0-9a-f]{12}", ia) and ia != ib
    assert jit.seed_identity(str(tmp_path / "none")) == "noseed"
    assert jit.seed_identity("") == "noseed"
    monkeypatch.setenv(jit.SEED_ENV, str(a))
    assert jit.cache_folder_name(7) == f"dmrdetector-numba-7-{ia}"
    # two builds never share a folder, and each holds its own kernels only
    got = {}
    for key, seed in (("a", a), ("b", b)):
        monkeypatch.setenv(jit.SEED_ENV, str(seed))
        monkeypatch.delenv("NUMBA_CACHE_DIR", raising=False)
        got[key] = Path(jit.set_default_cache_dir())
    assert got["a"] != got["b"] and got["a"].parent == got["b"].parent
    index = Path("fileio_abc") / "f.f-1.py312.nbi"
    assert (got["a"] / index).read_bytes() == (a / index).read_bytes()
    assert (got["b"] / index).read_bytes() == (b / index).read_bytes()
    assert (got["a"] / index).read_bytes() != (got["b"] / index).read_bytes()


def test_fallback_folder_goes_at_exit_and_the_user_folder_stays(tmp_path):
    seed = _seed(tmp_path / "seed")
    tmp = tmp_path / "tmp"
    tmp.mkdir(mode=0o700)  # the umask of the test machine must not matter
    env = {
        k: v
        for k, v in os.environ.items()
        if k not in ("NUMBA_CACHE_DIR", "DMRDETECTOR_NO_JIT")
    }
    env.update(PYTHONPATH=SRC,
               TMPDIR=str(tmp),
               DMRDETECTOR_NUMBA_SEED=str(seed))
    code = ("from dmrdetector.util import jit; "
            "print(jit.set_default_cache_dir())")

    def call():
        r = subprocess.run([sys.executable, "-c", code],
                           capture_output=True,
                           text=True,
                           env=env)
        assert r.returncode == 0, r.stderr
        return Path(r.stdout.strip().splitlines()[-1])

    kept = call()
    assert (kept / "fileio_abc" / "f.f-1.py312.nbi").is_file()
    assert call() == kept  # the per-user folder is reused and stays
    os.chmod(kept, 0o775)  # now it is not trusted: a folder for this run only
    scratch = call()
    assert scratch != kept and scratch.parent == tmp
    assert not scratch.exists()  # seeded copies went with the process
    assert kept.is_dir() and (kept / ".seeded").is_file()


def test_fallback_folder_is_removed_only_by_the_process_that_made_it(tmp_path):
    from dmrdetector.util import jit
    folder = tmp_path / "f"
    folder.mkdir()
    jit._remove_folder(str(folder), os.getpid() + 1)  # e.g. a forked worker
    assert folder.is_dir()
    jit._remove_folder(str(folder), os.getpid())
    assert not folder.exists()


def test_seed_copy_uses_unique_temporary_names(tmp_path, monkeypatch):
    from dmrdetector.util import jit
    seed, dest = _seed(tmp_path / "seed"), tmp_path / "cache"
    monkeypatch.setenv(jit.SEED_ENV, str(seed))
    sub = dest / "fileio_abc"
    sub.mkdir(parents=True, mode=0o700)
    # a file that another container with the same process id is writing
    foreign = sub / f"f.f-1.py312.nbi.tmp{os.getpid()}"
    foreign.write_bytes(b"foreign")
    assert jit.seed_cache(str(dest)) is True
    assert foreign.read_bytes() == b"foreign"
    assert (sub / "f.f-1.py312.nbi").read_bytes() == b"index seed 1"


def test_parallel_first_runs_seed_one_consistent_folder(tmp_path, monkeypatch):
    from dmrdetector.util import jit
    seed, dest = _seed(tmp_path / "seed"), tmp_path / "cache"
    dest.mkdir(mode=0o700)
    monkeypatch.setenv(jit.SEED_ENV, str(seed))
    results = []
    threads = [
        threading.Thread(
            target=lambda: results.append(jit.seed_cache(str(dest))))
        for _ in range(8)
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert len(results) == 8 and any(results)
    names = sorted(p.name for p in dest.rglob("*") if p.is_file())
    assert names == [".seeded", "f.f-1.py312.0.nbc",
                     "f.f-1.py312.nbi"]  # no temporary file is left
    assert (dest / "fileio_abc" /
            "f.f-1.py312.0.nbc").read_bytes() == b"data seed 1"


def test_seed_copy_writes_data_files_before_index_files(tmp_path, monkeypatch):
    from dmrdetector.util import jit
    seed, dest = _seed(tmp_path / "seed"), tmp_path / "cache"
    dest.mkdir(mode=0o700)
    monkeypatch.setenv(jit.SEED_ENV, str(seed))
    moved = []
    real = os.replace
    monkeypatch.setattr(
        os, "replace", lambda a, b:
        (moved.append(os.path.basename(b)), real(a, b))[1])
    assert jit.seed_cache(str(dest)) is True
    assert moved == ["f.f-1.py312.0.nbc", "f.f-1.py312.nbi", ".seeded"]
