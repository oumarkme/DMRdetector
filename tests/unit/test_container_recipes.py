"""Checks that env/Dockerfile and env/dmrdetector.def describe one image.

They must agree on everything that decides versions: base image, install
command, libdeflate, environment and label. The lock must be fully hashed.
pyproject.toml (the pip route) must agree with the lock and the source.
"""
import fnmatch
import os
from pathlib import Path
import re
import shlex
import subprocess
import tomllib

from dmrdetector import cli


ROOT = Path(__file__).resolve().parents[2]
ENV = ROOT / "env"
DOCKERFILE = (ENV / "Dockerfile").read_text()
DEF = (ENV / "dmrdetector.def").read_text()
LOCK = (ENV / "requirements.lock").read_text()
INSTALL = "--require-hashes --only-binary=:all:"


def _joined(text):
    return text.replace("\\\n", " ")


def _version():
    return re.search(r'__version__ = "([^"]+)"', (ROOT / "src" / "dmrdetector" /
                                                  "__init__.py").read_text())[1]


def test_base_image_pinned_by_digest_and_identical_in_both_files():
    df = re.findall(r"(?m)^FROM (\S+)", DOCKERFILE)
    de = re.findall(r"(?m)^From: (\S+)", DEF)
    assert len(df) == len(de) == 2 and len(set(df + de)) == 1
    assert re.search(r"@sha256:[0-9a-f]{64}$", df[0])


def test_python_of_the_base_tag_is_the_python_of_the_lock():
    py = re.match(r"python:(\d+\.\d+)\.",
                  re.findall(r"(?m)^FROM (\S+)", DOCKERFILE)[0])[1]
    assert py == "3.12"
    assert "Python 3.12" in LOCK.splitlines(
    )[0] and "--python-version 3.12" in LOCK
    assert f'python{py}/site-packages' in DOCKERFILE
    assert f'python{py}/site-packages' in DEF


def test_install_command_is_the_same_hash_locked_pip_call():
    for text in (_joined(DOCKERFILE), _joined(DEF)):
        calls = re.findall(r"pip --python /opt/venv/bin/python install ([^&]*)",
                           text)
        assert len(calls) == 1
        assert INSTALL in calls[0] and "-r /build/requirements.lock" in calls[0]
    assert "env/requirements.lock /build/requirements.lock" in DOCKERFILE
    assert "env/requirements.lock /build/requirements.lock" in DEF
    assert sorted(
        p.name for p in ENV.glob("requirements*")) == ["requirements.lock"]


def test_every_lock_entry_has_a_hash():
    text = _joined(LOCK)
    reqs = [
        ln for ln in text.splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    assert reqs
    for ln in reqs:
        assert re.match(r"^[A-Za-z0-9._-]+==\S+ ", ln), ln
        assert "--hash=sha256:" in ln, ln
    for name in ("numpy", "scipy", "pandas", "pysam", "numba", "llvmlite"):
        assert re.search(rf"(?m)^{name}==", LOCK), name
    assert "pytest" not in LOCK  # the tester installs pytest on top


def test_libdeflate_pin_is_in_both_files_and_checked_against_a_sha256():
    pins = []
    for text in (DOCKERFILE, DEF):
        v = re.search(r'LIBDEFLATE_DEB_VERSION[= "]+([0-9][^" \n]*)', text)[1]
        h = re.search(r'LIBDEFLATE_DEB_SHA256[= "]+([0-9a-f]{64})', text)[1]
        assert 'apt-get download libdeflate0=$LIBDEFLATE_DEB_VERSION' in text
        assert 'sha256sum -c' in text
        # where the code looks first (sys.prefix/lib)
        assert "/opt/venv/lib/libdeflate.so.0" in text
        pins.append((v, h))
    assert pins[0] == pins[1]
    assert pins[0][0] in (
        ENV /
        "smoke_test.py").read_text()  # the smoke test checks the pinned version


def _dockerfile_env():
    out = {}
    for line in _joined(DOCKERFILE).splitlines():
        if line.startswith("ENV "):
            for tok in shlex.split(line[4:]):
                k, _, v = tok.partition("=")
                out[k] = v
    return out


def _def_env():
    out = {}
    block = re.search(r"(?ms)^%environment\n(.*?)^%", DEF + "\n%")[1]
    for line in block.splitlines():
        m = re.match(r'\s*export (\w+)="(.*)"$', line)
        if m:
            v = re.sub(r"^\$\{\w+:-(.*)\}$", r"\1",
                       m[2])  # ${K:-default}: overridable at run time
            out[m[1]] = v
    return out


def test_environment_lines_are_the_same():
    final = {
        k: v
        for k, v in _dockerfile_env().items()
        if k not in ("PIP_NO_CACHE_DIR", "PIP_DISABLE_PIP_VERSION_CHECK",
                     "LIBDEFLATE_DEB_VERSION", "LIBDEFLATE_DEB_SHA256")
    }
    assert final == _def_env()
    assert final["PYTHONPATH"] == "" and final["NUMBA_NUM_THREADS"] == "1"
    # glibc picks FMA or plain libm variants from the CPU; the image fixes the
    # plain ones (env/smoke_test.py checks that the setting takes effect)
    assert final["GLIBC_TUNABLES"] == "glibc.cpu.hwcaps=-FMA,-FMA4"
    # the recipes leave NUMBA_CACHE_DIR to the caller: when it is unset, the
    # program picks a private per-user folder (util.jit.set_default_cache_dir)
    assert "NUMBA_CACHE_DIR" not in final
    assert "NUMBA_CACHE_DIR" not in DOCKERFILE + DEF


def test_smoke_test_and_entrypoint_are_the_same():
    cmd = "dmrdetector --version && python /opt/dmrdetector/smoke_test.py"
    assert f"RUN {cmd}" in DOCKERFILE
    strict = cmd.replace(
        "&& python",
        '&& DMRDETECTOR_REQUIRE_CACHE="{{ REQUIRE_CACHE }}" python')
    assert re.search(r"(?m)^%test\n\s*" + re.escape(strict) + r"\s*$", DEF)
    assert 'ENTRYPOINT ["dmrdetector"]' in DOCKERFILE
    assert 'exec dmrdetector' in DEF


def test_source_launcher_and_pth_are_installed_the_same_way():
    for text in (_joined(DOCKERFILE), _joined(DEF)):
        assert '/opt/dmrdetector/src' in text
        assert 'env/launcher.py /opt/dmrdetector/launcher.py' in text
        assert (
            "echo /opt/dmrdetector/src > "
            "/opt/venv/lib/python3.12/site-packages/dmrdetector.pth") in text
        assert ("install -m 755 /opt/dmrdetector/launcher.py "
                "/opt/venv/bin/dmrdetector") in text
        assert ("compileall -q --invalidation-mode unchecked-hash /opt/venv "
                "/opt/dmrdetector/src") in text
    assert (ROOT / "src" / "dmrdetector" / "data" / "assets" /
            "par_hg38.bed").is_file()
    launcher = (ENV / "launcher.py").read_text()
    assert launcher.startswith(
        "#!/opt/venv/bin/python\n"
    ) and "from dmrdetector.cli import main" in launcher


def test_version_has_one_source():
    """Checks that the version of __init__.py is the one everywhere.

    The recipes, the README, CITATION.cff and the top CHANGELOG heading must
    match it.
    """
    ver = _version()
    assert re.fullmatch(r"\d+\.\d+\.\d+", ver)
    assert f'org.opencontainers.image.version="{ver}"' in DOCKERFILE
    assert f"org.opencontainers.image.version {ver}" in DEF
    readme = (ROOT / "README.md").read_text()
    assert f'badge/version-{ver}-' in readme
    assert f'# dmrdetector {ver}' in readme
    assert f'dmrdetector:{ver}' in readme
    assert f'version: "{ver}"' in (ROOT / "CITATION.cff").read_text()
    released = re.findall(r"(?m)^## \[(\d[^\]]*)\]",
                          (ROOT / "CHANGELOG.md").read_text())
    assert released[0] == ver


def test_no_build_scripts_and_no_test_or_linter_config_files():
    # pyproject.toml is allowed (the pip route); it is checked below
    for name in ("pytest.ini", "setup.cfg", "tox.ini", "setup.py", "ruff.toml",
                 ".ruff.toml", ".pylintrc", ".flake8", "env/build.py",
                 "env/requirements.in", "env/requirements-test.in",
                 "env/requirements-test.lock", "env/build-requirements.in",
                 "env/build-requirements.lock"):
        assert not (ROOT / name).exists(), name
    assert (ROOT / "tests" /
            "conftest.py").is_file()  # puts src/ on sys.path for pytest


def _lock_versions():
    """Returns {name: version} of the top-level packages in the lock."""
    top = re.search(r"(?m)^# Top-level packages: ([^(]+)\(", LOCK)[1]
    names = [n.strip() for n in top.split(",")]
    return {n: re.search(rf"(?m)^{n}==(\S+) ", LOCK)[1] for n in names}


def test_pyproject_is_the_pip_route_of_the_same_package_and_versions():
    """Checks pyproject.toml against the lock, the source and the data files.

    The recipes do not build it, so this test is what keeps it from drifting:
    one version source, dependencies that admit the locked versions, the
    console script, every packaged data file, and no tool tables other than
    the three under [tool.setuptools].
    """
    meta = tomllib.loads((ROOT / "pyproject.toml").read_text())
    project, tool = meta["project"], meta["tool"]
    # no linter, formatter or pytest configuration (project rule)
    assert set(tool) == {"setuptools"}
    assert set(tool["setuptools"]) <= {"dynamic", "packages", "package-data"}
    assert meta["build-system"]["build-backend"] == "setuptools.build_meta"
    # one version source
    assert "version" not in project and project["dynamic"] == ["version"]
    assert tool["setuptools"]["dynamic"] == {
        "version": {
            "attr": "dmrdetector.__version__"
        }
    }
    assert project["name"] == "dmrdetector"
    # only the Python minor of the base image (and of the lock) was run, so
    # the bound admits that minor and no other
    base = re.findall(r"(?m)^FROM (\S+)", DOCKERFILE)[0]
    major, minor = re.match(r"python:(\d+)\.(\d+)\.", base).groups()
    assert project["requires-python"] == (
        f">={major}.{minor},<{major}.{int(minor) + 1}")
    assert project["requires-python"] == ">=3.12,<3.13"
    assert (f"Programming Language :: Python :: {major}.{minor}"
            in project["classifiers"])
    assert project["license"] == "MIT" and (ROOT / "LICENSE").is_file()
    # dependencies: the top-level packages of the lock, ranges that admit it
    locked = _lock_versions()
    assert set(locked) == {"numpy", "scipy", "pandas", "pysam", "numba"}
    specs = {}
    for dep in project["dependencies"]:
        m = re.fullmatch(r"([A-Za-z0-9_.-]+)~=(\d+(?:\.\d+)+)", dep)
        assert m, f"{dep}: use a compatible-release range (name~=X.Y.Z)"
        specs[m[1]] = tuple(int(x) for x in m[2].split("."))
    assert set(specs) == set(locked)
    for name, lower in specs.items():
        have = tuple(int(x) for x in locked[name].split("."))
        # ~=X.Y.Z means >=X.Y.Z and ==X.Y.*
        assert have >= lower and have[:len(lower) - 1] == lower[:-1], name
    # command
    assert project["scripts"] == {"dmrdetector": "dmrdetector.cli:main"}
    assert callable(cli.main)
    # data files: every non-Python file of the package is listed
    assert tool["setuptools"]["packages"]["find"]["where"] == ["src"]
    patterns = tool["setuptools"]["package-data"]["dmrdetector"]
    pkg = ROOT / "src" / "dmrdetector"
    data = [
        p.relative_to(pkg).as_posix()
        for p in pkg.rglob("*")
        if p.is_file() and p.suffix != ".py" and "__pycache__" not in p.parts
    ]
    assert (pkg / "data" / "assets" / "par_hg38.bed").is_file()

    def packaged(rel):
        # setuptools globs per path segment: '*' does not match '/', which
        # fnmatch does
        return any(
            fnmatch.fnmatchcase(rel, pat) and rel.count("/") == pat.count("/")
            for pat in patterns)

    assert not packaged("data/assets/subfolder/file.bed")
    for rel in data:
        assert packaged(rel), rel


def test_kernel_cache_is_compiled_in_the_build_stage_for_a_generic_cpu():
    """Checks that both recipes bake the numba cache the same way.

    The cache is compiled for the CPU name that the final image sets, so it is
    valid on every x86-64 host; the final stage only copies it, and the
    program looks for it where it is copied to.
    """
    from dmrdetector.util import jit
    assert _dockerfile_env()["NUMBA_CPU_NAME"] == "generic"
    assert _def_env()["NUMBA_CPU_NAME"] == "generic"
    run = "bash /build/warm_cache.sh /opt/dmrdetector/numba-cache /build/smoke"
    docker_build, docker_final = DOCKERFILE.split(" AS final")
    def_build, def_final = DEF.split("Stage: final")
    assert ("ENV NUMBA_CPU_NAME=generic\nARG DMRDETECTOR_REQUIRE_CACHE=0\n"
            "RUN " + run in docker_build)
    assert ('export NUMBA_CPU_NAME="generic"\n'
            '  export DMRDETECTOR_REQUIRE_CACHE="{{ REQUIRE_CACHE }}"\n  ' + run
            in def_build)
    for text in (docker_build, def_build):
        assert "warm_cache.sh /build/warm_cache.sh" in text
        for name in ("generate_data.py", "run_smoke.py", "warm_big.py"):
            assert name in text
    for text in (docker_final, def_final):
        assert "warm_cache" not in text and "generate_data" not in text
        assert "/opt/dmrdetector/numba-cache" in text
    assert "COPY --from=build /opt/dmrdetector/numba-cache " in docker_final
    assert ("/opt/dmrdetector/numba-cache /opt/dmrdetector/numba-cache"
            in def_final)
    assert jit.SEED_DIR == "/opt/dmrdetector/numba-cache"
    script = (ENV / "warm_cache.sh").read_text()
    assert script.startswith("#!/bin/bash\n")
    assert "DMRDETECTOR_NUMBA_SEED=" in script


def test_release_builds_can_require_the_kernel_cache_in_both_recipes():
    """Checks the strict option: off by default, one switch per recipe.

    With it, the warm step and the smoke test fail instead of leaving the
    cache out with a warning. The README and CONTRIBUTING.md give the command.
    """
    # Docker: the ARG is declared in the stage that uses it (an ARG does not
    # carry over a FROM) and defaults to 0
    stages = DOCKERFILE.split("AS final")
    assert all("ARG DMRDETECTOR_REQUIRE_CACHE=0" in t for t in stages)
    # Apptainer: %arguments in every stage (Apptainer needs it per stage)
    def_stages = DEF.split("Stage: final")
    assert all("%arguments\n  REQUIRE_CACHE=0\n" in t for t in def_stages)
    for name in ("README.md", "CONTRIBUTING.md"):
        text = (ROOT / name).read_text()
        assert "--build-arg REQUIRE_CACHE=1" in text, name
    smoke = (ENV / "smoke_test.py").read_text()
    assert 'os.environ.get("DMRDETECTOR_REQUIRE_CACHE") == "1"' in smoke
    assert "assert not REQUIRE_CACHE" in smoke


def _run_warm_cache(tmp_path, require):
    cache = tmp_path / "cache"
    cache.mkdir(exist_ok=True)
    (cache / "left.nbi").write_text("x")
    smoke = tmp_path / "smoke"  # no generate_data.py: the first step fails
    smoke.mkdir(exist_ok=True)
    env = {k: v for k, v in os.environ.items() if k != "NUMBA_CACHE_DIR"}
    if require:
        env["DMRDETECTOR_REQUIRE_CACHE"] = "1"
    else:
        env.pop("DMRDETECTOR_REQUIRE_CACHE", None)
    r = subprocess.run(
        ["bash", str(ENV / "warm_cache.sh"),
         str(cache), str(smoke)],
        capture_output=True,
        text=True,
        env=env)
    return r, cache


def test_warm_step_fails_only_when_the_cache_is_required(tmp_path):
    r, cache = _run_warm_cache(tmp_path, require=False)
    assert r.returncode == 0 and "WARNING: kernel cache not built" in r.stderr
    assert list(cache.iterdir()) == []  # no half-made cache is kept
    r, cache = _run_warm_cache(tmp_path, require=True)
    assert r.returncode == 1
    assert "ERROR: kernel cache not built (DMRDETECTOR_REQUIRE_CACHE=1)" in (
        r.stderr)
    assert list(cache.iterdir()) == []


def test_seed_file_names_the_version_the_cpu_and_a_hash_of_the_cache():
    script = (ENV / "warm_cache.sh").read_text()
    assert "manifest sha256" in script and "sha256sum" in script
    assert "no kernel was compiled" in script  # an empty cache is a failure
