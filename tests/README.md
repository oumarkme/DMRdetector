# Tests

| Directory | Content |
|---|---|
| `unit/` | pytest unit, fuzz, integration, thread-invariance, container-recipe and privacy tests. They use synthetic data only. `helpers.py` writes the synthetic inputs; `oracles/` holds two frozen reference implementations that some tests compare against. |
| `smoke/` | `run_smoke.py`: `build` on synthetic controls, then `run` on synthetic cases, with checks on every output file. `generate_data.py` writes the demo data set of the README (8 controls, 2 cases, three BED formats, a small GTF, relaxed thresholds), reusing the data model of `run_smoke.py`. |

Run inside any Python 3.12 environment that has the packages of `env/requirements.lock`, plus `pip install pytest`, from the repository root (pytest is not part of the image). The project has no `pytest.ini`, `setup.cfg` or `tox.ini`, and `pyproject.toml` holds no pytest configuration; `tests/conftest.py` puts `src/` on `sys.path`, so no installation is needed and the tests import the checkout, also when a copy of the package is installed in the same environment.

```bash
pytest tests/unit                                # compiled kernels
DMRDETECTOR_NO_JIT=1 pytest tests/unit           # plain-Python kernels; tests of the compiled paths are skipped
DMRDETECTOR_NO_LIBDEFLATE=1 pytest tests/unit    # zlib instead of libdeflate
python tests/smoke/run_smoke.py --threads 8
```

If a user-level pytest configuration interferes, add `-c /dev/null --rootdir . -p no:cacheprovider`. The test data folder is `tests/data` or `DMRDETECTOR_TEST_DATA`; the unit tests that need the real GTF read `Homo_sapiens.GRCh38.86.gtf.gz` there, or the file named by `DMRDETECTOR_TEST_GTF`.

Through the image, which has the locked packages but no pytest (this is the route we tested; the sources in the image are `/opt/dmrdetector/src`). Run from the repository root with the image built there (`dmrdetector.sif`):

```bash
pip install --target pytest_lib pytest      # pytest is pure Python; any pip will do
export DMRDETECTOR_TEST_GTF=/path/to/Homo_sapiens.GRCh38.86.gtf.gz   # needed by the smoke test and some unit tests
apptainer exec --bind "$PWD" --bind "$(dirname "$DMRDETECTOR_TEST_GTF")" --env PYTHONPATH="$PWD/pytest_lib" \
  dmrdetector.sif python -m pytest -c /dev/null --rootdir . -p no:cacheprovider tests/unit
apptainer exec --bind "$PWD" --bind "$(dirname "$DMRDETECTOR_TEST_GTF")" --env DMRDETECTOR_SRC=/opt/dmrdetector/src \
  --env DMRDETECTOR_TEST_GTF dmrdetector.sif python tests/smoke/run_smoke.py --threads 8
```

Apptainer passes `DMRDETECTOR_TEST_GTF` into the container when it is set in the shell. Add `--env DMRDETECTOR_NO_JIT=1` to repeat the unit tests with plain-Python kernels.

`settings.toml` in the repository root is generated, never edited by hand: `PYTHONPATH=src python -m dmrdetector config init --out settings.toml --force`. A unit test fails when it differs from `settings.default_toml()`.

`DMRDETECTOR_SRC` sets the source directory that the subprocess tests put on `PYTHONPATH` (default `src`). `test_container_recipes.py` checks that the Dockerfile and the Apptainer definition agree on the base image digest, the install command, the libdeflate pin, the environment and the smoke test, that every entry of `env/requirements.lock` has a hash, that the version is the same in the source, both recipes, the README and the CHANGELOG, and that `pyproject.toml` (the pip route) takes its version from the source, has dependency ranges that admit the locked versions, lists every data file of the package and holds no tool table besides the three under `[tool.setuptools]`. `test_pools.py` kills a worker of each pool and expects a RunError within a timeout, and `test_build_logs.py` runs real `build` calls to check that the logs hold no folder of the samples and that a resumed build repeats the warnings. `test_privacy.py` fails when a file that git would publish holds a lab storage path or a record number; it is skipped without git. `test_docs_commands.py` fails when a `dmrdetector` command, flag or settings key in a Markdown file does not exist, or when a link or a repository path in the docs is broken.

## Style

The code follows the [Google Python style guide](https://google.github.io/styleguide/pyguide.html): 80 columns, Google-style docstrings, sorted absolute imports of modules. The repository holds no formatter or linter configuration. We formatted with `yapf --style="{based_on_style: google, column_limit: 80, indent_width: 4}"` and checked with `ruff check --isolated --line-length 80 --select E,W,F,N,I,D`, leaving out the docstring rules that the Google convention disables. Test functions need no docstring. `tests/unit/oracles/` is frozen and not formatted.
