# Contributing

Bug reports and pull requests are welcome. For a bug, open an issue with the output of `dmrdetector --version`, the command line, and the `ERROR:` line or the log (`<sample_id>.log`). Remove sample IDs and paths that identify patients before you post anything.

## Run the tests

The tests run in any Python 3.12 environment that has the packages of `env/requirements.lock` plus `pytest`; pytest is not part of the image. From the repository root:

```bash
pytest tests/unit                                # compiled kernels
DMRDETECTOR_NO_JIT=1 pytest tests/unit           # plain-Python kernels
python tests/smoke/run_smoke.py --threads 8      # build and run on synthetic data
```

The smoke test needs an Ensembl GRCh38 GTF: run `dmrdetector get-genes --release 86 --outdir gtf` and set `DMRDETECTOR_TEST_GTF` to the file. [tests/README.md](tests/README.md) describes the suites and how to run them through the image.

Before a pull request, `tests/unit/test_privacy.py` must pass: no lab storage path or record number may appear in a tracked file. `tests/unit/test_docs_commands.py` checks that every `dmrdetector` command, flag and settings key in the Markdown files exists, so change the docs together with the code.

## Style

Python follows the [Google Python style guide](https://google.github.io/styleguide/pyguide.html): 80 columns, Google-style docstrings, sorted absolute imports of modules. Shell code follows the Google shell style guide and Markdown the Google Markdown style guide. The repository has no formatter or linter configuration; the commands we used are in [tests/README.md](tests/README.md). Do not change an algorithm, a default or an output file in a patch release; describe such a change in the pull request and in the CHANGELOG.

## Settings, versions and the lock file

- `settings.toml` is generated: `PYTHONPATH=src python -m dmrdetector config init --out settings.toml --force`. A unit test fails when it differs from `settings.default_toml()`.
- The version has one source, `src/dmrdetector/__init__.py`, which `pyproject.toml` reads. Both recipes in `env/`, the README badge and the top heading of the CHANGELOG must carry the same number (`tests/unit/test_container_recipes.py` checks this).
- `env/requirements.lock` is the only package list. A change of versions is a new release. To regenerate it, run the two commands in its header (they need [uv](https://docs.astral.sh/uv/)); keep the seven header lines. On 2026-10-07, uv 0.12.23 resolved the same versions and hashes as the committed file; only the `# via` comments differ between uv versions. Rebuild the image and run the full test set afterwards.
- `env/Dockerfile` and `env/dmrdetector.def` are written by hand and must stay in step; the test above compares them.
- `pyproject.toml` is the pip route. Its dependency ranges are the minor series of the locked versions (`name~=X.Y.Z`); change them together with the lock. The container recipes do not use it, so `test_pyproject_is_the_pip_route_of_the_same_package_and_versions` is what keeps it in step with the lock, the version source and the data files. It holds no linter or pytest configuration; add no `[tool.*]` table other than the three under `[tool.setuptools]`.

## Release

Update the CHANGELOG ([Keep a Changelog](https://keepachangelog.com/) format), bump the version, build the image of the release with the kernel cache required, run the tests through it, commit, and tag `vX.Y.Z` (annotated).

```bash
git archive vX.Y.Z | tar -x -C release_src && cd release_src
apptainer build --fakeroot --build-arg REQUIRE_CACHE=1 dmrdetector.sif env/dmrdetector.def
# Docker: docker build -f env/Dockerfile --build-arg DMRDETECTOR_REQUIRE_CACHE=1 -t dmrdetector:X.Y.Z .
```

With `REQUIRE_CACHE=1` the build stops when `env/warm_cache.sh` cannot compile the kernels, and the smoke test (`env/smoke_test.py`, run at the end of the build) fails unless the kernels load from the cache of the image and none is compiled again. Without it, a missing cache is a warning in the build log and the state `absent` in the smoke test, and the image compiles on first use. `env/warm_cache.sh` needs about 16 GB of memory.
