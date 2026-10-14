# Changelog

All notable changes to this project are documented in this file. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses [Semantic Versioning](https://semver.org/).

## [3.2.0] - 2026-10-13

Minor release. The package can be installed with pip. Algorithm, parameters, defaults, options and output files are unchanged, and so are `env/launcher.py` and `env/requirements.lock`; the image recipes change only in the version label and in the image tag of the comments and examples. The image is not built from `pyproject.toml`.

### Added
- README: a note at the top names v2.2 (git tag `v2.2-manuscript`) as the version used in the manuscript.
- `pyproject.toml` (PEP 621, setuptools backend, `src/` layout). In a Python 3.12 virtual environment, `pip install .` and `pip install git+<repository URL>` install the package, the files of `src/dmrdetector/data/assets/` and the command `dmrdetector` (`dmrdetector.cli:main`). The version is read from `src/dmrdetector/__init__.py`, which stays its only source. `requires-python` is `>=3.12,<3.13`: Python is limited to 3.12 because only 3.12 was run. Lift the bound after the three unit-test modes (compiled kernels, `DMRDETECTOR_NO_JIT=1`, `DMRDETECTOR_NO_LIBDEFLATE=1`) and the output gate have run on 3.13. The dependencies are the five top-level packages of `env/requirements.lock` in compatible-release ranges: `numpy~=2.5.3`, `scipy~=1.18.1`, `pandas~=3.0.6`, `pysam~=0.24.1` and `numba~=0.68.0`. The extra `test` installs pytest. The file has no linter, formatter or pytest configuration; its only `[tool.*]` tables are the three under `[tool.setuptools]` (`dynamic`, `packages.find` and `package-data`), which setuptools needs to read the version and to find the package and its data files.
- README, "Install": a pip route above the container route, with the commands that install exactly the locked versions (`pip install --require-hashes --only-binary=:all: -r env/requirements.lock`, then `pip install --no-deps .`) and a statement of what only the image guarantees.
- Unit test `test_pyproject_is_the_pip_route_of_the_same_package_and_versions`: the version is dynamic and read from the source, the dependencies are those of the lock and their ranges admit the locked versions, the console script is `dmrdetector.cli:main`, every non-Python file of the package is listed as package data, `requires-python` is the Python minor of the base image with the next minor excluded, and `pyproject.toml` has no tool table besides setuptools.
- `doc/dependencies.md`: how to set `GLIBC_TUNABLES` outside the image, and why the package does not set it.

### Changed
- `pyproject.toml` returns after its removal in 3.0.2. The image still has one package list, `env/requirements.lock`, and one settings file, `settings.toml`; the ranges in `pyproject.toml` have to change together with the lock (`CONTRIBUTING.md`).
- `tests/unit/test_container_recipes.py` no longer requires that `pyproject.toml` is absent. It still requires that `pytest.ini`, `setup.cfg`, `tox.ini`, `setup.py` and the removed `env/` build files are absent, and now also that `ruff.toml`, `.ruff.toml`, `.pylintrc` and `.flake8` are.
- `.gitignore`: `/dist/` and `*.egg-info/`. A non-editable `pip install .` in the checkout writes `build/` and `src/dmrdetector.egg-info/` there.
- `tests/README.md`, `tests/conftest.py` and `CONTRIBUTING.md` describe the new file.
- Source and test compaction. No output file, log line, option, default or `--help` text changes. The arithmetic of the numerical result path is untouched (kernel, FFT smoothing, `rstats`, `numparse`, BGZF reader and writer, JIT kernels, aggregation); in the raw atlas writer only the removed option's branches went. Duplicated code was merged. There is one `bed.warn` instead of three `_stderr`, and one `runtime.run_command` instead of two error wrappers in `build` and `run`. `runtime.temp_output` takes the temporary name, the cleanup and the removal of leftovers of eight steps; `runtime.commit` and `os.replace` stay as each step had them. `runtime.write_atomic` writes the JSON files, the stamps, the logs and `settings.toml`. `atlas_io.unread_path` resolves the atlas path of the CpG island and ICR panels, and `steps/target.py` has one empty-track table. `runtime.pool_map` is built on `runtime.fork_pool`.
- Removed code that no command used. The `return_arrays` option of `atlas.build_atlas` and its helper were called by one test only. The `force` parameter of `runtime.guarded` and the `force` keys of the task dictionaries were never true. `Sample.label` always returned `case`; the raw DMR table names keep it. `bed.bgzip_writer` had one caller, and `runtime._init_worker` had one. The test `test_raw_atlas_return_arrays_independent_of_window` and the arrays block of `test_build_fast_equals_legacy` went with the option. The byte comparison of `test_raw_atlas_bytes_do_not_depend_on_window_or_budget` covers the same property. Two tests were added: `test_temp_output_removes_what_was_not_moved` and `test_unread_path_names_the_file_only_while_it_is_not_loaded`.
- `runtime.write_atomic` (JSON files, stamps, logs and `settings.toml`) removes its temporary file and re-raises when the write or the rename fails. In 3.1.0 each of these writers left its temporary file in the folder in that case. New test: `test_write_atomic_failure_leaves_no_temporary_file_and_the_old_file`.
- Three edge cases differ from 3.1.0 without changing any output file. `atlas build --tmpdir ""` now indexes in its own temporary folder, which it removes, instead of leaving one behind. `config init` writes through `.tmp<pid>.<name>` in the target folder instead of `<name>.tmp<pid>`. A worker that dies while `pool_map` submits tasks now gives the RunError of `fork_pool` instead of `BrokenProcessPool`.
- `doc/dependencies.md` names functions in the table of floating-point operations instead of line numbers, which change with every edit.

### Release checks
The installation checks, the output comparison with 3.1.0, the timings and what was not verified are in [doc/validation.md](doc/validation.md), section "Verification of 3.2.0".

## [3.1.0] - 2026-10-09

Minor release. Its aim is that one version writes the same data files on every x86-64 CPU. Smoothed values change in the last printed digits compared with 3.0.x. The largest absolute differences in `REF_METH` and `REF_VAR` are 6.8e-12 and 1.8e-12 on the autosomal atlas of the 69 manuscript controls (4.3e-12 and 6.4e-13 on chrY), and 1.7e-12 and 7.5e-13 on the 12-control set. Calls, regions, genes, directions, `N`, `COV`, `PHI`, `RAW_METH` and the panels do not change, so existing references stay valid and rebuilding them is optional. A reference built by 3.0.x and one built by 3.1.0 from the same controls differ only in these digits.

What the change achieves: with numpy's AVX512 loops switched off, with its AVX2 and FMA3 loops switched off as well, with three OpenBLAS core types, with numba kernels for Haswell and Skylake-AVX512, and with all of these together, the 12-control set (161 data files) and the demo set (215 data files) are identical in content to the files of the unchanged Ice Lake server. In 3.0.10, the numpy and glibc settings changed 22 of the 161 files and 20 to 33 of the 215. These are emulations on one CPU model; a real AMD machine and a machine without AVX2 were not available (`doc/validation.md`).

### Changed
- The Gaussian kernel of the smoothing is computed with `decimal` arithmetic (50 digits, checked against 80) and rounded once to float64 (`rstats.exp_rounded`), once per sigma. Before, `np.exp` returned other last bits on a CPU with AVX512F than on one without: 44 of 901 entries for sigma 150 (about 5%). The new values equal `exp` of glibc 2.41 and of R 4.4.3 in all entries for sigma up to 500; for sigma 1000 and 2000 they differ from glibc and R in 4 entries each (one unit in the last place), where the correctly rounded value is used. The argument is computed as in v2.2, `-(x^2) / (2 sigma^2)` in float64.
- The two spectra are multiplied with separate float64 multiplications and additions (`fftsmooth.multiply_complex`, in blocks of 2^20 values). `a *= b` of numpy 2.5.3 fuses them on a CPU with AVX2 and FMA3 and does not on another, which changed the smoothed atlas by up to 2.9e-12 when the AVX2 loops were switched off.
- Both container recipes set `GLIBC_TUNABLES=glibc.cpu.hwcaps=-FMA,-FMA4`. glibc picks the FMA or the plain variant of `sin`, `cos`, `exp`, `log` and `pow` from the CPU, and `scipy.fft` (twiddle factors) and `scipy.stats.t.sf` call them. The two variants differ in the last bit of 0.05% to 0.08% of the arguments, which changes the smoothed atlas by up to 2.9e-12 and table cells by up to 4e-12 relative. The image now uses the plain variant on every host. An image without the setting gave different files on a CPU with FMA and AVX2 than on one without them. A native installation outside the image is not covered.
- `rstats.r_round` takes the binary exponent from `np.frexp` (R uses `logb`) and the constant `log10(2)` from `math`; `10^digits` is a rounded `Fraction`. `np.log2` can round up to the next integer just below a power of 2 and numpy may evaluate it with an AVX512 loop. The rounding of the table columns is unchanged for values below 2^46.
- `numparse.P10_EXACT` is built from Python integers instead of `10.0**k`.

### Added
- Unit tests: the kernel against bit patterns printed by R 4.4 (`sprintf("%a")`), the kernel without `np.exp` and `math.exp`, `exp_rounded` for special values and one subnormal, `multiply_complex` against the scalar formula in several block sizes, `r_round` for large values, and the `GLIBC_TUNABLES` line in both recipes.
- `env/smoke_test.py` checks that the glibc setting works: for two arguments, `math.exp` and `math.sin` must return the values of the plain variant.
- `doc/dependencies.md`, section "Floating-point arithmetic and the CPU": the audit of the operations on the result path. `doc/validation.md`, section "CPU independence": the emulation settings and their results.

### Release checks
Tests, kernel and CPU-emulation results, the comparison with 3.0.10, the rebuild of the manuscript atlas and the timings are in [doc/validation.md](doc/validation.md) (sections "CPU independence", "Verification of 3.1.0" and "Release checks of 3.0.8 and 3.1.0").

## [3.0.10] - 2026-10-09

Patch release after a second independent review of 3.0.6 to 3.0.9. For valid input, the data files of `build` and `run` are byte-identical to those of 3.0.9 (see "Release checks of 3.0.10"). These behaviors change:

- A `type` value in the samplesheet other than a control (`CTRL`, `control`) or `CASE` stops `build` again, as in v2.2. 3.0.8 dropped such rows.
- `build` and `run` print new warnings: a VCF without a sample column or without a heterozygous SNP call, and a sex label that disagrees with the number of heterozygous SNPs on chrX outside the PARs. They change no output file.
- A resumed `build` repeats the preprocessing warnings of the steps it skips.
- `build.log` and `reference.json` no longer show the folders of the sample files, the imprinting gene list or the PAR file when `log_input_paths = "name"` (the default).
- A worker process that the kernel kills in the atlas, tabix or smoothing steps ends the program with an error. Before, the program could wait for ever.
- The per-user numba cache folder has a new name, and a folder that others could have written to is not used.
- A release build can require the compiled kernels (`--build-arg REQUIRE_CACHE=1`).

### Fixed
- `build.log` listed the folders of the sample files in the messages and in the `[error]` block of a failed build, because only the samplesheet was known to the redaction. Every methylation BED, VCF and haplotype BED of the samplesheet, and the folder that holds each, is now mapped to its name. In `reference.json`, the options `icr_gene_list` and `par_bed` are recorded by file name (`log_input_paths = "full"` keeps absolute paths).
- A `type` column with a value other than `CTRL`, `control` or `CASE` (any case), including an empty one, raises an error that lists the distinct values with the first row of each and names the fix. v2.2 stopped there (`type must be CASE or CTRL`); 3.0.8 dropped the row. The types of all rows are checked before the files of any row. This reverses the 3.0.8 decision on purpose, and the docstring "as in v2.2" is now accurate.
- The raw atlas on the tabix route (`DMRDETECTOR_NO_JIT=1` or a file the fast reader rejects), the smoothing on the pandas route and the fork branch of the atlas reader used `multiprocessing.Pool`, which waits for ever for a worker that the out-of-memory killer ended. They use `concurrent.futures.ProcessPoolExecutor` with the fork context (`runtime.fork_pool`), as `pool_map` does, and a killed worker raises `a worker process was killed, probably out of memory; use fewer --threads`. Task order, chunks and the numbers are unchanged. `ProcessPoolExecutor` cannot replace a worker after each task (`max_tasks_per_child` is not available with fork), so the pandas smoothing runs the FFT of each chromosome in its own short-lived process (`fftsmooth.smooth_in_child`, as the fast route does), which frees the plan cache of scipy.fft in the same way.
- A resumed `build` dropped the warnings that preprocessing logs inside the step (empty sample, no CpG left, chrY, sex, contigs): the progress log is rewritten at the start, and a skipped step logged nothing. The step now stores its WARNING and NOTICE lines in its stamp, and a run that skips the step logs them again. The stamp key of the step also holds the sex label, the scope setting and the PAR file, so a resume with a corrected label runs the checks again instead of repeating the old warning. The warning about records with `END - START != 1` now starts with `WARNING:`, which puts it into the `[warnings]` section of the log and into the replay.
- The numba cache folder of a user (`dmrdetector-numba-<uid>`) was repaired with `chmod` and then trusted even when it was group or world writable. Numba unpickles cache files, so such a folder allows code execution as the user. A folder is now used only when the effective user owns it, it is not a link and it has no group or world write bit, in a temporary directory that is not writable by others without the sticky bit. Otherwise the process makes a new folder with `tempfile.mkdtemp` (mode 0700) in the temporary directory, or in `$XDG_CACHE_HOME` (`~/.cache`) when others could rename a folder out of the temporary directory; when neither is safe, `NUMBA_CACHE_DIR` stays unset and numba uses its own default. A folder that others may only read is closed to 0700 as before.
- That fallback folder (and the kernels copied into it) stayed in the temporary directory. The process that made it now removes it at exit (`atexit`); a forked worker does not, and the per-user folder is never removed.

### Added
- Warnings in `preprocess`, both routes: `the VCF has no sample column (sites-only?)` and `the VCF has no heterozygous SNP call` (the heterozygous-SNP filter then removes nothing).
- Sex label check from the other side. The lab convention for atlas controls (`doc/sex_chromosomes.md`) is one cutoff on the heterozygous SNVs of non-PAR chrX: at least 40,000 for a female, fewer for a male. `build` and `run` count the heterozygous SNP sites outside the PARs from the sites that `preprocess` already holds (no second pass over the VCF) and warn when a sample labeled M has 40,000 or more or a sample labeled F has fewer. The check is skipped for a VCF with no or several sample columns and when the run does not split by scope. The count differs slightly from the lab script (multi-allelic sites count, FILTER is ignored). Nothing else changes.
- `DMRDETECTOR_REQUIRE_CACHE=1` for release builds: `env/warm_cache.sh` exits with 1 when a step fails or no kernel was compiled, and `env/smoke_test.py` fails when the baked cache is absent. Apptainer: `--build-arg REQUIRE_CACHE=1` (declared in `%arguments` of both stages); Docker: `--build-arg DMRDETECTOR_REQUIRE_CACHE=1`. Without it, a missing cache stays a warning in the build log and the state `absent` in the smoke test.
- The `SEED` file of the kernel cache carries a hash of the cache files (`manifest sha256`).
- Tests: `test_pools.py` (a worker killed with SIGKILL in each pool gives a RunError within a timeout; no `multiprocessing.Pool` in `src`), `test_build_logs.py` (real `build` runs), and tests for the samplesheet types, the cache folder, the VCF and sex warnings and both recipes.

### Changed
- The per-user cache folder is named `dmrdetector-numba-<euid>-<id>`, where `<id>` is the first 12 hexadecimal digits of the SHA-256 of the `SEED` file of the image (version, CPU name and a hash of the kernel files, written when the image is built; reading the 150-byte file costs nothing measurable). Versions and builds of the image never share a folder. Without a seed the name ends in `-noseed`.
- Temporary files of the cache copy come from `tempfile.mkstemp` in the destination folder (exclusive creation), not from `<name>.tmp<pid>`: a process id is not unique across containers that share a folder. The order is unchanged: data files first, then index files, then `.seeded`, each moved into place with `os.replace`.
- `env/smoke_test.py` seeds a temporary cache folder from the baked seed and runs six synthetic samples through `preprocess`, the raw atlas, the smoothing and the CpG island detection, so that every module with baked entries calls a kernel. It asserts that kernels were loaded from the cache, that none was compiled and that no cache file was written. It needs 2.5 s.

### Release checks of 3.0.10
- Image. Built with `apptainer build --fakeroot --build-arg REQUIRE_CACHE=1` from `git archive` of the commit before the tag (3 min 19 s to 4 min; 168,267,776 bytes). The `%test` section ran the new smoke test with the cache required: 18 kernels loaded from the cache in 7 modules, none compiled, no cache file written (2.5 s). The smoke test was also run against variants of the 3.0.9 image: with the cache folder removed it prints `absent` and goes on, and with `DMRDETECTOR_REQUIRE_CACHE=1` it fails; with a changed source file of a baked module (`cpg_islands.py`) it fails with `a kernel was compiled and saved`; with an unreadable seed it fails with `the kernels were not copied`. Apptainer 1.4.5 accepts `%arguments` in both stages and substitutes `{{ REQUIRE_CACHE }}` in `%post` and `%test`; without the build argument the default is 0.
- Comparison with 3.0.9 (the 3.0.9 image ran the same inputs, 36 threads for the real-data set, 8 for the demo). The md5 sums of all 161 data files of the real-data set (12-control `build` and 4 `run` calls) and of all 215 data files of the demo set equal those of 3.0.9. The logs differ in the version, time stamps, work folder names, and, in the demo only, in one new warning for each female sample (`sex F, but the VCF has 0 heterozygous SNPs on chrX outside the PARs`): the synthetic VCFs have SNPs on chr22 only. The log of the real-data `build` has no line that differs after removing time stamps and paths. The same comparison on the working tree through the 3.0.9 image gave the same result.
- Plain-Python route (`DMRDETECTOR_NO_JIT=1`), the route that uses the tabix raw atlas and the pandas smoothing, which changed from `multiprocessing.Pool` to `ProcessPoolExecutor`. A demo `build` gave 13 files (atlases, indexes, panels) with the same md5 of the decompressed content as 3.0.9. A spy on `runtime.fork_pool` showed 8 calls for the raw atlas and 4 for the smoothing in that build (the demo has one chromosome per scope, so each smoothing pool had one worker). The smoothing of 4 real chromosomes (chr19 to chr22, 2,871,130 sites) with 4 workers gave the same md5 as 3.0.9, in 37.8 s against 38.0 s and with a peak resident set size of the main process of 2.86 GB against 2.84 GB. A complete plain-Python `build` of the 12 real controls was not run. In the compiled demo `build`, the fork branch of the atlas reader ran 4 times.
- A killed worker. A script that ends one worker of the pandas smoothing with SIGKILL did not finish within 40 s on 3.0.9; on 3.0.10 the same script ends with `smoothing: a worker process was killed, probably out of memory; use fewer --threads`. `tests/unit/test_pools.py` does this for each pool, with a timeout.
- Sex label check on real data. For the 14 VCFs of the real-data set (7 male and 7 female samples by label), the count of heterozygous SNP sites outside the PARs was 8,273 to 19,245 for males and 69,442 to 74,165 for females; the scan of the 14 files took 12 s. The log of the 12-control `build` has no sex warning. The ranges of the lab convention are 6,532 to 19,199 and 60,362 to 83,844 (`doc/sex_chromosomes.md`).
- numpy dispatch on the same server (3.0.9 image, `NPY_DISABLE_CPU_FEATURES` for the AVX512 loops of numpy 2.5.3), see `doc/validation.md`: the output is not byte-identical. `np.exp` in the Gaussian kernel changes by one unit in the last place in 44 of 901 entries (sigma 150). In the real-data set, 22 of 161 data files differ (the smoothed atlases with their indexes, and the tables with their temporary copies). Positions, `N`, `COV`, `PHI`, `RAW_METH` and the five panels are identical. `REF_METH` differs in 44.6% and `REF_VAR` in 61.1% of the autosomal cells (largest absolute difference 2.9e-12 and 7.3e-13); all cells pass the atlas criterion of `doc/validation.md`. In the tables, the non-float cells are equal text and the float cells differ by at most 4.2e-12 relative. Byte identity across `NUMBA_CPU_NAME` values, as stated for 3.0.9, holds on a CPU where numpy uses AVX512F; byte identity on other CPUs is not established.
- Tests, through the image with the checked-out sources. `pytest tests/unit` gives 430 passed and 2 skipped with the compiled kernels and 411 passed and 21 skipped with `DMRDETECTOR_NO_JIT=1` (the tests of the compiled paths are skipped). With `DMRDETECTOR_NO_LIBDEFLATE=1` the result is 430 passed and 2 skipped. `tests/smoke/run_smoke.py --threads 8`: 96 of 96 checks. The privacy, documentation-command and container-recipe tests run on the host Python without the image: 49 passed. The thread-invariance tests (2, 4 and 36 threads) are part of the suites.
- Docker. The Docker recipe was never built: no Docker daemon was usable on the build server. Its `ARG` lines, and any effect of Docker on file modification times, are not tested. The Apptainer image carries all the timing and byte-identity evidence of 3.0.9 and 3.0.10.
- Not checked. Hosts whose CPU lacks AVX512F, and other glibc or numpy builds; a complete plain-Python `run` or `build` on real data; the peak memory of the child processes of the pandas smoothing (only the main process was measured); a worker pool that is stopped by a signal other than SIGKILL (a fallback cache folder then stays in the temporary directory).

## [3.0.9] - 2026-10-09

Patch release: the image carries the compiled numerical kernels, and the kernels are compiled for a generic x86-64 CPU. No result changed: the data files of `build` and `run` (atlases, panels, DMR tables, target-region tables) are byte-identical to those of 3.0.8 on the demo data and on 12 real controls (see "Release checks of 3.0.9"). The gain is in the first call of a fresh container: on the demo data `run` took 3.9 s instead of 11.6 s and `build` 62 s instead of 80 s. The release also records the verification of 3.0.8 against the manuscript run in `doc/validation.md`.

### Changed
- The recipes compile the kernels when the image is built. `env/warm_cache.sh` runs the README demo (three BED formats, both sexes, a target-region report with haplotype tracks) and one large file per format (`env/warm_big.py`: 200,000 CpGs, several 4 MiB pieces, so that the kernels also see read-only buffers) in the build stage with `NUMBA_CACHE_DIR=/opt/dmrdetector/numba-cache`; the final stage copies the folder (3.1 MB, 71 files; the image grew by 1.0 MB). A step that fails (for example on a build machine with less memory than `build` budgets) leaves the cache out with a WARNING, and the image then compiles on first use as before. `env/smoke_test.py` prints whether the cache is present.
- On its first call, the program copies the kernels of the image into the cache folder in use: `NUMBA_CACHE_DIR` when the caller set it, otherwise the private per-user folder of 3.0.7. Numba cannot use the read-only folder of the image, because it tests that a cache folder is writable. The copy happens once per folder and image (the file `.seeded` names the image), files get mode 0600, files are moved into place from a temporary name so that parallel first runs do not read half-written files, and any error leaves numba to compile. `DMRDETECTOR_NUMBA_SEED` names another seed folder; empty skips the copy. `DMRDETECTOR_NO_JIT=1` is unchanged.
- Both recipes set `NUMBA_CPU_NAME=generic` (the Apptainer recipe as `${NUMBA_CPU_NAME:-generic}`, so a user can override it). Numba keys the cache on the CPU name and features, so a cache compiled for the build host would miss on other CPU models; the generic cache is valid on every x86-64 host. The kernels use SSE2 only, which cost 4% in `build` and 3% in `run` (medians of alternating runs on a shared server).
- `doc/usage.md`, `doc/dependencies.md` and the README describe the cache and the CPU setting.

### Documentation
- `doc/validation.md`: the section "Verification of 3.0.8" holds the rebuild of the manuscript atlas and the `run` of all 75 samples of the manuscript run. Atlas: same sites, `N` and `COV`; `REF_METH` and `REF_VAR` within tolerance; `PHI` and `RAW_METH` equal to the 15th printed digit, with a largest relative difference of 1.01e-14; panels byte-identical. Tables: text-identical calls; float cells differ by up to 1.3e-13 relative against the manuscript atlas, and by up to 8.5e-11 against the rebuilt atlas, where three tail p-values (below 1e-54) differ by up to 2.4e-10. The table of compared items gets a row for the four CPU settings.
- The README has a section "Validation" with these figures. `doc/ALGORITHM.md` no longer says that the kept behaviors "reproduce" the manuscript numbers without a qualifier: the calls are equal and the floats agree within tolerance.

### Release checks of 3.0.9
- Cache key. In numba 0.68.0 the index of a function is keyed on the file name of the function, the numba version, the SHA-256 of the source file (not its modification time), and for each compiled signature on the target triple, the CPU name, the CPU features, and hashes of the bytecode and of the closure. The cache folder is named from the directory of the source file, so the cache is valid at `/opt/dmrdetector/src` only. After setting every source file to another modification time and deleting every `__pycache__` in a writable copy of the image, the first `run` loaded 14 compiled entries and saved none.
- Misses on real data. With the seed copied into a fresh cache folder for each job, the 12-control `build` loaded 193 entries and saved none, and each of the four `run` scenarios (one with a target-region report and haplotype tracks) saved none. The first baked cache (demo data only) missed one signature on real data (`methbed._parse` with a read-only buffer), which is why `warm_big.py` exists. Demo runs on plain text and plain gzip inputs, with 2 and 8 threads, also saved none.
- First-call time in a fresh container. Apptainer 1.4.5, non-root user, read-only image, 8 threads, demo data, 3 repeats, 3.0.8 first and 3.0.9 second. `build` took 79.7 to 80.2 s and 61.5 to 63.9 s. The first `run` took 11.6 to 11.8 s and 3.8 to 4.1 s. A second `run` in the same cache folder took 3.9 to 4.1 s with both. Docker was not tried: no daemon was usable on the server.
- CPU settings (3.0.8 image). The md5 sums of all data files equal for the host CPU (Ice Lake server), Haswell with AVX2 and FMA, Skylake-AVX512, and generic: 12-control `build` plus four `run` calls (161 data files), and the demo set (215 data files). Logs differ in time stamps and durations only. Numba 0.68.0 does not use SVML here (`icc_rt` is not installed), and fast-math is off.
- Speed of generic against host kernels. 3.0.8 image, warm cache, alternating order, 36 threads, 3 rounds unless stated.
  - `build` of 12 controls: 102.0 s and 97.8 s (+4.3%).
  - `run` of one sample with the target-region report: 35.4 s and 34.2 s (+3.3%, 5 rounds).
  - `run` without a sex: 24.3 s and 26.3 s. `run` against the manuscript reference: 33.0 s and 34.3 s.
- Tests. `pytest tests/unit` gives 398 passed and 2 skipped with the compiled kernels, 379 passed and 21 skipped with `DMRDETECTOR_NO_JIT=1`, and 398 passed and 2 skipped with `DMRDETECTOR_NO_LIBDEFLATE=1`; five tests are new (the cache seed, and the two recipes). `tests/smoke/run_smoke.py --threads 8`: 95 of 95 checks. The privacy, documentation-command and container-recipe tests are part of the suite.
- Image of the release. Built with `apptainer build --fakeroot` from `git archive` of the final commit (168,255,488 bytes, 220 s including the kernel compilation); the checks below ran on it, in the default environment (kernels from the image, `NUMBA_CPU_NAME=generic`). A rebuild from the tag has the same files and another sha256.
- Comparison with 3.0.8 (the 3.0.8 image ran the same inputs). The md5 sums of all 161 data files of the real-data set (12-control `build` and 4 `run` calls, one of them against the manuscript reference) and of all 215 data files of the demo set (three BED formats, both sexes, target-region report with haplotype tracks, a RefSeq-contig run that stops) equal those of 3.0.8. The build logs and progress logs differ in time stamps, durations and the version; after removing those fields they are equal. The README Quick start, run verbatim in a fresh folder, gave the same md5 for the DMR table (7 DMRs), the target-region table, the atlas and the two panels, and the one WARNING that the README describes.
- Misses of the baked cache on the final image: 0 (12-control `build`: 193 entries loaded; each of the 4 `run` scenarios started from a freshly seeded folder: 12 to 21 loaded).

## [3.0.8] - 2026-10-08

Patch release after a methodology fidelity audit that compared 3.0.6 with v2.2 (git tag `v2.2-manuscript`) and with the manuscript text. The release restores one v2.2 behavior, makes every pandas reader exact, adds two warnings, one NOTICE and one record in `reference.json`, and corrects statements in the documents. The algorithm, the defaults and the options are unchanged. For valid input, the data files of `build` and `run` are byte-identical to those of 3.0.7 (checked, see `doc/validation.md`), with these exceptions: a samplesheet with a `type` column builds a different reference when it has rows other than `CTRL`; and a few `P_VALUE`, `P_FDR` and `DIFF` cells of the annotated DMR table and of the target-region table change in the last digits (at most 1.1e-12 relative), because the pandas parser of the annotation and of the target report is now exact. The atlases, the panels and the raw DMR tables are unchanged.

### Fixed
- `build` read every row of a samplesheet with a `type` column as a control (3.0.0 to 3.0.7, with a WARNING). v2.2 read only the `CTRL` rows. `build` now uses the rows whose type is `CTRL` or `control` (any case) and drops the others with a NOTICE that gives the number of rows dropped. A row with an empty or unknown type is dropped as well; v2.2 stopped with an error on a type other than `CASE` or `CTRL`. Dropped rows are not checked for files or duplicate IDs. `build` stops when no row is a control. `reference.json` lists the kept controls only, and `samplesheet_md5` stays the md5 of the whole file. A sheet without a `type` column behaves as before, and the 5-column samplesheets of the tests are unchanged.
- The pandas readers used the default pandas float parser. It is not correctly rounded, and it drops the last digits of some decimals: it reads `0.000535864718615665` as `0.0005358647186156`, and its `high` option does the same. The readers serve `DMRDETECTOR_NO_JIT=1`, files that the byte-level reader rejects, indexed files read with one thread, the target report, and the annotation step. Three effects followed. The annotated DMR tables of 3.0.7 had one or two digits fewer than the raw DMR tables in some `P_VALUE` and `P_FDR` cells. The text of `MEAN_REF`, `MEAN_DIFF`, `P_VALUE` and `P_FDR` from the pandas route differed from the text of the byte-level route (3.8e-12 relative on the synthetic data of the tests; 40% of the `P_VALUE` cells differed in text). A `build` under `DMRDETECTOR_NO_JIT=1` gave a smoothed atlas that differed from the compiled one in the demo data. Every pandas read of decimals in `src` now uses `float_precision="round_trip"`: the atlas, the raw atlas, the sample and atlas rows of `call-dmr`, the ICR atlas windows, the target-region tracks, the raw DMR tables and the methylation BEDs. A unit test fails when a new `pd.read_csv` of decimals lacks the option.
- `settings.toml` and `--help` said "at most" for the island CORE and SHORE thresholds. The code uses "below" (`cpg_islands.py`). The help of the DMR difference thresholds said "minimum" where the code needs a value above the threshold, and the help of `dmr_merge_distance` said "closer than" where a gap of exactly 150 bp merges. The wording follows the code; no behavior changed.
- `doc/target_region_methods.md` said that the target report requires more than 5 or more than 30 CpGs and that DMPK Region-A gets no `DIRECTION`. The code (`target.py`) applies no minimum CpG count; a region needs 2 CpGs for the t-test. The document now matches the code.

### Added
- `build` records the md5 of the imprinting gene list that it used in `reference.json` (`icr_gene_list_md5`). `run` computes the md5 of its list (the `icr_gene_list` setting or the bundled file) and prints `WARNING: the ICR gene list differs from the one used to build the reference (GENE labels may differ)` when the reference records another md5. A reference without the key (v2.2 layout, Zenodo record, `build` before 3.0.8) is not checked. The list gives the `GENE` labels of the ICR DMRs; the ICR panel is read from the reference. No output changes when the lists are equal.
- `preprocess` warns once per sample when a kept methylation record has `END - START` other than 1. DMRdetector selects the records of a panel region by `START`, and `tabix -R` of v2.2 selects them by overlap, so the calls near such records can differ. The check uses the arrays already in memory. A CpG record is one base long in the supported formats.
- `doc/ALGORITHM.md` lists three v2.2 behaviors that the manuscript text omits, with numbers from the manuscript run. The CpG at `END_EXTEND` of an extended island is never read: this holds for 25,206 of 25,719 islands, and 513 are covered by another query line. A CpG can enter a region m x m x k times; in the manuscript panels this is 8 for 2,323 CpGs and 27 for 226 ICR CpGs (the audit found up to 8 in the DMR tables). Repeated ICR DMR rows come from overlapping panel rows with different coordinates (57 extra rows in 52 of 75 tables), while rows with identical coordinates form one group with doubled counts.
- `doc/validation.md` states what was compared with which tolerance and environment, what was not compared, and has a placeholder section "Verification of 3.0.8" for the rebuild of the 69-control manuscript atlas.
- Unit tests for each item above, and smoke checks for the `type` column and for the gene-list record.

### Changed
- README, `doc/usage.md`, `doc/output.md`, `doc/ALGORITHM.md`, `doc/target_region_methods.md` and `doc/validation.md` describe the items above. The README states that the default of `--threads` is half of the available threads (at least 2); this agrees with `resolve_threads`.
- The text of the 3.0.5 entry now says that "byte-identical to 3.0.4" holds for the files that we compared in one environment.

### Release checks
The comparison with 3.0.7, the plain-Python route, the timings and the tests are in [doc/validation.md](doc/validation.md) (section "Release checks of 3.0.8 and 3.1.0").

## [3.0.7] - 2026-10-08

Patch release after an independent code review of 3.0.0 to 3.0.6. No result changed. For valid input the data files of `build` and `run` (atlases, panels, DMR tables, target-region tables) are byte-identical to those of 3.0.6; we compared the demo data (three BED formats, a female and a male case, all four scopes, with a target-region report) and real data (one `build` from 12 controls, and 4 `run` calls on 3 samples, one of them against the v2.2-layout manuscript reference) with md5. `reference.json` and the logs differ in the version, the time stamp, the file names (see Changed) and the new lines.

### Fixed
- `build` could delete a folder it did not create. It removed `<outdir>/work`, or `<tmp-dir>/dmrdetector-build`, after a success even when the user had put their own files there, and it wrote its progress log into such a folder. The work folder now carries a marker file `.dmrdetector-build` (version, outdir, process ID, time) and is locked with `flock`. An existing folder is used only when it holds the marker of a build of the same `--outdir`, or is empty; otherwise `build` stops with an error that names the folder and suggests `--tmp-dir`. The folder is removed only when its marker matches. Two builds cannot use one folder: the second stops. A file system that cannot lock gives a warning.
- With `--tmp-dir` the work folder is named `dmrdetector-build-<hash of the resolved outdir>`, so builds of different references share a `--tmp-dir` without sharing a folder.
- The warning that `build` overwrites a reference now also fires when only `ref_atlas_*` files are present, and not on the resume of an interrupted build.
- A failure of the target-region report no longer discards the finished DMR calls. `run` publishes `<id>.dmr.tsv` first. If the report then fails, the log and the screen say `WARNING: target-region report failed: ...`, the log status reads `FAILED (the DMR table was published; the target-region report failed)`, and the exit status is 1. A `<id>.target_region.tsv` of an earlier run stays in the folder, and the warning says it is not from this run.
- The haplotype tracks and the stand-alone `target` command skipped the checks that `run` makes on the combined BED. Every track is now checked for modkit rows of other modification codes or of both strands, and for records beyond the hg38 contig length. An unsorted track, or one with a contig in two blocks, stops with `cannot index the track` instead of `OSError: building of index failed`. The `target` command prints `ERROR:` and exits with 1 on these errors.
- The target-region report ignored `input_format`; it now uses the setting for the columns of every track. A `.tbi` that is older than its track is no longer trusted; the report reads such a track in memory or indexes it in the temporary folder. The same age check applies to the atlas.
- `ensure_tabix` named its link in the scratch folder by the base name only, so files with one name in different folders shared a link and an index. Each file now has a folder named by a hash of its path.
- `atlas` without `scratch` left a `dmrdetector_atlas_*` folder in the temporary directory; it is now removed.
- Stamps of finished steps include the program version, so a kept work folder of another release is not resumed.
- The memory budget is the smaller of 8 GB per thread and the memory the process may use (MemAvailable and the cgroup memory limit). A `--threads` value above the available threads is lowered, with a warning. A worker killed for lack of memory gives `a worker process was killed, probably out of memory; use fewer --threads` instead of a traceback. When the machine limits the budget, the `estimated peak memory` error says so instead of asking for more threads.
- The cgroup CPU quota now takes the smallest quota along the parent groups and reads the right group of a cgroup v1 hierarchy. Before, it used the first group of the file and looked below it.
- The Dockerfile and the Apptainer definition set `NUMBA_CACHE_DIR=/tmp/numba_cache`, one folder for all users of a machine; the cache cannot be written when another user owns it. They no longer set it. The program uses `dmrdetector-numba-<uid>` in the temporary directory (mode 0700, never a link or a folder of another user) when the variable is unset.
- `build` and `run` now remove the `.tmp<pid>.<name>` files that a killed process left in the output folder. A file goes only when no process with that ID exists and the file is at least one hour old. Nothing else is touched.

### Changed
- `<id>.log` and `build.log` show each path by its last component, and each input file also by its size, unless `log_input_paths = "full"`. This covers `[inputs]`, `[parameters]`, the command line, the settings file and the messages. The progress log in a kept work folder keeps full paths.
- `reference.json` stores the file name of the samplesheet and of the settings file instead of their absolute paths. The sample IDs and the sex of the controls stay, because they document the reference; `doc/output.md` says so. The new `log_input_paths` is not written to it.
- `tests/unit/test_privacy.py` also rejects study sample IDs, and its examples are built from parts that are not real IDs.
- The tests of `call_dmr` compare the text of the table from the byte-level reader and from the pandas reader. The text is equal except in `MEAN_REF`, `MEAN_DIFF`, `P_VALUE` and `P_FDR`. There the pandas reader is up to 1 ulp off per token, which gives 3.8e-12 relative in the p-values of the synthetic data. These four columns are compared with a relative tolerance of 1e-10 instead of 1e-6. The pandas route is only the fallback for `DMRDETECTOR_NO_JIT=1` and unusual files.

### Added
- Setting `log_input_paths` (`name` or `full`, default `name`) in `[paths]`.
- Flag `--input-format` of the `target` command.
- Warnings, none of which changes an output file:
  - no CpG is left after filtering a sample;
  - no CpG of a sample lies on chr1 to chr22 (for example RefSeq contig names), or a sex scope of the sample is empty;
  - the VCF names none of the BED contigs (`1` against `chr1`), or has several sample columns;
  - a sample labeled F has at least 50,000 CpGs on non-PAR chrY (`doc/sex_chromosomes.md` gives the rule and its limits; a male is not checked);
  - `run` replaces `<id>.dmr.tsv`, `<id>.log` or `<id>.target_region.tsv`, or leaves an old target-region report in place;
  - `--threads` or the memory budget was lowered.
- Error: `run` stops with `0 sample CpGs match the reference atlas` when no CpG of the sample matches the autosomal atlas. The message gives the counts and the likely causes. In a sex scope it is a warning. Before, such a run finished with exit status 0 and a DMR table that held no autosomal DMR (on the demo data with RefSeq contig names: the header and one chrX row).
- Unit tests for each item above.

## [3.0.6] - 2026-10-08

Patch release: the target-region report reads more input formats. The algorithm, the defaults and the numerical output are unchanged. For input that already worked, the target-region table and the other files of `run` are identical to those of 3.0.5.

### Fixed
- `run --target-region` stopped with `OSError: building of index for ... failed` when the combined or a haplotype BED was a pb-CpG-tools file (it starts with a `chrom` header line). The target report reads these files by region with tabix. `tabix -p bed` does not skip a `chrom` header line, and the in-memory reader used for BGZF files without an index rejected the line as well. The report now drops the leading header lines (`#`, `chrom`, blank) before it reads or indexes a track. Runs without `--target-region` were not affected.
- The target report also failed on other tracks that `build` and the rest of `run` accept: gzip that is not BGZF, uncompressed files, and modkit lines with spaces between the columns from the 10th on (also when the file has a tabix index). The report now rewrites such a track as a BGZF copy with tabs in the scratch folder and indexes the copy. The input file is not changed.
- Two unindexed tracks with the same file name in different folders (for example `hap1/S1.bed.gz` and `hap2/S1.bed.gz`) shared one index in the scratch folder, so the second was read from the first one's index. Each track now has a folder of its own.
- A methylation BED with only header lines or no lines no longer stops the target report; every region gets `NUM_CPG` 0.
- New unit tests of the target report: 5-column, modkit and pb-CpG-tools tracks (BGZF, gzip, plain), with and without haplotype BEDs, with mixed formats, on both read routes. The smoke test runs the same cases.

## [3.0.5] - 2026-10-07

Patch release: one crash fixed, one convenience value, clearer messages and documentation. The algorithm, the defaults and the numerical output are unchanged. For valid input, the files of `build` and `run` that we compared were byte-identical to those of 3.0.4. The comparison covered the demo data and 14 real samples and ran in one environment (one set of Python and package versions on one server). It does not show that the files are identical in another environment: a reference built twice from the same controls with other package versions differed in the 10th digit of `REF_METH` and `REF_VAR` (`doc/validation.md`). `reference.json` and the logs differ in the version, the time stamp and the new log lines.

### Added
- `--target-region builtin` (and `target_region = "builtin"` in `settings.toml`) selects the target BED packaged with the program. Inside the image the file exists only at a path that a user cannot pass from outside it. A path still works as before; a file named `builtin` in the working directory must be written `./builtin`.
- The first log lines of `build` and `run` give the settings file in use (`settings: ` and its absolute path, or `settings: built-in defaults`), and a NOTICE lists the keys of a file that differ from the packaged defaults.
- New warnings; nothing else changes. An empty methylation BED gives `no CpG records`. A `type` column in the `build` samplesheet is ignored, and every row is a control. `build` into a folder that already holds a reference overwrites its files; a resume of an interrupted `build` does not warn. A target-region report in which no region has a usable CpG gives a WARNING.
- `build --input /dev/null` (and any other input that is not a regular file) reports `not a regular file`.
- Tests for each of these, and `tests/unit/test_privacy.py` now also rejects the lab host and data set names, the scratch path and institutional addresses.

### Fixed
- `run --target-region BED` stopped with `ValueError: You are trying to merge on object and float64 columns for key 'CHROM'` and wrote no DMR table when no target region had a CpG that is also in the atlas (for example a BED for contigs that the sample lacks). The cause was an empty atlas frame with float64 columns. The report is now written with every region listed (`NUM_CPG` 0, `NA` in `P_FDR` and `DIRECTION`) and one WARNING line.
- `--set dmr_fdr=0.05` now suggests `dmr_fdr_threshold` (the nearest key was missed by the similarity cutoff).
- A wrong command name no longer lists the hidden stage commands.
- The help text of the hidden commands `atlas` and `preprocess` no longer names v2.2.

### Changed
- The `file not found` error for a relative path in a samplesheet adds that such a path counts from the working directory, when that is not the folder of the CSV. The behavior is unchanged.
- The output section of `doc/sex_chromosomes.md` describes the flat v3 output.
- `doc/output.md` and the README state the unit of the target-region `DIFF` (a fraction, not a percentage).
- The README and `settings.toml` separate the two values of 5: `cpg_min_samples_per_site` is a number of controls, `icr_min_sites` a number of CpGs.
- The README states that the default `--threads` is half of the available threads and should be set on shared servers. The Quick start uses `--target-region builtin`. The workaround row for the merge error is gone from `doc/usage.md`.
- Private identifiers (host name, data set name) are replaced by neutral wording in comments, tests, `.gitignore` and this file; the source of the DMPK regions is named as the authors' analysis.
- Docstrings for the public functions that lacked them; the ASCII codes in `fileio/methbed.py` have names (the generated code is the same); the unused `noqa` comment is removed.

## [3.0.4] - 2026-10-07

Patch release: documentation, repository hygiene and error wording for the first public release. The algorithm, the defaults, the options, the commands and the output files are unchanged (the unit and smoke tests pass, and the output of the README commands is identical for the three BED formats).

### Added
- `doc/usage.md` (mounting data, Docker, numba cache, a SLURM example, input checks, settings file, memory, and a table of error messages with fixes) and `doc/dependencies.md` (package licenses, a vulnerability scan, platform notes).
- `CITATION.cff`, `CONTRIBUTING.md`, `.gitattributes`.
- `tests/unit/test_docs_commands.py`: every `dmrdetector` command, flag and settings key in the Markdown files must exist, and every link and repository path must resolve.
- The README now states the supported platforms (Linux x86_64 and WSL2 only), builds without `--fakeroot`, what to do when `--fakeroot` is refused, and the minimum number of controls.

### Changed
- README rewritten: first screen, worked example on the synthetic demo data, supported inputs, a sample output table, a reference section (including the layout of the Zenodo record and how to build a reference from your own controls), measured run time and memory of `run`, a troubleshooting table and a FAQ. Statements about v2.2 times and removed files are dropped.
- `tests/smoke/generate_data.py` writes the demo data set of the README (8 controls, 2 cases, three BED formats, `controls.csv`, a demo `settings.toml`, a small GTF and a target-region BED) and reuses the data model of `run_smoke.py`. The old generator wrote the cases into the control samplesheet and an unusable reference.
- Error wording (no change of behavior or exit status): a missing file inside a container adds a hint to mount the folder; a truncated or invalid gzip input names the cause and `gzip -t`; `get-genes` prints one `ERROR:` line instead of a traceback when the network or the release fails; the empty-atlas and low-memory atlas errors name the fix.
- `doc/output.md` lists the `SAMPLE_TYPE` column of the target-region table, and the header of `settings.toml` no longer names `nextflow.config`.
- `doc/target_region_methods.md`, `doc/sex_chromosomes.md` and the target-region evidence table no longer cite lab scripts or paths; the earlier internal version is named as unpublished.
- `tests/unit/test_privacy.py` also rejects the lab's NAS mount names.

### Removed
- The stale statements of the README about the v2.2 timing table and a reference "extracted" from the Zenodo record.

### Image
- Two builds of the 3.0.3 recipe on one machine (one with and one without `--fakeroot`) gave files of the same size (167,186,432 bytes) and different sha256, so the image is not bit-reproducible. Record the sha256 of the file you publish.


## [3.0.3] - 2026-10-07

Patch release: internal refactor. There is no change in results, options or commands. The output files of `build` and `run` (all sex scopes and the target-region report) are byte-identical to those of 3.0.2 on the same inputs. The `--help` text and the error messages are identical.

### Changed
- Modules are renamed and merged. `io/` is now `fileio/` (`bed`, `bgzf`, `numparse`, `methbed`, `atlas_io`, `vcf`); `stages/` is now `steps/`; the code of `pipeline/` is in `build.py`, `run.py` and `runtime.py`; `config.py`, `defaults.py` and `inputcheck.py` are merged into `settings.py` and `runtime.py`. `cli.py` holds one command table and builds a sub-parser only for the command that is run.
- The functions of the pure-Python atlas path keep their behavior and have new names (`_build_atlas_tabix`, `_smooth_atlas_pool`). `DMRDETECTOR_NO_JIT=1` and the pandas and pysam fallbacks are unchanged.
- The code is formatted to the Google Python style guide (80 columns, Google docstrings with Args, Returns and Raises, absolute sorted imports, type annotations on public functions). The repository holds no formatter or linter configuration; the commands used are in `tests/README.md`.
- The deprecated `./dmrdetector.toml` is still read with the same warning; its removal is postponed.
- `doc/ALGORITHM.md` has a section on target-region statistics and implementation notes (fast readers, number parsing and formatting, fallbacks, output files).

### Removed
- Unused code: `bed.fetch_region_lines`, `tabix_contigs`, `read_bed5_lines`, `rstats.r_sd` and the class `panel_reader._Out`.
- The helper scripts `tests/unit/dump_tbi.py`, `gen_synth_io.py` and `synth_call_dmr.py`; their functions are in `tests/unit/helpers.py`. Duplicate tests are removed (274 collected tests became 257).

### Image
- The image of this release, built on a 96-core server with `apptainer build --fakeroot dmrdetector.sif env/dmrdetector.def` (Apptainer 1.4.5) from the committed tree, is 167,186,432 bytes with sha256 `edc2c0b7cdc974ecfa920d38b9be9f6cd9c03157e4a67bcdf0d69922de8c8666`. The package versions are those of 3.0.2 (Python 3.12.15, numpy 2.5.3, scipy 1.18.1, pandas 3.0.6, pysam 0.24.1, numba 0.68.0, llvmlite 0.50.0, libdeflate 1.23-2). The in-build test and the smoke test (91 of 91 checks) pass through the image.



## [3.0.2] - 2026-10-07

Patch release: one settings file, one environment file and no `pyproject.toml`. Algorithm, parameters, defaults and output files are unchanged.

### Changed
- The options file is now `settings.toml`. It lists every option that has a default, uncommented, with its default value and a one-line comment, in the sections of 3.0.1. Options without a default (`reference`, `gtf`, `target_region`, `tmp_dir`, `icr_gene_list`, `par_bed`, `threads`) stay as commented placeholders. A header states the precedence (flag > `--set` > file > default) and that deleting a line restores the default. `dmrdetector config init` writes `./settings.toml`; `build` and `run` read `./settings.toml` when it exists, or the file given by `--config`.
- A copy of the file is in the repository root. `config.default_toml()` is the only generator, and `tests/unit/test_run.py` fails when the root copy differs from its output. The packaged copy (`data/config/dmrdetector.toml`) is removed, and `data.config_file()` with it.
- `./dmrdetector.toml` (the name in 3.0.0 and 3.0.1) is still read in this release. `build`, `run` and `config show` print "dmrdetector.toml is deprecated; rename it to settings.toml" on stderr, and `build` and `run` also write it to the log. When both files exist, `settings.toml` is used and the warning says so. The fallback will be removed in the next release.
- `pyproject.toml` is removed. The image no longer builds a wheel: both recipes copy `src/` to `/opt/dmrdetector/src` (package data included), make it importable through `/opt/venv/lib/python3.12/site-packages/dmrdetector.pth`, and install the `dmrdetector` command as the launcher `env/launcher.py` in `/opt/venv/bin`.
- `env/` is reduced to the recipes and one environment file. `env/Dockerfile` and `env/dmrdetector.def` are written by hand. The generator and `env/build.py` are removed, with `render`, `stage`, `lock` and `clean`. `env/requirements.lock` is the only package list: exact versions and sha256 of every wheel for Python 3.12, with the five top-level packages and the regeneration command in its header. `requirements.in`, `requirements-test.in/.lock` and `build-requirements.in/.lock` are removed. The build is `apptainer build --fakeroot dmrdetector.sif env/dmrdetector.def` or `docker build -f env/Dockerfile -t dmrdetector:3.0.2 .`. pip comes from the digest-pinned base image; pytest is installed by the tester (`pip install pytest`) and is not in the image. `tests/unit/test_container_recipes.py` replaces `test_container_defs.py` and fails when the two recipes differ in base image digest, install command, libdeflate pin, environment, smoke test or version label.
- The log of every `build` and `run` has a `versions:` line with Python, numpy, scipy, pandas, pysam and numba.
- The Python version (3.12) is carried by the digest-pinned base image and the target in the header of `env/requirements.lock`. The `requires-python = ">=3.12"` of 3.0.1 lived in `pyproject.toml` and no longer exists.
- The version has one source, `src/dmrdetector/__init__.py`; both recipes, the README and the top CHANGELOG heading must match, and `tests/unit/test_container_recipes.py` checks this.
- pytest has no configuration file: `tests/conftest.py` puts `src/` on `sys.path`, and `pytest tests/unit` is the documented command (`tests/README.md`).

### Image
- The image of this release, built on a 96-core server with `apptainer build --fakeroot dmrdetector.sif env/dmrdetector.def` (Apptainer 1.4.5), is 167,182,336 bytes with sha256 `b95310a25d95c69bca06521b15181a2777a63019d28877bff76b050dc1f551bb`. It holds Python 3.12.15, numpy 2.5.3, scipy 1.18.1, pandas 3.0.6, pysam 0.24.1, numba 0.68.0, llvmlite 0.50.0 and libdeflate 1.23-2. A SIF embeds timestamps, so a rebuild has another hash and the same package versions.


## [3.0.1] - 2026-10-07

Patch release: startup time, input checks and documentation. Algorithm, parameters, defaults and output files are unchanged; the tables of the modkit, pb-CpG-tools and 5-column test files are identical to those of 3.0.0.

### Changed
- `dmrdetector --version`, `--help`, `run --help`, `build --help` and `config` start in about 1.2 s instead of 8.4 s (same server, wall time including the container start). The stage commands (`call-dmr`, `atlas`, `get-genes`, ...) build their parsers only when they are used, and `hashlib` is imported where it is needed, so the main parser imports no pandas, scipy or numba. `--help` output and the hidden stage commands are unchanged; `get-genes` still loads its module (and pandas) when it runs.
- `requires-python` is `>=3.12`, the version of the image and the only one the tests cover.
- The VCF is documented as it is read: only heterozygous genotypes are used, and phasing is not required. The `--vcf` help text says so.

### Added
- Genome-build check: preprocessing stops with `ERROR: ... hg38/GRCh38 ... hg19/GRCh37 ... not supported` when a record ends more than 1,000 bp beyond the hg38 length of its contig (chr1 to chr22, X, Y). It uses the rows already parsed, so it adds no pass over the file. It cannot prove that a file is hg38.
- modkit check: preprocessing stops when the first 20,000 data lines of a modkit bedMethyl file hold 5hmC or 6mA rows (or another code at the same position as an `m` row), or rows of both strands. The reader takes columns 10 and 12 of every row, so such rows count a CpG twice. The message says to run `modkit pileup --cpg --combine-strands` and keep the 5mC rows. Files made that way give the same tables as before.
- README: "Preparing input" (hg38 only, modkit, pb-CpG-tools, VCF), the numba cache in Docker and Apptainer, the memory floor of `run`, and the difference between `get-genes --release latest` and release 86.
- Tests for lazy imports (`tests/unit/test_cli_startup.py`) and the input checks (`tests/unit/test_inputcheck.py`).

## [3.0.0] - 2026-10-07

Python rewrite of the R/Nextflow pipeline (v2.2). The algorithm, parameters and output formats are those of v2.2; every stage was compared with v2.2 outputs on real data (`doc/validation.md`).

### Changed
- One Python package (numpy, scipy, pandas, pysam, numba) in one container image replaces Nextflow, R, data.table, bcftools, bedtools and conda. The image is built from `env/Dockerfile` (Docker) or `env/dmrdetector.def` (Apptainer, generated from the Dockerfile; a test checks that they stay equivalent). The base image is pinned by digest and every Python package by version and sha256 hash. Installation is by container only (Linux x86_64 and WSL2); macOS is not supported.
- `dmrdetector build` makes a reference from a samplesheet of control samples (`sample_id,sex,meth_bed,vcf`). `dmrdetector run` tests one sample per call, with the files as options (`--meth_bed_combined`, `--vcf`, `--sex`, optional `--meth_bed_hap1` and `--meth_bed_hap2`, `--reference`, `--gtf`, `--target-region FILE.bed`). `run` also reads v2.2 results directories as a reference.
- Output of `run` is `<sample_id>.dmr.tsv`, `<sample_id>.target_region.tsv` (with `--target-region`) and `<sample_id>.log`, directly in `--outdir` (default `dmrdetector/`). Intermediate files are removed.
- Flags cover only what changes per run. All other options are in one commented TOML file (`dmrdetector config init`) with the parameter names of v2.2, or are set for one run with `--set KEY=VALUE`.
- `--threads N` is the only resource option; memory is budgeted as 8 GB per thread (minimum 2 threads, default 50% of the available threads). Results do not depend on the thread count.
- Gene annotation reads an Ensembl GRCh38 GTF (`--gtf`, or `dmrdetector get-genes`) instead of `EnsDb.Hsapiens.v86`. With release 86 the annotated DMR tables of v2.2 are reproduced; 658 `LRG_gene` entries of EnsDb are not in the GTF.
- The shipped imprinting gene list (`src/dmrdetector/data/assets/icr_gene_list.bed`) is the 124-row list of the manuscript run. `assets/README.md` documents every bundled file.
- The raw atlas reads the inputs window by window through their tabix indexes and sizes the windows from the memory budget; the output is byte-identical for every budget. The FFT of each chromosome runs in a short-lived child process, which keeps the peak memory of smoothing 12 samples at `--threads 2` at 13.1 GB instead of 25.2 GB.
- The smoothed atlas is built at one bandwidth (σ = 150 bp). The σ = 500 bp atlas of v2.2 was not used by any later stage.
- The heterozygous-SNP filter types variants as `bcftools view -v snps` does (equal-length alleles that differ at one base count as SNPs).

### Added
- Target-region report and `src/dmrdetector/data/assets/target_region.hg38.bed`: 21 hg38 regions with disease and direction, an evidence table, and the methods in `doc/target_region_methods.md`. Columns include `DIFF`, `P_FDR` and `DIRECTION`; haplotype tracks are optional.
- `--set KEY=VALUE`, `dmrdetector config init` and `config show`, and `dmrdetector get-genes`.
- `doc/ALGORITHM.md`: every step with its parameters, and the behaviors of v2.2 that the manuscript text does not state.

### Performance
Measured on the 14-sample subset of the manuscript cohort at 24 to 36 threads: `build` 130 to 146 s (v2.2: 1,882 s), `run` of one sample 29 to 36 s (six cases in v2.2: 377 s), peak memory of `build` 115 to 120 GB at 36 threads (v2.2: 189 GB for one task). Numerical kernels are compiled with numba, BGZF blocks are read and written in parallel, and the VCF scan is a byte-level parser.

### Kept from v2.2 on purpose
Behaviors that the manuscript text does not describe but that produced the published numbers: CpG island DMRs are 1-bp anchors (`END = START + 1`) clustered on the island core start; a CpG inside overlapping panel rows is counted once per row (the options `dmr_dedupe_sites = true` and `icr_v1_duplicates = false` count it once); CpG islands need more than 5 controls per site while ICRs need at least 5.

### Removed
`main.nf`, `nextflow.config`, `nextflow_schema.json`, `conf/`, `environment.yml`, the R asset builder, the σ = 500 atlas, `install.sh`, and the `--force`, `--dry-run` and `--input-format` flags (set `input_format` in the options file).

## Earlier versions: R/Nextflow pipeline, 0.1.0 to 2.2

These versions were a Nextflow pipeline with R (data.table, Bioconductor). Version 2.2 produced the manuscript results; 3.0.0 reproduces its algorithm and parameters (git tag `v2.2-manuscript`). The method-relevant changes, newest first:

- 2.2 (2026-04-28): PAR1 and PAR2 are folded into the autosomal scope, because they are diploid in both sexes. Reads on chrY-PAR are translated to chrX-PAR coordinates (chrY-PAR2 by +98,813,480), so a locus counts once. The chrX and chrY scopes are strictly non-PAR.
- 2.2: `CALL_DMR` and `DETECT_ICR_REGIONS` read the atlas and the sample BED with tabix only at the panel regions (about 1% to 2% of hg38), which cut atlas I/O 50 to 200 times.
- 2.2: atlas files are headerless `.bed.gz` instead of `.tsv.gz`. The no-op `REF_METH` band filter in `CALL_DMR` (always 0 to 1) is removed.
- 2.1 (2026-04-28): sex is a required samplesheet column (M or F). Inference of sex from the VCF and its parameters (`sex_het_threshold_male`, `sex_het_threshold_female`, `sex_override_tsv`) are removed.
- 2.1: gene annotation uses `EnsDb.Hsapiens.v86` instead of a GTF (biotype priority, then overlap width and gene length).
- 2.1: the reference methylation range of ICRs widens from 0.40 to 0.60 to 0.35 to 0.65 (50% +/- 15%). The `PANEL_SCOPE` column is dropped from the DMR table.
- 2.1: scope-specific atlases and panels can be supplied to skip atlas building (`--ref_atlas_male_X`, `--ref_atlas_female_X`, `--ref_atlas_Y` and the matching `--cpg_panel_*`).
- 2.0 (2026-04-27): CpG island and ICR panels run by default; `icr_panel`, `icr_diff_threshold` and `icr_min_sites` replace `panel2`.
- 2.0: ICRs are hemi-methylated regions inside imprinting gene candidate intervals (`icr_gene_list`, `DETECT_ICR_REGIONS`). hg38 only. The bundled assets are `cpg_island.bed`, `icr_gene_list.bed`, `icr_region.bed` and `par_hg38.bed`.
- 2.0: female chrX uses a hemi-methylated CpG island model (`sex_f_ref_meth_min`, `sex_f_ref_meth_max`) with relaxed thresholds (`sex_f_diff_threshold`, `sex_f_min_sites`).
- 1.3 (2026-04-27): four scopes: autosomal, male_X, female_X and Y. PAR stays with chrX and chrY. CpG islands come from the male-derived scopes only; female chrX is tested against the female_X atlas with the male_X panel.
- 1.3: samples are routed by sex (male: autosomal, male_X, Y; female: autosomal, female_X; unknown sex: autosomal only). ICR analysis is autosomal only.
- 1.2 (2026-04-27): sex-stratified chrX and chrY analysis on by default (`include_sex_chromosomes`). Sex is inferred from the number of heterozygous SNPs on non-PAR chrX.
- 1.2: female chrX uses an imprinting-style rule on `REF_METH` (the X-inactivation band) with `sex_f_diff_threshold` 0.25 and `sex_f_min_sites` 5. The DMR table gets a `PANEL_SCOPE` column.
- 1.2: fixes: a phantom NA row for a scope without DMRs, and chromosome names are normalized (`1` and `chr1`).
- 1.1 (2026-04-21): release for the manuscript. `include_sex_chromosomes` option (default false). Samplesheet validation (`sample_id,type,meth_bed,vcf`; `type` CASE or CTRL; no duplicate IDs).
- 1.1: `DETECT_DMR` is split into `PARSE_GTF`, `CALL_DMR` (one-sample t-test, BH FDR) and `ANNOTATE_DMRS` (nearest-gene annotation).
- 1.0 (2026-03-12): first public release. DMR calling runs on CASE samples only. An Ensembl GRCh38 GTF is downloaded when `--gtf` is not given. `--tmp_dir` is added.
- 0.1.0 (2025-03-12): build the reference atlas, detect CpG islands, call DMRs. Input formats are detected (5-column BED, modkit bedMethyl with 18 columns, pb-CpG-tools with 9 columns).
- 0.1.0: FFT Gaussian smoothing at two bandwidths, empirical Bayes shrinkage of the dispersion, an optional second panel (`--panel2`), and prebuilt atlas and panel inputs (`--ref_atlas`, `--cpg_panel`).
