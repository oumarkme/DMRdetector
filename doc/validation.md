# Validation against v2.2

The comparison scripts need the lab's private sequencing data and the v2.2 reference outputs, so they are not part of this repository. This page lists what was compared, with which tolerance, in which environment, and what was not compared (October 2026).

## Reference outputs and environment

The R/Nextflow pipeline (v2.2, git tag `v2.2-manuscript`, code unchanged) ran in its container with R 4.4, data.table and `EnsDb.Hsapiens.v86` on three inputs:

- 14 whole-genome samples (12 controls, 2 patients) in atlas-building mode, all four scopes (autosomal, male chrX, female chrX, chrY).
- 6 patients against the 69-control manuscript atlas that v2.2 built, in skip mode (atlas supplied).
- The same two runs with the 124-row imprinting gene list of the manuscript run, which reproduces the manuscript's 446-line imprinted control region (ICR) panel.

The Python implementation ran on the same inputs on one Linux server, with 12 to 36 threads. Part of the runs used the image (Python 3.12.15, numpy 2.5.3, scipy 1.18.1, pandas 3.0.6, pysam 0.24.1, numba 0.68.0) and part a host environment. The environment of each comparison was not recorded. Gene annotation used the Ensembl release 86 gene transfer format (GTF) file. The comparisons were made on the code of 3.0.0 to 3.0.7 during development and were not repeated for every release; 3.0.8 has its own check at the end of this page. The CHANGELOG states which files were compared for each patch release.

## What was compared

A float cell is a cell whose text has a decimal point or an exponent. For tables, a float cell passes when the difference is at most 1e-10 times the larger value, and every other cell must be equal text. For the atlas, `REF_METH` and `REF_VAR` pass when |a - b| is at most 1e-12 + 1e-10 |b|. The maximum differences observed are given next to these criteria.

| Stage | Compared | Result |
|---|---|---|
| Preprocessing | 80 files (filtered, autosomal, chrX, chrY of 20 runs) | Identical bytes after decompression. |
| Raw and smoothed atlas | 12 controls only: 4 scopes, 27,968,997 autosomal sites | Position, `N` and `COV` identical. `PHI` and `RAW_METH` within 1.0e-14 relative, measured on the development code of 3.0.8 with its exact float parser (the earlier pandas parser gave up to 3.4e-14). `REF_METH` and `REF_VAR` pass the criterion above; about 1,000 autosomal cells near zero pass only through the absolute term. |
| CpG island panels (4 scopes) and ICR panel | 12-control build, and the ICR panel built from the manuscript atlas file | Byte-identical (15 of 15 panel files in the last comparison). |
| Annotated differentially methylated region (DMR) tables, 14-sample run and 6 patients | Same regions, genes, directions and `NUM_CPG` | Maximum relative difference 4.1e-12 (14-sample run) and 5.1e-13 (patients, 121-row gene list). |
| Annotated DMR tables of the manuscript run | 6 patients, 124-row gene list, v3 `run` against the files of `manuscript_run` | All non-float cells equal text; float cells within 6.7e-13 relative. The tables are not byte-identical: some `P_VALUE` and `P_FDR` cells differ in the last printed digits. Regions, genes and directions are identical. |
| Gene symbols | 296 symbols | 296 of 296 equal to `EnsDb.Hsapiens.v86` with the release 86 GTF. EnsDb has 658 `LRG_gene` entries that the GTF lacks; 39 of 3,948 intervals placed at their edges differ, none of the real DMRs. |
| Target-region report | 144 rows of four real samples against the earlier R implementation | Coverage, CpG counts, methylation, reference mean and SD equal. The `DIFF`, `P_FDR` and `DIRECTION` columns are new and have no reference. |
| Thread count | 2, 4 and 36 threads | Atlases, indexes, panels and tables identical; compressed raw tables differ only in the gzip header time stamp. |
| Container and host | Two builds of one reference from the same controls: 3.0.0 under a host Python 3.13 environment, 3.0.6 through the image | Same sites, same panels. `N`, `COV`, `PHI` and `RAW_METH` identical; `REF_METH` and `REF_VAR` differ by up to 1.7e-10 and 1.3e-10 relative where the value exceeds 0.001. The cause (Python and package versions, or code between the two builds) was not isolated. Two 3.0.6 builds through the image (30 and 36 threads) were identical. |
| `NUMBA_CPU_NAME` (3.0.8 image) | The same image with the numba kernels compiled for the host CPU (Ice Lake server), for Haswell (AVX2, FMA), for Skylake-AVX512, and for a generic x86-64 CPU (`NUMBA_CPU_NAME=generic`); the 12-control `build` with 4 `run` calls, and the demo data (three BED formats, both sexes, target-region report) | Every data file has the same md5 in the four settings (161 files in the real-data set, 215 in the demo set); the logs differ in time stamps and durations. The settings were forced by name and feature list on one CPU model, and numpy dispatched to its AVX512F code in all of them. This shows that the compiled kernels do not depend on the CPU name; it does not show that the output is the same on other CPUs (next row). |
| numpy SIMD dispatch, 3.0.10 (see "CPU independence" below) | The 3.0.10 image with `NPY_DISABLE_CPU_FEATURES` set so that numpy 2.5.3 does not use its AVX512 loops; 12-control `build` with 4 `run` calls, and the demo data | Not identical: `np.exp` in `gaussian_kernel` returned other values for 44 of 901 entries (sigma 150) and 22 of the 161 data files differed. 3.1.0 removed this cause. |

## What was not compared

- Before 3.0.8, v3 had not rebuilt the 69-control manuscript atlas: the atlas comparison in the table covers 12 controls, and the v3 `run` results for the six patients use the atlas files that v2.2 built. The 3.0.8 check at the end of this page rebuilds the manuscript atlas and runs all 75 samples.
- The image and a host installation give the same files only when their Python and package versions agree. The row "Container and host" of the table shows a difference in the 10th digit of `REF_METH` and `REF_VAR` between two environments.
- Docker: the recipe is written and its commands were inspected, but no image was built on a Docker host.
- Windows Subsystem for Linux (WSL).
- A real host with another CPU. Every setting of the section "CPU independence" was forced on one Ice Lake server (Intel Xeon Gold 6342). We had no AMD machine and no machine without AVX2. The switches emulate the selections that numpy, glibc, OpenBLAS and numba make from the CPU; they do not reproduce a CPU, and a library that selects code by another route would not be affected by them. The bytes of the compressed files (the gzip stream of the atlases and raw tables) were not compared, because the gate compares the decompressed content; libdeflate and zlib choose their instructions at run time.
- The native (non-container) installation: the glibc setting is part of the image only.
- A 130-sample atlas at the minimum memory budget.
- Timings on a quiet machine.

## Verification of 3.0.8

We rebuilt the manuscript atlas with 3.0.8 (image sha256 `592eaf68...`, built-in defaults, 36 threads) from the control rows of the manuscript samplesheet, with the sex labels of that sheet. The build took 285 s, and the peak resident set size summed over the process tree was 111 GB. The result was compared with the manuscript atlas files, which have the md5 sums of the Zenodo record 10.5281/zenodo.18978867. The panels in the record have the md5 sums of the panels rebuilt here.

| Item | Result |
|---|---|
| Atlas sites (autosomal 28,684,487; male_X 1,274,830; female_X 1,288,169; Y 173,591) | `CHROM`, `START`, `END`, `N` and `COV` identical in every row |
| `REF_METH`, `REF_VAR` | All cells pass the criterion above. Largest absolute differences in the autosomal scope: 8.7e-12 and 1.5e-12. Cells near zero pass only through the absolute term |
| `PHI`, `RAW_METH` | Equal text in 95 to 98% of the cells of every scope. The other cells differ in the 15th significant digit, the last one printed: largest absolute difference 1.1e-15, largest relative difference 1.01e-14. `RAW_METH` differs by at most 1.1 units of that digit and `PHI` by at most 3.1 units (9 autosomal cells above 2 units). 213 of 62.8 million cells exceed 1e-14 relative; none exceeds 1e-10 |
| CpG island panels (four scopes: 25,719, 942, 2,980 and 27 islands) and ICR panel (445 regions) | Byte-identical (5 of 5 files) |

The shipped imprinting gene list (124 rows, md5 `1b9166e5...`) gives the manuscript ICR panel, so we did not build with another list. The list that the manuscript run used is not stored with its files; the identical panel is the evidence that it is the shipped one.

`run` was executed for all 75 samples of the manuscript run (controls and patients, samplesheet sex labels, Ensembl release 86 GTF, no target region), once against the manuscript atlas and once against the rebuilt atlas. In both runs the 2,904 rows of the annotated tables have text-identical `CHROM`, `START`, `END`, `GENE`, `PANEL`, `NUM_CPG`, `NUM_READS`, `SAMPLE_METH (%)`, `REF_METH (%)`, `DIFF (%)` and `DIRECTION`. No table is byte-identical, because `P_VALUE` and `P_FDR` differ in the last printed digits. Against the manuscript atlas all 75 tables pass at 1e-10 relative (maximum 1.3e-13). Against the rebuilt atlas 72 tables pass (maximum 8.5e-11). In the other 3 tables one row differs in `P_VALUE` and `P_FDR` by 1.0e-10 to 2.4e-10 relative; these p-values are below 1e-54, and the rows and their direction are the same. Against the manuscript atlas, the same code and samples give at most 1.3e-13 for every table, so the difference comes from the rebuilt atlas. A likely cause is that a p-value this far in the tail amplifies the difference of one unit in the 15th digit of the atlas. We did not test this explanation.

The number of DMR rows per scope (autosomal, male_X, female_X, Y) and sample group equals that of the manuscript run, including the chrX rows of the controls and the absence of chrY rows. Two controls have a samplesheet sex that disagrees with their heterozygous-SNP count on chrX outside the pseudoautosomal regions; 3.0.8 flags one of them in `build` and `run`. We used the sheet labels, as the manuscript run did, so both controls are in the wrong sex scope of the atlas. Run with the sex corrected (atlas unchanged), these two controls have 36 and 0 chrX rows instead of 57 and 170. The atlases were not rebuilt with corrected labels.

Not verified in this check:

- the raw atlas and the preprocessed files of the manuscript run, which are not stored, so only the smoothed atlas was compared;
- the cause of the three p-value differences above 1e-10;
- the scope of each DMR row, which we derived from the coordinates and the sheet sex and did not read from the per-scope tables of `run`;
- the plain-Python route, CPU models other than the one used, and timings on an idle machine.

## CPU independence

Goal of 3.1.0: the same version writes the same data files on every x86-64 CPU. We could not test other CPUs. We forced the run-time CPU selections of numpy, glibc, OpenBLAS and numba on one Ice Lake server (Intel Xeon Gold 6342, 2026-10-09; numpy 2.5.3 finds the dispatch targets X86_V3, X86_V4 and AVX512_ICL there, and its baseline is X86_V2) and compared the data files with those of the unchanged server. The tests used the 12-control `build` with 4 `run` calls (161 data files) and the demo data (215 data files); the md5 sums of the decompressed content were compared. The section "Floating-point arithmetic and the CPU" of [dependencies.md](dependencies.md) lists the operations that we audited.

Kernel. For sigma 20, 40, 60, 100, 150, 300 and 500, every entry of the Gaussian kernel of 3.1.0 (7,027 entries in all) equals the value of an independent integer-arithmetic reference (Taylor series at scale 2^-300), the value of `math.exp` of glibc 2.41, and the value of `exp` in R 4.4.3 (glibc 2.35), bit for bit. `np.exp` of numpy 2.5.3 on the Ice Lake server differs from them in 12, 22, 22, 32, 44, 82 and 154 entries (in the order of the sigma values). For sigma 1000 and 2000, the reference differs from glibc 2.41 and R in 4 of 6,001 and 4 of 12,001 entries (two values of |x| each, one unit in the last place); the correctly rounded value is used.

| Setting (environment variable) | CPU selection that it emulates | 3.0.10: data files that differ from its default run (12 controls of 161; demo of 215) | 3.1.0: same comparison |
|---|---|---|---|
| `NPY_DISABLE_CPU_FEATURES="AVX512_ICL X86_V4"` | CPU with AVX2 and FMA3 but no AVX512 (numpy uses its X86_V3 loops) | 22; 33 | 0; 0 |
| `NPY_DISABLE_CPU_FEATURES="AVX512_ICL X86_V4 X86_V3"` | CPU without AVX2 and FMA3 (numpy uses its baseline loops) | 22; 27 | 0; 0 |
| `OPENBLAS_CORETYPE=Prescott` (3.1.0 also `Haswell` and `SkylakeX`) | OpenBLAS kernels of other CPU generations; the program does not call BLAS | 0; 0 | 0; 0 (all three) |
| `GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA,-FMA4,-AVX512F` | glibc without FMA, AVX2 and AVX512F: plain variants of `exp`, `log`, `pow`, `sin`, `cos` | 22; 20 | 0; 0 (3.1.0 sets the tunable itself, see below) |
| `NUMBA_CPU_NAME=haswell` and `skylake-avx512`, each with a feature list | numba kernels compiled with AVX2 and FMA, or with AVX512 | 0; 0 | 0; 0 |
| X86_V3 off, `NUMBA_CPU_NAME=haswell` and `OPENBLAS_CORETYPE=Prescott` together | CPU without AVX2 and FMA3, kernels for another CPU, other BLAS | not run | 0; 0 |
| 3.1.0 without `GLIBC_TUNABLES` | host with the FMA variants of glibc, as in an image without the setting | not applicable | 22; 4 |

Each setting was checked with a positive control: probes of `np.exp`, `np.log` and `np.log10` (numpy settings), of `math.exp`, `log`, `pow`, `sin`, `cos`, `atan` and `lgamma` (glibc setting) and of the numpy complex product (X86_V3 setting) that change under the setting; `OPENBLAS_VERBOSE=2`, which printed the core `Haswell`, `Katmai` (for Prescott) and `SkylakeX`; and the assembly of a test function, which has `vfmadd` and `ymm` registers under Haswell. The tunable names with the suffix `_Usable` (as in glibc before 2.33) have no effect in glibc 2.41. `-FMA` alone and `-AVX2` alone each select the plain variants.

What changed between 3.0.10 and 3.1.0, in the order in which the screen found it:

1. The kernel (`np.exp`). With AVX512 off, 44 of 901 entries (sigma 150) differed by one unit in the last place.
2. The product of the spectra (`fv *= kf`). With X86_V3 off, numpy rounds the complex product another way (separate instead of fused multiplication and addition), and the smoothed atlas differed although the kernel and the FFT were equal.
3. The twiddle factors of `scipy.fft` and the distribution functions of `scipy.stats`. With the plain glibc variants, the forward FFT of 1,000,000 and of 3,000,000 values and `t.sf` for 5 of the 7 degrees of freedom tried (2, 5, 12, 37, 74) differ in the last bit. Their effect on the smoothed atlas is as large as that of the kernel (up to 2.9e-12 in `REF_METH`). We found no cheap way to compute these functions ourselves, so the image sets `GLIBC_TUNABLES=glibc.cpu.hwcaps=-FMA,-FMA4`. `env/smoke_test.py` checks that the setting takes effect, with two arguments for which the FMA and the plain `exp` and `sin` differ.

Difference between 3.0.10 and 3.1.0 on the 69 manuscript controls (both rebuilt from the sample sheet, default settings): sites, `N`, `COV`, `PHI` and `RAW_METH` are identical in all four scopes. `REF_METH` differs in 46.9% of the autosomal cells (largest absolute difference 6.8e-12), 39.2% (chrX male, 7.3e-14), 35.8% (chrX female, 1.3e-13) and 90.9% (chrY, 4.3e-12); `REF_VAR` differs in 63.4% (1.8e-12), 61.0%, 46.3% and 87.8% (6.4e-13). All cells pass the atlas criterion.

Difference between 3.0.10 and 3.1.0 on the 12-control set, with default settings (`build` and four `run` calls, 161 data files): 139 files are identical in content, among them the five panels. The other 22 are the 4 smoothed atlases, their indexes, and the annotated and target-region tables with their temporary copies. In the atlases, positions, `N`, `COV`, `PHI` and `RAW_METH` are identical. `REF_METH` differs in 46.4% of the autosomal cells (largest absolute difference 1.7e-12) and `REF_VAR` in 62.7% (7.5e-13); in the sex scopes the largest absolute differences are 9.7e-13 (chrY) and below 2e-14 (chrX). All cells pass the atlas criterion above. In the tables, the non-float cells are equal text, and the float cells (`P_VALUE`, `P_FDR`, `DIFF`) differ by at most 6.2e-12 relative. The demo set shows the same pattern (largest relative difference 6.3e-14).

Not shown: that a CPU of another vendor behaves like these settings; that the compressed bytes of the output files are equal; and the effect of a user setting `GLIBC_TUNABLES` on purpose. A variable exported on the host is replaced by the image's value. A value given with `--env` or `APPTAINERENV_GLIBC_TUNABLES` (Apptainer) or `docker run -e` replaces the setting of the image, and the output then depends on the CPU again (22 of 161 data files differ in the FMA case).


## Verification of 3.1.0

We repeated the 3.0.8 check with 3.1.0 (image built from the release commit with the kernel cache required, built-in defaults, 36 threads) on the same inputs: the 69 control rows of the manuscript samplesheet, with its sex labels, and the 75 samples of the manuscript run. The build took 262 s (3.0.8: 285 s), and the peak resident set size summed over the process tree was 108.6 GB (3.0.8: 111.0 GB). The tools are those of the 3.0.8 check, with a count of text-identical float cells added.

| Item | 3.1.0 | 3.0.8 |
|---|---|---|
| Atlas sites, `CHROM`, `START`, `END`, `N`, `COV` | identical to the manuscript atlas in every row of the four scopes | same |
| Panels (four CpG island files and the ICR file) | byte-identical (5 of 5) | same |
| `REF_METH`, autosomal: cells equal as text; largest absolute difference | 11,669,449 of 28,684,487 (40.7%); 1.25e-11 | 11,676,625 (40.7%); 8.7e-12 |
| `REF_VAR`, autosomal: cells equal as text; largest absolute difference | 6,802,692 (23.7%); 2.68e-12 | 6,803,534 (23.7%); 1.5e-12 |
| `REF_METH`, `REF_VAR`, chrX scopes: largest absolute difference | 9.1e-14 and 1.9e-14 | 1.7e-13 and 3.3e-14 |
| `REF_METH`, `REF_VAR`, chrY: largest absolute difference | 5.2e-12 and 4.1e-13 | 4.4e-12 and 9.5e-13 |
| Cells outside \|a - b\| <= 1e-12 + 1e-10 \|b\| | 0 | 0 |
| `PHI`, `RAW_METH` | the same files as 3.0.8: 213 of 62.8 million cells exceed 1e-14 relative (largest 1.01e-14); none exceeds 1e-10 | same |
| `run`, 75 samples, manuscript atlas: tables with only float differences; largest relative difference; float cells equal as text | 75 of 75; 1.33e-13; 10,501 of 13,572 (the differences are in 1,536 `P_VALUE` and 1,535 `P_FDR` cells) | 75 of 75; 1.33e-13; same counts |
| `run`, 75 samples, rebuilt atlas: tables with only float differences at 1e-10; float cells equal as text | 71 of 75; 8,472 of 13,572 (2,547 `P_VALUE` and 2,553 `P_FDR` cells differ) | 72 of 75; 2,462 and 2,478 differing cells in the 72 |
| Same, tables above 1e-10 relative | 4: 1.1e-10 (p 2e-55), 1.8e-10 (4e-58), 1.0e-10 (2e-55) and 1.2e-9 (2e-218) | 3, largest 2.4e-10 |

In all 75 tables of both runs, the 2,904 rows have text-identical `CHROM`, `START`, `END`, `GENE`, `PANEL`, `NUM_CPG`, `NUM_READS`, `SAMPLE_METH (%)`, `REF_METH (%)`, `DIFF (%)` and `DIRECTION`. Only `P_VALUE` and `P_FDR` differ, in the last printed digits.

The rebuilt atlas does not agree more closely with the manuscript atlas than the 3.0.8 atlas did. The counts of equal cells and the largest differences changed by amounts of either sign, `PHI` and `RAW_METH` did not change, and one table left the list of tables above 1e-10 while two entered it. We conclude that the Gaussian kernel was not the main source of the difference from R: R computes the convolution with its own FFT (Singleton's mixed-radix algorithm), `scipy.fft` uses pocketfft, and each rounds differently at the 1e-12 level. The difference is amplified in the p-value of a region whose p-value is 1e-55 or smaller, because the spread of the per-site differences is tiny; we did not test this explanation. With the manuscript atlas as reference, the same code gives the 3.0.8 numbers, so the t-test code is not the source.

Not verified in this check: the cause of the four tail p-value differences; the raw atlas and the preprocessed files of the manuscript run (not stored); and a run of the plain-Python route on the 69 controls.

## Release checks of 3.0.8 and 3.1.0

The CHANGELOG used to hold these checks. The checks are kept here with the other validation results; the checks of 3.0.9 and 3.0.10 are in the CHANGELOG. Atlas and `run` results against the manuscript run are in the sections "Verification of 3.0.8" and "Verification of 3.1.0".

### 3.0.8

The same image ran 3.0.7 and the new `src` on two data sets. The demo set has three BED formats, a female and a male case, all four scopes, and a target-region report with haplotype tracks. The real-data set is one `build` from 12 controls and 4 `run` calls on 3 samples, one of them against the v2.2-layout manuscript reference.

- The md5 of every data file equals that of 3.0.7, except the annotated DMR tables of the four real `run` calls and the target-region tables of two of them. In the annotated tables 2 to 17 cells differ (`P_VALUE`, `P_FDR`; at most 7.8e-13 relative); in the target-region tables 2 to 4 cells (`P_FDR`, `DIFF`; at most 1.1e-12 relative). The cause is the exact float parser that 3.0.8 uses in the pandas readers. The raw DMR tables, the atlases with their indexes, the panels and all demo outputs are identical. `reference.json` differs in the version, the time stamp and the new key `icr_gene_list_md5`.
- For the sample that ran against the manuscript reference, 22 cells differ in text from the table of the manuscript run (25 with 3.0.7); the largest relative difference is 8.3e-14.
- Plain-Python route (`DMRDETECTOR_NO_JIT=1`): the demo `run` gives the same DMR and target-region tables as the compiled run (equal text), and the demo `build` gives atlases and panels that are equal after decompression (those of 3.0.7 differed). The `.tbi` files differ in both versions, because the plain-Python writer lays out the BGZF blocks differently. On real data we ran the stages `call-dmr` (autosomal and sex scopes), `annotate` and `target` on the preprocessed files of the compiled run for 2 samples (one male, one female); every table is equal text. A complete `run` under `DMRDETECTOR_NO_JIT=1` on real data did not finish within 2 hours, so the preprocessing of real data is not compared between the two routes.
- Speed, medians of 3 alternating runs after a warm-up on a shared server, 3.0.7 first and 3.0.8 second: demo `build` 63.9 s and 61.9 s; demo `run` 3.9 s and 4.1 s; `run` of one real sample 31.8 s and 32.1 s, and 39.8 s and 39.8 s with the target-region report and haplotype tracks; `build` of 12 controls at 36 threads 102.1 s and 106.0 s (+3.8%, the largest change).
- Tests: `pytest tests/unit` gives 393 passed and 2 skipped with the compiled kernels, 374 and 21 with `DMRDETECTOR_NO_JIT=1` (the tests of the compiled paths are skipped), and 393 and 2 with `DMRDETECTOR_NO_LIBDEFLATE=1`. `tests/smoke/run_smoke.py --threads 8`: 95 of 95 checks. The privacy, documentation-command and container-recipe tests (45 tests) pass.
- Not checked: a complete plain-Python `run` on real data, and the warning for records that are not one base long on a real modkit file (it did not fire on pb-CpG-tools or MethBat files).

### 3.1.0

- Image of the release commit, built on a Linux x86-64 server with `REQUIRE_CACHE=1`: `%test` passed with 18 kernels loaded from the baked cache in 7 modules and none compiled; `--version` prints 3.1.0 and `GLIBC_TUNABLES` is set.
- Unit tests on the final tree: 436 passed and 2 skipped with the compiled kernels, 417 and 21 with `DMRDETECTOR_NO_JIT=1`, 436 and 2 with `DMRDETECTOR_NO_LIBDEFLATE=1`; smoke test 96 of 96 checks. An earlier run of the compiled suite had one failure, `test_smooth_chromosome_matches_scipy_reference`, whose reference used the fused complex product of numpy; it now uses separate operations like the program.
- Threads: `build` and the four `run` calls with 8 and with 36 threads give the same 161 data files, and the demo set the same 215 (the logs differ in the thread count).
- README Quick start, run verbatim in a fresh folder: `build` 70 s, `run` 4.2 s. Panels and target-region table are identical to those of 3.0.10; the DMR table differs only in the last digits of one `P_FDR` (5.7500391984875e-292 against 5.75003919848776e-292).
- `REF_VAR`, 3.0.10 against 3.1.0 on the 69 manuscript controls: the largest absolute difference in chrX is 2.5e-14.
- Not verified, in addition to the list under "What was not compared": sigma values above 500 on real data.

## Verification of 3.2.0

Version 3.2.0 adds `pyproject.toml` (the pip route) and compacts the source and the tests. Algorithm, options and output files are meant to stay as they were in 3.1.0, so the checks compare files with those of 3.1.0.

### Release image and tests (commit of the tag)

- Image built from `git archive` of the release commit with `apptainer build --fakeroot --build-arg REQUIRE_CACHE=1`: 168,271,872 bytes, 519 s on a loaded server, sha256 `f4a6d27b3c0232be5f6f66a09c0dabf15a2a8d42a25ccdd4d51494d21622bd5b`. `--version` prints `dmrdetector 3.2.0`. The `%test` section loaded 18 kernels from the baked cache in 7 modules and compiled none.
- Unit tests through this image: 438 passed and 2 skipped with the compiled kernels, 420 and 20 with `DMRDETECTOR_NO_JIT=1`, 438 and 2 with `DMRDETECTOR_NO_LIBDEFLATE=1`. The 2 skips are the privacy tests, which need git; with git on the checkout the privacy, documentation-command and container-recipe tests give 49 passed. Smoke test: 96 of 96 checks.

### Output comparison with 3.1.0

The commit of the tag (3.2.0 image) and the tag `v3.1.0` (3.1.0 image) ran with at most 36 threads. The e2e set is `build` on 12 controls and four `run` scenarios on real data (238 files). The demo set is the README demo with three runs (232 files). Each set was also run with a second thread count and repeated three times warm.

The acceptance gate found no difference in the data files of either set. It exits with status 1 only because the image hash and the version string differ, as expected. An earlier run of the same gate on the tree after the compaction (3.1.0 image) reported ALL EQUAL. The gate compares the md5 of the decompressed content of each data file, and the text files that hold times, paths and process ids after those fields are removed. In the two sets, 162 and 169 files are identical in content to those of 3.1.0, and 17 and 15 text files are equal after the folders and times are normalized. The other files are resume stamps, whose names agree. All 45 timed jobs finished with status 0, and the `--help` text of every command is identical.

Time to run on a shared server, 3.2.0 divided by 3.1.0 (median of three warm repeats, load average 5 to 19): `build` of the e2e set 0.97, `run` scenarios 0.97 to 1.09, demo `build` 0.99 and demo `run` 0.97 to 1.00. In an interleaved comparison (3.1.0, 3.2.0, three pairs) the e2e `build` took 107.9 and 108.8 s (ratio 1.01) and the single-thread e2e `run` 35.6 and 39.5 s (ratio 1.11, higher in all three pairs). The same code run twice gave ratios up to 1.08 for that `run` scenario, so a difference of up to 10% is not separated from the noise of the server. No edited line runs once per site or per window. A timing on a quiet server has not been made.

### Installation with pip

These checks ran on a copy of the tree before the version change, so the package still carried the version 3.1.0.

- Install into a new virtual environment with Python 3.12.13 and an empty pip cache. pip resolved numpy 2.5.4 (the lock has 2.5.3) and the locked versions of scipy, pandas, pysam, numba, llvmlite, python-dateutil and six. It built `dmrdetector-3.1.0-py3-none-any.whl` (173,966 bytes) with setuptools 84.0.0 in an isolated build environment. The install took 247 s, almost all of it the download of 131 MB of wheels at 0.5 to 1.6 MB/s. With a warm pip cache, `pip install git+file://<copy of the tree>` took 17 s. Installing the lock with hashes and then the package with `--no-deps` took 38 s plus 3 s and gave the locked versions exactly. An editable install worked as well.
- In the non-editable, editable, locked and git-URL environments, `dmrdetector --version` printed `dmrdetector 3.1.0`, `python -m dmrdetector --version` agreed, `run --help` and `get-genes --help` worked, and `dmrdetector.data.asset()` found the four data files (`par_hg38.bed`, `icr_gene_list.bed`, `target_region.hg38.bed`, `target_region.sources.tsv`). `--target-region builtin` resolved the packaged BED in the demo run.
- README Quick start in the non-editable, the editable and the locked environment (8 threads, empty kernel cache that `build` filled): `build`, `run --target-region builtin` and `run --target-region target_region.bed` finished with exit status 0. The same demo run through the image of 3.1.0 (Python 3.12.15, numpy 2.5.3) is the comparison. Of 62 compared files, 61 were identical in content, and only `created` in `reference.json` differed. The autosomal atlas and its index were also identical byte for byte. The result was the same with and without `GLIBC_TUNABLES=glibc.cpu.hwcaps=-FMA,-FMA4` and `NUMBA_CPU_NAME=generic` exported. The host was one Intel Xeon (Ice Lake) with glibc 2.35, FMA, AVX2 and AVX512F.
- Resolution only (`pip install --dry-run`) with Python 3.13.12 gave the same versions. Python 3.13 was not run: no unit test, smoke test or output comparison. `requires-python` is `>=3.12,<3.13` for that reason.

Not verified: an AMD CPU, a CPU without AVX2, and a host without libdeflate (the program then uses zlib). The 215-file demo set was not run in a pip environment. Installation from a remote URL was tested with a `file://` URL. macOS and ARM Linux stay unsupported.
