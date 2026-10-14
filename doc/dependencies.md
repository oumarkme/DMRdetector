# Dependencies, licenses and platform

DMRdetector is released under the MIT license. The image installs the packages below from `env/requirements.lock` (exact versions and wheel sha256 hashes, Python 3.12, Linux x86_64). Licenses were read from the package metadata inside the image built from release 3.0.4.

## Python packages

| Package | Version | License (metadata) | Role |
|---|---|---|---|
| numpy | 2.5.3 | BSD-3-Clause and others (bundled OpenBLAS, libgfortran) | arrays |
| scipy | 1.18.1 | BSD-3-Clause (bundled OpenBLAS, libgfortran) | FFT, t distribution |
| pandas | 3.0.6 | BSD-3-Clause | tables |
| pysam | 0.24.1 | MIT (bundles htslib, bzip2, xz) | VCF and tabix access |
| numba | 0.68.0 | BSD-2-Clause | compiled kernels |
| llvmlite | 0.50.0 | BSD-2-Clause and Apache-2.0 with LLVM exception | numba backend |
| python-dateutil | 2.9.0.post0 | BSD-3-Clause or Apache-2.0 | pandas dependency |
| six | 1.17.0 | MIT | python-dateutil dependency |

Other components of the image: Python 3.12.15 (PSF license), the Debian trixie base layer (GPL, LGPL and other licenses of the Debian packages), and libdeflate 1.23-2 (MIT).

## License compatibility

All Python packages use permissive licenses (BSD, MIT, Apache-2.0) that are compatible with MIT source code. The numpy and scipy wheels bundle libgfortran and libquadmath (GPL-3.0 with the GCC Runtime Library Exception, and LGPL-2.1 or later). The libraries are dynamically linked and unmodified. The Debian base layer carries GPL packages. We distribute no modified copy of any of these. A person who redistributes the built image must follow the licenses of its components; the MIT license of this repository covers only its own source. We did not audit every bundled shared library one by one.

## Known vulnerabilities

`pip-audit 2.10.1 -r env/requirements.lock --require-hashes --disable-pip` reported "No known vulnerabilities found" on 2026-10-07 (PyPI advisory data of that day). The scan covers the Python packages only, not the Debian base layer. Re-run it before each release.

## Floating-point arithmetic and the CPU

Two installations of one release should write the same output files, whatever the CPU. We checked every floating-point operation on the result path of release 3.1.0 (October 2026). The result of an operation is fixed when the IEEE 754 standard defines it: addition, subtraction, multiplication, division, square root, comparison, floor and rounding to an integer. It is also fixed when our own code fixes the order of the operations. It depends on the host when a library chooses between code paths at run time from the CPU and the paths round differently. The grep covered `np.exp`, `np.log`, `np.power` and `**` with float exponents, trigonometric functions, `math`, `np.dot`, `@`, `np.linalg`, `np.sum`, `np.mean`, `np.var`, `np.cumsum`, `np.convolve`, `scipy.special`, `scipy.stats`, `scipy.fft`, the numba kernels and the pandas aggregations.

| Place (3.1.0) | Operation | Class | Handling |
|---|---|---|---|
| `steps/atlas.py` `_half_kernel` | `exp` of the Gaussian kernel | was CPU-dependent (numpy AVX512F loops, 1 ulp in about 5% of the entries) | decimal arithmetic, rounded once (`rstats.exp_rounded`) |
| `util/fftsmooth.py` `smooth_core`, inner function `inv` | product of two complex spectra | was CPU-dependent (numpy fuses multiply and add with FMA3) | separate multiplications and additions (`multiply_complex`) |
| `util/fftsmooth.py` `smooth_core`, inner functions `fwd_kernel`, `fwd` and `inv` | `scipy.fft.rfft` and `irfft` | libm-dependent: the twiddle factors use glibc `sin` and `cos`, which have an FMA and a plain variant | the image sets `GLIBC_TUNABLES`, so all hosts use the plain variant |
| `util/rstats.py` `t_test_one_sample_p` | `scipy.stats.t.sf` | libm-dependent (Boost calls `exp`, `log`, `pow`) | same setting |
| `util/rstats.py` `t_test_one_sample_p` | `mean`, `var` of at most a few hundred values | fixed-order reduction of numpy (pairwise sum) | no difference in 45 vector lengths under all numpy settings; R accumulates in long double (not reproduced) |
| `steps/atlas.py` `aggregate_sites`, `steps/call_dmr.py` `_nanmean_by` | `np.add.reduceat` on `longdouble` | fixed-order reduction, x87 80-bit arithmetic of every x86-64 CPU | none needed |
| `steps/atlas.py` `_summarize` | `np.median` | exact (selection, one mean of two values) | none needed |
| `util/rstats.py` `p_adjust_bh`, `r_round` | `p.adjust`, `round` | IEEE-exact; `10.0**digits` and `np.log2` were libm and SIMD calls | `Fraction` and `np.frexp` (R uses `logb`) |
| `fileio/numparse.py` `P10_EXACT` (module level) | powers of 10 | IEEE-exact; `10.0**k` was a libm call | `float(10**k)` |
| `fileio/numparse.py` `_digits15` | `np.log10` in a numba kernel | libm-dependent, but only an estimate that the loop below corrects | no effect on the output |
| all `@jit` kernels | float arithmetic in numba | IEEE-exact: `fastmath` is off and LLVM does not fuse | `NUMBA_CPU_NAME=generic`; Haswell and Skylake-AVX512 give the same files |
| `steps/target.py` `_summarise` | pandas `groupby().mean()`, `np.sqrt` | fixed algorithm of pandas (compensated sum), IEEE-exact square root | none needed |
| `np.argsort`, `np.lexsort`, `sort_values` | ordering of ties | every call is stable or has several keys; `np.sort` and `np.unique` return values only | none needed |
| `np.cumsum`, `np.add.reduceat` on integers | integer sums | exact | none needed |
| BLAS and LAPACK | none | `np.dot`, `@`, `np.linalg` and `np.convolve` do not occur in `src` | `OPENBLAS_CORETYPE` was set to three values anyway |

The kernel entries of sigma 150 and 500 equal those of R 4.4 (glibc 2.35) in all 901 and 3,001 entries. For sigma 1000 and 2000, two values of |x| in each kernel differ between the correctly rounded value and glibc 2.41 or R by one unit in the last place. The default sigma is 150.

What the table does not cover. The `scipy.fft` and `scipy.stats` rows rely on the glibc tunable, which only the image sets. A native installation uses the glibc of the host. Its `sin`, `cos`, `exp`, `log` and `pow` take the FMA path on a CPU with FMA3 and AVX2, and the plain path on an older CPU. The two paths differ in the last bit for 0.05% to 0.08% of random arguments of `exp`, `sin` and `pow` (300,000 arguments tried, glibc 2.41). In the smoothed atlas this gave differences of up to 2.9e-12 (absolute, `REF_METH`) and of up to 4.0e-12 relative in table cells, with the same calls. The image of 3.1.0 on an FMA-capable CPU without the tunable gave these differences from the image with it. A pip installation can export `GLIBC_TUNABLES=glibc.cpu.hwcaps=-FMA,-FMA4` before the command starts. We checked that this selects the plain variant for the two arguments of `env/smoke_test.py` (glibc 2.35 and glibc 2.41, FMA-capable CPU). Setting `os.environ` inside Python has no effect on the process, because glibc reads the tunable at start, so the package does not set it. On one host (Intel Xeon, glibc 2.35), a pip installation ran the demo of the README Quick start with and without that setting (and `NUMBA_CPU_NAME=generic`). In both runs the 61 data files were identical in content to those of the image; only the creation time in `reference.json` differed. The demo is small, so this does not show identity for other data.

## Platform notes

- Architecture: x86_64 only. The lock file targets `x86_64-manylinux_2_28`, and the base image digest is for `linux/amd64`. ARM (including Apple Silicon) and macOS are not supported.
- Host glibc and kernel: the image carries its own glibc (Debian trixie). The host needs only Apptainer or Docker and a Linux kernel they support.
- CPU: the image fixes every CPU-dependent choice that we found on the result path. The numba kernels are compiled for a generic x86-64 CPU (`NUMBA_CPU_NAME=generic`: SSE2, no AVX, no FMA). The Gaussian kernel and the spectrum product do not use the CPU-specific loops of numpy. `GLIBC_TUNABLES=glibc.cpu.hwcaps=-FMA,-FMA4` makes glibc use its plain libm variants on every host; do not override it with `--env`, `APPTAINERENV_GLIBC_TUNABLES` or `docker run -e`. The evidence and its limits are in the next section and in [validation.md](validation.md). A native installation outside the image is not covered.
- Python, numpy and numba: numba 0.68.0 declares `numpy>=1.22,<2.6` and `llvmlite>=0.50,<0.51`; the lock has numpy 2.5.3 and llvmlite 0.50.0. Because `uv` resolved all versions together and the smoke test and unit tests of this release pass with them, do not replace single packages in the image.
