# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-03-12

### Added

- `--tmp_dir` parameter for local temp/work directory (faster I/O when inputs are on NAS)
- Auto-download GTF from Ensembl when `--gtf` is not provided
- Show help text when `--input` is not provided

### Changed

- Work directory (`workDir`) defaults to `outdir/work` instead of Nextflow default (`./work`)
- Process `TMPDIR` defaults to `outdir/tmp` when `--tmp_dir` is not set
- Removed `nextflow.enable.dsl = 2` (DSL2 is default since Nextflow 22.12)
- Replaced `check_max` helper function with inline resource capping in base.config
- Removed timeline, trace, and DAG reports (kept execution report only)
- Updated license copyright holder
- Added license information to README

### Removed

- `gs4` profile (machine-specific paths)

## [0.1.0] - 2025-03-12

### Added

- Full pipeline: Build reference atlas, detect CpG islands, detect DMRs
- Auto-detection of input BED formats (5-col, modkit bedMethyl 18-col, pb-CpG-tools 9-col)
- FFT Gaussian smoothing at two bandwidths (narrow + wide)
- Empirical Bayes shrinkage for dispersion estimation
- Optional second panel (`--panel2`) for targeted region analysis
- Pre-built atlas support (`--ref_atlas`, `--cpg_panel`) to skip Stages 1-2
- Cross-platform resource auto-detection (Linux, macOS, WSL2)
- RAM-proportional CPU capping to prevent OOM on memory-limited systems
- Apptainer and Docker container support
- GitHub Actions CI for automatic container builds
- Configurable bind mounts via `--apptainer_bind`
