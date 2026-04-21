# Changelog

All notable changes to this project will be documented in this file.


## [1.1] - 2026-04-21

Publication-ready release accompanying the manuscript *"Quick detection of
differentially methylated regions for CpG islands and imprinted CpGs based on
an ONT methylation reference atlas."*

### Added

- `--include_sex_chromosomes` parameter (default: `false`). Autosomes (1-22)
  are always analysed; chrX/chrY are opt-in because they require special
  statistical treatment (ploidy, X-inactivation) not modelled here.
- Samplesheet validation at launch: required headers (`sample_id,type,meth_bed,vcf`),
  `type` ∈ {CASE, CTRL}, no duplicate `sample_id`s. Clear row-level error messages.
- `test/smoke/` end-to-end smoke test with synthetic chr22 data (3 CTRL + 1 CASE),
  exercising default, `--include_sex_chromosomes true`, and malformed-samplesheet
  scenarios.

### Changed

- Split the monolithic `DETECT_DMR` process into three focused processes:
  - `PARSE_GTF` — runs once, parses the Ensembl GTF into a slim gene-annotation
    table (avoids re-parsing the GTF for every CASE sample).
  - `CALL_DMR` — per-sample DMR calling (sample/atlas merge, panel overlap,
    t-test, BH-FDR). Shared `finalise_dmrs()` helper removes duplicated
    post-merge statistics between CpG and Panel2 branches.
  - `ANNOTATE_DMRS` — per-sample nearest-gene annotation and final formatting.
  Each process's R block now fits on one screen; `main.nf` dropped from 830 to
  ~560 lines.
- Container pinned from `:latest` to `:v1.0` for reproducible publication runs.
- Container image reference normalised to lowercase (`ghcr.io/oumarkme/dmrdetector`)
  to satisfy strict Apptainer reference-format validation.

### Removed

- `bioconductor-biomart` from the container — never imported by any script.
  Shrinks the image and simplifies the dependency surface.


## [1.0] - 2026-03-12

First public release.

### Added

- `--tmp_dir` parameter for local temp/work directory (faster I/O when inputs are on NAS)
- Auto-download latest GRCh38 GTF from Ensembl when `--gtf` is not provided
- Show help text when `--input` is not provided
- Apptainer cache directory automatically placed under `tmp_dir` (or `outdir`)
- CASE-only filtering for DMR detection stage (CTRL samples skip Stage 3)
- Conda environment setup instructions in README
- Professional README with badges, pipeline flow diagram, and collapsible parameter sections
- Author and citation information

### Changed

- Work directory (`workDir`) defaults to `outdir/work` instead of Nextflow default (`./work`)
- Process `TMPDIR` defaults to `outdir/tmp` when `--tmp_dir` is not set
- Replaced `check_max` helper function with inline resource capping in base.config
- Removed deprecated `nextflow.enable.dsl = 2` (DSL2 is default since Nextflow 22.12)
- Fixed implicit closure parameters and unused catch variables for Nextflow 25.x compatibility
- Fully qualified `nextflow.util.MemoryUnit` and `nextflow.util.Duration` for broad version compatibility
- Execution report output moved to `results/nf_report.html`
- Pipeline made fully portable (no machine-specific paths)

### Removed

- `gs4` profile (machine-specific paths)
- Timeline, trace, and DAG reports (kept execution report only)

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
