# Changelog

All notable changes to this project will be documented in this file.


## [Unreleased]

## [2.0] - 2026-04-27

### Changed

- CpG island and ICR panels now run by default. The former `panel2` interface
  is replaced by `icr_panel`, `icr_diff_threshold`, and `icr_min_sites`.
- Final DMR output remains one file per CASE sample:
  `results/dmr/dmr_{sample}_case.tsv`, with `PANEL` values `CpG` or `ICR`.
- Sex chromosome analysis now follows the manuscript strategy:
  autosomal chr1-22, male non-PAR chrX, female non-PAR chrX, and male chrY.
- Female chrX now uses a hemi-methylated CpG island model
  (`sex_f_ref_meth_min` to `sex_f_ref_meth_max`) and relaxed DMR thresholds
  (`sex_f_diff_threshold`, `sex_f_min_sites`).
- DMR-only runs with `--ref_atlas` no longer require `--cpg_panel`; CpG islands
  are detected from the provided atlas unless a panel override is supplied.
- ICRs are now defined as hemi-methylated regions within imprinting gene
  candidate intervals. `icr_gene_list` is the candidate BED source, while
  `icr_panel` must contain ICR regions rather than whole gene bodies.
- DMRdetector v2.0 is hg38-only; GRCh37 PAR selection has been removed.

### Added

- Bundled default BED assets under `assets/`: `cpg_island.bed`,
  `icr_gene_list.bed`, `icr_region.bed`, and `par_hg38.bed`.
- `DETECT_ICR_REGIONS` process to build ICRs from an autosomal reference atlas
  and imprinting gene candidate BED.

## [1.3] - 2026-04-27

### Changed

- Redesigned sex-chromosome analysis around component scopes:
  `autosomal`, `male_X`, `female_X`, and `Y`.
- PAR now stays with chrX/chrY in methylation analysis. PAR coordinates are
  still used for VCF-based sex inference, but PAR rows are no longer merged
  into the autosomal atlas.
- CpG island identification now uses male-derived atlas scopes only:
  autosomal, `male_X`, and `Y`. Female chrX DMRs reuse the `male_X` CpG
  panel while comparing against the `female_X` atlas.
- CASE DMR routing is sex-specific: male CASE samples run
  `autosomal + male_X + Y`, female CASE samples run
  `autosomal + female_X`, and unknown-sex CASE samples run autosomal only.
- ICR/panel2 analysis is restricted to autosomes and excludes X/Y rows.
- `sex_f_*` parameters are retained for compatibility but are inactive in the
  v1.3 DMR path.

## [1.2] - 2026-04-27

### Added

- Sex-stratified chrX/chrY analysis, enabled by default
  (`include_sex_chromosomes = true`; pass `--include_sex_chromosomes false`
  to opt out). The pipeline now:
  - Infers sample sex from the VCF by counting non-PAR chrX heterozygous
    SNPs (`INFER_SEX_FROM_VCF`); per-sample reports and a combined
    `pipeline_info/sex_inference_summary.tsv` are published.
  - Splits each methylation BED into autosomal + PAR / non-PAR chrX /
    non-PAR chrY streams (`SPLIT_PAR_AUTOSOMAL`) using built-in PAR
    coordinates for GRCh38 or GRCh37 (`assets/par_GRCh38.bed`,
    `assets/par_GRCh37.bed`; override with `--par_bed`).
  - Builds three reference atlases in one pass (scope-tagged channels into
    a single `BUILD_ATLAS` + `SMOOTH_ATLAS`): autosomal + PAR, male non-PAR
    chrX+chrY, and female non-PAR chrX. Outputs are suffixed with
    `_sex_M` / `_sex_F`.
  - Detects CpG panels per scope, with a dedicated
    `DETECT_CPG_ISLANDS_SEX_F` that applies an imprinting-style rule
    (REF_METH within the XCI band) for female chrX.
  - Calls DMRs three times per CASE sample with scope-specific thresholds
    (autosomal / sex_M use the standard CpG-panel thresholds; sex_F uses
    `sex_f_diff_threshold` 0.25, `sex_f_min_sites` 5, and a REF_METH gate
    of `[sex_f_ref_meth_min, sex_f_ref_meth_max]`). Results from all scopes
    are merged by `ANNOTATE_DMRS` and tagged with a `PANEL_SCOPE` column.
- New parameters: `genome_build`, `par_bed`, `sex_atlas_min_controls`,
  `sex_het_threshold_male`, `sex_het_threshold_female`, `sex_override_tsv`,
  `sex_f_ref_meth_min`, `sex_f_ref_meth_max`, `sex_f_diff_threshold`,
  `sex_f_min_sites`.
- Smoke test now covers `--include_sex_chromosomes true` via
  `test/smoke/sex_override.tsv`, asserting all three atlases exist,
  PAR/non-PAR separation is correct, the three CpG panels are produced,
  the DMR TSV carries `PANEL_SCOPE`, and `sex_inference_summary.tsv`
  reports all samples.

### Changed

- `include_sex_chromosomes` default flipped to `true`. Existing cohorts that
  lack VCFs, have fewer than `sex_atlas_min_controls` male/female controls,
  or should stay on the legacy autosomes-only path must now pass
  `--include_sex_chromosomes false` explicitly.
- `BUILD_ATLAS`, `SMOOTH_ATLAS`, and `DETECT_CPG_ISLANDS` now accept a
  scope-tagged tuple (`autosomal` / `sex_M` / `sex_F`) so a single process
  definition serves both autosomes-only and sex-stratified modes.
- `ANNOTATE_DMRS` now groups raw per-scope DMR files by sample and emits
  the `PANEL_SCOPE` column only when sex branches are active (autosomes-only
  output remains byte-identical with pre-feature releases).
- Default Apptainer execution now uses the GitHub Actions-published SIF image
  (`oras://ghcr.io/oumarkme/dmrdetector-apptainer:sha-a00657d`); Docker
  profile execution uses `ghcr.io/oumarkme/dmrdetector:sha-a00657d`.

### Fixed

- `CALL_DMR` emitted a phantom NA row when a scope produced zero DMRs,
  which later broke `ANNOTATE_DMRS` (`foverlaps` rejected the NA start).
  Empty-result files are now written as header-only with the canonical
  schema.
- GTF and raw DMR chromosome names are normalized before annotation so both
  Ensembl-style (`1`, `2`, `X`, `Y`) and UCSC-style (`chr1`, `chr2`, `chrX`,
  `chrY`) inputs annotate correctly. Final DMR output uses UCSC-style names.


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
