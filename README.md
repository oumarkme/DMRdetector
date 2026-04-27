<p align="center">
  <h1 align="center">DMRdetector</h1>
  <p align="center">
    ONT Methylation Analysis Pipeline<br>
    <em>Detect differentially methylated regions from Oxford Nanopore (or PacBio) whole-genome sequencing data</em>
  </p>
</p>

<p align="center">
  <a href="CHANGELOG.md"><img src="https://img.shields.io/badge/version-1.3-blue.svg" alt="Version"></a>
  <a href="https://www.nextflow.io/"><img src="https://img.shields.io/badge/Nextflow-%E2%89%A523.04.0-brightgreen.svg" alt="Nextflow"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://ghcr.io/oumarkme/DMRdetector"><img src="https://img.shields.io/badge/Container-GHCR-blue.svg" alt="Container"></a>
  <a href="https://apptainer.org/"><img src="https://img.shields.io/badge/Apptainer-supported-orange.svg" alt="Apptainer"></a>
  <a href="https://github.com/oumarkme/DMRdetector/actions"><img src="https://github.com/oumarkme/DMRdetector/actions/workflows/docker.yml/badge.svg" alt="CI"></a>
  <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20WSL2%20%7C%20macOS-lightgrey.svg" alt="Platform">
  <img src="https://img.shields.io/badge/DOI-pending-grey.svg" alt="DOI">
</p>

---

## Overview

**DMRdetector** builds a population-level methylation reference atlas from control samples and detects differentially methylated regions (DMRs) in case samples. It supports an optional second panel for targeted region analysis (e.g., imprinted control regions).

### Pipeline Stages

```
CTRL + CASE samples
  │
  ├─ PREPROCESS ──────────── All samples (filter HET SNPs, coverage filter)
  │
  ├─ Stage 1: BUILD ATLAS ── CTRL only ─── aggregate, dispersion, FFT smoothing
  │
  ├─ Stage 2: DETECT CpG ─── CTRL only ─── identify CORE/SHORE, cluster, extend
  │
  └─ Stage 3: DETECT DMRs ── CASE only
        ├─ PARSE_GTF       (once)      gene annotation table
        ├─ CALL_DMR        (per case)  sample/atlas merge, t-test, FDR
        └─ ANNOTATE_DMRS   (per case)  nearest-gene overlap, final TSV
```

Autosomes (1-22) are always analysed. Sex chromosomes (X, Y) are handled by
a sex-stratified workflow that runs by **default**: the pipeline infers
sample sex from the VCF and builds component references for autosomes,
male chrX, female chrX, and chrY. PAR stays with X/Y in v1.3. Pass
`--include_sex_chromosomes false` to restrict the run to autosomes only.
See [Sex-chromosome analysis](#sex-chromosome-analysis) below.

### Usage Scenarios

| Scenario | Input | What runs |
|----------|-------|-----------|
| **Full pipeline** | CTRL + CASE samples | Stages 1 → 2 → 3 |
| **Build atlas only** | CTRL samples only | Stages 1 → 2 |
| **DMR detection only** | CASE samples + `--ref_atlas` + `--cpg_panel` | Stage 3 only |

---

## Getting Started

### 1. Set up the environment

We recommend using [Conda](https://docs.conda.io/) to manage Nextflow and Apptainer:

```bash
conda create -n DMRdetector -c conda-forge -c bioconda nextflow apptainer
conda activate DMRdetector
```

### 2. Download the pipeline

```bash
git clone https://github.com/oumarkme/DMRdetector.git
cd DMRdetector
```

### 3. Download pre-built reference atlas (optional)

A pre-built reference atlas is available on Zenodo: [10.5281/zenodo.18978867](https://doi.org/10.5281/zenodo.18978867). You can use it to skip Stage 1–2 and run DMR detection directly (see [Advanced Usage](#skip-atlas-building-pre-built-references)).

### 4. Run

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results
```

The Apptainer image is pulled automatically from [GHCR](https://ghcr.io/oumarkme/DMRdetector) on first run. The cached image is stored under `tmp_dir/apptainer/` (or `outdir/apptainer/` if `--tmp_dir` is not set). Docker users can run the same pipeline with `-profile docker`.

> **Tip:** Use `--tmp_dir /local/fast/disk` when your input files are on network storage for significantly faster I/O.

---

## Input

### Samplesheet

A CSV file with the following columns:

```csv
sample_id,type,meth_bed,vcf
ONT001,CASE,/path/to/ONT001.combined.bed.gz,/path/to/ONT001.vcf
ONT004,CTRL,/path/to/ONT004.combined.bed.gz,/path/to/ONT004.vcf
```

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier |
| `type` | `CASE` or `CTRL` |
| `meth_bed` | Path to methylation BED file (gzipped) |
| `vcf` | Path to VCF file (used to exclude heterozygous SNP sites) |

### Methylation BED Formats

The pipeline auto-detects the input format based on column count:

| Format | Columns | COVERAGE | N_MOD |
|--------|---------|----------|-------|
| **5-column** (default) | CHROM, START, END, COVERAGE, N_MOD | col 4 | col 5 |
| **modkit pileup** | 18 columns (bedMethyl) | col 10 | col 12 |
| **pb-CpG-tools** | 7-9 columns (with header) | col 6 | col 7 |

### GTF Annotation

If `--gtf` is not provided, the **latest GRCh38 GTF is auto-downloaded** from Ensembl. To use a specific version:

```bash
# GRCh38 (hg38)
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

# GRCh37 (hg19)
wget https://ftp.ensembl.org/pub/grch37/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
```

Then pass with `--gtf Homo_sapiens.GRCh38.115.gtf.gz`. Use the GTF version matching your alignment reference.

### Optional Panel2

A headerless 3-column BED file (`CHROM, START, END`) for a second panel analysis (e.g., imprinted control regions). Gene annotation is performed via GTF overlap.

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --panel2 my_regions.bed \
    --outdir results
```

---

## Output

```
results/
├── atlas/
│   ├── ref_atlas_150.tsv.gz          # Autosomal chr1-22 reference atlas (narrow)
│   ├── ref_atlas_500.tsv.gz          # Autosomal chr1-22 reference atlas (wide)
│   ├── ref_atlas_150_male_X.tsv.gz   # Male CTRL chrX atlas
│   ├── ref_atlas_500_male_X.tsv.gz
│   ├── ref_atlas_150_female_X.tsv.gz # Female CTRL chrX atlas
│   ├── ref_atlas_500_female_X.tsv.gz
│   ├── ref_atlas_150_Y.tsv.gz        # Male CTRL chrY atlas
│   └── ref_atlas_500_Y.tsv.gz
├── cpg_islands/
│   ├── cpg_signal.bed                # Autosomal CpG panel
│   ├── cpg_signal_male_X.bed         # Male chrX-derived CpG panel
│   └── cpg_signal_Y.bed              # Male chrY-derived CpG panel
├── dmr/
│   └── dmr_{sample}_case.tsv     # DMR results per case sample (with PANEL_SCOPE column by default)
├── pipeline_info/
│   ├── sex_inference/            # Per-sample VCF-based sex inference reports
│   └── sex_inference_summary.tsv # Combined sex inference summary
└── nf_report.html                # Nextflow execution report
```

### DMR Output Columns

| Column | Description |
|--------|-------------|
| CHROM, START, END | Genomic coordinates |
| GENE | Nearest gene symbol |
| PANEL | `CpG` or `Panel2` |
| NUM_CPG | Number of CpG sites in region |
| NUM_READS | Mean read coverage |
| SAMPLE_METH (%) | Sample methylation level |
| REF_METH (%) | Reference methylation level |
| DIFF (%) | Difference (sample - reference) |
| P_VALUE | T-test p-value |
| P_FDR | Benjamini-Hochberg adjusted p-value |
| DIRECTION | `hyper` or `hypo` |

A `PANEL_SCOPE` column (`autosomal` / `male_X` / `female_X` / `Y`) identifies
which reference branch produced each row. It is omitted when
`--include_sex_chromosomes false` is passed.

---

## Parameters

Default values are defined in [`nextflow.config`](nextflow.config). Override on the command line:

```bash
nextflow run main.nf --input samplesheet.csv --dmr_fdr_threshold 0.05
```

<details>
<summary><strong>Reference Atlas (Stage 1)</strong></summary>

| Parameter | Default | Description |
|-----------|---------|-------------|
| `atlas_min_samples_per_site` | 3 | Minimum samples per CpG site |
| `atlas_min_read_coverage` | 5 | Minimum read coverage |
| `atlas_shrinkage_weight` | 0.1 | Empirical Bayes shrinkage weight |
| `smooth_sigma_narrow` | 150 | Narrow Gaussian bandwidth (bp) |
| `smooth_sigma_wide` | 500 | Wide Gaussian bandwidth (bp) |

</details>

<details>
<summary><strong>CpG Island Detection (Stage 2)</strong></summary>

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cpg_min_samples_per_site` | 5 | Minimum samples for pre-filtering |
| `cpg_core_max_methylation` | 0.1 | CORE site methylation threshold |
| `cpg_core_max_dispersion` | 0.05 | CORE site dispersion threshold |
| `cpg_shore_max_methylation` | 0.3 | SHORE site methylation threshold |
| `cpg_max_gap_distance` | 150 | Max gap for clustering (bp) |
| `cpg_min_island_length` | 150 | Minimum island length (bp) |

</details>

<details>
<summary><strong>DMR Detection (Stage 3)</strong></summary>

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dmr_min_coverage` | 10 | Minimum coverage for case samples |
| `dmr_cpg_diff_threshold` | 0.40 | \|DIFF\| threshold for CpG panel |
| `dmr_merge_distance` | 150 | Max gap for merging sites (bp) |
| `dmr_fdr_threshold` | 0.01 | FDR significance threshold |
| `dmr_cpg_min_sites` | 30 | Min CpG sites for CpG panel |
| `panel2_diff_threshold` | 0.25 | \|DIFF\| threshold for panel2 |
| `panel2_min_sites` | 5 | Min CpG sites for panel2 |

</details>

<details>
<summary><strong>Sex Chromosomes (default)</strong></summary>

Active only when `--include_sex_chromosomes true`. See
[Sex-chromosome analysis](#sex-chromosome-analysis) for the full workflow.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `include_sex_chromosomes` | true | Enable sex-stratified chrX/chrY analysis (default). Set to `false` for autosomes only |
| `genome_build` | `GRCh38` | PAR coordinate set (`GRCh38` or `GRCh37`) |
| `par_bed` | null | Override built-in PAR BED (4-col: CHROM START END NAME) |
| `sex_atlas_min_controls` | 3 | Minimum male/female CTRL count required for sex-specific CASE comparisons |
| `sex_het_threshold_male` | 100 | Non-PAR chrX HET SNP count below → `M` |
| `sex_het_threshold_female` | 1000 | Non-PAR chrX HET SNP count above → `F` (between → `U`, autosomes only) |
| `sex_override_tsv` | null | Optional 2-col TSV (`sample_id`, `sex`) to override inference |
| `sex_f_ref_meth_min` | 0.35 | Deprecated compatibility parameter; inactive in v1.3 |
| `sex_f_ref_meth_max` | 0.65 | Deprecated compatibility parameter; inactive in v1.3 |
| `sex_f_diff_threshold` | 0.25 | Deprecated compatibility parameter; inactive in v1.3 |
| `sex_f_min_sites` | 5 | Deprecated compatibility parameter; inactive in v1.3 |

</details>

<details>
<summary><strong>Resource Management</strong></summary>

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_cpus` | auto | Maximum total CPUs |
| `max_memory` | auto | Maximum total memory |
| `max_time` | 168.h | Maximum time per process |
| `tmp_dir` | null | Local temp directory for work/tmp files |
| `apptainer_bind` | null | Comma-separated host paths to bind-mount |

Resources are auto-detected at startup (80% of system RAM; CPUs capped by RAM/4GB ratio):

| System | Cores | RAM | max_cpus | max_memory |
|--------|-------|-----|----------|------------|
| Laptop | 8 | 16 GB | 3 | 12 GB |
| Workstation | 16 | 64 GB | 12 | 51 GB |
| Server | 96 | 2 TB | 76 | 1.6 TB |

Override: `--max_cpus 4 --max_memory '16.GB'`

</details>

---

## Container and Bind Mounts

### Apptainer (Linux / WSL2) — default

Default image: `oras://ghcr.io/oumarkme/dmrdetector-apptainer:sha-a00657d`

If input files are outside the working directory, bind-mount host paths:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --apptainer_bind '/data,/mnt/nas' \
    --outdir results
```

### Docker (macOS)

Docker profile image: `ghcr.io/oumarkme/dmrdetector:sha-a00657d`

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --outdir results
```

Docker Desktop shares `/Users` and `/tmp` by default. Add other paths in Docker Desktop > Settings > Resources > File Sharing.

---

## Advanced Usage

### Skip atlas building (pre-built references)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --ref_atlas results/atlas/ref_atlas_150.tsv.gz \
    --cpg_panel results/cpg_islands/cpg_signal.bed \
    --outdir results
```

### Resume after failure

Nextflow caches completed tasks. Add `-resume` to continue where you left off:

```bash
nextflow run main.nf [same options] -resume
```

### Profiles

| Profile | Description |
|---------|-------------|
| `test` | Reduced resources (4 CPUs, 16 GB) for CI/testing |
| `docker` | Use Docker instead of Apptainer (required for macOS) |

---

## Sex-chromosome analysis

Sex-stratified analysis of chrX/chrY runs by default. Pass
`--include_sex_chromosomes false` to disable it (see
[Autosomes-only mode](#autosomes-only-mode) below). In the default path the
pipeline:

1. **Infers sex from the VCF** (`INFER_SEX_FROM_VCF`) by counting heterozygous
   SNPs on non-PAR chrX. Males are hemizygous (HET ≈ 0); females are diploid
   (HET in the thousands). Counts below `sex_het_threshold_male` → `M`,
   above `sex_het_threshold_female` → `F`, between → `U` (autosomes only).
   An audit row per sample is written to
   `pipeline_info/sex_inference_summary.tsv`.
2. **Splits each methylation BED** (`SPLIT_SEX_COMPONENTS`) into component
   scopes: autosomal chr1-22, full chrX, and full chrY. PAR stays with X/Y;
   the built-in PAR coordinates are used only for VCF-based sex inference.
3. **Builds four reference atlases** (`BUILD_ATLAS` + `SMOOTH_ATLAS`):
   - **autosomal** — chr1-22 from all CTRLs, written as
     `ref_atlas_{150,500}.tsv.gz`.
   - **male_X** — full chrX from male CTRLs.
   - **female_X** — full chrX from female CTRLs.
   - **Y** — full chrY from male CTRLs.
   If a CASE requires a sex-specific comparison and too few matching CTRLs are
   available, the pipeline exits with a clear error.
4. **Detects CpG panels from male-derived atlases only**: autosomal, `male_X`,
   and `Y` use the standard CORE (REF_METH < 10%) + SHORE
   (REF_METH < 30%) rules. Female chrX DMRs reuse the `male_X` CpG panel but
   compare methylation against the `female_X` atlas.
5. **Calls DMRs against sex-appropriate component atlases**:

   | Scope | `diff` threshold | `min_sites` | REF_METH gate |
   |---|---|---|---|
   | `autosomal` | `dmr_cpg_diff_threshold` (0.40) | `dmr_cpg_min_sites` (30) | none |
   | `male_X` | `dmr_cpg_diff_threshold` (0.40) | `dmr_cpg_min_sites` (30) | none |
   | `female_X` | `dmr_cpg_diff_threshold` (0.40) | `dmr_cpg_min_sites` (30) | none |
   | `Y` | `dmr_cpg_diff_threshold` (0.40) | `dmr_cpg_min_sites` (30) | none |

   Male CASE samples run `autosomal + male_X + Y`; female CASE samples run
   `autosomal + female_X`; unknown-sex CASE samples run autosomal only with a
   warning. Results are merged per sample by `ANNOTATE_DMRS`, tagged with a
   `PANEL_SCOPE` column, and published as a single `dmr_{sample}_case.tsv`.

### Overriding inferred sex

If a sample has an aneuploid karyotype (Turner 45,X, Klinefelter 47,XXY,
etc.) or clinical records that supersede the VCF inference, pass a 2-column
TSV:

```tsv
sample_id	sex
ONT042	F
ONT098	M
```

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --include_sex_chromosomes true \
    --sex_override_tsv sex_calls.tsv \
    --outdir results
```

Overrides are applied after inference; the audit summary records the source
(`inferred` vs. `override`).

### Autosomes-only mode

Pass `--include_sex_chromosomes false` to opt out of the sex-stratified
workflow. Only the autosomal branch runs, no sex inference is performed,
and outputs match the pre-feature layout: `ref_atlas_{150,500}.tsv.gz`,
`cpg_signal.bed`, and `dmr_{sample}_case.tsv` without a `PANEL_SCOPE`
column.

---

## Citation

If you use this pipeline in your research, please cite:

> Ou JH, Tsai M, et al. A map of CpG islands based on an ONT methylation reference atlas and a framework for detection of differentially methylated regions. *(manuscript in preparation)*

## Authors

Jen-Hsiang (Mark) Ou\*, Melody Tsai\*, Ni-Chung Lee†, Wuh-Liang Hwu†

\* Equal contribution · † Corresponding authors

**Pipeline maintainer:** Jen-Hsiang (Mark) Ou — [oujenhsiang@ntu.edu.tw](mailto:oujenhsiang@ntu.edu.tw)

## License

This project is licensed under the [MIT License](LICENSE).

Copyright (c) 2025 Jen-Hsiang (Mark) Ou
