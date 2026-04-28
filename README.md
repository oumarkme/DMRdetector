<p align="center">
  <h1 align="center">DMRdetector</h1>
  <p align="center">
    ONT Methylation Analysis Pipeline<br>
    <em>Detect differentially methylated regions from Oxford Nanopore (or PacBio) whole-genome sequencing data</em>
  </p>
</p>

<p align="center">
  <a href="CHANGELOG.md"><img src="https://img.shields.io/badge/version-2.0-blue.svg" alt="Version"></a>
  <a href="https://www.nextflow.io/"><img src="https://img.shields.io/badge/Nextflow-%E2%89%A523.04.0-brightgreen.svg" alt="Nextflow"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://docs.conda.io/"><img src="https://img.shields.io/badge/Conda-supported-44A833.svg" alt="Conda"></a>
  <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20WSL2%20%7C%20macOS-lightgrey.svg" alt="Platform">
  <img src="https://img.shields.io/badge/DOI-pending-grey.svg" alt="DOI">
</p>

---

## Overview

**DMRdetector** builds a population-level methylation reference atlas from control samples and detects differentially methylated regions (DMRs) in case samples. CpG island and imprinted control region (ICR) panels are both analysed by default.

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
        ├─ CALL_DMR        (per case)  sample/atlas merge, t-test, FDR
        └─ ANNOTATE_DMRS   (per case)  EnsDb.Hsapiens.v86 gene overlap, final TSV
```

Autosomes (1-22) are always analysed. Sex chromosomes (X, Y) are handled by
a sex-stratified workflow that runs by **default**: sex (`M`/`F`) is read
from the samplesheet, and the pipeline builds component references for
autosomes, male non-PAR chrX, female non-PAR chrX, and chrY. Pass
`--include_sex_chromosomes false` to restrict the run to autosomes only.
See [Sex-chromosome analysis](#sex-chromosome-analysis) below.

### Usage Scenarios

| Scenario | Input | What runs |
|----------|-------|-----------|
| **Full pipeline** | CTRL + CASE samples | Stages 1 → 2 → 3 |
| **Build atlas only** | CTRL samples only | Stages 1 → 2 |
| **DMR detection only** | CASE samples + `--ref_atlas` (optional `--cpg_panel`) | Stage 3 only |

---

## Getting Started

### 1. Clone and set up the conda environment

The bundled `environment.yml` installs **everything** you need (Nextflow,
R + Bioconductor packages, bcftools, bedtools, htslib) into a single
conda environment. No Docker or Apptainer required.

```bash
git clone https://github.com/oumarkme/DMRdetector.git
cd DMRdetector
conda env create -f environment.yml
conda activate dmrdetector
```

### 2. Download pre-built reference atlas (optional)

A pre-built reference atlas is available on Zenodo: [10.5281/zenodo.18978867](https://doi.org/10.5281/zenodo.18978867). You can use it to skip Stage 1–2 and run DMR detection directly (see [Advanced Usage](#skip-atlas-building-pre-built-references)).

### 3. Run

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results
```

> **Tip:** Use `--tmp_dir /local/fast/disk` when your input files are on network storage for significantly faster I/O.

---

## Input

### Samplesheet

A CSV file with the following columns:

```csv
sample_id,type,sex,meth_bed,vcf
ONT001,CASE,M,/path/to/ONT001.combined.bed.gz,/path/to/ONT001.vcf
ONT004,CTRL,F,/path/to/ONT004.combined.bed.gz,/path/to/ONT004.vcf
```

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier |
| `type` | `CASE` or `CTRL` |
| `sex` | `M` or `F` (required; routes chrX/chrY to the correct component scope) |
| `meth_bed` | Path to methylation BED file (gzipped) |
| `vcf` | Path to VCF file (used to exclude heterozygous SNP sites) |

### Methylation BED Formats

The pipeline auto-detects the input format based on column count:

| Format | Columns | COVERAGE | N_MOD |
|--------|---------|----------|-------|
| **5-column** (default) | CHROM, START, END, COVERAGE, N_MOD | col 4 | col 5 |
| **modkit pileup** | 18 columns (bedMethyl) | col 10 | col 12 |
| **pb-CpG-tools** | 7-9 columns (with header) | col 6 | col 7 |

### Gene Annotation

DMRs are annotated with the nearest gene using the `EnsDb.Hsapiens.v86`
Bioconductor package (GRCh38, Ensembl release 86) bundled with the conda
environment. No external GTF download is needed; per-region selection
prefers `protein_coding` over `lncRNA` over pseudogenes, with overlap
width and gene length as tie-breakers.

### ICR Panel

ICR analysis runs by default. When a reference atlas is available, DMRdetector defines ICRs as hemi-methylated regions within imprinting gene candidate intervals from [`assets/icr_gene_list.bed`](assets/icr_gene_list.bed). To override the ICR panel, pass an ICR BED such as [`assets/icr_region.bed`](assets/icr_region.bed); do not pass whole imprinting gene body coordinates as `--icr_panel`.

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --icr_panel my_icr_regions.bed \
    --outdir results
```

---

## Output

```
results/
├── atlas/
│   ├── ref_atlas_150.tsv.gz          # Autosomal chr1-22 reference atlas (narrow)
│   ├── ref_atlas_500.tsv.gz          # Autosomal chr1-22 reference atlas (wide)
│   ├── ref_atlas_150_male_X.tsv.gz   # Male CTRL non-PAR chrX atlas
│   ├── ref_atlas_500_male_X.tsv.gz
│   ├── ref_atlas_150_female_X.tsv.gz # Female CTRL non-PAR chrX atlas
│   ├── ref_atlas_500_female_X.tsv.gz
│   ├── ref_atlas_150_Y.tsv.gz        # Male CTRL chrY atlas
│   └── ref_atlas_500_Y.tsv.gz
├── cpg_islands/
│   ├── cpg_signal.bed                # Autosomal CpG panel
│   ├── cpg_signal_male_X.bed         # Male non-PAR chrX CpG panel
│   ├── cpg_signal_female_X.bed       # Female non-PAR chrX hemi-methylated panel
│   └── cpg_signal_Y.bed              # Male chrY CpG panel
├── icr_regions/
│   └── icr_region.bed                # ICRs defined from imprinting genes + atlas
├── dmr/
│   └── dmr_{sample}_case.tsv     # DMR results per case sample (with PANEL_SCOPE column by default)
└── nf_report.html                # Nextflow execution report
```

### DMR Output Columns

| Column | Description |
|--------|-------------|
| CHROM, START, END | Genomic coordinates |
| GENE | Nearest gene symbol |
| PANEL | `CpG` or `ICR` |
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
| `icr_panel` | null | Advanced ICR BED override; must contain ICR regions, not whole gene bodies |
| `icr_gene_list` | `assets/icr_gene_list.bed` | Imprinting gene candidate BED used to define ICRs from an atlas |
| `icr_ref_meth_min` | 0.40 | Lower REF_METH bound for hemi-methylated ICR definition |
| `icr_ref_meth_max` | 0.60 | Upper REF_METH bound for hemi-methylated ICR definition |
| `icr_diff_threshold` | 0.25 | \|DIFF\| threshold for ICR panel |
| `icr_min_sites` | 5 | Min CpG sites for ICR panel |

</details>

<details>
<summary><strong>Sex Chromosomes (default)</strong></summary>

Active only when `--include_sex_chromosomes true`. See
[Sex-chromosome analysis](#sex-chromosome-analysis) for the full workflow.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `include_sex_chromosomes` | true | Enable sex-stratified chrX/chrY analysis (default). Set to `false` for autosomes only |
| `par_bed` | null | Override built-in hg38 PAR BED (4-col: CHROM START END NAME) |
| `sex_atlas_min_controls` | 3 | Minimum male/female CTRL count required for sex-specific CASE comparisons |
| `sex_f_ref_meth_min` | 0.40 | Female_X hemi-methylated CpG island lower REF_METH bound |
| `sex_f_ref_meth_max` | 0.60 | Female_X hemi-methylated CpG island upper REF_METH bound |
| `sex_f_diff_threshold` | 0.25 | \|DIFF\| threshold for female_X DMRs |
| `sex_f_min_sites` | 5 | Min CpG sites for female_X DMRs |

</details>

<details>
<summary><strong>Resource Management</strong></summary>

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_cpus` | auto | Maximum total CPUs |
| `max_memory` | auto | Maximum total memory |
| `max_time` | 168.h | Maximum time per process |
| `tmp_dir` | null | Local temp directory for work/tmp files |

Resources are auto-detected at startup (80% of system RAM; CPUs capped by RAM/4GB ratio):

| System | Cores | RAM | max_cpus | max_memory |
|--------|-------|-----|----------|------------|
| Laptop | 8 | 16 GB | 3 | 12 GB |
| Workstation | 16 | 64 GB | 12 | 51 GB |
| Server | 96 | 2 TB | 76 | 1.6 TB |

Override: `--max_cpus 4 --max_memory '16.GB'`

</details>

---

## Advanced Usage

### Skip atlas building (pre-built references)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --ref_atlas results/atlas/ref_atlas_150.tsv.gz \
    --outdir results
```

When `--ref_atlas` is provided without `--cpg_panel`, DMRdetector detects a
CpG island panel from that atlas. Pass `--cpg_panel` only when you want to
override this with a pre-built BED such as [`assets/cpg_island.bed`](assets/cpg_island.bed).

### Resume after failure

Nextflow caches completed tasks. Add `-resume` to continue where you left off:

```bash
nextflow run main.nf [same options] -resume
```

### Profiles

| Profile | Description |
|---------|-------------|
| `test` | Reduced resources (4 CPUs, 16 GB) for CI/testing |

---

## Sex-chromosome analysis

Sex-stratified analysis of chrX/chrY runs by default. Pass
`--include_sex_chromosomes false` to disable it (see
[Autosomes-only mode](#autosomes-only-mode) below). Sex (`M` or `F`) is a
required column in the input samplesheet — the pipeline does not infer sex.
In the default path the pipeline:

1. **Reads sex per sample** from the samplesheet `sex` column (M/F). Validation
   at launch rejects missing or non-M/F values.
2. **Splits each methylation BED** (`SPLIT_SEX_COMPONENTS`) into component
   scopes: autosomal chr1-22, non-PAR chrX, and chrY. Built-in PAR coordinates
   are used for non-PAR chrX/Y splitting.
3. **Builds four reference atlases** (`BUILD_ATLAS` + `SMOOTH_ATLAS`):
   - **autosomal** — chr1-22 from all CTRLs, written as
     `ref_atlas_{150,500}.tsv.gz`.
   - **male_X** — non-PAR chrX from male CTRLs.
   - **female_X** — non-PAR chrX from female CTRLs.
   - **Y** — full chrY from male CTRLs.
   If a CASE requires a sex-specific comparison and too few matching CTRLs are
   available, the pipeline exits with a clear error.
4. **Detects CpG panels per scope**: autosomal, `male_X`, and `Y` use the
   standard CORE (REF_METH < 10%) + SHORE (REF_METH < 30%) rules. `female_X`
   uses a hemi-methylated island rule (`sex_f_ref_meth_min` to
   `sex_f_ref_meth_max`, default 40-60%) to model X-inactivation.
5. **Defines ICRs from the autosomal atlas** by selecting 40-60%
   hemi-methylated regions inside imprinting gene candidates from
   `assets/icr_gene_list.bed`. If `--icr_panel` is provided, that ICR BED is
   used directly.
6. **Calls DMRs against sex-appropriate component atlases**:

   | Scope | `diff` threshold | `min_sites` | REF_METH gate |
   |---|---|---|---|
   | `autosomal` | `dmr_cpg_diff_threshold` (0.40) | `dmr_cpg_min_sites` (30) | none |
   | `male_X` | `dmr_cpg_diff_threshold` (0.40) | `dmr_cpg_min_sites` (30) | none |
   | `female_X` | `sex_f_diff_threshold` (0.25) | `sex_f_min_sites` (5) | none |
   | `Y` | `dmr_cpg_diff_threshold` (0.40) | `dmr_cpg_min_sites` (30) | none |

   CpG and ICR DMRs are both called by default. ICR calling is autosomal-only.
   Male CASE samples run `autosomal + male_X + Y`; female CASE samples run
   `autosomal + female_X`. Results are merged per sample by `ANNOTATE_DMRS`,
   tagged with `PANEL` and `PANEL_SCOPE`, and published as a single
   `dmr_{sample}_case.tsv`.

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
