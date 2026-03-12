# nf-ontmeth

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Container](https://img.shields.io/badge/Container-GHCR-blue.svg)](https://ghcr.io/oumarkme/nf-ontmeth)
[![CI](https://github.com/oumarkme/nf-ontmeth/actions/workflows/docker.yml/badge.svg)](https://github.com/oumarkme/nf-ontmeth/actions)

ONT Methylation Analysis Pipeline — detects differentially methylated regions (DMRs) from Oxford Nanopore Technologies (or PacBio) whole-genome sequencing data. Supports an optional user-provided second panel for targeted region analysis.

## Pipeline Architecture

| Stage | Process | Description |
|-------|---------|-------------|
| 1 | **Build Atlas** | Aggregate control samples, estimate dispersion, FFT Gaussian smoothing |
| 2 | **Detect CpG Islands** | Identify CORE/SHORE sites, cluster, filter, extend |
| 3 | **Detect DMRs** | Compare each sample vs reference within CpG panel + optional panel2 |

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04
- [Apptainer](https://apptainer.org/) (Linux/WSL) or [Docker](https://www.docker.com/) (macOS)
- **OS:** Linux, WSL2, or macOS

## Installation

```bash
# Clone the repository
git clone https://github.com/oumarkme/nf-ontmeth.git
cd nf-ontmeth
```

The container image is hosted on [GitHub Container Registry](https://ghcr.io/oumarkme/nf-ontmeth). Apptainer and Docker pull it automatically on first run — no manual build needed.

The image is auto-built via GitHub Actions whenever the `Dockerfile` changes on `main` or a new release is published.

### Building the container locally (optional)

```bash
docker build -t ghcr.io/oumarkme/nf-ontmeth:latest .
# Convert to Apptainer SIF (optional, for offline use)
apptainer build nf-ontmeth.sif docker://ghcr.io/oumarkme/nf-ontmeth:latest
```

## Quick Start

### Full pipeline (build atlas + detect DMRs)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results
```

The GTF annotation is auto-downloaded from Ensembl (latest GRCh38 release). To use a specific GTF file:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --gtf Homo_sapiens.GRCh38.115.gtf.gz \
    --outdir results
```

### With optional second panel

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --panel2 my_regions.bed \
    --outdir results
```

### Build atlas only (no CASE samples in samplesheet)

```bash
nextflow run main.nf \
    --input controls_only.csv \
    --outdir results
```

### DMR detection only (pre-built atlas)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --ref_atlas results/atlas/ref_atlas_150.tsv.gz \
    --cpg_panel results/cpg_islands/cpg_signal.bed \
    --outdir results
```

## Container and Bind Mounts

The pipeline runs inside a container (Apptainer or Docker). If your input files live outside the working directory, the container needs access to those paths.

### Apptainer (Linux / WSL2) — default

Use `--apptainer_bind` to bind-mount host directories into the container. Paths that don't exist are silently skipped.

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --gtf annotation.gtf.gz \
    --apptainer_bind '/data,/mnt/nas,/scratch' \
    --outdir results
```

### Docker (macOS)

Use the `docker` profile. Docker Desktop shares `/Users` and `/tmp` by default. Add other paths in Docker Desktop > Settings > Resources > File Sharing.

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --gtf annotation.gtf.gz \
    --outdir results
```

## Input: Samplesheet

CSV file with columns:

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

## Input BED Formats (auto-detected)

The pipeline auto-detects the input BED format based on column count:

| Format | Columns | COVERAGE | N_MOD |
|--------|---------|----------|-------|
| **5-column** (default) | CHROM, START, END, COVERAGE, N_MOD | col 4 | col 5 |
| **modkit pileup** (bedMethyl) | 18 columns | col 10 (Nvalid_cov) | col 12 (Nmod) |
| **pb-CpG-tools** | 7-9 columns (with header) | col 6 (cov) | col 7 (est_mod_count) |

## Optional Panel2

A headerless 3-column BED file (CHROM, START, END) specifying genomic regions for a second panel analysis. Gene annotation is performed via GTF overlap (same as CpG panel).

## Reference Files

The pipeline uses an Ensembl GTF annotation file for gene annotation. If `--gtf` is not provided, the **latest GRCh38 GTF is auto-downloaded** from Ensembl.

To use a specific version or GRCh37, download manually:

```bash
# GRCh38 (hg38) — Ensembl release 115
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

# GRCh37 (hg19) — Ensembl release 115
wget https://ftp.ensembl.org/pub/grch37/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
```

Then pass with `--gtf Homo_sapiens.GRCh38.115.gtf.gz`. Use the GTF version matching your alignment reference genome.

## Parameters

### Reference Atlas (Stage 1)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `atlas_min_samples_per_site` | 3 | Minimum samples per CpG site |
| `atlas_min_read_coverage` | 5 | Minimum read coverage |
| `atlas_shrinkage_weight` | 0.1 | Empirical Bayes shrinkage weight |
| `smooth_sigma_narrow` | 150 | Narrow Gaussian bandwidth (bp) |
| `smooth_sigma_wide` | 500 | Wide Gaussian bandwidth (bp) |

### CpG Island Detection (Stage 2)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cpg_min_samples_per_site` | 5 | Minimum samples for pre-filtering |
| `cpg_core_max_methylation` | 0.1 | CORE site REF_METH threshold |
| `cpg_core_max_dispersion` | 0.05 | CORE site PHI threshold |
| `cpg_shore_max_methylation` | 0.3 | SHORE site REF_METH threshold |
| `cpg_max_gap_distance` | 150 | Max gap for clustering (bp) |
| `cpg_min_island_length` | 150 | Minimum island length (bp) |

### DMR Detection (Stage 3)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dmr_min_coverage` | 10 | Minimum coverage for samples |
| `dmr_cpg_diff_threshold` | 0.40 | \|DIFF\| threshold for CpG panel |
| `dmr_merge_distance` | 150 | Max gap for merging sites (bp) |
| `dmr_fdr_threshold` | 0.01 | FDR significance threshold |
| `dmr_cpg_min_sites` | 30 | Min CpG sites for CpG panel |
| `panel2_diff_threshold` | 0.25 | \|DIFF\| threshold for panel2 |
| `panel2_min_sites` | 5 | Min CpG sites for panel2 |

### Resource Management

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_cpus` | auto | Maximum total CPUs (see below) |
| `max_memory` | auto | Maximum total memory (see below) |
| `max_time` | 168.h | Maximum time per process |
| `apptainer_bind` | null | Comma-separated host paths to bind-mount |
| `tmp_dir` | null | Custom temporary directory |

**Auto-detection:** Resources are auto-detected at startup:
- **Memory:** 80% of system RAM
- **CPUs:** minimum of (80% of cores) and (total RAM / 4 GB)

This prevents out-of-memory on RAM-constrained systems. For example:

| System | Cores | RAM | max_cpus | max_memory |
|--------|-------|-----|----------|------------|
| Laptop | 8 | 16 GB | 3 | 12 GB |
| Workstation | 16 | 64 GB | 12 | 51 GB |
| Server | 96 | 2 TB | 76 | 1.6 TB |

Override with `--max_cpus` and `--max_memory`:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --gtf annotation.gtf.gz \
    --max_cpus 4 \
    --max_memory '16.GB' \
    --outdir results
```

## Output

```
results/
├── atlas/
│   ├── ref_atlas_150.tsv.gz    # Smoothed atlas (narrow)
│   └── ref_atlas_500.tsv.gz    # Smoothed atlas (wide)
├── cpg_islands/
│   └── cpg_signal.bed          # CpG island panel
├── dmr/
│   └── dmr_{sample}_{type}.tsv # DMR results per sample
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.svg
```

### DMR Output Format

| Column | Description |
|--------|-------------|
| CHROM | Chromosome |
| START | Region start |
| END | Region end |
| GENE | Gene symbol |
| PANEL | CpG or Panel2 |
| NUM_CPG | Number of CpG sites |
| NUM_READS | Mean read coverage |
| SAMPLE_METH (%) | Sample methylation percentage |
| REF_METH (%) | Reference methylation percentage |
| DIFF (%) | Difference (sample - reference) |
| P_VALUE | T-test p-value |
| P_FDR | BH-adjusted p-value |
| DIRECTION | hyper or hypo |

## Profiles

| Profile | Description |
|---------|-------------|
| `test` | Reduced resources (4 CPUs, 16 GB) for CI/testing |
| `docker` | Use Docker instead of Apptainer (required for macOS) |

## Resuming Failed Runs

Nextflow caches completed tasks. To resume after a failure or interruption:

```bash
nextflow run main.nf [same options] -resume
```

## Technology

- **Nextflow** (DSL2) — workflow orchestration
- **R** (data.table) — all statistical computation (FFT smoothing, atlas building, DMR detection)
- **bash** — preprocessing (bcftools/bedtools)
- **Apptainer/Docker** — containerized execution

## Citation

If you use this pipeline in your research, please cite:

> [Your manuscript citation here]

## License

This project is licensed under the [MIT License](LICENSE).

Copyright (c) 2025 Mark Ou
