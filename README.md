<p align="center">
  <picture>
    <img alt="DMRdetector" src="https://img.shields.io/badge/DMRdetector-v2.2-1f6feb?style=for-the-badge">
  </picture>
</p>

<h1 align="center">DMRdetector</h1>

<p align="center">
  <strong>Quick detection of differentially methylated regions for CpG islands<br>
  and imprinted CpGs based on an ONT methylation reference atlas.</strong>
</p>

<p align="center">
  <em>A three-stage Nextflow pipeline for Oxford Nanopore (ONT) whole-genome
  sequencing data. Builds a population methylation reference atlas from
  controls, defines CpG islands and imprinted control regions (ICRs)
  directly from the atlas, and calls per-patient DMRs with
  locus-type-specific thresholds.</em>
</p>

<p align="center">
  <a href="CHANGELOG.md"><img src="https://img.shields.io/badge/version-2.2-1f6feb.svg" alt="Version"></a>
  <a href="https://www.nextflow.io/"><img src="https://img.shields.io/badge/Nextflow-%E2%89%A524.10-2bb673.svg" alt="Nextflow"></a>
  <a href="https://docs.conda.io/"><img src="https://img.shields.io/badge/Conda-bundled-44a833.svg" alt="Conda"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-eab308.svg" alt="License"></a>
  <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20WSL2%20%7C%20macOS-6b7280.svg" alt="Platform">
</p>

---

## At a glance

|  | DMRdetector v2.2 |
|---|---|
| **Input** | ONT methylation BEDs (5-col, modkit, or pb-CpG-tools) + phased VCFs |
| **Atlas method** | Beta-Binomial dispersion + empirical Bayes shrinkage + FFT-based Gaussian smoothing |
| **Atlas format** | Tabix-indexable bgzipped BED (`ref_atlas_*.bed.gz`); panel-restricted region queries cut atlas I/O 50–100× |
| **CpG islands** | CORE/SHORE classification from population methylation (no GC-content threshold) |
| **ICRs** | Hemi-methylated sites (50% ± 10%) within 122 imprinted genes ± 5 kb |
| **DMR test** | One-sample t-test on per-CpG `DIFF`, BH FDR per panel |
| **Sex chromosomes** | Ploidy-aware partition: autosomal (chr1–22 + PAR pooled); non-PAR male_X / female_X / Y atlases; female_X uses hemi-methylated (40–60%) X-inactivation model |
| **Annotation** | `EnsDb.Hsapiens.v86` (offline GRCh38; biotype-prioritised) |
| **Runtime** | Single conda env, ~5–15 min to install, no Docker/Apptainer needed |

---

## Quick start

```bash
git clone https://github.com/oumarkme/DMRdetector.git
cd DMRdetector

bash install.sh                    # creates the `dmrdetector` conda env
conda activate dmrdetector
bash bin/check_install.sh          # all green? you're ready

nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results
```

Outputs land under `results/`. Subsequent runs can reuse the atlas and
panels — see [Reusing references](#reusing-references-for-faster-reruns).

> **No control cohort?** Skip the atlas build entirely and call DMRs against
> the published 63-genome reference on Zenodo
> ([10.5281/zenodo.18978867](https://doi.org/10.5281/zenodo.18978867)) —
> see [Download the published reference atlas](#download-the-published-reference-atlas-zenodo).

> **No conda yet?** Install [Miniforge](https://github.com/conda-forge/miniforge#install)
> (recommended — bundles `mamba` and uses `conda-forge` by default), then
> rerun `bash install.sh`.

---

## Table of contents

1. [Why DMRdetector](#why-dmrdetector)
2. [Pipeline overview](#pipeline-overview)
3. [Installation](#installation)
4. [Input format](#input-format)
5. [Running the pipeline](#running-the-pipeline)
6. [Output](#output)
7. [Reusing references for faster reruns](#reusing-references-for-faster-reruns)
8. [Sex-chromosome handling](#sex-chromosome-handling)
9. [Methodology in detail](#methodology-in-detail)
10. [Parameters reference](#parameters-reference)
11. [Troubleshooting](#troubleshooting)
12. [Citation](#citation)
13. [License](#license)

---

## Why DMRdetector

DNA methylation aberrations cause **imprinting disorders** (Beckwith-Wiedemann,
Silver-Russell, Prader-Willi, Angelman), **episignature syndromes**, and
**cancer**. ONT long-read sequencing detects 5-methylcytosine in native DNA
without bisulfite conversion, but the standard ONT analysis stops at
per-sample methylation calls — there is no convenient way to detect DMRs
against a long-read population reference.

DMRdetector closes that gap. It combines three capabilities that no other
publicly available tool offers in one workflow:

1. **Population-level reference atlas construction** — coverage-weighted
   per-CpG aggregation across controls, with Beta-Binomial dispersion
   estimation, empirical Bayes shrinkage, and FFT-based Gaussian kernel
   smoothing.
2. **De novo CpG island and ICR identification** — feature boundaries are
   anchored in observed methylation rather than sequence composition,
   capturing functional regions that classical CG-content rules miss and
   excluding sequence-defined islands that are methylated in whole blood.
3. **Single-sample DMR detection with locus-type-specific thresholds** —
   CpG islands and ICRs are tested with biology-appropriate cutoffs in a
   single per-CASE TSV.

The same pipeline scales from a single-patient diagnostic run (with a
pre-built atlas) to a full cohort atlas build.

---

## Pipeline overview

```
Preprocessing (all samples)
    Heterozygous-SNP filter (bcftools + bedtools)
    Minimum read-coverage filter
    Auto-format BED → 5-col uniform representation
        │
        ├──► CTRL stream
        │       │
        │       Stage 1 — Reference atlas
        │           ├─ Per-site aggregation (μ, COV, N̄, S²)
        │           ├─ Beta-Binomial dispersion estimation
        │           ├─ Empirical Bayes shrinkage  (φ = (1−w)·φ_raw + w·φ_trend)
        │           └─ FFT-based Gaussian smoothing (σ = 150 + 500 bp, parallel)
        │       ⇒ ref_atlas_{150,500}.bed.gz (bgzipped, tabix-indexable)
        │
        │       Stage 2 — Feature panels
        │           ├─ CORE  (REF_METH < 10%, φ < 0.05)
        │           ├─ SHORE (REF_METH < 30%, ≤ 150 bp from CORE)
        │           ├─ Cluster, length & density filter, ±150 bp boundary extension
        │           └─ ICRs  (50% ± 10% in 122 imprinted genes ± 5 kb)
        │       ⇒ cpg_signal.bed, icr_region.bed
        │
        └──► CASE stream — Stage 3 — Targeted DMRs
                ├─ Sample-vs-atlas merge → per-CpG DIFF
                ├─ Threshold + cluster (gap 150 bp)
                ├─ One-sample t-test on per-CpG DIFFs, BH FDR per panel
                └─ Annotate with nearest gene (EnsDb.Hsapiens.v86)
            ⇒ dmr_<sample>_case.tsv
```

When `--include_sex_chromosomes true` (default), Stages 1 and 2 run per scope
(`autosomal`, `male_X`, `female_X`, `Y`). The autosomal scope pools chr1–22
plus PAR1 + PAR2 (PAR is biologically diploid in both sexes; reads on chrY-PAR
are translated to chrX-PAR coordinates). The sex-chromosome scopes are strictly
**non-PAR**: `male_X` and `female_X` are non-PAR chrX, `Y` is non-PAR chrY.
The `female_X` atlas uses a hemi-methylated CpG-island model (REF_METH 40–60%)
to capture X-inactivation; the other scopes use the autosomal hypomethylated
model.

---

## Installation

DMRdetector ships as a single conda environment. The bundled `environment.yml`
includes Nextflow itself plus every binary and R/Bioconductor package the
pipeline invokes — no Docker, no Apptainer, no Ensembl downloads.

```bash
bash install.sh                    # create
bash install.sh --update           # refresh after pulling pipeline updates
bash install.sh --force            # remove and recreate
```

`install.sh` auto-detects `mamba`, `micromamba`, or `conda` and uses the
fastest solver available. After install:

```bash
conda activate dmrdetector
bash bin/check_install.sh          # verifies every CLI + R package
```

### Bundled dependencies

| Tool | Min version | Purpose |
|------|-------------|---------|
| `nextflow`                          | 24.10 | Workflow runner |
| `openjdk`                           | 17    | Java runtime for Nextflow |
| `r-base`                            | 4.4   | R interpreter |
| `r-data.table`                      | 1.17  | Atlas / DMR table operations |
| `r-r.utils`                         | 2.13  | gzip support for `data.table::fwrite` |
| `r-readxl`                          | 1.4   | xlsx samplesheet helpers |
| `bioconductor-ensembldb`            | 2.30  | Gene annotation queries |
| `bioconductor-ensdb.hsapiens.v86`   | 2.99  | GRCh38 gene database (offline) |
| `bioconductor-genomicranges`        | 1.58  | Range overlap |
| `bcftools`                          | 1.20  | VCF parsing |
| `bedtools`                          | 2.31  | Interval intersection |
| `htslib`                            | 1.20  | bgzip / tabix |

### System requirements

- **OS:** Linux, WSL2, or x86_64 macOS.
  (`EnsDb.Hsapiens.v86` has no native Apple Silicon build; use Rosetta on M-series Macs.)
- **RAM:** ≥ 16 GB for ~10-sample cohorts; ≥ 64 GB for full hg38 autosomal
  smoothing of ~30 M CpG sites.
- **Disk:** ≥ 30 GB for atlas TSVs and Nextflow work cache (per cohort).

---

## Input format

The samplesheet is a 5-column CSV with this exact header:

```csv
sample_id,type,sex,meth_bed,vcf
ONT001,CASE,M,/path/to/ONT001.combined.bed.gz,/path/to/ONT001.phased.vcf
ONT004,CTRL,F,/path/to/ONT004.combined.bed.gz,/path/to/ONT004.phased.vcf
```

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | string | Unique identifier; appears in output filenames. |
| `type`      | enum   | `CASE` or `CTRL`. CTRLs build the atlas; CASEs are tested. |
| `sex`       | enum   | `M` or `F`. Routes chrX/chrY beds to the correct scope. |
| `meth_bed`  | path   | Methylation BED (gzipped). Auto-detects format. |
| `vcf`       | path   | Phased VCF used to mask heterozygous SNP positions. |

The header, missing/duplicate IDs, and `type`/`sex` values are validated
before any process runs — bad samplesheets fail in under a second.

### Methylation BED formats (auto-detected)

| Format | Cols | Coverage column | Modified-read column |
|--------|------|-----------------|----------------------|
| 5-col plain                    | 5    | col 4 | col 5 |
| modkit `pileup` (bedMethyl)    | 18   | col 10 (`Nvalid_cov`) | col 12 (`Nmod`) |
| pb-CpG-tools `combined`        | 7–9 (with header) | col 6 | col 7 |

All formats are normalised to a uniform 5-column representation during
preprocessing.

---

## Running the pipeline

### Default — build atlas and run DMRs

```bash
nextflow run main.nf \
    --input  samplesheet.csv \
    --outdir results
```

Both autosomal and sex-stratified analyses run. Sex chromosomes are on by
default; pass `--include_sex_chromosomes false` for chr1–22 only.

### Common flags

```bash
nextflow run main.nf \
    --input    samplesheet.csv \
    --outdir   results \
    --tmp_dir  /raid/nf/dmrdetector \   # local scratch (NAS-friendly)
    --max_cpus 16 \                      # cap concurrency
    --max_memory '64.GB' \
    -ansi-log false \                    # readable in piped logs
    -with-report results/nf_report.html
```

### Resume after failure

Nextflow caches every completed task. Add `-resume` to restart only failed
and downstream tasks:

```bash
nextflow run main.nf --input samplesheet.csv --outdir results -resume
```

### Test profile

```bash
nextflow run main.nf -profile test --input samplesheet.csv --outdir results
```

Caps `max_cpus=4`, `max_memory=16.GB`, `max_time=1.h` — useful for CI or
quick local checks.

---

## Output

```
results/
├── atlas/
│   ├── ref_atlas_150.bed.gz             # Autosomal atlas (σ=150, primary)
│   ├── ref_atlas_500.bed.gz             # Autosomal atlas (σ=500, regional)
│   ├── ref_atlas_150_male_X.bed.gz      # Male non-PAR chrX
│   ├── ref_atlas_150_female_X.bed.gz    # Female non-PAR chrX (hemi-meth)
│   ├── ref_atlas_150_Y.bed.gz           # Male chrY
│   └── (matching σ=500 files)
├── cpg_islands/
│   ├── cpg_signal.bed                   # Autosomal CpG panel
│   ├── cpg_signal_male_X.bed
│   ├── cpg_signal_female_X.bed
│   └── cpg_signal_Y.bed
├── icr_regions/
│   └── icr_region.bed                   # ICR panel (autosomal)
├── dmr/
│   └── dmr_<sample_id>_case.tsv         # Per-CASE annotated DMRs
└── nf_report.html                       # Nextflow execution report
```

### Reference atlas columns (`ref_atlas_*.bed.gz`)

The atlas is a headerless bgzipped BED with 9 columns; column names are
supplied by downstream consumers. Compatible with `tabix -R panel.bed`
for fast random access (`-p bed`).

`CHROM, START, END, N, COV, PHI, RAW_METH, REF_METH, REF_VAR`

- `N`        — number of contributing control samples at the site.
- `COV`      — total read coverage pooled across those controls.
- `PHI`      — shrunken inter-individual dispersion (empirical Bayes).
- `RAW_METH` — coverage-weighted mean methylation **before** smoothing (0–1).
- `REF_METH` — coverage-weighted, Gaussian-smoothed mean methylation (0–1).
- `REF_VAR`  — expected biological variance under Beta-Binomial
               (`REF_METH·(1−REF_METH)·PHI`).

### DMR TSV columns (`dmr_<sample>_case.tsv`)

| Column | Description |
|--------|-------------|
| `CHROM, START, END`  | Region coordinates (chr-prefixed, hg38). |
| `GENE`               | Nearest gene from `EnsDb.Hsapiens.v86` (biotype-prioritised). |
| `PANEL`              | `CpG` or `ICR`. |
| `NUM_CPG`            | Number of CpG sites contributing to the region. |
| `NUM_READS`          | Mean per-CpG read coverage. |
| `SAMPLE_METH (%)`    | Mean methylation in this CASE. |
| `REF_METH (%)`       | Mean methylation in the reference atlas. |
| `DIFF (%)`           | `SAMPLE_METH − REF_METH`. |
| `P_VALUE`            | One-sample t-test on per-CpG `DIFF` within the region. |
| `P_FDR`              | Benjamini-Hochberg adjusted, per panel. |
| `DIRECTION`          | `hyper` (DIFF > 0) or `hypo` (DIFF < 0). |

### Default locus-type-specific thresholds

| Panel | `\|DIFF\|` | `NUM_CPG` | `FDR` |
|---|---|---|---|
| **CpG island** | > 40 % | > 30 | < 0.01 |
| **ICR**        | > 25 % | > 5  | < 0.01 |
| **female_X** (hemi-methylated CpG model) | > 25 % | > 5 | < 0.01 |

ICR thresholds are deliberately relaxed: imprinted loci are hemi-methylated
(baseline ≈ 50%), so allele-specific methylation changes of smaller magnitude
are biologically meaningful.

---

## Reusing references for faster reruns

Stage 1 (atlas) and Stage 2 (CpG / ICR panels) are the slow steps. Once
they exist, **subsequent runs can skip them** and only recompute the
per-CASE DMR calls.

### Download the published reference atlas (Zenodo)

Most diagnostic users have **no in-house control cohort** — and don't need
one. A ready-to-use reference set built from **69 ONT control whole-blood
genomes** is published on Zenodo and drops straight into the reuse flags:

> 📦 **[10.5281/zenodo.18978867](https://doi.org/10.5281/zenodo.18978867)**
> — pre-built σ = 150 bp and σ = 500 bp atlases · 28,624,561 autosomal CpG
> sites · 25,719 CpG islands · 445 ICRs across 122 imprinted genes (hg38).
> MIT-licensed, free to redistribute.

Download the record and unpack it into a `reference/` directory. The archive
mirrors the pipeline's own output layout, so the filenames match the
`--ref_atlas` / `--cpg_panel` / `--icr_panel` flags one-to-one:

```
reference/
├── ref_atlas_150.bed.gz            ref_atlas_150_male_X.bed.gz
├── ref_atlas_150_female_X.bed.gz   ref_atlas_150_Y.bed.gz
├── cpg_signal.bed                  cpg_signal_male_X.bed
├── cpg_signal_female_X.bed         cpg_signal_Y.bed
└── icr_region.bed
```

Then run a **CASE-only** samplesheet against it — Stages 1 and 2 are bypassed
and only per-CASE DMRs are computed:

```bash
nextflow run main.nf \
    --input  cases.csv \
    --outdir results \
    --ref_atlas          reference/ref_atlas_150.bed.gz \
    --ref_atlas_male_X   reference/ref_atlas_150_male_X.bed.gz \
    --ref_atlas_female_X reference/ref_atlas_150_female_X.bed.gz \
    --ref_atlas_Y        reference/ref_atlas_150_Y.bed.gz \
    --cpg_panel          reference/cpg_signal.bed \
    --cpg_panel_male_X   reference/cpg_signal_male_X.bed \
    --cpg_panel_female_X reference/cpg_signal_female_X.bed \
    --cpg_panel_Y        reference/cpg_signal_Y.bed \
    --icr_panel          reference/icr_region.bed
```

The published atlas is hg38/GRCh38 only. To call DMRs against your **own**
freshly built references instead, use the same flags with the paths emitted
under `results/` — see below.

### Sex-stratified reuse (recommended)

```bash
nextflow run main.nf \
    --input  samplesheet.csv \
    --outdir results \
    --ref_atlas          results/atlas/ref_atlas_150.bed.gz \
    --ref_atlas_male_X   results/atlas/ref_atlas_150_male_X.bed.gz \
    --ref_atlas_female_X results/atlas/ref_atlas_150_female_X.bed.gz \
    --ref_atlas_Y        results/atlas/ref_atlas_150_Y.bed.gz \
    --cpg_panel          results/cpg_islands/cpg_signal.bed \
    --cpg_panel_male_X   results/cpg_islands/cpg_signal_male_X.bed \
    --cpg_panel_female_X results/cpg_islands/cpg_signal_female_X.bed \
    --cpg_panel_Y        results/cpg_islands/cpg_signal_Y.bed \
    --icr_panel          results/icr_regions/icr_region.bed
```

When all three sex-specific atlases are supplied alongside the autosomal
one, Stages 1 and 2 are bypassed completely while CASE samples are still
sex-stratified into autosomal/male_X/female_X/Y for DMR calling.

### Autosomal-only reuse

```bash
nextflow run main.nf \
    --input  samplesheet.csv \
    --outdir results \
    --ref_atlas results/atlas/ref_atlas_150.bed.gz \
    --cpg_panel results/cpg_islands/cpg_signal.bed \
    --icr_panel results/icr_regions/icr_region.bed \
    --include_sex_chromosomes false
```

> **Tip:** No control cohort of your own? Use the published 63-genome atlas —
> see [Download the published reference atlas](#download-the-published-reference-atlas-zenodo).

---

## Sex-chromosome handling

DMRdetector partitions the genome into four ploidy-aware compartments to
handle X-inactivation, male hemizygosity, and the diploid behaviour of PAR.
With `--include_sex_chromosomes true` (default):

1. Each preprocessed BED is split using bundled hg38 PAR coordinates
   (`assets/par_hg38.bed`) into:
   - **autosomal** — chr1–22 PLUS PAR1 + PAR2 (PAR is biologically
     diploid in both sexes; reads on chrY-PAR are translated to chrX-PAR
     coordinates so the same locus is counted once).
   - **non-PAR chrX** (split downstream into male_X / female_X by sex).
   - **non-PAR chrY** (male-only).
2. CTRL components feed four reference atlases:
   - **autosomal** — all CTRLs (M and F), pooled PAR included.
     Standard CORE/SHORE hypomethylated model.
   - **male_X**    — male CTRLs only; autosomal-style hypomethylated
                    (REF_METH < 10 %).
   - **female_X**  — female CTRLs only; **hemi-methylated** model
                    (REF_METH 40–60 %) capturing X-inactivation.
   - **Y**         — male CTRLs only (non-PAR chrY); exploratory.
3. CASE samples are routed by the `sex` column:
   - **Male CASE**   → autosomal + male_X + Y.
   - **Female CASE** → autosomal + female_X (no Y).
4. ICR analysis remains autosomal-only by design (the curated imprinting
   gene list spans chr1–22).

The minimum same-sex CTRLs required per scope is `--sex_atlas_min_controls`
(default 3); cohorts below threshold for one scope skip that scope and emit
a warning.

To restrict the run to autosomes only:

```bash
nextflow run main.nf --input samplesheet.csv --outdir results \
    --include_sex_chromosomes false
```

---

## Methodology in detail

The descriptions below mirror the *Materials and Methods* section of the
manuscript. The pipeline supports hg38 / GRCh38 only.

### Preprocessing

Heterozygous SNP positions are extracted from each sample's VCF using
`bcftools` and converted to BED; CpG sites overlapping any heterozygous SNP
are removed using `bedtools` to **eliminate allele-specific methylation
artifacts caused by C-to-T polymorphisms** at or near CpG dinucleotides.
A per-sample read-coverage filter is applied — `> 5×` for atlas building
and `> 10×` for DMR detection. The pipeline uses the lower threshold for
a single preprocessing pass that satisfies both stages.

### Stage 1 — Reference atlas construction

1. **Per-site aggregation.** For each CpG site: number of contributing
   samples (`N`), total coverage (`COV`), total methylated read count,
   sample variance of per-sample methylation ratios (`S²`), mean per-sample
   coverage (`N̄ = COV / N`). Sites with fewer than three contributing
   samples are excluded. The aggregate reference methylation
   `μ = methylated_reads / COV` weighs each sample proportionally to its
   sequencing depth.
2. **Beta-Binomial dispersion estimation** quantifies inter-individual
   methylation variability via a moment-based estimator for the
   overdispersion parameter φ; negative estimates are zeroed and values
   above 1 are capped (Feng et al., 2014).
3. **Empirical Bayes shrinkage** stabilises low-coverage estimates:
   `φ = (1 − w)·φ_raw + w·φ_trend`, with `w = 0.1` and `φ_trend` =
   genome-wide median of positive `φ_raw`. Modest weight preserves
   site-specific signal while dampening noise (Love et al., 2014;
   Wu et al., 2013; Feng et al., 2014).
4. **FFT-based Gaussian kernel smoothing.** A 1-D Gaussian kernel
   `K(x) = exp(−x²/2σ²)` (truncated at `R = |3σ|`, > 99.7% of mass) is
   convolved per chromosome with dense numerator (`μ × COV`) and
   denominator (`COV`) vectors via FFT; smoothed
   `REF_METH = numerator / denominator` (coverage-weighted Nadaraya–Watson
   regression — avoids the bias of directly smoothing methylation ratios
   at sites with heterogeneous coverage). Two bandwidths run in parallel:
   `σ = 150 bp` (CpG-island scale, used downstream) and `σ = 500 bp`
   (regional trends). FFT cuts the smoother from O(nk) to O(n log n).

### Stage 2a — De novo CpG island identification

Sites are first filtered to retain only those with **more than five
contributing control samples** (a stricter threshold than the
three-sample minimum used during atlas construction), then classified:

- **CORE**  — `REF_METH < 10 %` and shrunken dispersion `φ < 0.05`.
              Positions consistently and stably unmethylated across the cohort.
- **SHORE** — `REF_METH < 30 %`, within 150 bp of the nearest CORE site,
              not themselves CORE.

CORE sites are clustered (new cluster when the gap > 150 bp), filtered for
minimum island length ≥ 150 bp and minimum CORE density of 1 site / 100 bp,
then extended ±150 bp to incorporate flanking SHORE sites. The two-tier
classification draws a conceptual parallel with MethylSeekR's UMR framework
(Burger et al., 2013); the difference is that CORE/SHORE here are derived
from **population-level** statistics, providing cohort-consensus boundaries
rather than sample-specific segmentations.

### Stage 2b — Imprinted control regions (ICRs)

Imprinted loci are hemi-methylated (≈ 50%) and are therefore excluded from
the hypomethylated CpG island set. A separate panel is built: a curated
list of **122 imprinted genes** from geneimprint.com is intersected with
Ensembl GRCh38 coordinates and extended ±5 kb. Within each region, CpG
sites with **REF_METH 50% ± 10%** (`[0.40, 0.60]`) and inter-site distances
< 150 bp are grouped into ICRs.

### Stage 3 — Targeted DMR detection

Preprocessed CpG sites from each CASE are joined with the reference atlas
by an inner merge on coordinates; `DIFF = SAMPLE_METH − REF_METH` is
computed per site. Sites below the panel-specific `|DIFF|` cutoff are
removed; remaining sites are clustered into contiguous regions (gap
threshold 150 bp). Each region is tested with a **one-sample t-test** on
the per-CpG `DIFF` values within that region; regions with fewer than two
sites get NA p-values. **Benjamini-Hochberg FDR is applied per panel
independently.** The single-sample testing approach follows
Rezwan et al. (2015) and Grolaux et al. (2022).

### Gene annotation

Each DMR is annotated with the nearest gene symbol from the offline
`EnsDb.Hsapiens.v86` Bioconductor package. When multiple genes overlap a
DMR, biotype priority is applied (protein_coding → IG/TR genes →
lncRNA/antisense → miRNA → sno/snRNA → rRNA → misc_RNA →
processed_transcript → pseudogenes), with overlap-width and gene-length
tie-breaks. Genes lacking a `gene_name` attribute fall back to their
Ensembl `gene_id`.

---

## Parameters reference

<details>
<summary><strong>Input / output</strong></summary>

| Parameter   | Default     | Description |
|-------------|-------------|-------------|
| `input`     | required    | Samplesheet CSV. |
| `outdir`    | `./results` | Output directory. |
| `tmp_dir`   | null        | Local scratch (recommended for NAS-backed inputs). |

</details>

<details>
<summary><strong>Pre-built references (skip Stages 1 / 2)</strong></summary>

| Parameter            | Description |
|----------------------|-------------|
| `ref_atlas`          | Autosomal atlas (`.bed.gz` preferred; legacy `.tsv.gz` accepted). Skips Stage 1 autosomal. |
| `ref_atlas_male_X`   | Male non-PAR chrX atlas. |
| `ref_atlas_female_X` | Female non-PAR chrX atlas. |
| `ref_atlas_Y`        | chrY atlas. |
| `cpg_panel`          | Autosomal CpG island BED. |
| `cpg_panel_male_X`   | male_X CpG panel. |
| `cpg_panel_female_X` | female_X CpG panel. |
| `cpg_panel_Y`        | Y CpG panel. |
| `icr_panel`          | ICR BED. |
| `icr_gene_list`      | Imprinting candidate gene list (BED4: CHROM/START/END/GENE). Defaults to `assets/icr_gene_list.bed`. |

</details>

<details>
<summary><strong>Stage 1 — Atlas</strong></summary>

| Parameter                    | Default | Description |
|------------------------------|---------|-------------|
| `atlas_min_samples_per_site` | 3       | Min CTRLs covering a site to keep it. |
| `atlas_min_read_coverage`    | 5       | Per-sample coverage filter during preprocessing. |
| `atlas_shrinkage_weight`     | 0.1     | Linear shrinkage `w` toward genome-wide φ trend. |
| `smooth_sigma_narrow`        | 150     | Primary FFT-Gaussian σ (CpG-island scale). |
| `smooth_sigma_wide`          | 500     | Secondary σ (regional trend). |

</details>

<details>
<summary><strong>Stage 2 — CpG islands</strong></summary>

| Parameter                   | Default | Description |
|-----------------------------|---------|-------------|
| `cpg_min_samples_per_site`  | 5       | Stricter sample threshold for island calling. |
| `cpg_core_max_methylation`  | 0.1     | CORE upper REF_METH bound. |
| `cpg_core_max_dispersion`   | 0.05    | CORE upper φ bound. |
| `cpg_shore_max_methylation` | 0.3     | SHORE upper REF_METH bound. |
| `cpg_max_gap_distance`      | 150     | Max gap (bp) before a new cluster. |
| `cpg_min_island_length`     | 150     | Minimum island length (bp). |

</details>

<details>
<summary><strong>Stage 2 — ICR panel</strong></summary>

| Parameter            | Default | Description |
|----------------------|---------|-------------|
| `icr_ref_meth_min`   | 0.40    | Lower REF_METH bound (50% − 10%). |
| `icr_ref_meth_max`   | 0.60    | Upper REF_METH bound (50% + 10%). |
| `icr_diff_threshold` | 0.25    | DMR `|DIFF|` cutoff for ICR panel. |
| `icr_min_sites`      | 5       | Min CpG count per ICR DMR. |

</details>

<details>
<summary><strong>Stage 3 — DMR detection</strong></summary>

| Parameter                | Default | Description |
|--------------------------|---------|-------------|
| `dmr_min_coverage`       | 10      | Per-sample coverage threshold for CpGs entering DMR test. |
| `dmr_cpg_diff_threshold` | 0.40    | DMR `|DIFF|` cutoff for CpG panel. |
| `dmr_merge_distance`     | 150     | Max gap (bp) before clustering breaks. |
| `dmr_fdr_threshold`      | 0.01    | Benjamini-Hochberg cutoff. |
| `dmr_cpg_min_sites`      | 30      | Min CpG count per CpG-panel DMR. |

</details>

<details>
<summary><strong>Sex chromosomes</strong></summary>

| Parameter                  | Default     | Description |
|----------------------------|-------------|-------------|
| `include_sex_chromosomes`  | true        | Master toggle. |
| `par_bed`                  | bundled hg38 | Override PAR coordinates. |
| `sex_atlas_min_controls`   | 3           | Min same-sex CTRLs needed for a sex-specific atlas. |
| `sex_f_ref_meth_min`       | 0.40        | Female-X hemi-methylated island lower bound. |
| `sex_f_ref_meth_max`       | 0.60        | Female-X upper bound. |
| `sex_f_diff_threshold`     | 0.25        | DMR cutoff on female_X. |
| `sex_f_min_sites`          | 5           | Min CpGs per female_X DMR. |

</details>

<details>
<summary><strong>Resources</strong></summary>

`max_cpus` and `max_memory` auto-detect to 80 % of the host (CPU capped to
1 core per 4 GB RAM). Override with `--max_cpus N --max_memory '32.GB'`.

</details>

---

## Troubleshooting

> **`samplesheet ... is missing required column(s): sex`**
> The samplesheet must include all 5 columns including `sex`. See
> [Input format](#input-format).

> **`include_sex_chromosomes=true but --ref_atlas is provided without --ref_atlas_male_X / _female_X / _Y`**
> Pass all three sex-specific atlases or accept the autosomal-only fallback.
> See [Reusing references](#reusing-references-for-faster-reruns).

> **Autosomal `SMOOTH_ATLAS` runs slowly**
> Full hg38 covers ~30 M CpG sites; FFT smoothing per chromosome needs
> ≥ 30 GB RAM and several minutes per σ. Use `--max_cpus` to control
> parallelism and `--tmp_dir` on a fast local disk.

> **`Cannot find script file: main.nf`**
> Run `nextflow run` from the repo root, or pass the full path:
> `nextflow run /path/to/DMRdetector/main.nf …`.

> **`-resume` does not pick up a hand-edited atlas**
> Nextflow caches by input file hash. If you edit a pre-built reference in
> place, also remove the work cache directory (`rm -rf <tmp_dir>/work`).

---

## Citation

If you use DMRdetector, please cite:

> Ou J-H\*, Tsai M-JM\*, Hsu R-H, Chen H-A, Wu C-S, Chao M-Y, Hsu C,
> Chien Y-H, Lee N-C, Hwu W-L. **Detecting differentially methylated regions
> at CpG islands and imprinted loci with an ONT methylation atlas.**
> Department of Medical Genetics, National Taiwan University Hospital,
> Taipei, Taiwan. (Manuscript in preparation, 2026.)
> *\*These authors contributed equally.*

The pipeline and pre-built atlases are freely available at
<https://github.com/oumarkme/DMRdetector> under the MIT licence. A
pre-built reference atlas constructed from 69 ONT control whole-blood
genomes (28,624,561 autosomal CpG sites; 25,719 CpG islands; 445 ICRs
across 122 imprinted genes) is archived at Zenodo:
[10.5281/zenodo.18978867](https://doi.org/10.5281/zenodo.18978867).

### Selected method references

- Burger L. *et al.* MethylSeekR. *Nucleic Acids Research* (2013).
- Feng H. *et al.* DSS / Beta-Binomial dispersion. *Bioinformatics* (2014).
- Gardiner-Garden M. & Frommer M. CpG island sequence definition.
  *J. Mol. Biol.* (1987).
- Grolaux R. *et al.* Single-patient DMR detection. *Clinical Epigenetics* (2022).
- Hansen K.D. *et al.* BSmooth. *Genome Biology* (2012).
- Love M.I. *et al.* DESeq2 dispersion shrinkage. *Genome Biology* (2014).
- Peters T.J. *et al.* DMRcate. *Genome Biology* (2015).
- Rezwan F.I. *et al.* Single-case methylation diagnosis.
  *Clinical Epigenetics* (2015).
- Sigurpalsdottir B. *et al.* deCODE 7,179-sample ONT methylation atlas.
  *Nature Genetics* (2024).

---

## License

[MIT](LICENSE) © Jen-Hsiang Ou, Department of Medical Genetics, National Taiwan University Hospital.
