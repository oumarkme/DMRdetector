# Sex-Chromosome Analysis in DMRdetector

> **v1.3 note:** This document records the original v1.2 sex-chromosome
> implementation. It has been superseded by the v1.3 component-scope model:
> `autosomal`, `male_X`, `female_X`, and `Y`. In v1.3 PAR remains with X/Y
> during methylation analysis, CpG island detection uses male-derived scopes,
> female chrX DMRs compare against the `female_X` atlas, and ICR/panel2
> analysis is autosomal only.

## Problem

Before this change, DMRdetector analysed autosomes only. Sex chromosomes
(chrX / chrY) were filtered out because they require statistical treatment
that the autosomal code path cannot model:

- **Ploidy asymmetry.** Males have 1 copy of chrX and 1 copy of chrY;
  females have 2 copies of chrX and no chrY. Pooling both sexes in one
  reference atlas produces meaningless per-site means.
- **X-chromosome inactivation (XCI).** In females, one of the two X
  chromosomes is silenced allele-by-allele across cells. On a
  combined-strand methylation BED (summed across both alleles), X-linked
  promoters appear at ~50% methylation — not because the site is truly
  hemi-methylated, but because roughly half the cells contribute a
  methylated read and half contribute an unmethylated read. Autosome-style
  thresholds (CORE < 10%, SHORE < 30%) mis-classify every such site.
- **Pseudoautosomal regions (PAR).** PAR1 and PAR2 on chrX / chrY are
  diploid in both sexes and recombine like autosomes. They should be
  analysed with the autosomal pool, not the sex-specific pool.

Clinically this matters: MECP2 duplication syndrome, FMR1 CGG-expansion
carriers, ATRX-associated skewed XCI, and several X-linked intellectual
disabilities present as altered methylation on non-PAR chrX. None of these
were detectable with the previous pipeline.

## Solution overview

A sex-stratified workflow is now wired into the pipeline and runs by
default. For every input methylation BED the pipeline:

1. **Infers sample sex from the VCF** by counting heterozygous SNPs on
   non-PAR chrX. Males are hemizygous (HET ≈ 0); females are diploid
   (HET in the thousands). Counts between the two thresholds are flagged
   `U` (unknown — sample is still analysed on autosomes).
2. **Splits each BED into three streams**: autosomal + PAR (diploid,
   pooled), non-PAR chrX (sex-stratified), non-PAR chrY (males only).
3. **Builds three reference atlases** in parallel: autosomal + PAR (all
   controls), sex\_M (male controls, non-PAR X + Y), sex\_F (female
   controls, non-PAR X).
4. **Detects CpG panels per scope**, with a dedicated imprinting-style
   rule for female non-PAR chrX that targets XCI-bimodal sites
   (`REF_METH ∈ [0.35, 0.65]`).
5. **Calls DMRs three times per case sample** with scope-specific
   thresholds, then merges the results into a single annotated TSV.

### Biological mapping

| Region | Biology | Treatment |
| --- | --- | --- |
| chr1–22 | Diploid, sex-agnostic | Pooled autosomal atlas (unchanged) |
| PAR1 / PAR2 | Diploid in both sexes | Treated as autosomal |
| Non-PAR chrX (M) | Hemizygous, promoters ~0% methylated | Male-only atlas, autosome-like thresholds |
| Non-PAR chrX (F) | XCI mosaic, combined ~50% methylated | Female-only atlas, imprinting-style thresholds |
| Non-PAR chrY | Males only | Male-only atlas, autosome-like thresholds; females skip |

The "female non-PAR chrX = imprinting-style" mapping is motivated by the
XCI signature: it mirrors the combined-strand bimodal profile of imprinted
control regions (one parental allele methylated, one unmethylated). The
same detection logic used for ICRs therefore applies directly.

### Pipeline diagram

```
samplesheet.csv  ─▶  PREPROCESS_SAMPLE
                           │
                           ├─▶ INFER_SEX_FROM_VCF (non-PAR chrX HET count)
                           │   └─ meta.sex ∈ {M, F, U}
                           ▼
                    SPLIT_PAR_AUTOSOMAL
                    ├─ autosomal + PAR
                    ├─ non-PAR chrX
                    └─ non-PAR chrY  (M only)
                           │
       ┌───────────────────┼────────────────────┐
       ▼                   ▼                    ▼
  autosomal            sex_M                sex_F
  (all CTRL)           (M CTRL)             (F CTRL)
  BUILD_ATLAS          BUILD_ATLAS          BUILD_ATLAS
  SMOOTH_ATLAS         SMOOTH_ATLAS         SMOOTH_ATLAS
  DETECT_CPG_ISLANDS   DETECT_CPG_ISLANDS   DETECT_CPG_ISLANDS_SEX_F
  CALL_DMR (auto)      CALL_DMR (M cases)   CALL_DMR (F cases, ref∈[0.35,0.65])
       │                   │                    │
       └──────────── ANNOTATE_DMRS ─────────────┘
                     dmr_{ID}_{type}.tsv
                     (PANEL_SCOPE column)
```

## Changes at a glance

### New processes

- **`INFER_SEX_FROM_VCF`** — per-sample sex call from non-PAR chrX HET
  SNP density; writes `pipeline_info/sex_inference/<id>_sex_inference.tsv`.
- **`SPLIT_PAR_AUTOSOMAL`** — bedtools-backed split of each methylation
  BED into autosomal + PAR, non-PAR chrX, non-PAR chrY.
- **`DETECT_CPG_ISLANDS_SEX_F`** — XCI-aware CpG-island detection for
  female non-PAR chrX.

### Refactored processes

- **`BUILD_ATLAS`**, **`SMOOTH_ATLAS`**, **`DETECT_CPG_ISLANDS`** now
  accept a scope-tagged tuple (`autosomal` / `sex_M` / `sex_F`). One
  process definition serves both the autosomes-only and sex-stratified
  paths; the channel layer picks which scopes to instantiate.
- **`CALL_DMR`** now takes a scope, a REF\_METH gate, and per-scope
  threshold arguments. It runs up to three times per case sample.
- **`ANNOTATE_DMRS`** now groups raw per-scope DMR files by sample
  (`groupTuple()`) and merges them into one annotated TSV. A
  `PANEL_SCOPE` column is added only when sex branches are active, so
  autosomes-only output stays byte-identical with the pre-feature release.

### New parameters

All live under the **Chromosome Filter** schema section.

| Parameter | Default | Purpose |
| --- | --- | --- |
| `include_sex_chromosomes` | **`true`** | Master toggle (flipped from `false` in this change) |
| `genome_build` | `GRCh38` | Selects built-in PAR coordinates (`GRCh38` or `GRCh37`) |
| `par_bed` | `null` | Override built-in PAR BED (4-col: CHROM START END NAME) |
| `sex_atlas_min_controls` | `3` | Skip a sex branch if fewer M or F controls are available |
| `sex_het_threshold_male` | `100` | Non-PAR chrX HET count below this → `M` |
| `sex_het_threshold_female` | `1000` | Non-PAR chrX HET count above this → `F` (between → `U`) |
| `sex_override_tsv` | `null` | Optional 2-col TSV (`sample_id`, `sex`) to override inference |
| `sex_f_ref_meth_min` | `0.35` | Lower bound of the XCI-bimodal band for female chrX |
| `sex_f_ref_meth_max` | `0.65` | Upper bound of the XCI-bimodal band |
| `sex_f_diff_threshold` | `0.25` | \|DIFF\| threshold for female non-PAR chrX DMRs (looser than autosomal 0.40) |
| `sex_f_min_sites` | `5` | Minimum CpG sites per female non-PAR chrX DMR (smaller than autosomal 30) |

### New assets

- `assets/par_GRCh38.bed` — PAR1/PAR2 coordinates (both `X/Y` and
  `chrX/chrY` contig naming).
- `assets/par_GRCh37.bed` — same for GRCh37/hg19 cohorts.

### New outputs (with the default-on toggle)

```
results/
├── atlas/
│   ├── ref_atlas_{150,500}.tsv.gz          # autosomal + PAR (same filenames as before)
│   ├── ref_atlas_{150,500}_sex_M.tsv.gz    # male non-PAR chrX + chrY
│   └── ref_atlas_{150,500}_sex_F.tsv.gz    # female non-PAR chrX
├── cpg_islands/
│   ├── cpg_signal.bed                      # autosomal + PAR CpG panel
│   ├── cpg_signal_sex_M.bed                # male non-PAR CpG panel
│   └── cpg_signal_sex_F.bed                # female non-PAR chrX (imprinting-style)
├── dmr/
│   └── dmr_{sample}_case.tsv               # extra PANEL_SCOPE column
└── pipeline_info/
    ├── sex_inference/                      # per-sample VCF-based sex reports
    └── sex_inference_summary.tsv           # combined summary (sample_id, HET count, source, call)
```

### DMR-table additions

- **`PANEL_SCOPE`** column — one of `autosomal`, `sex_M`, `sex_F`. Only
  present when `--include_sex_chromosomes true` (the default). Omitted
  when the user opts out with `--include_sex_chromosomes false`, which
  keeps the schema byte-identical to the pre-feature release.

## Default-flipped

The pipeline previously defaulted to autosomes-only
(`include_sex_chromosomes = false`). It now defaults to the
sex-stratified workflow (`true`). To stay on the old path, pass
`--include_sex_chromosomes false` explicitly.

Cohorts that should opt out:

- Samples without VCFs (sex inference requires the VCF).
- Cohorts with fewer than `sex_atlas_min_controls` male or female
  controls — sex branches would be skipped anyway with a warning, but
  opting out avoids the extra inference step.
- Downstream consumers that depend on the legacy filename set or on DMR
  TSVs without a `PANEL_SCOPE` column.

## Overriding the inferred sex

For aneuploid karyotypes (Turner 45,X, Klinefelter 47,XXY, triple-X,
etc.), or when clinical records should supersede the VCF, pass a 2-column
TSV:

```tsv
sample_id	sex
ONT042	F
ONT098	M
```

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --sex_override_tsv sex_calls.tsv \
    --outdir results
```

Overrides are applied after inference; the audit summary records the
source (`inferred` vs. `override`).

## Bug fixed in the same change

`CALL_DMR` used to emit a phantom NA row when a scope produced zero DMRs.
`rbindlist(list(empty_dt, empty_dt), fill = TRUE)` returns a 0×0
data.table; the subsequent `DT[, col := NA]` loop then created one row
(because `NA` is length 1 and there were no existing rows to recycle
against). The next stage, `ANNOTATE_DMRS`, called `foverlaps()` on that
row and died with "NA values in data.table x 'start' column". Empty-result
scopes now write a header-only file with the canonical typed schema, so
downstream `fread` returns a 0-row table and `foverlaps` is short-circuited.

## Validation

The `test/smoke/run_smoke.sh` end-to-end test now covers both paths:

| Scenario | Exercises | Assertions |
| --- | --- | --- |
| `default` (runs with `--include_sex_chromosomes false`) | Autosomes-only opt-out path; byte-identical output shape | 3 |
| `withsex` (runs under the new default + `sex_override.tsv`) | All three atlases, PAR / non-PAR separation, three CpG panels, PANEL\_SCOPE column, sex-inference summary | 7 |
| `bad_sheet` | Samplesheet validator rejects a missing required column | 1 |

All eleven assertions pass on the GS4 test host with synthetic chr22 +
chrX + chrY data (3 CTRL + 1 CASE) and `sex_atlas_min_controls=1` /
`atlas_min_samples_per_site=1` to accommodate the tiny fixture.

## Files touched

- `main.nf` — new processes, refactored processes, rewritten workflow
  block, channel-level sex routing.
- `nextflow.config` — new parameter block, default flip.
- `nextflow_schema.json` — `chromosome_filter` section expanded.
- `conf/modules.config` — publish rules for the new processes.
- `assets/par_GRCh38.bed`, `assets/par_GRCh37.bed` — new.
- `test/smoke/run_smoke.sh`, `test/smoke/sex_override.tsv` — extended /
  added.
- `README.md`, `CHANGELOG.md` — user-facing documentation updates.

## Scope left for later

- **Phased methylation.** When `modkit --phased` produces
  `hp1.bedmethyl` / `hp2.bedmethyl`, female chrX can be called
  allele-specifically instead of combined-strand. Not wired in yet.
- **Aneuploid karyotypes.** The current classifier assumes 46,XX or
  46,XY. Turner, Klinefelter, and mosaic karyotypes are flagged as `U`
  and routed through autosomes only; a manual override via
  `--sex_override_tsv` is the current workaround.
- **XCI escape-gene masking.** A curated list of genes that escape
  X-inactivation could be subtracted from the sex\_F panel to improve
  specificity. Not implemented.
- **Non-human genomes.** PAR coordinates ship for GRCh38 and GRCh37
  only; other builds require `--par_bed`.
