# Output files

## `dmrdetector run`

All files go directly into one folder (default `dmrdetector/`). Every file name starts with the sample ID, so samples can share the folder.

```
dmrdetector/
  <sample_id>.dmr.tsv              annotated DMRs
  <sample_id>.target_region.tsv    only with --target-region
  <sample_id>.log                  command, parameters, timings, warnings (file names, not paths)
```

Intermediate files (preprocessed BEDs, raw DMR tables) are written to a temporary directory and removed at the end of the run, also after an error.

The log lists each input file by its name and size, and shows every other path (also in the command line and the settings) by its last component, so that the log can be shared. Set `log_input_paths = "full"` (`--set log_input_paths=full`) to write absolute paths. The progress log in a kept temporary folder (`--keep-tmp`) always has full paths. If the target-region report fails, `<sample_id>.dmr.tsv` is still written, the log says `status: FAILED (the DMR table was published; the target-region report failed)`, and the exit status is 1.

### DMR table

| Column | Meaning |
|---|---|
| `CHROM, START, END` | Region (hg38). CpG island DMRs are 1-bp anchors (`END = START + 1`), as in v2.2. |
| `GENE` | Gene symbol from the Ensembl annotation (overlap, with biotype priority). |
| `PANEL` | `CpG` or `ICR`. |
| `NUM_CPG`, `NUM_READS` | CpGs in the region and their mean read coverage. |
| `SAMPLE_METH (%)`, `REF_METH (%)`, `DIFF (%)` | Mean methylation of the sample, of the atlas, and their difference. |
| `P_VALUE`, `P_FDR` | One-sample t-test on the per-CpG differences; Benjamini-Hochberg FDR within each panel. |
| `DIRECTION` | `hyper` if the mean difference is above 0, `hypo` if below. |

Default thresholds (all adjustable in the settings file):

| Panel | Absolute difference | CpGs | FDR |
|---|---|---|---|
| CpG island | above 0.40 | above 30 | below 0.01 |
| ICR | above 0.25 | above 5 | below 0.01 |
| female chrX (hemi-methylated model) | above 0.25 | above 5 | below 0.01 |

### Target-region table

One row per region and track: merged reads, plus `hap1` and `hap2` when `--meth_bed_hap1` and `--meth_bed_hap2` are given.

| Column | Meaning |
|---|---|
| `SAMPLE_ID`, `SAMPLE_TYPE` | The sample ID and its type (`CASE` for a `run` sample). |
| `CHROM, START, END, REGION_NAME, GENE_SYMBOL, DISEASE_NAME` | The BED row. |
| `HAPLOTYPE` | `merged`, `hap1` or `hap2` (phasing labels, not parent of origin). |
| `NUM_CPG`, `COVERAGE` | CpGs that pass the filters and their mean coverage. |
| `SAMPLE_METH (%)`, `REF_METH (%)`, `REF_SD (%)` | Sample, atlas and between-individual SD. |
| `DIFF`, `P_FDR`, `DIRECTION` | Mean per-CpG difference as a fraction (sample minus atlas, from -1 to 1; not a percentage, unlike `DIFF (%)` of the DMR table), FDR across the regions of the sample, and `hyper`, `hypo` or `normal`. Empty for `hap1` and `hap2` rows. |

A region without a usable CpG (`NUM_CPG` 0, for example when the BED names contigs that the sample lacks) has `NA` in the methylation columns, `P_FDR` and `DIRECTION`; `run` prints a WARNING when this holds for every region of the sample. `DIRECTION` is the state of the sample against the atlas. A region is called when `P_FDR` is below `dmr_fdr_threshold` and the absolute difference exceeds the ICR threshold (regions overlapping the ICR panel or the imprinted gene list) or the CpG island threshold (all other regions). There is no minimum CpG count; a region with fewer than two usable CpGs gets `NA`. Read `DIRECTION` together with the direction in `DISEASE_NAME`: for the Prader-Willi/Angelman center, `hyper` fits Prader-Willi and `hypo` fits Angelman. The atlas is built from merged reads, so it is the wrong reference for a single haplotype; the `hap1` and `hap2` rows show the raw haplotype methylation only. Methods and references: [target_region_methods.md](target_region_methods.md).

## `dmrdetector build`

```
reference/
  ref_atlas_150.bed.gz[.tbi]             autosomal atlas (sigma 150 bp)
  ref_atlas_150_{male_X,female_X,Y}.bed.gz[.tbi]
  cpg_signal.bed, cpg_signal_{male_X,female_X,Y}.bed    CpG island panels
  icr_region.bed                         ICR panel (autosomal)
  reference.json                         version, parameters, sample IDs, counts
  build.log                              like the run log: file names, not paths
```

`reference.json` holds the controls (sample ID and sex, which document how the reference was made), the name and md5 of the samplesheet but not its path (the md5 is that of the whole file; a `type` column keeps only the control rows in the list of controls), the md5 of the imprinting gene list that `build` used (`icr_gene_list_md5`, checked by `run`), the name of the settings file, and the parameters (a path option such as `par_bed` or `icr_gene_list` is recorded by its file name, unless `log_input_paths = "full"`). Remove or change the sample IDs before you share it when they identify people.

Atlas columns (headerless, bgzip-compressed, tabix-indexed BED): `CHROM, START, END, N, COV, PHI, RAW_METH, REF_METH, REF_VAR`.

| Column | Meaning |
|---|---|
| `N`, `COV` | Control samples and pooled read coverage at the site. |
| `PHI` | Shrunken inter-individual dispersion. |
| `RAW_METH` | Coverage-weighted mean methylation before smoothing (0 to 1). |
| `REF_METH` | Smoothed reference methylation (0 to 1). |
| `REF_VAR` | `REF_METH × (1 − REF_METH) × PHI`. |
