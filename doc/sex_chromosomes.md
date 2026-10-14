# Sex chromosomes in DMRdetector v3

This note describes how `build` and `run` treat chrX and chrY. Every statement below was checked against the v3 code (`src/dmrdetector/build.py`, `run.py`, `runtime.py`, `settings.py`, `steps/cpg_islands.py`, `steps/call_dmr.py`, `steps/preprocess.py`) and, for the numbers, against the settings file. The biological rationale is the one given in the manuscript; it was not re-examined against the literature for this note.

## Why separate scopes

A pooled reference of both sexes has no meaning on chrX and chrY. Males carry one copy of each, females carry two copies of chrX and none of chrY. In females one X is inactivated, so CpG island promoters on chrX read as intermediate methylation (about 50%) on a combined-strand BED, and the autosomal island criterion (`REF_METH` below 0.10) misses them. The pseudoautosomal regions (PAR1 and PAR2) are diploid in both sexes and behave like autosomal sequence.

## Scopes

The genome is split into four scopes during preprocessing (`preprocess split`) and each scope has its own atlas and CpG island panel.

| Scope | Region | Controls used for the atlas | CpG island model | DMR thresholds |
|---|---|---|---|---|
| `autosomal` | chr1 to chr22, PAR1 and PAR2 | all CTRL samples | CORE: `REF_METH` < 0.10 and `PHI` < 0.05; SHORE: `REF_METH` < 0.30 within 150 bp of a CORE site | absolute difference > 0.40, more than 30 CpGs, FDR < 0.01 |
| `male_X` | non-PAR chrX | male CTRL samples | as autosomal | as autosomal |
| `female_X` | non-PAR chrX | female CTRL samples | CORE: 0.40 <= `REF_METH` <= 0.60; no `PHI` filter, no SHORE extension, no density rule (length of at least 150 bp only) | absolute difference > 0.25, more than 5 CpGs, FDR < 0.01 |
| `Y` | non-PAR chrY | male CTRL samples | as autosomal | as autosomal |

The ICR panel and the ICR DMR calls (absolute difference > 0.25, more than 5 CpGs) use the autosomal scope only.

## Pseudoautosomal regions

The bundled `par_hg38.bed` holds PAR1 and PAR2 of GRCh38 on chrX and chrY (with and without the `chr` prefix). Rows of chrX and of chrY inside a PAR go to the autosomal part. Rows of chrY PAR2 are moved by +98,813,480 bp onto the chrX PAR2 coordinates; PAR1 has the same coordinates on both chromosomes. Only non-PAR chrX and non-PAR chrY go to the sex scopes. The genome build is hg38, and the `par_bed` setting replaces the file.

## Sex of a sample

The `sex` column of the samplesheet (`M` or `F`) decides the routing. v3 does not infer sex from the VCF. `build` and `run` stop with an error that names the samples when a sex-specific scope is needed and a sample has no `M` or `F`.

* A male sample is tested in `autosomal`, `male_X` and `Y`.
* A female sample is tested in `autosomal` and `female_X`.
* `build` makes a sex-specific scope only when the samplesheet has at least `sex_atlas_min_controls` (default 3) CTRL samples of that sex; otherwise it prints a warning and writes no reference for that scope. `run` uses the sex scopes only when the reference has an atlas and a CpG panel for all three of `male_X`, `female_X` and `Y`; otherwise it prints a warning and calls DMRs in the autosomal scope only.
* `--no-include-sex-chromosomes` (`include_sex_chromosomes = false`) restricts the analysis to chr1 to chr22 and skips the split.

## Check of the sex

`build` and `run` compare the `sex` label of each sample with two markers from its own files. They only warn; the label still decides the scopes.

* **chrY CpGs.** A sample labeled F that keeps at least 50,000 CpGs on non-PAR chrY after the filters gets a warning (`CHRY_MALE_MIN_CPG` in `steps/preprocess.py`). The cutoff is our rule for calling a male in whole-blood ONT controls: in 143 candidate controls, males had 68,408 to 171,117 chrY CpG records and females 4,928 to 8,570. The coverage filter removes CpGs and never adds them, so a female sample above the cutoff has at least that many records before filtering. This rule is one-sided. A male with few CpGs on chrY is not reported, because the filtered count is lower than the raw count that the cutoff was set on and we have no cutoff for it.
* **chrX heterozygous SNPs.** The lab convention for selecting atlas controls is one cutoff on the number of heterozygous SNVs on non-PAR chrX: at least 40,000 for a female, fewer for a male (`FEMALE_CHRX_HET_MIN`). In the same 143 candidates, females had 60,362 to 83,844 and males 6,532 to 19,199. A sample labeled M with 40,000 or more, or labeled F with fewer, gets a warning. The count uses the heterozygous SNP sites that the preprocess step already reads from the VCF, outside the PARs of `par_bed`, on the contigs `chrX` or `X`; it adds no pass over the file. It differs slightly from the count of the lab script, which reads the first sample only and counts single-base REF and ALT only. Here a site counts when some ALT is a SNP and some sample is heterozygous, so multi-allelic sites count too. The check is skipped when the VCF has no sample column or more than one, and when the run does not split by scope (`include_sex_chromosomes = false`). A VCF with few chrX sites for a technical reason (for example targeted sequencing) makes every female label warn.

The ranges come from one technology and one tissue; check them before you rely on a warning for other data. The warnings are in the `[warnings]` section of the log.

## Output

`build` writes `ref_atlas_150.bed.gz` (plus `.tbi`) and `cpg_signal.bed` for the autosomal scope, and the same files with the suffixes `_male_X`, `_female_X` and `_Y` for the sex scopes, directly into its `--outdir`; the number in the atlas name is `smooth_sigma`. `run` writes the DMRs of all scopes of the sample into one table, `<sample_id>.dmr.tsv`, in the flat output folder; the table has no scope column. The per-scope raw DMR tables are temporary files that `run` removes at the end. The files and columns are listed in [output.md](output.md).

## Not handled

* Karyotypes other than 46,XX and 46,XY (for example 45,X or 47,XXY): the scope follows the `sex` entry of the samplesheet. The only checks against the data are the two warnings above.
* Skewed X inactivation, and genes that escape X inactivation: neither is modeled or masked.
* Assemblies other than hg38.
