# Algorithm

This document describes what DMRdetector computes, in the order of the pipeline, with the parameter names of the settings file (`settings.toml`) and the default values. It follows the methods of the manuscript (Ou et al., in preparation) and the code of v2.2 (the R/Nextflow pipeline). Where the code does something the manuscript does not state, the behavior is listed under "Behaviors kept from v2.2".

## Overview

1. Preprocess every sample: remove CpGs on heterozygous SNPs, apply a coverage filter, split the genome into four scopes.
2. Build a reference atlas from control samples: per-site statistics, dispersion, shrinkage, and Gaussian smoothing.
3. Derive two region panels from the atlas: CpG islands and imprinted control regions (ICRs).
4. Test each sample against the atlas inside the two panels and call DMRs.
5. Annotate each DMR with a gene.

Steps 1 to 3 are `dmrdetector build`. Steps 1, 4 and 5 are `dmrdetector run`. The genome build is hg38 (GRCh38).

## 1. Preprocessing

**Input.** A methylation BED with coverage and modified-read counts (plain 5-column, modkit bedMethyl with columns 10 and 12, or pb-CpG-tools with columns 6 and 7) and a VCF.

**Heterozygous SNPs.** The VCF need not be phased: a record counts when some sample has two or more different called alleles in the genotype (`0/1` and `0|1` alike). A CpG is removed when its position falls on a heterozygous single-nucleotide variant of the sample. Variants are typed as `bcftools view -v snps -g het` types them. A record counts as a SNP when REF and ALT, after trimming the common prefix and suffix, differ at one base. This includes equal-length alleles such as `CGTGTGT` and `TGTGTGT`. Only the position of the CpG C is tested; a SNP at the G of the CpG is not removed.

**Input checks.** Preprocessing stops in two cases. The first is a modkit file with 5hmC or 6mA rows next to the 5mC rows, or with one row per strand; each CpG would count twice. The second is a contig with a record that ends beyond its hg38 length. Neither check changes the computation for correct input. See "Preparing input" in the README.

**Coverage.** A CpG is kept when its read coverage is above `atlas_min_read_coverage` (5, controls in `build`) or `dmr_min_coverage` (10, tested samples in `run`).

**Scopes.** The genome is split into four scopes by the sex of the sample:

| Scope | Content | Samples |
|---|---|---|
| autosomal | chr1 to 22 and the pseudoautosomal regions (chrY PAR2 is shifted by +98,813,480 onto chrX PAR2) | all |
| male_X | non-PAR chrX | males |
| female_X | non-PAR chrX | females |
| Y | non-PAR chrY | males |

A sex-specific scope is built only when at least `sex_atlas_min_controls` (3) controls of that sex are available.

## 2. Reference atlas

**Control samples.** `build` reads a samplesheet with the columns `sample_id`, `sex`, `meth_bed` and `vcf`. A sheet without a `type` column lists controls only. In a sheet with a `type` column (the v2.2 samplesheet), only the rows of type `CTRL` or `control` (any case) are controls, as in v2.2. Rows of type `CASE` (any case) are dropped without reading their files, and a NOTICE gives the number of dropped rows. Any other value, an empty one included, stops `build` with an error, as in v2.2. It also stops when no row is a control.

**Per-site aggregation.** For each CpG the pipeline records, over the control samples, the number of samples N, the pooled coverage COV, the pooled methylated reads, and the sample variance S² of the per-sample methylation ratios. The mean coverage is N̄ = COV / N. Sites with N below `atlas_min_samples_per_site` (3) are dropped. The reference methylation is μ = methylated reads / COV, so each sample weighs in proportion to its depth.

**Dispersion.** The Beta-Binomial overdispersion φ is estimated by moments:

    φ_raw = (S² · N̄ − μ(1 − μ)) / (μ(1 − μ) · (N̄ − 1))

Negative values are set to 0 and values above 1 are set to 1. The estimate is 0 when μ(1 − μ)(N̄ − 1) is at most 1e-9 (sites at 0% or 100%).

**Shrinkage.** φ = (1 − w) · φ_raw + w · φ_trend, with w = `atlas_shrinkage_weight` (0.1) and φ_trend the median of the positive φ_raw of the scope (0.01 if there is none).

**Smoothing.** The atlas is smoothed per chromosome with a Gaussian kernel exp(−x² / 2σ²), σ = `smooth_sigma` (150 bp), truncated at radius 3σ. The numerator (μ · COV) and the denominator (COV) are placed on dense position-indexed vectors, convolved with the kernel by FFT, and divided:

    REF_METH = smooth(μ · COV) / smooth(COV)
    REF_VAR  = REF_METH · (1 − REF_METH) · φ

The kernel entries exp(−x² / 2σ²), x = −3σ to 3σ, are computed once per σ with 50-digit decimal arithmetic and rounded once to the nearest float64 (`rstats.exp_rounded`), not with `np.exp`. The argument −x² / 2σ² is evaluated in float64 as in v2.2. numpy and the C library return other last bits for about 5% of the entries depending on the CPU (AVX512 loops, FMA), so the library call would make the atlas depend on the host. For σ = 150 and σ = 500, the kernel equals the one that R 4.4 computes in all 901 and 3,001 entries. The spectra are multiplied with separate float64 multiplications and additions (`fftsmooth.multiply_complex`), because numpy fuses them on a CPU with FMA3.

The FFT length is the smallest integer above the span plus the kernel length with prime factors 2, 3 and 5. The convolution output is centered on the input, so no site is removed; sites near the chromosome ends are normalized by the ratio. Smoothing the numerator and the denominator separately (a Nadaraya-Watson estimator) avoids the bias of smoothing ratios at sites with different coverage.

The atlas columns are `CHROM, START, END, N, COV, PHI, RAW_METH, REF_METH, REF_VAR`.

## 3. Region panels

### CpG islands

Atlas sites with N above `cpg_min_samples_per_site` (5) are classified:

- CORE: `REF_METH` below `cpg_core_max_methylation` (0.10) and PHI below `cpg_core_max_dispersion` (0.05).
- SHORE: `REF_METH` below `cpg_shore_max_methylation` (0.30), within 150 bp of the nearest CORE site, and not CORE.

CORE sites are clustered: a gap above `cpg_max_gap_distance` (150 bp) starts a new cluster. A cluster is kept as an island when its length is at least `cpg_min_island_length` (150 bp) and it has at least one CORE site per 100 bp. The island is extended to the outermost SHORE sites within 150 bp of its CORE boundaries. The panel records the core (START, END), the extended boundaries (START_EXTEND, END_EXTEND), the length and the number of CORE sites.

In female chrX, the hemi-methylated band (`sex_f_ref_meth_min` to `sex_f_ref_meth_max`, 0.40 to 0.60) replaces the CORE rule; the model has no dispersion criterion, no density filter and no SHORE extension. Male chrX and chrY use the autosomal rule.

### Imprinted control regions

`assets/icr_gene_list.bed` lists the imprinted gene regions (from geneimprint.com; the 5 kb extension on each side is part of the file). Atlas sites with N at least `cpg_min_samples_per_site` (5) and `REF_METH` between `icr_ref_meth_min` and `icr_ref_meth_max` (0.40 to 0.60) are grouped per gene: a gap above 150 bp starts a new group. A group with at least `icr_min_sites` (5) CpGs is an ICR. ICRs exist in the autosomal scope only.

## 4. DMR calling

The sample is joined to the atlas on position, and for each CpG DIFF = SAMPLE_METH − REF_METH. CpGs are tested inside each panel separately:

| Panel | CpGs kept | Cluster | Region kept when |
|---|---|---|---|
| CpG island (extended boundaries) | \|DIFF\| above `dmr_cpg_diff_threshold` (0.40) | all kept CpGs of one island | P_FDR below `dmr_fdr_threshold` (0.01) and NUM_CPG above `dmr_cpg_min_sites` (30) |
| ICR | \|DIFF\| above `icr_diff_threshold` (0.25) | CpGs of one ICR, gap above `dmr_merge_distance` (150 bp) starts a new region | P_FDR below 0.01 and NUM_CPG above `icr_min_sites` (5) |
| female_X islands | \|DIFF\| above `sex_f_diff_threshold` (0.25) | as CpG island | P_FDR below 0.01 and NUM_CPG above `sex_f_min_sites` (5) |

For each region the pipeline records NUM_CPG, the mean read coverage (NUM_READS), the mean sample and reference methylation, and the mean difference (DIFF). A one-sample, two-sided t-test compares the per-CpG differences with 0. Regions with fewer than two CpGs are skipped, and a near-constant vector gives NA. The p values are adjusted by the Benjamini-Hochberg method within one sample, scope and panel (`P_FDR`). `DIRECTION` is `hyper` when the mean difference is above 0 and `hypo` otherwise.

## 5. Annotation

Each DMR gets the symbol of an Ensembl gene that overlaps it (a DMR without an overlapping gene has an empty `GENE`). When more than one gene overlaps, biotype decides, in this order: protein-coding, immunoglobulin and T-cell receptor genes, lncRNA and antisense, miRNA, snoRNA and snRNA, rRNA, miscellaneous RNA, processed transcripts, pseudogenes. The largest overlap and then the longest gene break ties. ICR DMRs take the symbol of the ICR gene region with the largest overlap (ties: the smallest region in the gene list). The manuscript used Ensembl release 86 (`EnsDb.Hsapiens.v86`); release 86 reproduces its symbols, and other releases can differ for a few genes.

## Behaviors kept from v2.2

The manuscript text does not describe these, but they produced its numbers, and `dmrdetector` keeps them by default. On the 75 samples of the manuscript run the calls are the same; the floating-point cells agree within the tolerances in [validation.md](validation.md).

| Behavior | Effect | Switch |
|---|---|---|
| CpG island DMRs are clustered on the island core start, and the output row is a 1-bp anchor (`END = START + 1`) | One region per island, whatever the CpG spacing; the table gives no DMR span | none |
| A CpG is counted once per query line that contains it and once per panel row that overlaps it (htslib region queries return a record for each query line). A CpG in m query lines and k panel rows enters a region m x m x k times: m copies from the sample, m from the atlas, and one per panel row | NUM_CPG and the t-test sample size can be larger than the number of distinct CpGs. The factor is 1 for most CpGs. In the manuscript panels it is 8 for 608 island CpGs and 1,715 ICR CpGs, and 27 for 226 ICR CpGs (computed from the panel files). Factors up to 8 were found in the DMR tables of the manuscript run. The P value rests on the repeated CpGs, and the rules NUM_CPG above 30 (CpG islands) and above 5 (ICRs) can be met by fewer distinct CpGs than the numbers say | `dmr_dedupe_sites = true` (DMR calling), `icr_v1_duplicates = false` (ICR panel) |
| ICR panel rows with identical coordinates form one group, because DMRs are grouped by region coordinates, not by gene | The group holds each CpG once per row, so NUM_CPG and the t-test sample size double for two identical rows. The DMR table has one row for the group | none |
| ICR panel rows that overlap but differ in their coordinates (a gene and its antisense gene, for example) are separate regions. Each region is clustered on its CpGs, and two regions can end on the same first and last CpG | The same DMR (same CHROM, START and END) is written twice. The 75 DMR tables of the manuscript run hold 57 such extra rows in 52 tables; all 57 come from panel rows with different coordinates | none |
| The CpG at `END_EXTEND` of an extended island is not read. The query is the half-open interval [`START_EXTEND`, `END_EXTEND`) (BED, 0-based), and `END_EXTEND` is the 0-based start of the last CpG | The last CpG is missing from NUM_CPG and the t-test unless the query line of another panel row contains it. In the manuscript panel (25,719 islands) this holds for 25,206 islands; 513 are covered by another query line | none |
| Atlas sites need N above 5 for CpG islands but at least 5 for ICRs | Different cut-offs in the two panels | `cpg_min_samples_per_site` |
| Strict comparisons (`>`) on coverage, \|DIFF\| and NUM_CPG; non-strict on island length, density, the ICR band and the ICR CpG count | Thresholds read as "above" | none |

## Target regions

The target-region report is separate from DMR detection: see `doc/target_region_methods.md`.

### Target-region statistics

The report has one row per sample, region and track (merged, `hap1`, `hap2`). Samples are ordered CASE first, then CTRL (sample ID within type), and regions by chromosome (string order), START and END. Only the merged rows are tested, because the atlas is a merged mixture; `hap1` and `hap2` rows have NA in `DIFF`, `P_FDR` and `DIRECTION`.

- `DIFF` is the mean of the per-CpG differences `SAMPLE_METH - REF_METH`, a fraction, as `MEAN_DIFF` in the DMR table.
- `P_FDR` comes from a two-sided one-sample t-test (mu = 0) of the per-CpG differences, adjusted by Benjamini-Hochberg across the evaluable regions of one sample. Only evaluable regions enter the adjustment; the other rows get NA. A region is not evaluable when the t-test is undefined (fewer than 2 CpGs, or constant differences). There is no other minimum CpG count.
- `DIRECTION` is `hyper` when `P_FDR` is below `dmr_fdr_threshold` and `DIFF` is above the threshold, `hypo` when it is below the negative threshold, and `normal` otherwise. The threshold is `icr_diff_threshold` (0.25) when the region overlaps the ICR panel or the imprinted gene list, else `dmr_cpg_diff_threshold` (0.40). Without an ICR file the bundled `assets/icr_gene_list.bed` is used. Every evaluable region gets a value, whether or not `DISEASE_NAME` names a direction.

CpG filtering follows v2: coverage must be above `target_min_coverage` (merged, 8) or `target_hap_min_coverage` (hap1 and hap2, 4). CpGs on a heterozygous SNP of the sample VCF are removed when a VCF is given. The VCF is required for the haplotype tracks, because hap1 and hap2 BEDs never pass through preprocessing. Contig names are compared as exact strings. Regions on chrX, chrY or any other non-autosome are skipped with a message, because the reference is the autosomal atlas.

## Implementation notes

These notes describe how the fast readers and writers work. They change no result.

- Threads: BGZF blocks are inflated in a thread pool with libdeflate through ctypes (zlib when libdeflate is missing or `DMRDETECTOR_NO_LIBDEFLATE=1`). The numba kernels are `nogil` and use plain float64 and integer arithmetic, with no fastmath and no parallel reductions. Each output cell comes from one input line, so the result does not depend on the number of threads. `DMRDETECTOR_NO_JIT=1` runs the kernels as plain Python (slow; used by the tests to compare compiled and interpreted results).
- Numbers: A decimal of at most 15 significant digits and an exponent within +-22 is one correctly rounded IEEE operation (the mantissa times or divided by an exact power of ten). Up to 19 digits a double-double product is used. When its rounding is within 1e-9 of half a unit in the last place of a tie, the token goes to Python `float()`. No value is then off by one unit in the last place. The pandas fallbacks read decimals with `float_precision="round_trip"`, which is also correctly rounded, so both routes give the same doubles. The default pandas parser, and its `high` option, can be off in the last digits: it reads `0.000535864718615665` as `0.0005358647186156`. A token a kernel cannot read (more than 19 digits, `NA`, `nan`) sets a flag and the caller reads the file with pandas.
- Fallbacks: The byte-level methylation BED reader accepts ordinary files only: tab or space separated, the same number of fields on every data line, plain ASCII, digit-only START, END, coverage and N_mod. A `#` or blank line after the header, quotes, carriage returns, decimals, signs, exponents, a short or long line, or more than 65,535 contigs per chunk sends the whole file to the pandas code. The VCF scan falls back to the pysam loop on a carriage return, a non-ASCII REF or ALT, a short line, a bad POS, a genotype allele that is not a number or points past the ALT list, or an empty allele.
- Number text: Output numbers follow R `fwrite`: 15 significant digits, the shorter of fixed and scientific notation, fixed on ties. The Python formatter rounds exactly; `fwrite` computes the 15th digit with its own arithmetic and can differ in about 3 % of random ratios (relative size 1e-15).
- Output files: Writers deflate in a thread pool and report the virtual offset of every row start. The tabix index is written from those offsets, with no second pass over the file, and equals the index htslib builds from the finished file.
- FFT: The convolution of each chromosome runs in a short-lived child process, because scipy.fft keeps the plan of every transform length (about 8 bytes per padded base) until the process ends. The child runs the same function with the same arguments, so the numbers are identical; the parent computes `REF_VAR` from the returned `REF_METH`. `DMRDETECTOR_FFT_LOG=1` makes the child report its peak memory on stderr.

## Differences between the manuscript text and this description

The manuscript describes the atlas statistics and the two panels as above. Six differences matter to a reader who reimplements the method. The code does not trim the convolution output; it centers it. Annotation uses overlap, not the nearest gene. The island DMR clustering and the duplicate counting are not in the text. The CpG at the end of an extended island is not tested. The two panels use different N cut-offs. The φ_trend median is taken per scope.
