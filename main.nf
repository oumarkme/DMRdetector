#!/usr/bin/env nextflow

/*
 * DMRdetector: ONT Methylation Analysis Pipeline
 * ==============================================
 * Detects differentially methylated regions (DMRs) from Oxford Nanopore
 * Technologies (or PacBio) whole-genome sequencing data.
 *
 * Three stages:
 *   1. Build reference atlas from control samples
 *   2. Detect hypomethylated CpG islands
 *   3. Detect DMRs per sample vs reference (CpG island + ICR panels)
 *
 * Three usage scenarios (auto-detected):
 *   - Build atlas only:  --ref_atlas NOT set, samplesheet has only CTRL
 *   - DMR detection only: --ref_atlas provided → skip Stages 1-2
 *   - Full pipeline:      --ref_atlas NOT set, CASE samples exist
 *
 * Supported input BED formats (auto-detected):
 *   - 5-col: CHROM, START, END, COVERAGE, N_MOD
 *   - modkit pileup (bedMethyl): 18 cols → uses cols 10 (Nvalid_cov), 12 (Nmod)
 *   - pb-CpG-tools: 9 cols (with header) → uses cols 6 (cov), 7 (est_mod_count)
 *
 * Runs natively inside the bundled `dmrdetector` conda environment
 * (see environment.yml).
 */

// ============================================================================
// Helper functions
// ============================================================================

def show_help() {
    log.info """
    =================================================================
    DMRdetector v${workflow.manifest.version}
    =================================================================

    Usage:
      nextflow run main.nf --input samplesheet.csv [options]

    Required:
      --input           Samplesheet CSV (sample_id, type, sex, meth_bed, vcf)
                        type = CASE|CTRL; sex = M|F

    Optional:

    Panels:
      --icr_panel       ICR BED override (default: build from reference atlas)
      --icr_gene_list   Imprinting gene candidate BED (default: assets/icr_gene_list.bed)

    Optional (skip or override panel/reference generation):
      --ref_atlas       Pre-built reference atlas (tsv.gz)
      --cpg_panel       Pre-built CpG island panel (bed)

    Chromosome coverage:
      --include_sex_chromosomes
                        Include sex-specific chrX/chrY analysis [default: true].
                        Builds component atlases for autosomes, non-PAR
                        male_X, non-PAR female_X, and Y. Use false for
                        autosomes only.

    Output:
      --outdir          Output directory [default: ./results]

    Performance:
      --tmp_dir         Local temp directory for work/tmp files (recommended
                        when inputs are on network storage)

    Input BED formats (auto-detected):
      - 5-column:   CHROM, START, END, COVERAGE, N_MOD
      - modkit:     bedMethyl 18-col (uses cols 10, 12)
      - pb-CpG-tools: 9-col with header (uses cols 6, 7)

    Profiles:
      -profile test         Reduced resources for testing

    Full parameter list: nextflow run main.nf --help
    =================================================================
    """.stripIndent()
}

def validate_params() {
    if (params.help) {
        show_help()
        exit 0
    }
    if (!params.input) {
        show_help()
        error "ERROR: --input samplesheet is required"
    }
}

// ============================================================================
// Process: PREPROCESS_SAMPLE
// ============================================================================

process PREPROCESS_SAMPLE {
    label 'standard'
    tag "${meta.id}"

    input:
        tuple val(meta), path(meth_bed), path(vcf)
        val min_coverage

    output:
        tuple val(meta), path("${meta.id}_filtered.bed.gz"), optional: true

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Step 1: Compress and index VCF if needed
    if [[ "${vcf}" != *.vcf.gz ]]; then
        bgzip -c -@ ${task.cpus} "${vcf}" > temp.vcf.gz
        tabix -p vcf temp.vcf.gz
        VCF_FILE="temp.vcf.gz"
    else
        cp "${vcf}" temp.vcf.gz
        tabix -p vcf temp.vcf.gz
        VCF_FILE="temp.vcf.gz"
    fi

    # Step 2: Extract heterozygous SNP positions
    bcftools view --threads ${task.cpus} -v snps -g het "\${VCF_FILE}" \\
        | bcftools query -f '%CHROM\\t%POS0\\t%POS\\n' \\
        | sort -k1,1 -k2,2n \\
        > het_sites.bed

    # Step 3: Auto-detect BED format and extract CHROM, START, END, COV, N_MOD
    # Detect number of columns from first non-header data line
    NCOLS=\$(set +o pipefail; (zcat "${meth_bed}" 2>/dev/null || cat "${meth_bed}") \\
        | grep -v '^#' | grep -v '^chrom' | head -1 | awk '{print NF}'; set -o pipefail)

    if [ "\$NCOLS" -ge 18 ]; then
        # modkit bedMethyl format (18 cols): cols 10=Nvalid_cov, 12=Nmod
        echo "Detected modkit bedMethyl format (\${NCOLS} columns)" >&2
        bedtools intersect -v -a "${meth_bed}" -b het_sites.bed \\
            | awk -v min_cov=${min_coverage} '
                BEGIN{OFS="\\t"}
                /^#/{next} /^chrom/{next}
                { cov=\$10+0; nmod=\$12+0;
                  if (cov > min_cov) print \$1, \$2, \$3, cov, nmod }
              ' \\
            | bgzip -@ ${task.cpus} \\
            > ${meta.id}_filtered.bed.gz
    elif [ "\$NCOLS" -ge 7 ] && [ "\$NCOLS" -le 9 ]; then
        # pb-CpG-tools format (9 cols with header): cols 6=cov, 7=est_mod_count
        echo "Detected pb-CpG-tools format (\${NCOLS} columns)" >&2
        bedtools intersect -v -a "${meth_bed}" -b het_sites.bed \\
            | awk -v min_cov=${min_coverage} '
                BEGIN{OFS="\\t"}
                /^#/{next} /^chrom/{next}
                { cov=\$6+0; nmod=\$7+0;
                  if (cov > min_cov) print \$1, \$2, \$3, cov, nmod }
              ' \\
            | bgzip -@ ${task.cpus} \\
            > ${meta.id}_filtered.bed.gz
    else
        # Standard 5-column format: CHROM, START, END, COVERAGE, N_MOD
        echo "Detected 5-column format (\${NCOLS} columns)" >&2
        bedtools intersect -v -a "${meth_bed}" -b het_sites.bed \\
            | awk -v min_cov=${min_coverage} '
                BEGIN{OFS="\\t"}
                /^#/{next}
                \$4 > min_cov { print \$1, \$2, \$3, \$4, \$5 }
              ' \\
            | bgzip -@ ${task.cpus} \\
            > ${meta.id}_filtered.bed.gz
    fi

    rm -f temp.vcf.gz temp.vcf.gz.tbi het_sites.bed
    """
}

// ============================================================================
// Process: SPLIT_SEX_COMPONENTS
// ----------------------------------------------------------------------------
// Splits a preprocessed 5-col methylation BED into component files:
//   - autosomal: chr1-22 only
//   - chrX: non-PAR chrX only
//   - chrY: full chrY
// Only runs when params.include_sex_chromosomes is true.
// ============================================================================

process SPLIT_SEX_COMPONENTS {
    label 'light'
    tag "${meta.id}"

    input:
        tuple val(meta), path(bed)
        path par_bed

    output:
        tuple val(meta),
              path("${meta.id}_autosomal.bed.gz"),
              path("${meta.id}_chrX.bed.gz"),
              path("${meta.id}_chrY.bed.gz")

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Autosomal rows (chr1-22 with optional 'chr' prefix).
    zcat ${bed} \\
        | awk 'BEGIN{OFS="\\t"} \$1 ~ /^(chr)?([0-9]+)\$/' \\
        | sort -k1,1 -k2,2n \\
        | bgzip -@ ${task.cpus} \\
        > ${meta.id}_autosomal.bed.gz

    # Non-PAR chrX. The bundled PAR BED contains both chr-prefixed and
    # unprefixed contig names, so this works with either BED convention.
    zcat ${bed} \\
        | awk 'BEGIN{OFS="\\t"} \$1 ~ /^(chr)?X\$/' \\
        | bedtools intersect -v -a - -b ${par_bed} \\
        | sort -k1,1 -k2,2n \\
        | bgzip -@ ${task.cpus} \\
        > ${meta.id}_chrX.bed.gz

    # Full chrY. Y is handled as a male-specific exploratory component.
    zcat ${bed} \\
        | awk 'BEGIN{OFS="\\t"} \$1 ~ /^(chr)?Y\$/' \\
        | sort -k1,1 -k2,2n \\
        | bgzip -@ ${task.cpus} \\
        > ${meta.id}_chrY.bed.gz
    """
}

// ============================================================================
// Process: BUILD_ATLAS
// ============================================================================

process BUILD_ATLAS {
    label 'heavy'
    tag "${scope}"

    input:
        tuple val(scope), path(bed_files)

    output:
        tuple val(scope), path("raw_atlas.rds"), emit: atlas_rds

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(data.table))
    setDTthreads(${task.cpus})

    # 1. Discover input files
    files <- list.files(".", pattern = "\\\\.bed\\\\.gz\$", full.names = TRUE)
    files <- files[file.info(files)\$size > 50]
    if (length(files) == 0) stop("No valid input files found")
    message(sprintf("Found %d valid input files.", length(files)))

    # 2. Read all BED files
    dt_list <- lapply(files, function(f) {
        tryCatch({
            fread(f, select = 1:5,
                  col.names = c("CHROM", "START", "END", "COVERAGE", "N_METH"),
                  header = FALSE, nThread = 1)
        }, error = function(e) { message("Warning: Failed to read ", f); NULL })
    })
    dt_list <- Filter(Negate(is.null), dt_list)
    if (length(dt_list) == 0) stop("All file reads failed")
    dt <- rbindlist(dt_list); rm(dt_list); gc()
    message(sprintf("Loaded %d total CpG observations.", nrow(dt)))

    # 3. Compute per-sample methylation ratio
    dt[, ratio := N_METH / COVERAGE]

    # 4. Aggregate per CpG site
    message("Aggregating data per CpG site...")
    atlas <- dt[, .(
        N = .N,
        COV = sum(COVERAGE),
        sum_meth = sum(N_METH),
        sum_sq_r = sum((ratio - mean(ratio))^2),
        mean_cov = mean(COVERAGE)
    ), by = .(CHROM, START, END)]
    rm(dt); gc()

    # 5. Filter by min samples
    atlas <- atlas[N >= ${params.atlas_min_samples_per_site}]
    if (nrow(atlas) == 0) stop("No sites passed min_samples filter")
    message(sprintf("Retained %d CpG sites after filtering.", nrow(atlas)))

    # 6. Compute raw methylation and dispersion
    atlas[, RAW_METH := sum_meth / COV]
    atlas[, `:=`(S2 = sum_sq_r / (N - 1), N_bar = COV / N, mu_var = RAW_METH * (1 - RAW_METH))]
    atlas[, phi_num := S2 * N_bar - mu_var]
    atlas[, phi_denom := mu_var * (N_bar - 1)]
    atlas[, phi_raw := 0.0]
    atlas[phi_denom > 1e-9 & phi_num > 0, phi_raw := phi_num / phi_denom]
    atlas[phi_raw > 1, phi_raw := 1]

    # 7. Empirical Bayes shrinkage
    valid_phis <- atlas[phi_raw > 0, phi_raw]
    phi_trend <- if (length(valid_phis) == 0) 0.01 else median(valid_phis)
    message(sprintf("phi_trend (median): %.6f", phi_trend))
    w <- ${params.atlas_shrinkage_weight}
    atlas[, PHI := (1 - w) * phi_raw + w * phi_trend]

    # 8. Output
    final <- atlas[, .(CHROM, START, END, N, COV, PHI, RAW_METH)]
    setkey(final, CHROM, START)
    saveRDS(final, "raw_atlas.rds")
    message(sprintf("Saved %d sites to raw_atlas.rds", nrow(final)))
    """
}

// ============================================================================
// Process: SMOOTH_ATLAS
// ============================================================================

process SMOOTH_ATLAS {
    label 'heavy'
    tag "sigma_${sigma}_${scope}"

    input:
        tuple val(scope), path(atlas_rds), val(sigma)

    output:
        tuple val(scope), val(sigma), path("${out_file}"), emit: smoothed

    script:
    out_file = scope == 'autosomal' ? "ref_atlas_${sigma}.tsv.gz" : "ref_atlas_${sigma}_${scope}.tsv.gz"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(data.table)
        library(parallel)
        library(stats)
    })
    setDTthreads(1)  # single-threaded within parallel workers

    sigma <- ${sigma}
    threads <- ${task.cpus}
    message(sprintf("=== FFT Gaussian Smoothing (sigma = %d bp) ===", sigma))

    # Gaussian kernel
    get_gaussian_kernel <- function(sigma) {
        radius <- ceiling(3 * sigma)
        x <- seq(-radius, radius)
        exp(-(x^2) / (2 * sigma^2))
    }

    # FFT convolution
    fft_convolve <- function(data_vec, kernel) {
        n_data <- length(data_vec)
        n_kernel <- length(kernel)
        full_len <- n_data + n_kernel - 1
        fft_len <- nextn(full_len)
        data_padded <- c(data_vec, numeric(fft_len - n_data))
        kernel_padded <- c(kernel, numeric(fft_len - n_kernel))
        result_fft <- fft(data_padded) * fft(kernel_padded)
        result_conv <- Re(fft(result_fft, inverse = TRUE)) / fft_len
        result_conv[1:full_len]
    }

    # Process one chromosome
    process_chromosome <- function(chr, dt_all, kernel) {
        sub_dt <- dt_all[CHROM == chr]
        if (nrow(sub_dt) == 0) return(NULL)
        min_pos <- min(sub_dt\$START); max_pos <- max(sub_dt\$END)
        span <- max_pos - min_pos + 1
        vec_meth <- numeric(span); vec_cov <- numeric(span)
        indices <- sub_dt\$START - min_pos + 1
        vec_meth[indices] <- sub_dt\$RAW_METH * sub_dt\$COV
        vec_cov[indices] <- sub_dt\$COV
        smooth_meth <- fft_convolve(vec_meth, kernel)
        smooth_cov <- fft_convolve(vec_cov, kernel)
        k_len <- length(kernel)
        start_trim <- floor(k_len / 2)
        valid_range <- (start_trim + 1):(start_trim + span)
        smooth_meth <- smooth_meth[valid_range]
        smooth_cov <- smooth_cov[valid_range]
        smooth_vals <- numeric(span)
        valid_mask <- smooth_cov > 1e-9
        smooth_vals[valid_mask] <- smooth_meth[valid_mask] / smooth_cov[valid_mask]
        sub_dt[, REF_METH := smooth_vals[indices]]
        sub_dt[, REF_VAR := REF_METH * (1 - REF_METH) * PHI]
        sub_dt
    }

    # Load data
    dt <- readRDS("${atlas_rds}")
    include_sex <- ${params.include_sex_chromosomes ? 'TRUE' : 'FALSE'}
    normalize_primary_chrom <- function(x) {
        chrom <- as.character(x)
        base <- toupper(sub("^chr", "", chrom, ignore.case = TRUE))
        primary <- base %in% c(as.character(1:22), "X", "Y")
        chrom[primary] <- paste0("chr", base[primary])
        chrom
    }

    allowed_chroms <- paste0("chr", 1:22)
    if (include_sex) allowed_chroms <- c(allowed_chroms, "chrX", "chrY")
    dt[, CHROM := normalize_primary_chrom(CHROM)]
    dt <- dt[CHROM %in% allowed_chroms]
    if (nrow(dt) == 0) stop("No allowed chromosomes found (include_sex_chromosomes=", include_sex, ")")
    message(sprintf("Chromosome filter: include_sex_chromosomes=%s (%d chroms retained)",
                    include_sex, length(unique(dt\$CHROM))))
    message(sprintf("Loaded %d CpG sites.", nrow(dt)))
    dt[, `:=`(REF_METH = as.numeric(NA), REF_VAR = as.numeric(NA))]

    kernel <- get_gaussian_kernel(sigma)
    chroms <- unique(dt\$CHROM)
    use_threads <- min(threads, length(chroms))
    message(sprintf("Processing %d chromosomes using %d cores...", length(chroms), use_threads))

    result_list <- mclapply(chroms, function(chr) {
        tryCatch(process_chromosome(chr, dt, kernel), error = function(e) {
            message(sprintf("Error processing %s: %s", chr, e\$message)); NULL
        })
    }, mc.cores = use_threads, mc.preschedule = FALSE)

    result_list <- Filter(Negate(is.null), result_list)
    if (length(result_list) == 0) stop("All chromosomes failed")
    final_dt <- rbindlist(result_list)
    setkey(final_dt, CHROM, START)
    final_output <- final_dt[, .(CHROM, START, END, N, COV, PHI, RAW_METH, REF_METH, REF_VAR)]

    out_file <- "${out_file}"
    fwrite(final_output, out_file, sep = "\\t", nThread = threads)
    message(sprintf("Written %d sites to %s", nrow(final_output), out_file))
    """
}

// ============================================================================
// Process: DETECT_CPG_ISLANDS
// ============================================================================

process DETECT_CPG_ISLANDS {
    label 'medium'
    tag "${scope}"

    input:
        tuple val(scope), path(ref_atlas)

    output:
        tuple val(scope), path("${out_file}"), emit: cpg_bed

    script:
    out_file = scope == 'autosomal' ? "cpg_signal.bed" : "cpg_signal_${scope}.bed"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(data.table)
    })
    setDTthreads(${task.cpus})

    cat("=== Detect CpG Island Panel ===\\n")

    # Parameters
    min_samples  <- ${params.cpg_min_samples_per_site}
    core_max_meth <- ${params.cpg_core_max_methylation}
    core_max_phi  <- ${params.cpg_core_max_dispersion}
    shore_max_meth <- ${params.cpg_shore_max_methylation}
    female_min_meth <- ${params.sex_f_ref_meth_min}
    female_max_meth <- ${params.sex_f_ref_meth_max}
    max_gap       <- ${params.cpg_max_gap_distance}
    min_length    <- ${params.cpg_min_island_length}
    scope         <- "${scope}"

    # Empty output helper
    write_empty <- function() {
        fwrite(data.table(CHROM=character(), START=integer(), END=integer(),
                          START_EXTEND=integer(), END_EXTEND=integer(),
                          LENGTH=integer(), NUM_CPG=integer()),
               "${out_file}", sep="\\t", col.names=TRUE)
    }

    # Step 1: Load and filter atlas
    ref_atlas <- fread("${ref_atlas}")
    ref_atlas[, CHROM := as.character(CHROM)]
    cat("  Loaded", nrow(ref_atlas), "sites\\n")

    # Apply chromosome filter (handles the DMR-only case where a user-supplied
    # atlas may contain sex chromosomes that the current run wants to exclude)
    include_sex <- ${params.include_sex_chromosomes ? 'TRUE' : 'FALSE'}
    normalize_primary_chrom <- function(x) {
        chrom <- as.character(x)
        base <- toupper(sub("^chr", "", chrom, ignore.case = TRUE))
        primary <- base %in% c(as.character(1:22), "X", "Y")
        chrom[primary] <- paste0("chr", base[primary])
        chrom
    }

    allowed_chroms <- paste0("chr", 1:22)
    if (include_sex) allowed_chroms <- c(allowed_chroms, "chrX", "chrY")
    ref_atlas <- ref_atlas[CHROM %in% allowed_chroms]
    cat("  After chromosome filter (include_sex=", include_sex, "):", nrow(ref_atlas), "sites\\n")

    ref_atlas <- ref_atlas[N > min_samples]
    cat("  After N >", min_samples, "filter:", nrow(ref_atlas), "sites\\n")

    if (scope == "female_X") {
        cat("  Female_X hemi-methylated model: REF_METH between",
            female_min_meth, "and", female_max_meth, "\\n")
        ref_atlas[, CORE := REF_METH >= female_min_meth & REF_METH <= female_max_meth]
        n_core <- sum(ref_atlas\$CORE)
        cat("  Hemi-methylated sites:", n_core, "\\n")
        if (n_core == 0) { cat("No hemi-methylated sites. Empty output.\\n"); write_empty(); quit(save="no") }

        ref_atlas[CORE == TRUE, island_id := {
            d <- START - shift(START, fill = START[1])
            cumsum(d > max_gap)
        }, by = CHROM]

        island_stats <- ref_atlas[CORE == TRUE, .(
            START = min(START), END = max(START), NUM_CPG = .N
        ), by = .(CHROM, island_id)]
        island_stats[, LENGTH := END - START + 1]
        island_stats <- island_stats[LENGTH >= min_length]
        cat("  Islands after filtering:", nrow(island_stats), "\\n")
        if (nrow(island_stats) == 0) { write_empty(); quit(save="no") }

        island_stats[, `:=`(START_EXTEND = START, END_EXTEND = END)]
        result <- island_stats[, .(CHROM, START, END, START_EXTEND, END_EXTEND, LENGTH, NUM_CPG)]
        result <- result[order(CHROM, START)]
        fwrite(result, "${out_file}", sep = "\\t", col.names = TRUE)
        cat("Total islands:", nrow(result), "\\n")
        quit(save="no")
    }

    # Step 2: CORE and SHORE sites for autosomal, male_X, and Y scopes
    ref_atlas[, CORE := REF_METH < core_max_meth & PHI < core_max_phi]
    n_core <- sum(ref_atlas\$CORE)
    cat("  CORE sites:", n_core, "\\n")
    if (n_core == 0) { cat("No CORE sites. Empty output.\\n"); write_empty(); quit(save="no") }

    # SHORE: rolling join to nearest CORE
    cores <- ref_atlas[CORE == TRUE, .(CHROM, START, CORE_START = START)]
    setkey(cores, CHROM, START)
    setkey(ref_atlas, CHROM, START)
    ref_atlas[, DIST := cores[ref_atlas, abs(x.CORE_START - i.START), roll = "nearest"]]
    ref_atlas[, SHORE := REF_METH < shore_max_meth & DIST <= max_gap & CORE == FALSE]
    cat("  SHORE sites:", sum(ref_atlas\$SHORE), "\\n")

    # Step 3: Cluster CORE sites
    ref_atlas[CORE == TRUE, island_id := {
        d <- START - shift(START, fill = START[1])
        cumsum(d > max_gap)
    }, by = CHROM]

    # Step 4: Island stats and filter
    island_stats <- ref_atlas[CORE == TRUE, .(
        START = min(START), END = max(START), NUM_CPG = .N
    ), by = .(CHROM, island_id)]
    island_stats[, LENGTH := END - START + 1]
    island_stats <- island_stats[LENGTH >= min_length & (NUM_CPG * 100 >= LENGTH)]
    cat("  Islands after filtering:", nrow(island_stats), "\\n")
    if (nrow(island_stats) == 0) { write_empty(); quit(save="no") }

    # Step 5: Extend with SHORE sites
    shores <- ref_atlas[SHORE == TRUE, .(CHROM, START, END = START)]
    setkey(shores, CHROM, START, END)
    island_stats[, `:=`(s_win = START - max_gap, e_win = END + max_gap)]
    setkey(island_stats, CHROM, s_win, e_win)
    ov <- foverlaps(shores, island_stats, nomatch = NULL)
    if (nrow(ov) > 0) {
        ext <- ov[, .(S_MIN = min(i.START), E_MAX = max(i.START)), by = .(CHROM, island_id)]
        island_stats <- merge(island_stats, ext, by = c("CHROM", "island_id"), all.x = TRUE)
        island_stats[, `:=`(START_EXTEND = pmin(START, S_MIN, na.rm = TRUE),
                            END_EXTEND = pmax(END, E_MAX, na.rm = TRUE))]
    } else {
        island_stats[, `:=`(START_EXTEND = START, END_EXTEND = END)]
    }

    # Step 6: Output
    result <- island_stats[, .(CHROM, START, END, START_EXTEND, END_EXTEND, LENGTH, NUM_CPG)]
    result <- result[order(CHROM, START)]
    fwrite(result, "${out_file}", sep = "\\t", col.names = TRUE)
    cat("Total islands:", nrow(result), "\\n")
    """
}

// ============================================================================
// Process: DETECT_ICR_REGIONS
// ============================================================================

process DETECT_ICR_REGIONS {
    label 'medium'
    tag "${scope}"

    input:
        tuple val(scope), path(ref_atlas), path(icr_gene_list)

    output:
        tuple val(scope), path("icr_region.bed"), emit: icr_bed

    when:
        scope == 'autosomal'

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(data.table))
    setDTthreads(${task.cpus})

    cat("=== Detect Imprinted Control Regions (ICR) ===\\n")

    ref_lo <- ${params.icr_ref_meth_min}
    ref_hi <- ${params.icr_ref_meth_max}
    min_sites <- ${params.icr_min_sites}
    merge_distance <- ${params.dmr_merge_distance}

    normalize_primary_chrom <- function(x) {
        chrom <- as.character(x)
        base <- toupper(sub("^chr", "", chrom, ignore.case = TRUE))
        primary <- base %in% c(as.character(1:22), "X", "Y")
        chrom[primary] <- paste0("chr", base[primary])
        chrom
    }

    write_empty <- function() {
        fwrite(data.table(CHROM=character(), START=integer(), END=integer(),
                          GENE=character(), ICR_ID=character(), NUM_CPG=integer()),
               "icr_region.bed", sep="\\t", col.names=TRUE)
    }

    candidates <- fread("${icr_gene_list}")
    candidates <- candidates[, 1:min(4, ncol(candidates)), with = FALSE]
    setnames(candidates, c("CHROM", "START", "END", "GENE"))
    candidates[, CHROM := normalize_primary_chrom(CHROM)]
    candidates <- candidates[CHROM %in% paste0("chr", 1:22)]
    candidates[, `:=`(START = as.integer(START), END = as.integer(END))]
    if (nrow(candidates) == 0) {
        cat("No autosomal imprinting gene candidates. Empty output.\\n")
        write_empty()
        quit(save="no")
    }

    atlas <- fread(cmd = sprintf("zcat %s", "${ref_atlas}"),
                   select = c("CHROM", "START", "END", "N", "REF_METH"))
    atlas[, CHROM := normalize_primary_chrom(CHROM)]
    atlas <- atlas[CHROM %in% paste0("chr", 1:22)]
    atlas[, `:=`(START = as.integer(START), END = as.integer(END))]
    atlas <- atlas[N >= ${params.cpg_min_samples_per_site}]
    if (nrow(atlas) == 0) {
        cat("No autosomal atlas sites after filtering. Empty output.\\n")
        write_empty()
        quit(save="no")
    }

    setkey(candidates, CHROM, START, END)
    setkey(atlas, CHROM, START, END)
    hemi <- foverlaps(atlas, candidates, nomatch = NULL)[
        REF_METH >= ref_lo & REF_METH <= ref_hi,
        .(CHROM, START = i.START, GENE)
    ]
    cat("  Hemi-methylated ICR CpGs:", nrow(hemi), "\\n")
    if (nrow(hemi) == 0) {
        write_empty()
        quit(save="no")
    }

    setorder(hemi, GENE, CHROM, START)
    hemi[, cluster := {
        gaps <- c(Inf, diff(START))
        cumsum(gaps > merge_distance)
    }, by = .(GENE, CHROM)]

    icr <- hemi[, .(
        START = min(START),
        END = max(START) + 1L,
        NUM_CPG = .N
    ), by = .(GENE, CHROM, cluster)]
    icr <- icr[NUM_CPG >= min_sites]
    if (nrow(icr) == 0) {
        cat("No ICR passed min_sites. Empty output.\\n")
        write_empty()
        quit(save="no")
    }

    setorder(icr, CHROM, START, GENE)
    icr[, ICR_ID := paste(GENE, CHROM, sprintf("icr%02d", seq_len(.N)), sep = "_"),
        by = .(GENE, CHROM)]
    out <- icr[, .(CHROM, START, END, GENE, ICR_ID, NUM_CPG)]
    fwrite(out, "icr_region.bed", sep="\\t", col.names=TRUE)
    cat("Total ICRs:", nrow(out), "across", uniqueN(out\$GENE), "genes\\n")
    """
}

// Process: CALL_DMR
// ----------------------------------------------------------------------------
// Per-sample DMR calling: merges the sample with the reference atlas, finds
// over-threshold sites inside each panel, merges them, runs a t-test plus
// Benjamini-Hochberg FDR, and emits raw DMR records (unannotated, unformatted).
// ============================================================================

process CALL_DMR {
    label 'medium'
    tag "${meta.id}_${meta.type}_${panel_scope}"

    input:
        tuple val(meta), path(sample_bed), path(ref_atlas), path(cpg_panel),
              val(panel_scope), val(diff_threshold), val(min_sites),
              val(ref_meth_min), val(ref_meth_max)
        path icr_panel

    output:
        tuple val(meta), path("dmr_raw_${meta.id}_${meta.type}_${panel_scope}.tsv.gz"), optional: true

    script:
    def run_icr = panel_scope == 'autosomal'
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(data.table))
    setDTthreads(${task.cpus})

    # Parameters
    sample_id   <- "${meta.id}"
    sample_type <- "${meta.type}"
    panel_scope <- "${panel_scope}"
    cpg_diff_threshold    <- ${diff_threshold}
    cpg_min_sites         <- ${min_sites}
    ref_meth_min          <- ${ref_meth_min}
    ref_meth_max          <- ${ref_meth_max}
    icr_diff_threshold    <- ${params.icr_diff_threshold}
    merge_distance        <- ${params.dmr_merge_distance}
    fdr_threshold         <- ${params.dmr_fdr_threshold}
    icr_min_sites         <- ${params.icr_min_sites}
    run_icr               <- ${run_icr ? 'TRUE' : 'FALSE'}
    include_sex           <- ${params.include_sex_chromosomes ? 'TRUE' : 'FALSE'}
    output_file <- "dmr_raw_${meta.id}_${meta.type}_${panel_scope}.tsv.gz"

    normalize_primary_chrom <- function(x) {
        chrom <- as.character(x)
        base <- toupper(sub("^chr", "", chrom, ignore.case = TRUE))
        primary <- base %in% c(as.character(1:22), "X", "Y")
        chrom[primary] <- paste0("chr", base[primary])
        chrom
    }

    allowed_chroms <- paste0("chr", 1:22)
    if (include_sex) allowed_chroms <- c(allowed_chroms, "chrX", "chrY")

    message(sprintf("Processing sample: %s (%s) | scope=%s | diff>%.2f | ref_meth [%.2f, %.2f]",
                    sample_id, sample_type, panel_scope,
                    cpg_diff_threshold, ref_meth_min, ref_meth_max))

    # ---- Helper functions ----
    find_sites_in_regions <- function(sites, regions, start_col, end_col) {
        regions_dt <- copy(regions)
        setnames(regions_dt, c(start_col, end_col), c("REGION_START", "REGION_END"))
        regions_dt[, REGION_ID := .I]
        setkey(regions_dt, CHROM, REGION_START, REGION_END)
        sites_dt <- copy(sites)
        sites_dt[, END := START + 1L]
        setnames(sites_dt, "START", "START_SITE")
        setnames(sites_dt, "END", "END_SITE")
        result <- foverlaps(sites_dt, regions_dt,
                            by.x = c("CHROM", "START_SITE", "END_SITE"),
                            by.y = c("CHROM", "REGION_START", "REGION_END"),
                            nomatch = NULL)
        setnames(result, "START_SITE", "START")
        result[, END_SITE := NULL]
        result
    }

    merge_sites <- function(sites, max_gap) {
        if (nrow(sites) == 0) return(data.table())
        setorder(sites, CHROM, START)
        sites[, gap := c(Inf, diff(START)), by = CHROM]
        sites[, cluster := cumsum(gap > max_gap), by = CHROM]
        regions <- sites[, .(
            START = min(START), END = max(START) + 1L, NUM_CPG = .N,
            NUM_READS = mean(COVERAGE, na.rm = TRUE),
            MEAN_SAMPLE = mean(SAMPLE_METH, na.rm = TRUE),
            MEAN_REF = mean(REF_METH, na.rm = TRUE),
            MEAN_DIFF = mean(DIFF, na.rm = TRUE),
            DIFF_VALUES = list(DIFF)
        ), by = .(CHROM, cluster)]
        regions[, cluster := NULL]
        regions
    }

    run_ttest <- function(regions) {
        if (nrow(regions) == 0) return(regions)
        regions[, P_VALUE := sapply(DIFF_VALUES, function(diffs) {
            if (length(diffs) < 2) return(NA_real_)
            tryCatch(t.test(diffs, mu = 0)\$p.value, error = function(e) NA_real_)
        })]
        regions[, DIFF_VALUES := NULL]
        regions
    }

    # Shared post-merge stats: t-test, BH-FDR, final filter
    finalise_dmrs <- function(merged, fdr_threshold, min_sites) {
        merged <- run_ttest(merged)
        if (nrow(merged) == 0) return(data.table())
        if (any(!is.na(merged\$P_VALUE))) {
            merged[, P_FDR := p.adjust(P_VALUE, method = "BH")]
        } else {
            merged[, P_FDR := NA_real_]
        }
        merged[!is.na(P_FDR) & P_FDR < fdr_threshold & NUM_CPG > min_sites]
    }

    # ---- Load data (with chromosome filter) ----
    message("Loading sample data...")
    sample_data <- fread(cmd = sprintf("zcat %s", "${sample_bed}"),
                         col.names = c("CHROM", "START", "END", "COVERAGE", "N_METH"))
    sample_data[, CHROM := normalize_primary_chrom(CHROM)]
    sample_data <- sample_data[CHROM %in% allowed_chroms]
    sample_data[, SAMPLE_METH := N_METH / COVERAGE]

    message("Loading reference atlas...")
    ref_atlas <- fread(cmd = sprintf("zcat %s", "${ref_atlas}"),
                       select = c("CHROM", "START", "END", "REF_METH", "REF_VAR"))
    ref_atlas[, CHROM := normalize_primary_chrom(CHROM)]
    ref_atlas <- ref_atlas[CHROM %in% allowed_chroms]
    # Restrict atlas to the expected REF_METH band. Component scopes use
    # [0, 1] for DMR calling by default; female_X island selection happens
    # upstream in DETECT_CPG_ISLANDS.
    ref_atlas <- ref_atlas[REF_METH >= ref_meth_min & REF_METH <= ref_meth_max]
    message(sprintf("  Atlas after REF_METH filter: %d sites", nrow(ref_atlas)))

    message("Loading CpG panel...")
    cpg_panel <- fread("${cpg_panel}")
    cpg_panel[, CHROM := normalize_primary_chrom(CHROM)]
    cpg_panel <- cpg_panel[CHROM %in% allowed_chroms]
    # Ensure integer types for foverlaps (needed when panel is empty/header-only)
    for (col in c("START","END","START_EXTEND","END_EXTEND","LENGTH","NUM_CPG")) {
        if (col %in% names(cpg_panel)) cpg_panel[, (col) := as.integer(get(col))]
    }

    # ---- Merge sample and reference ----
    message("Merging sample and reference data...")
    setkey(sample_data, CHROM, START)
    setkey(ref_atlas, CHROM, START)
    merged <- merge(sample_data[, .(CHROM, START, COVERAGE, N_METH, SAMPLE_METH)],
                    ref_atlas[, .(CHROM, START, REF_METH, REF_VAR)],
                    by = c("CHROM", "START"))
    merged[, DIFF := SAMPLE_METH - REF_METH]
    message(sprintf("  Merged: %d sites", nrow(merged)))

    # ---- CpG panel: merge across regions within each chromosome ----
    message("Processing CpG panel...")
    cpg_sites    <- find_sites_in_regions(merged, cpg_panel, "START_EXTEND", "END_EXTEND")
    cpg_filtered <- cpg_sites[abs(DIFF) > cpg_diff_threshold]
    cpg_merged   <- merge_sites(cpg_filtered, merge_distance)
    cpg_dmrs     <- finalise_dmrs(cpg_merged, fdr_threshold, cpg_min_sites)
    message(sprintf("  CpG DMRs: %d (from %d sites)", nrow(cpg_dmrs), nrow(cpg_sites)))
    if (nrow(cpg_dmrs) > 0) {
        cpg_dmrs[, PANEL := "CpG"]
        cpg_dmrs[, DIRECTION := fifelse(MEAN_DIFF > 0, "hyper", "hypo")]
    }

    # ---- ICR panel: autosomal only; X/Y rows are ignored if present. ----
    if (run_icr) {
        message("Processing ICR panel...")
        icr_dt <- fread("${icr_panel}")
        icr_dt <- icr_dt[, 1:min(3, ncol(icr_dt)), with = FALSE]
        setnames(icr_dt, c("CHROM", "START", "END"))
        icr_dt[, CHROM := normalize_primary_chrom(CHROM)]
        panel_allowed_chroms <- paste0("chr", 1:22)
        icr_dt <- icr_dt[CHROM %in% panel_allowed_chroms]
        for (col in c("START", "END")) icr_dt[, (col) := as.integer(get(col))]
        icr_sites    <- find_sites_in_regions(merged, icr_dt, "START", "END")
        icr_filtered <- icr_sites[abs(DIFF) > icr_diff_threshold]

        if (nrow(icr_filtered) > 0) {
            icr_merged <- icr_filtered[, {
                dt <- copy(.SD); setorder(dt, START)
                dt[, gap := c(Inf, diff(START))]
                dt[, cluster := cumsum(gap > merge_distance)]
                dt[, .(START = min(START), END = max(START) + 1L, NUM_CPG = .N,
                       NUM_READS = mean(COVERAGE, na.rm = TRUE),
                       MEAN_SAMPLE = mean(SAMPLE_METH, na.rm = TRUE),
                       MEAN_REF = mean(REF_METH, na.rm = TRUE),
                       MEAN_DIFF = mean(DIFF, na.rm = TRUE),
                       DIFF_VALUES = list(DIFF)), by = cluster][, cluster := NULL]
            }, by = .(CHROM, REGION_START, REGION_END)]
            icr_merged[, c("REGION_START", "REGION_END") := NULL]
        } else {
            icr_merged <- data.table()
        }
        icr_dmrs <- finalise_dmrs(icr_merged, fdr_threshold, icr_min_sites)
        message(sprintf("  ICR DMRs: %d (from %d sites)", nrow(icr_dmrs), nrow(icr_sites)))
        if (nrow(icr_dmrs) > 0) {
            icr_dmrs[, PANEL := "ICR"]
            icr_dmrs[, DIRECTION := fifelse(MEAN_DIFF > 0, "hyper", "hypo")]
        }
    } else {
        icr_dmrs <- data.table()
    }

    # ---- Combine and emit raw (unannotated) DMRs ----
    raw_cols <- c("CHROM","START","END","PANEL","PANEL_SCOPE","NUM_CPG","NUM_READS",
                  "MEAN_SAMPLE","MEAN_REF","MEAN_DIFF","P_VALUE","P_FDR","DIRECTION")
    all_dmrs <- rbindlist(list(cpg_dmrs, icr_dmrs), use.names = TRUE, fill = TRUE)
    if (nrow(all_dmrs) > 0) {
        all_dmrs[, PANEL_SCOPE := panel_scope]
        for (col in raw_cols) if (!col %in% names(all_dmrs)) all_dmrs[, (col) := NA]
        all_dmrs <- all_dmrs[, ..raw_cols]
        setorder(all_dmrs, CHROM, START)
    } else {
        # Empty result: emit header-only file with the canonical schema. Assigning
        # NA via := on a 0x0 data.table would otherwise create a phantom row
        # (NA has length 1 and there are no existing rows to match against).
        all_dmrs <- data.table(
            CHROM = character(), START = integer(), END = integer(),
            PANEL = character(), PANEL_SCOPE = character(),
            NUM_CPG = integer(), NUM_READS = numeric(),
            MEAN_SAMPLE = numeric(), MEAN_REF = numeric(), MEAN_DIFF = numeric(),
            P_VALUE = numeric(), P_FDR = numeric(), DIRECTION = character())
    }

    message(sprintf("Total DMRs: %d (CpG: %d, ICR: %d)",
                    nrow(all_dmrs), nrow(cpg_dmrs), nrow(icr_dmrs)))
    fwrite(all_dmrs, output_file, sep = "\\t", quote = FALSE)
    """
}

// ============================================================================
// Process: ANNOTATE_DMRS
// ----------------------------------------------------------------------------
// Adds nearest-gene annotation (biotype-prioritised) and produces the final
// publishable TSV with percentage-scaled methylation columns.
// ============================================================================

process ANNOTATE_DMRS {
    label 'light'
    tag "${meta.id}_${meta.type}"

    input:
        tuple val(meta), path(dmr_raws)

    output:
        path "dmr_${meta.id}_${meta.type}.tsv", optional: true

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(data.table)
        library(ensembldb)
        library(EnsDb.Hsapiens.v86)
        library(GenomicRanges)
    })
    setDTthreads(${task.cpus})

    message(sprintf("Annotating DMRs for %s (%s) via EnsDb.Hsapiens.v86",
                    "${meta.id}", "${meta.type}"))

    # --- Gather raw DMR files ---
    raw_files <- list.files(".",
        pattern = "^dmr_raw_${meta.id}_${meta.type}_.*\\\\.tsv\\\\.gz\$",
        full.names = TRUE)
    if (length(raw_files) == 0) {
        stop("No raw DMR files found for ${meta.id}_${meta.type}")
    }
    message(sprintf("Found %d raw DMR file(s): %s",
                    length(raw_files),
                    paste(basename(raw_files), collapse = ", ")))

    dmrs <- rbindlist(lapply(raw_files, fread), use.names = TRUE, fill = TRUE)

    # --- Final output schema ---
    output_cols  <- c("CHROM","START","END","GENE","PANEL","NUM_CPG","NUM_READS",
                      "MEAN_SAMPLE","MEAN_REF","MEAN_DIFF","P_VALUE","P_FDR","DIRECTION")
    output_names <- c("CHROM","START","END","GENE","PANEL","NUM_CPG","NUM_READS",
                      "SAMPLE_METH (%)","REF_METH (%)","DIFF (%)","P_VALUE","P_FDR","DIRECTION")
    out_file <- "dmr_${meta.id}_${meta.type}.tsv"

    if (nrow(dmrs) == 0) {
        empty <- setNames(data.table::setDT(setNames(replicate(length(output_cols), character(), simplify=FALSE),
                                                     output_cols))[0], output_names)
        fwrite(empty, out_file, sep = "\\t", quote = FALSE)
        message("No DMRs to annotate; wrote empty TSV.")
        quit(save = "no", status = 0)
    }

    normalize_primary_chrom <- function(x) {
        chrom <- as.character(x)
        base <- toupper(sub("^chr", "", chrom, ignore.case = TRUE))
        primary <- base %in% c(as.character(1:22), "X", "Y")
        chrom[primary] <- paste0("chr", base[primary])
        chrom
    }
    dmrs[, CHROM := normalize_primary_chrom(CHROM)]

    # --- Cache EnsDb gene GRanges + biotype-priority rank ---
    .genes_gr <- genes(EnsDb.Hsapiens.v86,
                       columns = c("gene_id", "gene_name", "gene_biotype"),
                       return.type = "GRanges")
    .symbol_lookup <- ifelse(
        !is.na(.genes_gr\$gene_name) & nzchar(.genes_gr\$gene_name),
        .genes_gr\$gene_name,
        ifelse(!is.na(.genes_gr\$gene_id), .genes_gr\$gene_id, "")
    )
    .biotype_rank <- c(
        protein_coding = 1L,
        IG_C_gene = 2L, IG_D_gene = 2L, IG_J_gene = 2L, IG_V_gene = 2L,
        TR_C_gene = 2L, TR_D_gene = 2L, TR_J_gene = 2L, TR_V_gene = 2L,
        lncRNA = 3L, lincRNA = 3L, antisense = 3L,
        miRNA = 4L,
        snoRNA = 5L, snRNA = 5L,
        rRNA = 6L,
        Mt_tRNA = 7L, Mt_rRNA = 7L,
        misc_RNA = 8L,
        processed_transcript = 9L,
        pseudogene = 10L, processed_pseudogene = 10L, unprocessed_pseudogene = 10L,
        transcribed_processed_pseudogene = 10L, transcribed_unprocessed_pseudogene = 10L,
        IG_pseudogene = 10L, TR_pseudogene = 10L
    )
    .gene_rank <- {
        r <- .biotype_rank[.genes_gr\$gene_biotype]
        r[is.na(r)] <- 99L
        unname(r)
    }

    # --- Vectorised lookup: best gene per query (rank → overlap → length) ---
    q <- GRanges(as.character(dmrs\$CHROM),
                 IRanges(start = dmrs\$START, end = dmrs\$END))
    seqlevelsStyle(q) <- "Ensembl"
    hits <- findOverlaps(q, .genes_gr)

    gene_vec <- rep(NA_character_, length(q))
    if (length(hits) > 0L) {
        qh <- queryHits(hits); sh <- subjectHits(hits)
        ov_width <- width(pintersect(q[qh], .genes_gr[sh]))
        gene_len <- width(.genes_gr)[sh]
        rk       <- .gene_rank[sh]
        ord <- order(qh, rk, -ov_width, -gene_len)
        qh_o <- qh[ord]; sh_o <- sh[ord]
        best <- !duplicated(qh_o)
        gene_vec[qh_o[best]] <- .symbol_lookup[sh_o[best]]
    }
    dmrs[, GENE := gene_vec]

    # --- Final formatting + write ---
    for (col in output_cols) if (!col %in% names(dmrs)) dmrs[, (col) := NA]
    dmrs <- dmrs[, ..output_cols]
    dmrs[, NUM_READS   := round(NUM_READS)]
    dmrs[, MEAN_SAMPLE := round(MEAN_SAMPLE * 100, 1)]
    dmrs[, MEAN_REF    := round(MEAN_REF * 100, 1)]
    dmrs[, MEAN_DIFF   := round(MEAN_DIFF * 100, 1)]
    setnames(dmrs, output_cols, output_names)
    setorder(dmrs, CHROM, START)

    fwrite(dmrs, out_file, sep = "\\t", quote = FALSE)
    message(sprintf("Wrote %d annotated DMRs to %s", nrow(dmrs), out_file))
    """
}

// ============================================================================
// Functions
// ============================================================================

def parse_samplesheet(csv_file) {
    // Eagerly validate header, type values, sex values, and duplicate sample_ids
    // before the channel is consumed — catches user errors at launch instead of
    // deep in a process execution.
    def required_cols = ['sample_id', 'type', 'sex', 'meth_bed', 'vcf']
    def allowed_types = ['case', 'ctrl']
    def allowed_sexes = ['m', 'f']
    def sheet = file(csv_file, checkIfExists: true)
    def lines = sheet.readLines()
    if (lines.isEmpty()) {
        error "ERROR: samplesheet ${csv_file} is empty"
    }

    def header = lines[0].split(',').collect { it.trim() }
    def missing = required_cols - header
    if (missing) {
        error "ERROR: samplesheet ${csv_file} is missing required column(s): ${missing.join(', ')}. Expected header: ${required_cols.join(',')}"
    }

    def seen_ids = [] as Set
    lines.drop(1).eachWithIndex { line, idx ->
        if (!line?.trim()) return
        def cells = line.split(',').collect { it.trim() }
        def row   = [header, cells].transpose().collectEntries()
        def row_n = idx + 2
        if (!row.sample_id) {
            error "ERROR: samplesheet row ${row_n}: sample_id is empty"
        }
        if (!seen_ids.add(row.sample_id)) {
            error "ERROR: samplesheet row ${row_n}: duplicate sample_id '${row.sample_id}'"
        }
        if (!(row.type?.toLowerCase() in allowed_types)) {
            error "ERROR: samplesheet row ${row_n} (sample_id=${row.sample_id}): type must be CASE or CTRL, got '${row.type}'"
        }
        if (!(row.sex?.toLowerCase() in allowed_sexes)) {
            error "ERROR: samplesheet row ${row_n} (sample_id=${row.sample_id}): sex must be M or F, got '${row.sex}'"
        }
    }

    return channel.fromPath(csv_file, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:    row.sample_id,
                type:  row.type.toLowerCase(),
                sex:   row.sex.toUpperCase()
            ]
            def meth_bed = file(row.meth_bed, checkIfExists: true)
            def vcf      = file(row.vcf, checkIfExists: true)
            return [meta, meth_bed, vcf]
        }
}

// ============================================================================
// Workflow
// ============================================================================

workflow {
    validate_params()

    ch_input = parse_samplesheet(params.input)

    ch_icr_gene_list = channel.value(
        params.icr_gene_list
            ? file(params.icr_gene_list, checkIfExists: true)
            : file("${projectDir}/assets/icr_gene_list.bed", checkIfExists: true)
    )

    // Execution mode flags
    def build_atlas_mode = !params.ref_atlas
    def has_sex_atlases  = params.ref_atlas_male_X && params.ref_atlas_female_X && params.ref_atlas_Y
    def use_sex_branches = params.include_sex_chromosomes && (build_atlas_mode || has_sex_atlases)

    if (params.include_sex_chromosomes && !build_atlas_mode && !has_sex_atlases) {
        log.warn "include_sex_chromosomes=true but --ref_atlas is provided without " +
                 "--ref_atlas_male_X / _female_X / _Y. Falling back to autosomal-only DMR calling."
    }

    // Use the lower of the two coverage thresholds so a single
    // PREPROCESS_SAMPLE call satisfies both atlas and DMR stages.
    def min_cov = build_atlas_mode
        ? Math.min(params.atlas_min_read_coverage as int, params.dmr_min_coverage as int)
        : params.dmr_min_coverage

    ch_all_preprocessed = PREPROCESS_SAMPLE(ch_input, min_cov)

    // ====================================================================
    // Build per-scope CTRL bed channel (input to BUILD_ATLAS) and per-scope
    // CASE bed channel (input to CALL_DMR).
    //
    //   Autosomal-only (toggle off): 1 scope.
    //   Sex-stratified (toggle on):  component scopes = autosomal,
    //   male_X, female_X, and Y.
    // ====================================================================
    if (use_sex_branches) {
        def par_bed_file = params.par_bed
            ? file(params.par_bed, checkIfExists: true)
            : file("${projectDir}/assets/par_hg38.bed", checkIfExists: true)
        ch_par_bed = channel.value(par_bed_file)

        // Sex is supplied per sample in the samplesheet (M/F) and already
        // attached to meta by parse_samplesheet. Split each preprocessed BED
        // into chr1-22 / non-PAR chrX / chrY using meta.sex.
        ch_split = SPLIT_SEX_COMPONENTS(ch_all_preprocessed, ch_par_bed)

        // ---- Per-scope CTRL beds for atlas building ----
        ch_ctrl_autosomal = ch_split
            .filter { meta, _a, _x, _y -> meta.type == 'ctrl' }
            .map    { _meta, a, _x, _y -> a }
            .collect()
            .map    { beds -> ['autosomal', beds] }

        ch_ctrl_male_X = ch_split
            .filter { meta, _a, _x, _y -> meta.type == 'ctrl' && meta.sex == 'M' }
            .map    { _meta, _a, x, _y -> x }
            .collect()
            .filter { beds -> beds.size() >= (params.sex_atlas_min_controls as int) }
            .map    { beds -> ['male_X', beds] }

        ch_ctrl_female_X = ch_split
            .filter { meta, _a, _x, _y -> meta.type == 'ctrl' && meta.sex == 'F' }
            .map    { _meta, _a, x, _y -> x }
            .collect()
            .filter { beds -> beds.size() >= (params.sex_atlas_min_controls as int) }
            .map    { beds -> ['female_X', beds] }

        ch_ctrl_Y = ch_split
            .filter { meta, _a, _x, _y -> meta.type == 'ctrl' && meta.sex == 'M' }
            .map    { _meta, _a, _x, y -> y }
            .collect()
            .filter { beds -> beds.size() >= (params.sex_atlas_min_controls as int) }
            .map    { beds -> ['Y', beds] }

        ch_build_atlas_input = ch_ctrl_autosomal
            .mix(ch_ctrl_male_X)
            .mix(ch_ctrl_female_X)
            .mix(ch_ctrl_Y)

        // ---- Per-scope CASE beds for CALL_DMR ----
        ch_case_autosomal = ch_split
            .filter { meta, _a, _x, _y -> meta.type == 'case' }
            .map    { meta, a, _x, _y -> ['autosomal', meta, a] }

        ch_case_male_X = ch_split
            .filter { meta, _a, _x, _y -> meta.type == 'case' && meta.sex == 'M' }
            .map    { meta, _a, x, _y -> ['male_X', meta, x] }

        ch_case_Y = ch_split
            .filter { meta, _a, _x, _y -> meta.type == 'case' && meta.sex == 'M' }
            .map    { meta, _a, _x, y -> ['Y', meta, y] }

        ch_case_female_X = ch_split
            .filter { meta, _a, _x, _y -> meta.type == 'case' && meta.sex == 'F' }
            .map    { meta, _a, x, _y -> ['female_X', meta, x] }

        ch_case_scopes = ch_case_autosomal
            .mix(ch_case_male_X)
            .mix(ch_case_female_X)
            .mix(ch_case_Y)

    } else {
        // Autosomal-only: single 'autosomal' scope
        ch_build_atlas_input = ch_all_preprocessed
            .filter { meta, _bed -> meta.type == 'ctrl' }
            .map    { _meta, bed -> bed }
            .collect()
            .map    { beds -> ['autosomal', beds] }

        ch_case_scopes = ch_all_preprocessed
            .filter { meta, _bed -> meta.type == 'case' }
            .map    { meta, bed -> ['autosomal', meta, bed] }
    }

    // ====================================================================
    // Build / load per-scope ref atlas + CpG panel
    // ====================================================================
    if (build_atlas_mode) {
        ch_raw_atlas = BUILD_ATLAS(ch_build_atlas_input)

        ch_smooth_input = ch_raw_atlas.atlas_rds
            .combine(channel.of(params.smooth_sigma_narrow, params.smooth_sigma_wide))
        // (scope, rds, sigma)

        ch_smoothed = SMOOTH_ATLAS(ch_smooth_input)
        // ch_smoothed.smoothed: (scope, sigma, tsv.gz)

        ch_narrow_atlas_by_scope = ch_smoothed.smoothed
            .filter { _scope, sigma, _f -> (sigma as int) == (params.smooth_sigma_narrow as int) }
            .map    { scope, _sigma, f -> [scope, f] }

        // CpG panels are built per biological scope. female_X uses a
        // hemi-methylated island model inside DETECT_CPG_ISLANDS.
        ch_cpg_panels_by_scope = DETECT_CPG_ISLANDS(ch_narrow_atlas_by_scope).cpg_bed
        if (params.icr_panel) {
            ch_icr_panel = channel.value(file(params.icr_panel, checkIfExists: true))
        } else {
            ch_icr_panel = DETECT_ICR_REGIONS(
                ch_narrow_atlas_by_scope
                    .filter { scope, _f -> scope == 'autosomal' }
                    .combine(ch_icr_gene_list)
            ).icr_bed.map { _scope, f -> f }.collect().map { files -> files[0] }
        }
    } else {
        // Skip mode: user-provided references.
        // Build the per-scope atlas channel — autosomal only, or all 4 scopes
        // when sex-stratified pre-built atlases are supplied.
        def atlas_entries = [['autosomal', file(params.ref_atlas, checkIfExists: true)]]
        if (use_sex_branches) {
            atlas_entries += [
                ['male_X',   file(params.ref_atlas_male_X,   checkIfExists: true)],
                ['female_X', file(params.ref_atlas_female_X, checkIfExists: true)],
                ['Y',        file(params.ref_atlas_Y,        checkIfExists: true)]
            ]
        }
        ch_narrow_atlas_by_scope = channel.of(*atlas_entries)

        // Per-scope CpG panel: prefer user-provided; fall back to detection.
        def has_sex_cpg = params.cpg_panel_male_X && params.cpg_panel_female_X && params.cpg_panel_Y
        if (params.cpg_panel && (!use_sex_branches || has_sex_cpg)) {
            def cpg_entries = [['autosomal', file(params.cpg_panel, checkIfExists: true)]]
            if (use_sex_branches) {
                cpg_entries += [
                    ['male_X',   file(params.cpg_panel_male_X,   checkIfExists: true)],
                    ['female_X', file(params.cpg_panel_female_X, checkIfExists: true)],
                    ['Y',        file(params.cpg_panel_Y,        checkIfExists: true)]
                ]
            }
            ch_cpg_panels_by_scope = channel.of(*cpg_entries)
        } else {
            ch_cpg_panels_by_scope = DETECT_CPG_ISLANDS(ch_narrow_atlas_by_scope).cpg_bed
        }

        if (params.icr_panel) {
            ch_icr_panel = channel.value(file(params.icr_panel, checkIfExists: true))
        } else {
            ch_icr_panel = DETECT_ICR_REGIONS(
                ch_narrow_atlas_by_scope
                    .filter { scope, _f -> scope == 'autosomal' }
                    .combine(ch_icr_gene_list)
            ).icr_bed.map { _scope, f -> f }.collect().map { files -> files[0] }
        }
    }

    // (scope, atlas, panel) — one entry per scope
    ch_refs_by_scope = ch_narrow_atlas_by_scope.join(ch_cpg_panels_by_scope)

    // ====================================================================
    // Assemble CALL_DMR inputs and run
    // ====================================================================
    ch_dmr_input = ch_case_scopes
        .combine(ch_refs_by_scope, by: 0)
        // (scope, meta, sample_bed, atlas, panel)
        .map { scope, meta, bed, atlas, panel ->
            def is_female_x = scope == 'female_X'
            def diff_thresh = is_female_x ? params.sex_f_diff_threshold : params.dmr_cpg_diff_threshold
            def min_sites   = is_female_x ? params.sex_f_min_sites : params.dmr_cpg_min_sites
            def ref_min     = 0.0
            def ref_max     = 1.0
            [meta, bed, atlas, panel, scope, diff_thresh, min_sites, ref_min, ref_max]
        }

    ch_dmr_raw = CALL_DMR(ch_dmr_input, ch_icr_panel)

    // Group raw component-scope DMR files per sample, then annotate in one pass.
    // meta is stable across scopes for a given sample, so groupTuple on the
    // first slot (meta) correctly merges autosomal + X/Y records.
    ch_dmr_grouped = ch_dmr_raw.groupTuple()

    ANNOTATE_DMRS(ch_dmr_grouped)
}
