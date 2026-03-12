#!/usr/bin/env nextflow

/*
 * nf-ontmeth: ONT Methylation Analysis Pipeline
 * ==============================================
 * Detects differentially methylated regions (DMRs) from Oxford Nanopore
 * Technologies (or PacBio) whole-genome sequencing data.
 *
 * Three stages:
 *   1. Build reference atlas from control samples
 *   2. Detect hypomethylated CpG islands
 *   3. Detect DMRs per sample vs reference (CpG panel + optional panel2)
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
 * Runs via Apptainer containers on Linux, WSL2, or macOS.
 */

// ============================================================================
// Helper functions
// ============================================================================

def show_help() {
    log.info """
    =================================================================
    nf-ontmeth v${workflow.manifest.version}
    =================================================================

    Usage:
      nextflow run main.nf --input samplesheet.csv [options]

    Required:
      --input           Samplesheet CSV (sample_id, type, meth_bed, vcf)

    Optional:
      --gtf             Ensembl GTF annotation file (auto-downloaded if not provided)

    Optional panels:
      --panel2          Optional second panel BED (headerless 3-col: CHROM, START, END)

    Optional (skip Stages 1-2):
      --ref_atlas       Pre-built reference atlas (tsv.gz)
      --cpg_panel       Pre-built CpG panel (bed)

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
// Process: BUILD_ATLAS
// ============================================================================

process BUILD_ATLAS {
    label 'heavy'

    input:
        path bed_files

    output:
        path "raw_atlas.rds", emit: atlas_rds

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
    tag "sigma_${sigma}"

    input:
        path atlas_rds
        val sigma

    output:
        path "ref_atlas_${sigma}.tsv.gz", emit: smoothed

    script:
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
    standard_chroms <- c(as.character(1:22), "X", "Y", paste0("chr", c(1:22, "X", "Y")))
    dt <- dt[CHROM %in% standard_chroms]
    if (nrow(dt) == 0) stop("No standard chromosomes found")
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

    out_file <- "ref_atlas_${sigma}.tsv.gz"
    fwrite(final_output, out_file, sep = "\\t", nThread = threads)
    message(sprintf("Written %d sites to %s", nrow(final_output), out_file))
    """
}

// ============================================================================
// Process: DETECT_CPG_ISLANDS
// ============================================================================

process DETECT_CPG_ISLANDS {
    label 'medium'

    input:
        path ref_atlas

    output:
        path "cpg_signal.bed", emit: cpg_bed

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(data.table)
    })
    setDTthreads(${task.cpus})

    cat("=== Detect Hypomethylated CpG Islands ===\\n")

    # Parameters
    min_samples  <- ${params.cpg_min_samples_per_site}
    core_max_meth <- ${params.cpg_core_max_methylation}
    core_max_phi  <- ${params.cpg_core_max_dispersion}
    shore_max_meth <- ${params.cpg_shore_max_methylation}
    max_gap       <- ${params.cpg_max_gap_distance}
    min_length    <- ${params.cpg_min_island_length}

    # Empty output helper
    write_empty <- function() {
        fwrite(data.table(CHROM=character(), START=integer(), END=integer(),
                          START_EXTEND=integer(), END_EXTEND=integer(),
                          LENGTH=integer(), NUM_CPG=integer()),
               "cpg_signal.bed", sep="\\t", col.names=TRUE)
    }

    # Step 1: Load and filter atlas
    ref_atlas <- fread("${ref_atlas}")
    cat("  Loaded", nrow(ref_atlas), "sites\\n")
    ref_atlas <- ref_atlas[N > min_samples]
    cat("  After N >", min_samples, "filter:", nrow(ref_atlas), "sites\\n")

    # Step 2: CORE and SHORE sites
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
    fwrite(result, "cpg_signal.bed", sep = "\\t", col.names = TRUE)
    cat("Total islands:", nrow(result), "\\n")
    """
}

// ============================================================================
// Process: DETECT_DMR
// ============================================================================

process DETECT_DMR {
    label 'medium'
    tag "${meta.id}_${meta.type}"

    input:
        tuple val(meta), path(sample_bed)
        path ref_atlas
        path cpg_panel
        path panel2
        path gtf_file

    output:
        path "dmr_${meta.id}_${meta.type}.tsv", optional: true

    script:
    def has_panel2 = panel2.name != 'NO_FILE'
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(data.table))
    setDTthreads(${task.cpus})

    # Parameters
    sample_id   <- "${meta.id}"
    sample_type <- "${meta.type}"
    cpg_diff_threshold <- ${params.dmr_cpg_diff_threshold}
    panel2_diff_threshold <- ${params.panel2_diff_threshold}
    merge_distance <- ${params.dmr_merge_distance}
    fdr_threshold <- ${params.dmr_fdr_threshold}
    cpg_min_sites <- ${params.dmr_cpg_min_sites}
    panel2_min_sites <- ${params.panel2_min_sites}
    has_panel2 <- ${has_panel2 ? 'TRUE' : 'FALSE'}
    output_file <- "dmr_${meta.id}_${meta.type}.tsv"

    message(sprintf("Processing sample: %s (%s)", sample_id, sample_type))

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

    annotate_genes <- function(regions, gtf) {
        if (nrow(regions) == 0) { regions[, GENE := character()]; return(regions) }
        gtf_dt <- copy(gtf)
        setkey(gtf_dt, CHROM, gene_start, gene_end)
        regions_dt <- copy(regions)
        regions_dt[, REGION_ID := .I]
        setnames(regions_dt, c("START", "END"), c("region_start", "region_end"))
        overlaps <- foverlaps(regions_dt, gtf_dt,
                              by.x = c("CHROM", "region_start", "region_end"),
                              by.y = c("CHROM", "gene_start", "gene_end"),
                              nomatch = NA)
        gene_annot <- overlaps[, .(GENE = {
            valid <- !is.na(gene_name)
            if (!any(valid)) NA_character_ else gene_name[valid][which.min(priority[valid])]
        }), by = REGION_ID]
        setnames(regions_dt, c("region_start", "region_end"), c("START", "END"))
        regions_dt <- merge(regions_dt, gene_annot, by = "REGION_ID", sort = FALSE)
        regions_dt[, REGION_ID := NULL]
        regions_dt
    }

    # ---- Load data ----
    message("Loading sample data...")
    sample_data <- fread(cmd = sprintf("zcat %s", "${sample_bed}"),
                         col.names = c("CHROM", "START", "END", "COVERAGE", "N_METH"))
    sample_data[, CHROM := as.character(CHROM)]
    sample_data[, SAMPLE_METH := N_METH / COVERAGE]

    message("Loading reference atlas...")
    ref_atlas <- fread(cmd = sprintf("zcat %s", "${ref_atlas}"),
                       select = c("CHROM", "START", "END", "REF_METH", "REF_VAR"))
    ref_atlas[, CHROM := as.character(CHROM)]

    message("Loading CpG panel...")
    cpg_panel <- fread("${cpg_panel}")
    cpg_panel[, CHROM := as.character(CHROM)]
    # Ensure integer types for foverlaps (needed when panel is empty/header-only)
    for (col in c("START","END","START_EXTEND","END_EXTEND","LENGTH","NUM_CPG")) {
        if (col %in% names(cpg_panel)) cpg_panel[, (col) := as.integer(get(col))]
    }

    message("Loading GTF file...")
    gtf <- fread(cmd = sprintf("zcat %s | grep -v '^#'", "${gtf_file}"),
                 sep = "\\t", header = FALSE,
                 col.names = c("CHROM","source","feature","start","end","score","strand","frame","attributes"))
    gtf <- gtf[feature == "gene"]
    gtf[, CHROM := as.character(CHROM)]
    gtf[, gene_name := sub('.*gene_name "([^"]+)".*', "\\\\1", attributes)]
    gtf[, gene_biotype := sub('.*gene_biotype "([^"]+)".*', "\\\\1", attributes)]
    gtf[, c("source","feature","score","strand","frame","attributes") := NULL]
    setnames(gtf, c("start","end"), c("gene_start","gene_end"))
    biotype_priority <- c("protein_coding"=1,"lncRNA"=2,"processed_pseudogene"=3)
    gtf[, priority := fifelse(gene_biotype %in% names(biotype_priority), biotype_priority[gene_biotype], 99L)]

    # ---- Merge sample and reference ----
    message("Merging sample and reference data...")
    setkey(sample_data, CHROM, START)
    setkey(ref_atlas, CHROM, START)
    merged <- merge(sample_data[, .(CHROM, START, COVERAGE, N_METH, SAMPLE_METH)],
                    ref_atlas[, .(CHROM, START, REF_METH, REF_VAR)],
                    by = c("CHROM", "START"))
    merged[, DIFF := SAMPLE_METH - REF_METH]
    message(sprintf("  Merged: %d sites", nrow(merged)))

    # ---- Process CpG panel ----
    message("Processing CpG panel...")
    cpg_sites <- find_sites_in_regions(merged, cpg_panel, "START_EXTEND", "END_EXTEND")
    message(sprintf("  Sites in CpG regions: %d", nrow(cpg_sites)))
    cpg_filtered <- cpg_sites[abs(DIFF) > cpg_diff_threshold]
    cpg_merged <- merge_sites(cpg_filtered, merge_distance)
    cpg_merged <- run_ttest(cpg_merged)
    if (nrow(cpg_merged) > 0) {
        if (any(!is.na(cpg_merged\$P_VALUE))) {
            cpg_merged[, P_FDR := p.adjust(P_VALUE, method = "BH")]
        } else {
            cpg_merged[, P_FDR := NA_real_]
        }
        cpg_dmrs <- cpg_merged[!is.na(P_FDR) & P_FDR < fdr_threshold & NUM_CPG > cpg_min_sites]
    } else {
        cpg_dmrs <- data.table()
    }
    message(sprintf("  CpG DMRs: %d", nrow(cpg_dmrs)))
    if (nrow(cpg_dmrs) > 0) {
        cpg_dmrs <- annotate_genes(cpg_dmrs, gtf)
        cpg_dmrs[, PANEL := "CpG"]
        cpg_dmrs[, DIRECTION := fifelse(MEAN_DIFF > 0, "hyper", "hypo")]
    }

    # ---- Process panel2 (optional) ----
    if (has_panel2) {
        message("Processing panel2...")
        panel2_dt <- fread("${panel2}", header = FALSE, col.names = c("CHROM", "START", "END"))
        panel2_dt[, CHROM := as.character(CHROM)]
        p2_sites <- find_sites_in_regions(merged, panel2_dt, "START", "END")
        message(sprintf("  Sites in panel2 regions: %d", nrow(p2_sites)))
        p2_filtered <- p2_sites[abs(DIFF) > panel2_diff_threshold]

        # Merge within each panel2 region separately
        if (nrow(p2_filtered) > 0) {
            p2_merged <- p2_filtered[, {
                dt <- copy(.SD); setorder(dt, START)
                dt[, gap := c(Inf, diff(START))]
                dt[, cluster := cumsum(gap > merge_distance)]
                dt[, .(START = min(START), END = max(START) + 1L, NUM_CPG = .N,
                       NUM_READS = mean(COVERAGE, na.rm=TRUE),
                       MEAN_SAMPLE = mean(SAMPLE_METH, na.rm=TRUE),
                       MEAN_REF = mean(REF_METH, na.rm=TRUE),
                       MEAN_DIFF = mean(DIFF, na.rm=TRUE),
                       DIFF_VALUES = list(DIFF)), by = cluster][, cluster := NULL]
            }, by = .(CHROM, REGION_START, REGION_END)]
            p2_merged[, c("REGION_START", "REGION_END") := NULL]
        } else {
            p2_merged <- data.table()
        }

        p2_merged <- run_ttest(p2_merged)
        if (nrow(p2_merged) > 0) {
            if (any(!is.na(p2_merged\$P_VALUE))) {
                p2_merged[, P_FDR := p.adjust(P_VALUE, method = "BH")]
            } else {
                p2_merged[, P_FDR := NA_real_]
            }
            p2_dmrs <- p2_merged[!is.na(P_FDR) & P_FDR < fdr_threshold & NUM_CPG > panel2_min_sites]
        } else {
            p2_dmrs <- data.table()
        }
        message(sprintf("  Panel2 DMRs: %d", nrow(p2_dmrs)))
        if (nrow(p2_dmrs) > 0) {
            p2_dmrs <- annotate_genes(p2_dmrs, gtf)
            p2_dmrs[, PANEL := "Panel2"]
            p2_dmrs[, DIRECTION := fifelse(MEAN_DIFF > 0, "hyper", "hypo")]
        }
    } else {
        p2_dmrs <- data.table()
    }

    # ---- Combine and output ----
    all_dmrs <- rbindlist(list(cpg_dmrs, p2_dmrs), use.names = TRUE, fill = TRUE)
    output_cols <- c("CHROM","START","END","GENE","PANEL","NUM_CPG","NUM_READS",
                     "MEAN_SAMPLE","MEAN_REF","MEAN_DIFF","P_VALUE","P_FDR","DIRECTION")
    output_names <- c("CHROM","START","END","GENE","PANEL","NUM_CPG","NUM_READS",
                      "SAMPLE_METH (%)","REF_METH (%)","DIFF (%)","P_VALUE","P_FDR","DIRECTION")
    for (col in output_cols) if (!col %in% names(all_dmrs)) all_dmrs[, (col) := NA]
    all_dmrs <- all_dmrs[, ..output_cols]
    if (nrow(all_dmrs) > 0) {
        all_dmrs[, NUM_READS := round(NUM_READS)]
        all_dmrs[, MEAN_SAMPLE := round(MEAN_SAMPLE * 100, 1)]
        all_dmrs[, MEAN_REF := round(MEAN_REF * 100, 1)]
        all_dmrs[, MEAN_DIFF := round(MEAN_DIFF * 100, 1)]
    }
    setnames(all_dmrs, output_cols, output_names)
    setorder(all_dmrs, CHROM, START)
    message(sprintf("Total DMRs: %d (CpG: %d, Panel2: %d)", nrow(all_dmrs), nrow(cpg_dmrs), nrow(p2_dmrs)))
    fwrite(all_dmrs, output_file, sep = "\\t", quote = FALSE)
    message("Output written to: ", output_file)
    """
}

// ============================================================================
// Process: DOWNLOAD_GTF
// ============================================================================

process DOWNLOAD_GTF {
    label 'light'

    output:
        path "*.gtf.gz", emit: gtf

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Fetch the latest Ensembl release number
    RELEASE=\$(wget -qO- 'https://rest.ensembl.org/info/data/?content-type=application/json' \
        | grep -oP '"releases":\\[\\K[0-9]+')
    echo "Latest Ensembl release: \${RELEASE}" >&2

    URL="https://ftp.ensembl.org/pub/release-\${RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.\${RELEASE}.gtf.gz"
    echo "Downloading: \${URL}" >&2
    wget -q "\${URL}"
    """
}

// ============================================================================
// Functions
// ============================================================================

def parse_samplesheet(csv_file) {
    return channel.fromPath(csv_file, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:    row.sample_id,
                type:  row.type.toLowerCase()
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
    // Validate parameters
    validate_params()

    // Parse samplesheet
    ch_input = parse_samplesheet(params.input)

    ch_all = ch_input

    // Optional panel2 channel (value channel — reused across samples)
    ch_panel2 = params.panel2
        ? channel.value(file(params.panel2, checkIfExists: true))
        : channel.value(file("${projectDir}/assets/NO_FILE"))

    // Determine mode
    def build_atlas_mode = !params.ref_atlas

    // Use the lower of the two coverage thresholds so a single
    // PREPROCESS_SAMPLE call satisfies both atlas and DMR stages.
    def min_cov = build_atlas_mode
        ? Math.min(params.atlas_min_read_coverage as int, params.dmr_min_coverage as int)
        : params.dmr_min_coverage

    // Preprocess ALL samples once
    ch_all_preprocessed = PREPROCESS_SAMPLE(ch_all, min_cov)

    // GTF: use provided file or auto-download from Ensembl
    if (params.gtf) {
        ch_gtf = channel.value(file(params.gtf, checkIfExists: true))
    } else {
        ch_gtf = DOWNLOAD_GTF().gtf
    }

    if (build_atlas_mode) {
        // ====================================================================
        // Stages 1-2: Build atlas + CpG panel
        // ====================================================================

        // Filter CTRL preprocessed beds for atlas building
        ch_ctrl_beds = ch_all_preprocessed
            .filter { meta, _bed -> meta.type == 'ctrl' }
            .map { _meta, bed -> bed }
            .collect()

        ch_raw_atlas = BUILD_ATLAS(ch_ctrl_beds)

        // Smooth at both bandwidths in a single process call
        ch_sigmas = channel.of(params.smooth_sigma_narrow, params.smooth_sigma_wide)
        ch_smoothed = SMOOTH_ATLAS(
            ch_raw_atlas.atlas_rds,
            ch_sigmas
        )

        // Stage 2: Detect CpG islands (narrow atlas only)
        ch_atlas_narrow = ch_smoothed.smoothed
            .filter { f -> f.name.contains("_${params.smooth_sigma_narrow}.") }

        ch_cpg = DETECT_CPG_ISLANDS(ch_atlas_narrow)

        // Set references for Stage 3 (value channels for reuse across samples)
        ch_ref_atlas = ch_atlas_narrow.collect().map { files -> files[0] }
        ch_cpg_panel = ch_cpg.cpg_bed.collect().map { files -> files[0] }

    } else {
        // ====================================================================
        // Skip Stages 1-2: Use provided references
        // ====================================================================

        ch_ref_atlas = channel.value(file(params.ref_atlas, checkIfExists: true))
        ch_cpg_panel = channel.value(file(params.cpg_panel, checkIfExists: true))
    }

    // ====================================================================
    // Stage 3: DMR detection
    // ====================================================================

    DETECT_DMR(
        ch_all_preprocessed,
        ch_ref_atlas,
        ch_cpg_panel,
        ch_panel2,
        ch_gtf
    )
}
