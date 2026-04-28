#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[[hit + 1]]
}

atlas_file <- get_arg("--atlas", "/lab_nas_home_user/Mark/ont-methylation/results/atlas/ref_atlas_150.tsv.gz")
icr_gene_file <- get_arg("--icr-gene-list", "/lab_nas_home_user/Mark/ont-methylation/results/icr/imprinted_gene.bed")
outdir <- get_arg("--outdir", "assets")
ref_lo <- as.numeric(get_arg("--ref-min", "0.40"))
ref_hi <- as.numeric(get_arg("--ref-max", "0.60"))
min_sites <- as.integer(get_arg("--min-sites", "5"))
merge_distance <- as.integer(get_arg("--merge-distance", "150"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

normalize_primary_chrom <- function(x) {
  chrom <- as.character(x)
  base <- toupper(sub("^chr", "", chrom, ignore.case = TRUE))
  primary <- base %in% c(as.character(1:22), "X", "Y")
  chrom[primary] <- paste0("chr", base[primary])
  chrom
}

message("Reading imprinting gene candidates: ", icr_gene_file)
genes <- fread(icr_gene_file)
genes <- genes[, 1:min(4, ncol(genes)), with = FALSE]
setnames(genes, c("CHROM", "START", "END", "GENE"))
genes[, CHROM := normalize_primary_chrom(CHROM)]
genes <- genes[CHROM %in% paste0("chr", 1:22)]
genes[, `:=`(START = as.integer(START), END = as.integer(END))]
setorder(genes, CHROM, START, END, GENE)
fwrite(genes, file.path(outdir, "icr_gene_list.bed"), sep = "\t", col.names = TRUE)
message("Wrote ", nrow(genes), " imprinting gene candidate intervals")

message("Reading atlas: ", atlas_file)
atlas <- fread(cmd = sprintf("zcat %s", shQuote(atlas_file)),
               select = c("CHROM", "START", "END", "N", "REF_METH"))
atlas[, CHROM := normalize_primary_chrom(CHROM)]
atlas <- atlas[CHROM %in% paste0("chr", 1:22)]
atlas[, `:=`(START = as.integer(START), END = as.integer(END))]

setkey(genes, CHROM, START, END)
setkey(atlas, CHROM, START, END)
hemi <- foverlaps(atlas, genes, nomatch = NULL)[
  REF_METH >= ref_lo & REF_METH <= ref_hi,
  .(CHROM, START = i.START, GENE)
]
message("Hemi-methylated CpGs in imprinting gene candidates: ", nrow(hemi))

if (nrow(hemi) == 0) {
  icr <- data.table(CHROM=character(), START=integer(), END=integer(),
                    GENE=character(), ICR_ID=character(), NUM_CPG=integer())
} else {
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
  setorder(icr, CHROM, START, GENE)
  icr[, ICR_ID := paste(GENE, CHROM, sprintf("icr%02d", seq_len(.N)), sep = "_"),
      by = .(GENE, CHROM)]
  icr <- icr[, .(CHROM, START, END, GENE, ICR_ID, NUM_CPG)]
}

fwrite(icr, file.path(outdir, "icr_region.bed"), sep = "\t", col.names = TRUE)
message("Wrote ", nrow(icr), " ICRs across ", uniqueN(icr$GENE), " genes")
