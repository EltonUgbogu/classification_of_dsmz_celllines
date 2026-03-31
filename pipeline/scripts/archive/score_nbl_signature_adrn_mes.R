#!/usr/bin/env Rscript
#
# score_nbl_signature_adrn_mes.R
#
# Purpose: Score neuroblastoma (NBL) samples along multiple biological axes:
#   1. Custom NBL neuronal differentiation signature (from ENSG gene list)
#   2. ADRN (adrenergic/noradrenergic) cell state score
#   3. MES (mesenchymal-like) cell state score
#   4. Delta score (MES - ADRN) representing position on the ADRN↔MES axis
#
# Background: NBL tumors exhibit cellular plasticity between ADRN and MES states,
# which affects differentiation, drug response, and prognosis. This script quantifies
# where each sample falls on these biological axes.
#
# Input: 
#   - counts.tsv: gene-by-sample count matrix (gene_id + samples)
#   - nbl_features.txt: NBL signature gene list (Ensembl IDs)
#
# Output:
#   - *_scores.tsv: per-sample scores for all axes
#   - *_summary.tsv: statistics including correlations and gene overlap counts
#   - *_plots.pdf: diagnostic scatter plots
#   - *_ADRN_markers_symbols.txt: ADRN marker gene list used
#   - *_MES_markers_symbols.txt: MES marker gene list used
#   - *.log: execution log with timestamps
#
# Method: DESeq2 VST transformation → gene-wise z-scores → mean score per module
#

# Check for required packages and stop early if missing
suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Missing DESeq2", call.=FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing ggplot2", call.=FALSE)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Missing AnnotationDbi", call.=FALSE)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Missing org.Hs.eg.db", call.=FALSE)
})

library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)

# =============================
# Configuration: File paths
# =============================
# Input: gene-by-sample count matrix
expr_file <- "/work/ugbogu/pipeline/results/unsupervised/nbl/deseq2_inputs/counts.tsv"

# Input: NBL signature gene list (one Ensembl gene ID per line)
sig_file  <- "/work/ugbogu/pipeline/results/unsupervised/nbl/deseq2_markers/markers/nbl_features.txt"

# Output directory and file prefix
out_dir    <- "/work/ugbogu/pipeline/results/unsupervised/nbl/gsea_ova/nbl_signature_axes"
out_prefix <- file.path(out_dir, "nbl_signature_adrn_mes")

# Option: write out the marker lists used for reproducibility
write_marker_lists <- TRUE

# =============================
# ADRN / MES State Marker Genes
# =============================
# Core "state marker" gene programs widely used in ADRN↔MES neuroblastoma research.
# These lists are kept modest (~15-17 genes each) for robustness across datasets.
# 
# ADRN markers: genes involved in adrenergic/noradrenergic neuron identity and
# catecholamine biosynthesis/transport. High ADRN score indicates differentiated,
# noradrenergic-like tumor cells.
ADRN_SYMBOLS <- c(
  "PHOX2A","PHOX2B",     # Transcription factors for noradrenergic differentiation
  "DBH","TH","DDC",       # Catecholamine biosynthesis enzymes
  "SLC6A2","SLC18A1",     # Norepinephrine transporter, vesicular monoamine transporter
  "CHGA","CHGB",          # Chromogranins (neuroendocrine markers)
  "GATA3","HAND2","ISL1", # Core regulatory transcription factors
  "ASCL1",                # Neuronal differentiation TF
  "ALK","CHRNA3"          # Additional ADRN-associated genes
)

# MES markers: genes involved in mesenchymal cell identity, extracellular matrix,
# cell motility, and epithelial-mesenchymal transition (EMT). High MES score
# indicates undifferentiated, invasive tumor cells.
MES_SYMBOLS <- c(
  "PRRX1",                           # Mesenchymal transcription factor
  "FN1","VIM",                       # Fibronectin, vimentin (ECM/cytoskeleton)
  "COL1A1","COL1A2","COL3A1","COL5A1", # Collagens (ECM components)
  "ITGA5","ITGB1",                   # Integrins (cell adhesion)
  "CDH2",                            # N-cadherin (mesenchymal adhesion)
  "SNAI2","TWIST1","ZEB1",          # EMT transcription factors
  "YAP1","WWTR1",                    # Hippo pathway effectors (proliferation)
  "FOSL1","JUN"                      # AP-1 transcription factors
)

# =============================
# Set up logging
# =============================
# Create output directory if needed
dir_create <- function(path) if (!dir.exists(path)) dir.create(path, recursive=TRUE)
dir_create(out_dir)

# Define log file path with timestamp
log_file <- paste0(out_prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")

# Open log file connection
log_con <- file(log_file, open = "wt")

# Custom logging function: writes to both console and log file
log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  message(msg)  # Print to console
  # Write to log file with timestamp
  writeLines(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg), con = log_con)
  flush(log_con)  # Force immediate write to disk
}

# Wrapper for warnings that also logs
log_warning <- function(...) {
  msg <- paste0("WARNING: ", ..., collapse = "")
  warning(msg, call. = FALSE)
  writeLines(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg), con = log_con)
  flush(log_con)
}

# Ensure log file is closed on exit (even if script errors)
on.exit({
  if (exists("log_con") && isOpen(log_con)) {
    log_msg("Script completed.")
    close(log_con)
    message("Log saved to: ", log_file)
  }
}, add = TRUE)

# Log script start
log_msg("========================================")
log_msg("NBL Signature ADRN/MES Analysis")
log_msg("========================================")
log_msg("Script started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_msg("Log file: ", log_file)
log_msg("")

# =============================
# Helper functions
# =============================

# Validate file exists, stop with error if not found
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    log_msg("ERROR: File not found: ", path)
    stop("File not found: ", path, call.=FALSE)
  }
}

# Read expression matrix from TSV file
# Handles two formats: (1) explicit gene_id column, (2) genes as first column
read_expr_tsv <- function(path) {
  x <- read.table(path, header=TRUE, sep="\t", quote="", comment.char="",
                  stringsAsFactors=FALSE, check.names=FALSE)
  if ("gene_id" %in% colnames(x)) {
    rownames(x) <- x$gene_id
    x$gene_id <- NULL
  } else {
    rownames(x) <- x[[1]]
    x[[1]] <- NULL
  }
  m <- as.matrix(x)
  mode(m) <- "numeric"
  m
}

# Remove Ensembl version suffix (e.g., ENSG00000012345.2 → ENSG00000012345)
strip_ensg_version <- function(x) sub("\\.[0-9]+$", "", x)

# Compute row-wise z-scores (standardize each gene across samples)
zscore_rows <- function(mat) {
  # gene-wise z-score: (x - mean) / sd for each row
  t(scale(t(mat)))
}

# Score a gene set using mean z-score approach
# Handles Ensembl IDs with or without version suffixes
score_gene_set_ensg <- function(expr_mat, ensg_vec) {
  # Match genes by stripping version suffixes
  rn0 <- strip_ensg_version(rownames(expr_mat))
  hit_rows <- rownames(expr_mat)[rn0 %in% ensg_vec]

  # Need at least 3 genes for meaningful score
  if (length(hit_rows) < 3) {
    return(list(score = rep(NA_real_, ncol(expr_mat)), genes_present = character(0)))
  }
  
  # Z-score genes, then average across genes for each sample
  ez <- zscore_rows(expr_mat[hit_rows, , drop=FALSE])
  list(score = colMeans(ez, na.rm=TRUE), genes_present = hit_rows)
}

# Convert gene symbols to Ensembl IDs using org.Hs.eg.db annotation
symbols_to_ensg <- function(symbols) {
  symbols <- unique(trimws(symbols))
  symbols <- symbols[nzchar(symbols)]
  
  # Map SYMBOL → ENSEMBL via Bioconductor annotation package
  m <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    columns = c("ENSEMBL"),
    keytype = "SYMBOL"
  )
  
  # Filter out missing mappings
  m <- m[!is.na(m$ENSEMBL) & nzchar(m$ENSEMBL), , drop=FALSE]
  unique(m$ENSEMBL)
}

# Compute Spearman correlation with safe error handling
spearman_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(list(rho=NA_real_, p=NA_real_, n=sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# =============================
# Load and validate input data
# =============================
log_msg("Checking input files...")
stop_if_missing(expr_file)
stop_if_missing(sig_file)
log_msg("All input files found.")
log_msg("")

log_msg("Reading count matrix: ", expr_file)
counts <- read_expr_tsv(expr_file)
samples <- colnames(counts)
log_msg("Expression matrix dimensions: ", nrow(counts), " genes x ", length(samples), " samples")

if (length(samples) < 3) {
  log_msg("ERROR: Expected >= 3 samples, found: ", length(samples))
  stop("Expected >= 3 samples, found: ", length(samples), call.=FALSE)
}
log_msg("")

log_msg("Reading NBL signature ENSG list: ", sig_file)
sig_ensg <- unique(trimws(readLines(sig_file, warn=FALSE)))
sig_ensg <- sig_ensg[nzchar(sig_ensg)]
log_msg("NBL signature genes loaded: ", length(sig_ensg))
log_msg("")

# =============================
# Transform count data (VST)
# =============================
log_msg("Applying DESeq2 VST transform...")
# Create minimal sample metadata (single condition)
coldata <- data.frame(row.names=samples, condition=factor(rep("nbl", length(samples))))

# Build DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=~1)

# Filter low-count genes for VST stability
keep <- rowSums(counts(dds)) >= 10
log_msg("Filtering low-count genes: ", sum(keep), " of ", length(keep), " genes retained (rowSum >= 10)")
dds <- dds[keep, ]

# Apply variance-stabilizing transformation
expr <- assay(vst(dds, blind=TRUE))
log_msg("VST transformation complete.")
log_msg("Transform method: DESeq2::vst (blind=TRUE; genes with rowSum>=10 kept)")
log_msg("")

# =============================
# Compute biological axis scores
# =============================
log_msg("Computing biological axis scores...")
log_msg("")

# AXIS 1: NBL neuronal signature
log_msg("AXIS 1: NBL neuronal signature")
sig_res <- score_gene_set_ensg(expr, sig_ensg)
sig_score <- sig_res$score
log_msg("  Signature genes provided: ", length(sig_ensg))
log_msg("  Signature genes matched in expr: ", length(sig_res$genes_present))
log_msg("  Signature score computed for ", length(sig_score), " samples")
log_msg("")

# AXIS 2: ADRN (adrenergic) state
log_msg("AXIS 2: ADRN (adrenergic/noradrenergic) state")
log_msg("  ADRN marker symbols provided: ", length(ADRN_SYMBOLS))

# Convert gene symbols to Ensembl IDs
adrn_ensg <- symbols_to_ensg(ADRN_SYMBOLS)
log_msg("  ADRN markers mapped to ENSG IDs: ", length(adrn_ensg))

# Compute ADRN score
adrn_res <- score_gene_set_ensg(expr, adrn_ensg)
adrn_score <- adrn_res$score
log_msg("  ADRN genes matched in expr: ", length(adrn_res$genes_present))
log_msg("  ADRN score computed for ", length(adrn_score), " samples")
log_msg("")

# AXIS 3: MES (mesenchymal) state
log_msg("AXIS 3: MES (mesenchymal-like) state")
log_msg("  MES marker symbols provided: ", length(MES_SYMBOLS))

# Convert gene symbols to Ensembl IDs
mes_ensg <- symbols_to_ensg(MES_SYMBOLS)
log_msg("  MES markers mapped to ENSG IDs: ", length(mes_ensg))

# Compute MES score
mes_res <- score_gene_set_ensg(expr, mes_ensg)
mes_score <- mes_res$score
log_msg("  MES genes matched in expr: ", length(mes_res$genes_present))
log_msg("  MES score computed for ", length(mes_score), " samples")
log_msg("")

# AXIS 4: Delta score (MES - ADRN)
# Positive delta = more MES-like, negative delta = more ADRN-like
log_msg("AXIS 4: Delta score (MES - ADRN)")
delta_score <- mes_score - adrn_score
log_msg("  Delta score computed for ", length(delta_score), " samples")
log_msg("  Interpretation: positive = MES-shifted, negative = ADRN-shifted")
log_msg("")

# =============================
# Statistical tests
# =============================
log_msg("Running statistical tests (Spearman correlations)...")
log_msg("")

# TEST 1: NBL signature vs ADRN state
log_msg("TEST 1: NBL signature vs ADRN score")
stat_sig_adrn <- spearman_safe(sig_score, adrn_score)
log_msg("  rho = ", round(stat_sig_adrn$rho, 3), ", p = ", format(stat_sig_adrn$p, digits = 3), 
        ", n = ", stat_sig_adrn$n)
log_msg("")

# TEST 2: NBL signature vs MES state
log_msg("TEST 2: NBL signature vs MES score")
stat_sig_mes <- spearman_safe(sig_score, mes_score)
log_msg("  rho = ", round(stat_sig_mes$rho, 3), ", p = ", format(stat_sig_mes$p, digits = 3), 
        ", n = ", stat_sig_mes$n)
log_msg("")

# TEST 3: NBL signature vs Delta (MES-ADRN axis)
log_msg("TEST 3: NBL signature vs Delta (MES - ADRN)")
stat_sig_delta <- spearman_safe(sig_score, delta_score)
log_msg("  rho = ", round(stat_sig_delta$rho, 3), ", p = ", format(stat_sig_delta$p, digits = 3), 
        ", n = ", stat_sig_delta$n)
log_msg("")

# TEST 4: ADRN vs MES (sanity check - should be negative)
log_msg("TEST 4: ADRN vs MES (sanity check - expect negative correlation)")
stat_adrn_mes <- spearman_safe(adrn_score, mes_score)
log_msg("  rho = ", round(stat_adrn_mes$rho, 3), ", p = ", format(stat_adrn_mes$p, digits = 3), 
        ", n = ", stat_adrn_mes$n)
log_msg("")

# =============================
# Write output files
# =============================
log_msg("Writing output files...")

# OUTPUT 1: Per-sample scores table
scores <- data.frame(
  sample = samples,
  nbl_signature = as.numeric(sig_score),    # NBL differentiation signature
  ADRN = as.numeric(adrn_score),            # Adrenergic state score
  MES  = as.numeric(mes_score),             # Mesenchymal state score
  MES_minus_ADRN = as.numeric(delta_score), # Position on ADRN↔MES axis
  stringsAsFactors = FALSE
)

out_scores <- paste0(out_prefix, "_scores.tsv")
write.table(scores, out_scores, sep="\t", quote=FALSE, row.names=FALSE)
log_msg("  Wrote scores table: ", out_scores)

# OUTPUT 2: Summary statistics and metadata
summary <- data.frame(
  key = c(
    "expr_file","sig_file","transform",
    "sig_genes_total","sig_genes_matched",
    "adrn_symbols_n","adrn_ensg_n","adrn_genes_matched",
    "mes_symbols_n","mes_ensg_n","mes_genes_matched",
    "spearman_sig_vs_ADRN_rho","spearman_sig_vs_ADRN_p","spearman_sig_vs_ADRN_n",
    "spearman_sig_vs_MES_rho","spearman_sig_vs_MES_p","spearman_sig_vs_MES_n",
    "spearman_sig_vs_DELTA_rho","spearman_sig_vs_DELTA_p","spearman_sig_vs_DELTA_n",
    "spearman_ADRN_vs_MES_rho","spearman_ADRN_vs_MES_p","spearman_ADRN_vs_MES_n"
  ),
  value = c(
    expr_file, sig_file, "DESeq2::vst(blind=TRUE; rowSum>=10 kept)",
    length(sig_ensg), length(sig_res$genes_present),
    length(ADRN_SYMBOLS), length(adrn_ensg), length(adrn_res$genes_present),
    length(MES_SYMBOLS), length(mes_ensg), length(mes_res$genes_present),
    stat_sig_adrn$rho, stat_sig_adrn$p, stat_sig_adrn$n,
    stat_sig_mes$rho,  stat_sig_mes$p,  stat_sig_mes$n,
    stat_sig_delta$rho,stat_sig_delta$p,stat_sig_delta$n,
    stat_adrn_mes$rho, stat_adrn_mes$p, stat_adrn_mes$n
  ),
  stringsAsFactors = FALSE
)

out_summary <- paste0(out_prefix, "_summary.tsv")
write.table(summary, out_summary, sep="\t", quote=FALSE, row.names=FALSE)
log_msg("  Wrote summary table: ", out_summary)

# OUTPUT 3: Write out marker gene lists for reproducibility (optional)
if (write_marker_lists) {
  adrn_markers_file <- paste0(out_prefix, "_ADRN_markers_symbols.txt")
  mes_markers_file <- paste0(out_prefix, "_MES_markers_symbols.txt")
  
  writeLines(ADRN_SYMBOLS, adrn_markers_file)
  writeLines(MES_SYMBOLS, mes_markers_file)
  
  log_msg("  Wrote ADRN marker list: ", adrn_markers_file)
  log_msg("  Wrote MES marker list: ", mes_markers_file)
}
log_msg("")

# =============================
# Generate diagnostic plots
# =============================
log_msg("Generating diagnostic plots...")

out_pdf <- paste0(out_prefix, "_plots.pdf")
pdf(out_pdf, width=7, height=6)

df <- scores

# Helper function to create scatter plots with correlation stats
p_scatter <- function(x, y, xlab, ylab, title, stat) {
  ggplot(df, aes(x=.data[[x]], y=.data[[y]])) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE) +  # Linear regression line
    theme_minimal() +
    labs(
      title = title,
      subtitle = sprintf("Spearman rho=%.3f, p=%.3g (n=%d)", stat$rho, stat$p, stat$n),
      x = xlab, y = ylab
    )
}

# PLOT 1: NBL signature vs ADRN state
print(p_scatter("ADRN", "nbl_signature",
                "ADRN score (z-mean)", "NBL signature score (z-mean)",
                "NBL signature vs ADRN state", stat_sig_adrn))
log_msg("  Plot 1: NBL signature vs ADRN (scatter)")

# PLOT 2: NBL signature vs MES state
print(p_scatter("MES", "nbl_signature",
                "MES score (z-mean)", "NBL signature score (z-mean)",
                "NBL signature vs MES state", stat_sig_mes))
log_msg("  Plot 2: NBL signature vs MES (scatter)")

# PLOT 3: NBL signature vs Delta (MES-ADRN axis)
print(p_scatter("MES_minus_ADRN", "nbl_signature",
                "MES − ADRN", "NBL signature score (z-mean)",
                "NBL signature vs MES−ADRN axis", stat_sig_delta))
log_msg("  Plot 3: NBL signature vs Delta (scatter)")

# PLOT 4: MES vs ADRN (sanity check - should show negative correlation)
print(p_scatter("MES", "ADRN",
                "MES score (z-mean)", "ADRN score (z-mean)",
                "MES vs ADRN (sanity check)", stat_adrn_mes))
log_msg("  Plot 4: MES vs ADRN sanity check (scatter)")

dev.off()
log_msg("  Wrote plots PDF: ", out_pdf)
log_msg("")

log_msg("========================================")
log_msg("Analysis complete!")
log_msg("========================================")
log_msg("Output files:")
log_msg("  - Scores: ", out_scores)
log_msg("  - Summary: ", out_summary)
log_msg("  - Plots: ", out_pdf)
if (write_marker_lists) {
  log_msg("  - ADRN markers: ", paste0(out_prefix, "_ADRN_markers_symbols.txt"))
  log_msg("  - MES markers: ", paste0(out_prefix, "_MES_markers_symbols.txt"))
}