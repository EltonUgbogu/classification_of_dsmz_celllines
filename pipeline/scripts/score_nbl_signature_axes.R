#!/usr/bin/env Rscript
#
# score_nbl_signature_axes.R
#
# Purpose: Analyze neuroblastoma (NBL) gene expression data along three biological axes:
#   1. NBL neuronal differentiation signature (from a custom gene list)
#   2. MYCN oncogene expression and its relationship to the signature
#   3. Cell-cycle/proliferation activity and its relationship to the signature
#
# Outputs:
#   - TSV file with per-sample scores for all three axes
#   - Summary TSV with correlation statistics and metadata
#   - PDF with diagnostic plots (signature vs MYCN, signature vs cell-cycle)
#   - Log file with all console output
#
# NOTE: Input is a count matrix. The script transforms counts using DESeq2's
# variance-stabilizing transformation (VST) only; DESeq2 is required.

# Suppress startup messages and check for optional dependencies
suppressPackageStartupMessages({
  ok_deseq2 <- requireNamespace("DESeq2", quietly = TRUE)
  ok_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
})

# =============================
# Configuration: File paths
# =============================
# Input: gene-by-sample count matrix (tab-separated)
expr_file <- "/work/ugbogu/pipeline/results/unsupervised/nbl/deseq2_inputs/counts.tsv"

# Input: NBL neuronal signature gene list (one Ensembl gene ID per line)
sig_file  <- "/work/ugbogu/pipeline/results/unsupervised/nbl/deseq2_markers/markers/nbl_features.txt"

# Cell-cycle gene list (one Ensembl gene ID per line) for proliferation proxy
cell_cycle_file <- "/work/ugbogu/pipeline/results/unsupervised/nbl/gsea_ova/cell_cycle_ensg.txt"

# Output directory and file prefix (all results written under out_dir)
out_dir <- "/work/ugbogu/pipeline/results/unsupervised/nbl/gsea_ova/nbl_signature_axes"
out_prefix <- "nbl_signature_axes"

# =============================
# Set up logging
# =============================
# Create output directory if it doesn't exist
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Define log file path
log_file <- file.path(out_dir, paste0(out_prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

# Open log file connection
log_con <- file(log_file, open = "wt")

# Custom logging function that writes to both console and file
log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  # Write to console
  message(msg)
  # Write to log file with timestamp
  writeLines(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg), con = log_con)
  flush(log_con)  # Ensure immediate write to disk
}

# Wrapper for warnings that also logs
log_warning <- function(...) {
  msg <- paste0("WARNING: ", ..., collapse = "")
  warning(msg, call. = FALSE)
  writeLines(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg), con = log_con)
  flush(log_con)
}

# Ensure log file is closed on exit (even if script errors out)
on.exit({
  if (exists("log_con") && isOpen(log_con)) {
    log_msg("Script completed.")
    close(log_con)
    message("Log saved to: ", log_file)
  }
}, add = TRUE)

# Log script start
log_msg("========================================")
log_msg("NBL Signature Axes Analysis")
log_msg("========================================")
log_msg("Script started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_msg("Log file: ", log_file)
log_msg("")

# =============================
# Helper functions
# =============================

# Validate that a file exists, otherwise stop with error
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    log_msg("ERROR: File not found: ", path)
    stop("File not found: ", path, call. = FALSE)
  }
}

# Read expression matrix from TSV file
# Handles two formats: (1) gene_id column present, (2) genes already as rownames
read_expr_tsv <- function(path) {
  x <- read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!("gene_id" %in% colnames(x))) {
    # Format 2: first column is already gene identifiers
    rownames(x) <- x[[1]]
    x[[1]] <- NULL
    return(as.matrix(x))
  }
  
  # Format 1: explicit gene_id column
  rownames(x) <- x$gene_id
  x$gene_id <- NULL
  as.matrix(x)
}

# Compute row-wise z-scores (center and scale each gene across samples)
zscore_rows <- function(mat) {
  # Transpose to make genes into columns, scale, then transpose back
  t(scale(t(mat)))
}

# Compute gene set score as mean z-score across a list of genes.
# Always returns a named vector of length ncol(expr_mat) to avoid downstream row mismatches.
score_gene_set <- function(expr_mat, genes) {
  genes <- intersect(rownames(expr_mat), genes)

  # Always return a vector aligned to samples (named, length = ncol)
  score <- rep(NA_real_, ncol(expr_mat))
  names(score) <- colnames(expr_mat)

  if (length(genes) < 3) {
    return(list(score = score, genes_present = genes))
  }

  ez <- t(scale(t(expr_mat[genes, , drop = FALSE])))
  score <- colMeans(ez, na.rm = TRUE)
  names(score) <- colnames(expr_mat)

  list(score = score, genes_present = genes)
}

# =============================
# Load input data
# =============================
log_msg("Checking input files...")
stop_if_missing(expr_file)
stop_if_missing(sig_file)
log_msg("All input files found.")
log_msg("")

log_msg("Reading expression matrix: ", expr_file)
counts <- read_expr_tsv(expr_file)

# Ensure all values are numeric (coerce if needed)
mode(counts) <- "numeric"

# Extract sample identifiers
samples <- colnames(counts)
log_msg("Expression matrix dimensions: ", nrow(counts), " genes x ", length(samples), " samples")
if (length(samples) < 3) {
  log_msg("ERROR: Expected >=3 samples, found ", length(samples))
  stop("Expected >=3 samples, found ", length(samples), call. = FALSE)
}
log_msg("")

# Load NBL signature gene list (one gene ID per line)
log_msg("Reading NBL signature gene list: ", sig_file)
sig <- unique(trimws(readLines(sig_file)))
sig <- sig[nzchar(sig)]  # Remove empty lines
log_msg("NBL signature genes loaded: ", length(sig))
log_msg("")

# =============================
# Transform count data (VST only)
# =============================
# For correlation analysis we use DESeq2's variance-stabilizing transformation.
if (!ok_deseq2) {
  log_msg("ERROR: DESeq2 is required but not found.")
  stop("DESeq2 is required. Install with: BiocManager::install('DESeq2')", call. = FALSE)
}

log_msg("Applying DESeq2 VST transform...")
suppressPackageStartupMessages(library(DESeq2))

# Create minimal sample metadata (all samples labeled as 'nbl')
coldata <- data.frame(row.names = samples, condition = factor(rep("nbl", length(samples))))

# Build DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = coldata,
                              design = ~ 1)  # No actual design (intercept-only)

# Filter out very low-count genes for VST stability
keep <- rowSums(counts(dds)) >= 10
log_msg("Filtering low-count genes: ", sum(keep), " of ", length(keep), " genes retained (rowSum >= 10)")
dds <- dds[keep, ]

# Apply variance-stabilizing transformation (blind=TRUE: ignore experimental design)
expr <- assay(vst(dds, blind = TRUE))
transform_used <- "DESeq2::vst (blind=TRUE; genes with rowSum>=10 kept)"
log_msg("VST transformation complete.")
log_msg("Transform method: ", transform_used)
log_msg("")

# =============================
# Compute biological axis scores
# =============================
log_msg("Computing biological axis scores...")
log_msg("")

# AXIS 1: NBL neuronal differentiation signature
log_msg("AXIS 1: NBL neuronal differentiation signature")
# Handle ENSG version suffixes: match signature list to expression rownames
expr_genes_nover <- sub("\\.[0-9]+$", "", rownames(expr))
sig_nover <- sub("\\.[0-9]+$", "", sig)
sig_present_rows <- rownames(expr)[expr_genes_nover %in% sig_nover]

log_msg("  Signature genes total: ", length(sig_nover))
log_msg("  Signature genes matched in expr: ", length(sig_present_rows))

sig_res <- score_gene_set(expr, sig_present_rows)
sig_score <- sig_res$score           # Per-sample scores
sig_present <- sig_res$genes_present # Genes from signature found in data
log_msg("  Signature score computed for ", length(sig_score), " samples")
log_msg("")

# AXIS 2: MYCN oncogene expression (version-safe ENSG matching)
log_msg("AXIS 2: MYCN oncogene expression")
mycn_gene <- "ENSG00000134323"  # unversioned MYCN ENSG

expr_genes_nover <- sub("\\.[0-9]+$", "", rownames(expr))
mycn_row <- rownames(expr)[expr_genes_nover == mycn_gene]

mycn_present <- length(mycn_row) == 1
mycn_expr <- if (mycn_present) as.numeric(expr[mycn_row, ]) else rep(NA_real_, ncol(expr))
names(mycn_expr) <- colnames(expr)

if (mycn_present) {
  log_msg("  MYCN gene found: ", mycn_row)
  log_msg("  MYCN expression range: ", round(min(mycn_expr), 2), " to ", round(max(mycn_expr), 2))
} else {
  log_warning("MYCN gene (", mycn_gene, ") not found in expression matrix")
}
log_msg("")

# AXIS 3: Cell-cycle/proliferation score (ENSG IDs from file)
log_msg("AXIS 3: Cell-cycle/proliferation score")
stop_if_missing(cell_cycle_file)
cell_cycle_ensg <- unique(trimws(readLines(cell_cycle_file, warn = FALSE)))
cell_cycle_ensg <- cell_cycle_ensg[nzchar(cell_cycle_ensg)]
log_msg("  Cell-cycle genes loaded from file: ", length(cell_cycle_ensg))

# Match to expression matrix: rownames may have version suffix (e.g. ENSG00000012345.2)
expr_genes_nover <- sub("\\.[0-9]+$", "", rownames(expr))
cc_genes_present <- rownames(expr)[expr_genes_nover %in% cell_cycle_ensg]
log_msg("  Cell-cycle genes matched in expr: ", length(cc_genes_present))

if (length(cc_genes_present) < 3) {
  log_warning("Cell-cycle proxy not computed: fewer than 3 cell-cycle genes found in expression matrix. Check that ", cell_cycle_file, " contains ENSG IDs.")
}
cc_res <- score_gene_set(expr, cc_genes_present)
cc_score <- cc_res$score
log_msg("  Cell-cycle score computed for ", sum(!is.na(cc_score)), " samples")
log_msg("")

# =============================
# Statistical tests
# =============================
log_msg("Running statistical tests...")
log_msg("")

# Helper: compute Spearman correlation with safe error handling
spearman_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# TEST 1: Correlation between NBL signature and MYCN expression
log_msg("TEST 1: NBL signature vs MYCN expression (Spearman correlation)")
stat_sig_mycn <- spearman_safe(sig_score, mycn_expr)
log_msg("  rho = ", round(stat_sig_mycn$rho, 3), ", p = ", format(stat_sig_mycn$p, digits = 3), ", n = ", stat_sig_mycn$n)
log_msg("")

# TEST 2: Correlation between NBL signature and cell-cycle score
log_msg("TEST 2: NBL signature vs cell-cycle score (Spearman correlation)")
stat_sig_cc   <- spearman_safe(sig_score, cc_score)
log_msg("  rho = ", round(stat_sig_cc$rho, 3), ", p = ", format(stat_sig_cc$p, digits = 3), ", n = ", stat_sig_cc$n)
log_msg("")

# TEST 3: Wilcoxon test comparing signature scores in MYCN-high vs MYCN-low samples
log_msg("TEST 3: NBL signature in MYCN-high vs MYCN-low (Wilcoxon rank-sum test)")
# Split samples into quartiles based on MYCN expression
mycn_group <- rep(NA_character_, length(mycn_expr))
names(mycn_group) <- names(mycn_expr)
wilcox_p <- NA_real_

if (mycn_present) {
  # Define high = top quartile, low = bottom quartile
  q <- quantile(mycn_expr, probs = c(0.25, 0.75), na.rm = TRUE)
  mycn_group[mycn_expr >= q[2]] <- "MYCN_high"
  mycn_group[mycn_expr <= q[1]] <- "MYCN_low"
  
  log_msg("  MYCN quartile thresholds: low <= ", round(q[1], 2), ", high >= ", round(q[2], 2))
  log_msg("  MYCN_low samples: ", sum(mycn_group == "MYCN_low", na.rm = TRUE))
  log_msg("  MYCN_high samples: ", sum(mycn_group == "MYCN_high", na.rm = TRUE))
  
  # Prepare data for Wilcoxon test
  df2 <- data.frame(sig = sig_score, group = mycn_group)
  df2 <- df2[!is.na(df2$group) & is.finite(df2$sig), , drop = FALSE]
  
  # Need at least 3 samples per group (6 total) for meaningful test
  if (nrow(df2) >= 6) {
    wilcox_p <- suppressWarnings(wilcox.test(sig ~ group, data = df2)$p.value)
    log_msg("  Wilcoxon p = ", format(wilcox_p, digits = 3))
  } else {
    log_warning("Insufficient samples for Wilcoxon test (need >= 6, found ", nrow(df2), ")")
  }
} else {
  log_msg("  Skipped (MYCN not present)")
}
log_msg("")

# =============================
# Write output files
# =============================
log_msg("Writing output files...")

# OUTPUT 1: Per-sample scores table
results <- data.frame(
  sample = names(sig_score),
  nbl_signature = as.numeric(sig_score),       # NBL differentiation score
  MYCN = as.numeric(mycn_expr),                # MYCN expression level
  MYCN_group = mycn_group,                     # MYCN_high / MYCN_low / NA
  cell_cycle_proxy = as.numeric(cc_score),     # Proliferation score
  stringsAsFactors = FALSE
)

out_tsv <- file.path(out_dir, paste0(out_prefix, "_scores.tsv"))
write.table(results, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("  Wrote scores table: ", out_tsv)

# OUTPUT 2: Summary metadata and statistics
meta <- data.frame(
  key = c("expr_file","sig_file","transform_used",
          "sig_genes_total","sig_genes_present",
          "mycn_present",
          "spearman_sig_vs_mycn_rho","spearman_sig_vs_mycn_p","spearman_sig_vs_mycn_n",
          "spearman_sig_vs_cellcycle_rho","spearman_sig_vs_cellcycle_p","spearman_sig_vs_cellcycle_n",
          "wilcox_sig_MYCN_high_vs_low_p"),
  value = c(expr_file, sig_file, transform_used,
            length(sig), length(sig_present),
            mycn_present,
            stat_sig_mycn$rho, stat_sig_mycn$p, stat_sig_mycn$n,
            stat_sig_cc$rho, stat_sig_cc$p, stat_sig_cc$n,
            wilcox_p),
  stringsAsFactors = FALSE
)

out_meta <- file.path(out_dir, paste0(out_prefix, "_summary.tsv"))
write.table(meta, file = out_meta, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("  Wrote summary table: ", out_meta)
log_msg("")

# =============================
# Generate diagnostic plots
# =============================
if (ok_ggplot) {
  log_msg("Generating diagnostic plots...")
  suppressPackageStartupMessages(library(ggplot2))

  out_pdf <- file.path(out_dir, paste0(out_prefix, "_plots.pdf"))
  pdf(out_pdf, width = 7, height = 6)

  # PLOT 1: Scatter plot of NBL signature vs MYCN expression
  dfp <- results
  if (mycn_present) {
    p1 <- ggplot(dfp, aes(x = MYCN, y = nbl_signature, label = sample)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +  # Linear regression line
      theme_minimal() +
      labs(
        title = "NBL signature vs MYCN expression",
        subtitle = sprintf("Spearman rho=%.3f, p=%.3g (n=%d)", 
                          stat_sig_mycn$rho, stat_sig_mycn$p, stat_sig_mycn$n),
        x = "MYCN (transformed expression)",
        y = "NBL signature score (mean z across genes)"
      )
    print(p1)
    log_msg("  Plot 1: NBL signature vs MYCN (scatter)")

    # PLOT 2: Boxplot comparing signature in MYCN-high vs MYCN-low samples
    dfb <- dfp[!is.na(dfp$MYCN_group), ]
    if (nrow(dfb) > 0) {
      p2 <- ggplot(dfb, aes(x = MYCN_group, y = nbl_signature)) +
        geom_boxplot() +
        geom_point() +  # Overlay individual sample points
        theme_minimal() +
        labs(
          title = "Signature by MYCN high/low (quartiles)",
          subtitle = if (is.finite(wilcox_p)) sprintf("Wilcoxon p=%.3g", wilcox_p) else "Wilcoxon p=NA",
          x = NULL, y = "NBL signature score"
        )
      print(p2)
      log_msg("  Plot 2: NBL signature by MYCN group (boxplot)")
    }
  } else {
    # MYCN gene not found in expression matrix
    plot.new()
    title("MYCN not found in expression matrix (ENSG00000134323)")
    log_msg("  Plot 1-2: Skipped (MYCN not found)")
  }

  # PLOT 3: Scatter plot of NBL signature vs cell-cycle proxy
  if (all(is.na(results$cell_cycle_proxy))) {
    # No cell-cycle genes were found
    plot.new()
    title("Cell-cycle proxy not computed (no matching genes found).")
    mtext("Provide CELL_CYCLE_genes.txt with ENSG IDs to enable this test.", side = 3)
    log_msg("  Plot 3: Skipped (cell-cycle genes not found)")
  } else {
    p3 <- ggplot(dfp, aes(x = cell_cycle_proxy, y = nbl_signature, label = sample)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +  # Linear regression line
      theme_minimal() +
      labs(
        title = "NBL signature vs cell-cycle proxy",
        subtitle = sprintf("Spearman rho=%.3f, p=%.3g (n=%d)", 
                          stat_sig_cc$rho, stat_sig_cc$p, stat_sig_cc$n),
        x = "Cell-cycle proxy score",
        y = "NBL signature score"
      )
    print(p3)
    log_msg("  Plot 3: NBL signature vs cell-cycle (scatter)")
  }

  dev.off()
  log_msg("  Wrote plots PDF: ", out_pdf)
} else {
  log_msg("ggplot2 not available -> skipping plots.")
}
log_msg("")

log_msg("========================================")
log_msg("Analysis complete!")
log_msg("========================================")