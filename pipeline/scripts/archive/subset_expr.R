#!/usr/bin/env Rscript
# =============================================================================
# subset_expr.R
# =============================================================================
#
# PURPOSE:
# Subsets a VST-normalized expression matrix to include only the genes
# specified in a gene list file. This is used in the feature selection
# pipeline to create method-specific expression matrices.
#
# INPUTS:
#   --vst_rds  : Path to VST-normalized expression matrix (.rds file)
#   --genes    : Path to text file with one gene ID per line
#   --out_rds  : Path to output subsetted expression matrix (.rds file)
#
# OUTPUTS:
#   - Subsetted expression matrix (genes × samples) saved as .rds
#
# USAGE:
#   Rscript subset_expr.R \
#     --vst_rds results/vst_joint.rds \
#     --genes results/featuresets/Variance/genes_top500.txt \
#     --out_rds results/featuresets/Variance/expr_submatrix.rds
#
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# -----------------------------------------------------------------------------
# COMMAND-LINE INTERFACE
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(
    c("--vst_rds"),
    type = "character",
    help = "Path to VST-normalized expression matrix (.rds file)"
  ),
  make_option(
    c("--genes"),
    type = "character",
    help = "Path to text file with one gene ID per line"
  ),
  make_option(
    c("--out_rds"),
    type = "character",
    help = "Path to output subsetted expression matrix (.rds file)"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$vst_rds) || is.null(opt$genes) || is.null(opt$out_rds)) {
  print_help(opt_parser)
  stop("All arguments (--vst_rds, --genes, --out_rds) are required")
}

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

message("[INFO] Loading VST matrix from: ", opt$vst_rds)
vst_matrix <- readRDS(opt$vst_rds)

if (!is.matrix(vst_matrix)) {
  stop("VST object must be a matrix")
}

if (is.null(rownames(vst_matrix))) {
  stop("Matrix must have gene rownames")
}

n_genes_orig <- nrow(vst_matrix)
n_samples <- ncol(vst_matrix)
message("[INFO] Original matrix: ", n_genes_orig, " genes × ", n_samples, " samples")

# -----------------------------------------------------------------------------
# LOAD GENE LIST
# -----------------------------------------------------------------------------

message("[INFO] Loading gene list from: ", opt$genes)
gene_list <- readLines(opt$genes)
gene_list <- gene_list[gene_list != ""]  # Remove empty lines
gene_list <- trimws(gene_list)  # Remove whitespace
gene_list <- unique(gene_list)  # Remove duplicates

message("[INFO] Gene list contains ", length(gene_list), " unique genes")

# -----------------------------------------------------------------------------
# SUBSET EXPRESSION MATRIX
# -----------------------------------------------------------------------------

# Find which genes from the list are present in the matrix
genes_found <- intersect(gene_list, rownames(vst_matrix))
genes_missing <- setdiff(gene_list, rownames(vst_matrix))

if (length(genes_found) == 0) {
  stop("No genes from the gene list were found in the VST matrix rownames")
}

if (length(genes_missing) > 0) {
  message("[WARNING] ", length(genes_missing), " genes from list not found in matrix:")
  if (length(genes_missing) <= 10) {
    message("[WARNING]   ", paste(genes_missing, collapse = ", "))
  } else {
    message("[WARNING]   ", paste(genes_missing[1:10], collapse = ", "), " ... (and ", 
            length(genes_missing) - 10, " more)")
  }
}

message("[INFO] Subsetting to ", length(genes_found), " genes")

# Subset the matrix (preserve sample order)
expr_subset <- vst_matrix[genes_found, , drop = FALSE]

message("[INFO] Subsetted matrix: ", nrow(expr_subset), " genes × ", ncol(expr_subset), " samples")

# -----------------------------------------------------------------------------
# SAVE OUTPUT
# -----------------------------------------------------------------------------

message("[INFO] Saving subsetted matrix to: ", opt$out_rds)

# Create output directory if it doesn't exist
out_dir <- dirname(opt$out_rds)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(expr_subset, opt$out_rds)

message("[INFO] Done! Subsetted expression matrix saved.")

