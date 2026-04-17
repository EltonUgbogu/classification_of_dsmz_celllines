#!/usr/bin/env Rscript
# =============================================================================
# compute_similarity_graph.R
# =============================================================================
#
# PURPOSE:
# Standalone helper to compute sample-sample similarity matrices and kNN
# graphs from an expression matrix (genes × samples). The main Snakefile uses
# the legacy compute_cell_line_similarity.R; this script is optional and not
# invoked automatically. Use it if you want to prototype alternative
# similarity graphs (e.g., on a specific feature-set expression matrix) and
# then feed those outputs into downstream exploration/visualisation steps.
#
# SAFE DEFAULTS FOR THIS PIPELINE:
# * Works for any profile (BRCA/RBL/NBL) as long as --expr_rds points to a
#   gene × sample matrix produced by the pipeline (e.g., a feature-set
#   expression matrix under results/feature_selection_unsupervised or a
#   tumour_neighbourhood input matrix).
# * The default Spearman correlation is robust for expression ranks; set
#   --similarity to "pearson" if you want linear correlation on log-scale
#   data.
# * Choose --k smaller than your sample count (k < n) to avoid kNN errors; a
#   common heuristic is k ≈ min(30, n/3).
#
# INPUTS:
#   --expr_rds    : Path to expression matrix (.rds file, genes × samples)
#   --similarity  : Similarity metric ("spearman" or "pearson", default: "spearman")
#   --k           : Number of nearest neighbors for kNN graph (default: 30)
#   --sim_out     : Path to output similarity matrix (.rds file)
#   --edges_out   : Path to output kNN edge list (.tsv file)
#
# OUTPUTS:
#   - Similarity matrix (samples × samples) saved as .rds
#   - kNN edge list (from, to, weight) saved as .tsv
#
# USAGE:
#   Rscript compute_similarity_graph.R \
#     --expr_rds results/featuresets/Variance/expr_submatrix.rds \
#     --similarity spearman \
#     --k 30 \
#     --sim_out results/featuresets/Variance/similarity_matrix.rds \
#     --edges_out results/featuresets/Variance/knn_edges.tsv
#
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
})

# -----------------------------------------------------------------------------
# COMMAND-LINE INTERFACE
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(
    c("--expr_rds"),
    type = "character",
    help = "Path to expression matrix (.rds file, genes × samples)"
  ),
  make_option(
    c("--similarity"),
    type = "character",
    default = "spearman",
    help = "Similarity metric: 'spearman' or 'pearson' (default: spearman)"
  ),
  make_option(
    c("--k"),
    type = "integer",
    default = 30,
    help = "Number of nearest neighbors for kNN graph (default: 30)"
  ),
  make_option(
    c("--sim_out"),
    type = "character",
    help = "Path to output similarity matrix (.rds file)"
  ),
  make_option(
    c("--edges_out"),
    type = "character",
    help = "Path to output kNN edge list (.tsv file)"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$expr_rds) || is.null(opt$sim_out) || is.null(opt$edges_out)) {
  print_help(opt_parser)
  stop("All arguments (--expr_rds, --sim_out, --edges_out) are required")
}

# Validate similarity method
if (!opt$similarity %in% c("spearman", "pearson")) {
  stop("--similarity must be 'spearman' or 'pearson'")
}

# Validate k
if (opt$k < 1) {
  stop("--k must be >= 1")
}

# -----------------------------------------------------------------------------
# LOAD EXPRESSION MATRIX
# -----------------------------------------------------------------------------

cat("=== Computing similarity graph from expression matrix ===\n")
cat("Expression matrix:\n  ", opt$expr_rds, "\n")
cat("Similarity method: ", opt$similarity, "\n")
cat("kNN parameter k: ", opt$k, "\n\n")

expr_matrix <- readRDS(opt$expr_rds)

if (!is.matrix(expr_matrix)) {
  stop("Expression object must be a matrix")
}

if (is.null(rownames(expr_matrix)) || is.null(colnames(expr_matrix))) {
  stop("Matrix must have both gene rownames and sample colnames")
}

n_genes <- nrow(expr_matrix)
n_samples <- ncol(expr_matrix)
cat("Expression matrix dimensions (genes × samples):\n")
cat("  ", n_genes, " genes × ", n_samples, " samples\n\n")

if (n_samples < opt$k + 1) {
  stop("Number of samples (", n_samples, ") must be >= k+1 (", opt$k + 1, ")")
}

# -----------------------------------------------------------------------------
# COMPUTE SAMPLE-SAMPLE SIMILARITY MATRIX
# -----------------------------------------------------------------------------
# Similarity between two samples is defined as the correlation of their
# expression vectors across all genes.
#
# The cor() function computes pairwise correlations.
# t(expr_matrix) transposes so that cor() computes correlations between
# columns (samples) rather than rows (genes).

cat("Computing ", opt$similarity, " correlation between samples...\n")

# Transpose so samples are columns (for cor() function)
expr_t <- t(expr_matrix)

# Compute correlation matrix (samples × samples)
sim_mat <- cor(expr_t, method = opt$similarity, use = "pairwise.complete.obs")

# Handle any NA values
if (any(is.na(sim_mat))) {
  cat("[WARNING] Found NA values in similarity matrix, replacing with 0\n")
  sim_mat[is.na(sim_mat)] <- 0
}

# Ensure symmetric (should already be symmetric, but enforce it)
sim_mat <- (sim_mat + t(sim_mat)) / 2

# Set diagonal to 1 (self-similarity)
diag(sim_mat) <- 1

cat("Similarity matrix dimensions (samples × samples):\n")
cat("  ", nrow(sim_mat), " × ", ncol(sim_mat), "\n\n")

# -----------------------------------------------------------------------------
# EXPORT SIMILARITY MATRIX
# -----------------------------------------------------------------------------

# Create output directory if it doesn't exist
dir.create(dirname(opt$sim_out), showWarnings = FALSE, recursive = TRUE)

saveRDS(sim_mat, opt$sim_out)
cat("Saved similarity matrix to:\n  ", opt$sim_out, "\n")

# -----------------------------------------------------------------------------
# BUILD kNN GRAPH
# -----------------------------------------------------------------------------
# For each sample, find its k nearest neighbors (excluding self).
# We use the similarity scores directly (higher = more similar).

cat("\nBuilding kNN graph with k=", opt$k, "...\n")

edges_list <- list()
edge_count <- 0

for (i in seq_len(n_samples)) {
  sample_i <- colnames(expr_matrix)[i]
  
  # Get similarities to all other samples
  similarities <- sim_mat[i, ]
  
  # Exclude self (set to -Inf so it won't be selected)
  similarities[i] <- -Inf
  
  # Find top k neighbors
  top_k_indices <- order(similarities, decreasing = TRUE)[1:min(opt$k, length(similarities) - 1)]
  top_k_similarities <- similarities[top_k_indices]
  
  # Create edges
  for (j_idx in seq_along(top_k_indices)) {
    j <- top_k_indices[j_idx]
    sample_j <- colnames(expr_matrix)[j]
    weight <- top_k_similarities[j_idx]
    
    # Only add edge if weight is valid (not -Inf)
    if (is.finite(weight)) {
      edge_count <- edge_count + 1
      edges_list[[edge_count]] <- data.frame(
        from = sample_i,
        to = sample_j,
        weight = weight,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine all edges into a data frame
if (length(edges_list) > 0) {
  edges_df <- do.call(rbind, edges_list)
} else {
  stop("No edges created in kNN graph")
}

cat("Created ", nrow(edges_df), " edges in kNN graph\n")

# -----------------------------------------------------------------------------
# EXPORT kNN EDGES
# -----------------------------------------------------------------------------

# Create output directory if it doesn't exist
dir.create(dirname(opt$edges_out), showWarnings = FALSE, recursive = TRUE)

readr::write_tsv(edges_df, opt$edges_out)
cat("Saved kNN edges to:\n  ", opt$edges_out, "\n")

# -----------------------------------------------------------------------------
# SUMMARY STATISTICS
# -----------------------------------------------------------------------------

cat("\n=== Similarity graph summary ===\n")
cat("Expression matrix: ", n_genes, " genes × ", n_samples, " samples\n")
cat("Similarity metric: ", opt$similarity, "\n")
cat("kNN parameter k: ", opt$k, "\n")
cat("Similarity matrix: ", nrow(sim_mat), " × ", ncol(sim_mat), "\n")
cat("Number of edges: ", nrow(edges_df), "\n")
cat("Average degree: ", sprintf("%.2f", nrow(edges_df) / n_samples), "\n")
cat("Similarity range: [", sprintf("%.3f", min(sim_mat)), ", ", 
    sprintf("%.3f", max(sim_mat)), "]\n")
cat("Edge weight range: [", sprintf("%.3f", min(edges_df$weight)), ", ", 
    sprintf("%.3f", max(edges_df$weight)), "]\n")
cat("\n=== Done ===\n")

