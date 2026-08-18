#!/usr/bin/env Rscript

# pca_hc_cell.R
# PCA + hierarchical clustering on CELL-only VST matrix
# Wrapper around run_hclust_kmeans(kind = "pca_hc_cell")

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(
    "--cell_rds",
    type = "character",
    help = "Path to CELL VST RDS (matrix or SummarizedExperiment-like).",
    metavar = "file"
  ),
  make_option(
    "--outdir",
    type = "character",
    help = "Output directory for clustering results.",
    metavar = "dir"
  ),
  make_option(
    "--base_dir",
    type = "character",
    default = "R",
    help = "Directory containing hclust_kmeans_utils.R [default: %default].",
    metavar = "dir"
  ),
  make_option(
    "--label",
    type = "character",
    default = "CELL",
    help = "Dataset label for CELL samples [default: %default].",
    metavar = "label"
  ),
  make_option(
    "--n_pcs",
    type = "integer",
    default = 20L,
    help = "Number of principal components to use [default: %default].",
    metavar = "integer"
  ),
  make_option(
    "--max_k",
    type = "integer",
    default = 8L,
    help = "Maximum K for hierarchical clustering [default: %default].",
    metavar = "integer"
  ),
  make_option(
    "--dist_method",
    type = "character",
    default = "euclidean",
    help = "Distance method for HC ('euclidean' or 'correlation') [default: %default].",
    metavar = "method"
  ),
  make_option(
    "--seed",
    type = "integer",
    default = 42L,
    help = "Random seed for clustering [default: %default].",
    metavar = "integer"
  ),
  make_option(
    "--cluster_rds",
    type = "character",
    help = "Path to output RDS file to save cluster object.",
    metavar = "file"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# ─────────────────────────────────────────────────────────────
# Basic argument checks
# ─────────────────────────────────────────────────────────────
if (is.null(opt$cell_rds)) {
  stop("--cell_rds is required.", call. = FALSE)
}
if (is.null(opt$outdir)) {
  stop("--outdir is required.", call. = FALSE)
}
if (is.null(opt$cluster_rds)) {
  stop("--cluster_rds is required.", call. = FALSE)
}

# ─────────────────────────────────────────────────────────────
# Load shared HC/k-means clustering utilities
# ─────────────────────────────────────────────────────────────
utils_path <- file.path(opt$base_dir, "hclust_kmeans_utils.R")
if (!file.exists(utils_path)) {
  stop("Cannot find hclust_kmeans_utils.R at: ", utils_path)
}
source(utils_path)

# ─────────────────────────────────────────────────────────────
# Run HC/k-means clustering: pca_hc_cell
# ─────────────────────────────────────────────────────────────
run_hclust_kmeans(
  kind              = "pca_hc_cell",
  cell_rds          = opt$cell_rds,
  tumour_rds        = NULL,
  outdir            = opt$outdir,
  n_pcs             = opt$n_pcs,
  max_k_hc          = opt$max_k,
  k_min             = 2L,              # not used for HC, but kept for API consistency
  k_max             = 8L,              # not used for HC, but kept for API consistency
  dist_method       = opt$dist_method,
  seed              = opt$seed,
  cluster_rds_path  = opt$cluster_rds,
  cell_label        = opt$label,
  tumour_label      = "TUMOUR"
)
