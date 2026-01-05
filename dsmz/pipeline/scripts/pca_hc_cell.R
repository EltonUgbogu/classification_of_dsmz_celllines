#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--cell_rds",    type = "character"),
  make_option("--outdir",      type = "character"),
  make_option("--base_dir",    type = "character", default = "R"),
  make_option("--label",       type = "character", default = "CELL"),
  make_option("--n_pcs",       type = "integer",   default = 20L),
  make_option("--max_k",       type = "integer",   default = 8L),
  make_option("--dist_method", type = "character", default = "euclidean"),
  make_option("--seed",        type = "integer",   default = 42L),
  make_option("--cluster_rds", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$cell_rds))    stop("--cell_rds is required.", call. = FALSE)
if (is.null(opt$outdir))      stop("--outdir is required.", call. = FALSE)
if (is.null(opt$cluster_rds)) stop("--cluster_rds is required.", call. = FALSE)

utils_path <- file.path(opt$base_dir, "agnostic_clustering_utils.R")
if (!file.exists(utils_path)) stop("Cannot find agnostic_clustering_utils.R at: ", utils_path)
source(utils_path)

run_agnostic_clustering(
  kind             = "pca_hc_cell",
  cell_rds         = opt$cell_rds,
  tumour_rds       = NULL,
  outdir           = opt$outdir,
  n_pcs            = opt$n_pcs,
  max_k_hc         = opt$max_k,
  k_min            = 2L,
  k_max            = 8L,
  dist_method      = opt$dist_method,
  seed             = opt$seed,
  cluster_rds_path = opt$cluster_rds,
  cell_label       = opt$label,
  tumour_label     = "TUMOUR"
)