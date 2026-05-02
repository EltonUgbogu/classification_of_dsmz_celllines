#!/usr/bin/env Rscript

# build_dsmz_tcga_pam50matrix.R
# Build integrated DSMZ+TCGA PAM50 expression matrix + mapping
# for tumour neighbourhood analysis (BRCA-only pipeline).

options(stringsAsFactors = FALSE)
set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(optparse)
  library(yaml)
})

# Source config helper (handles profile merging)
# Find lib_config.R relative to this script's location
# When run via Rscript, the script path is in commandArgs
script_path <- commandArgs(trailingOnly = FALSE)[4]
if (is.na(script_path) || script_path == "") {
  # Fallback: assume we're in the pipeline root
  lib_config_path <- "scripts/lib_config.R"
} else {
  script_dir <- dirname(normalizePath(script_path))
  lib_config_path <- file.path(script_dir, "lib_config.R")
  if (!file.exists(lib_config_path)) {
    # Fallback: try from current working directory
    lib_config_path <- "scripts/lib_config.R"
  }
}
if (!file.exists(lib_config_path)) {
  stop("Cannot find lib_config.R. Tried: ", 
       ifelse(exists("script_dir"), file.path(script_dir, "lib_config.R"), "N/A"), 
       " and scripts/lib_config.R")
}
source(lib_config_path)

# --------------------------------------------------------------------
# 0. CLI + config
# --------------------------------------------------------------------
option_list <- list(
  make_option(
    "--config",
    type = "character",
    default = "config/config.yaml",
    help = "Path to config.yaml [default: %default]"
  ),
  make_option(
    "--profile",
    type = "character",
    default = NULL,
    help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Using config file:", opt$config, "\n")
# Priority: opt$profile -> SNAKEMAKE_PROFILE env var -> config default
# Note: %||% is defined in lib_config.R, so we can use it after sourcing
profile_override <- if (!is.null(opt$profile)) opt$profile else Sys.getenv("SNAKEMAKE_PROFILE", unset = "")
if (identical(profile_override, "")) profile_override <- NULL
cfg <- read_profiled_config(opt$config, profile_override = profile_override)

unsup_root  <- cfg$paths$unsup_root
input_root  <- cfg$paths$tumour_nh_input_root

# Get PAM50 file paths with validation
tcga_file <- cfg$paths$tcga_brca_pam50_expr %||% cfg$paths$tcga_pam50_expr
dsmz_file <- cfg$paths$dsmz_bcc_pam50_expr %||% cfg$paths$dsmz_pam50_expr

stopifnot(
  "tcga_brca_pam50_expr or tcga_pam50_expr not found in config" = !is.null(tcga_file),
  "dsmz_bcc_pam50_expr or dsmz_pam50_expr not found in config" = !is.null(dsmz_file),
  "tcga file path must be a single string" = is.character(tcga_file) && length(tcga_file) == 1,
  "dsmz file path must be a single string" = is.character(dsmz_file) && length(dsmz_file) == 1,
  "tcga file does not exist" = file.exists(tcga_file),
  "dsmz file does not exist" = file.exists(dsmz_file)
)

tcga_mat <- readRDS(tcga_file)
dsmz_mat <- readRDS(dsmz_file)

# align genes
common_genes <- intersect(rownames(tcga_mat), rownames(dsmz_mat))
stopifnot(length(common_genes) >= 40)
tcga_mat <- tcga_mat[common_genes, , drop = FALSE]
dsmz_mat <- dsmz_mat[common_genes, , drop = FALSE]

# DSMZ metadata
meta <- read.csv(cfg$paths$dsmz_meta_csv, stringsAsFactors = FALSE)
stopifnot("sample_name" %in% colnames(meta))
cell_col <- intersect(c("DSMZ_Cell_line_norm", "Cell_line", "Cell_Line", "CellLine"), colnames(meta))[1]
stopifnot(!is.na(cell_col))

cell_line_of <- setNames(meta[[cell_col]], meta$sample_name)

dsmz_original_ids <- colnames(dsmz_mat)
tcga_ids          <- colnames(tcga_mat)

expr_mat <- rbind(
  t(dsmz_mat),  # rownames = original DSMZ sample IDs
  t(tcga_mat)   # rownames = TCGA barcodes
)

dataset_vec <- ifelse(rownames(expr_mat) %in% dsmz_original_ids, "DSMZ", "TCGA")
names(dataset_vec) <- rownames(expr_mat)

# Get output paths from config (profile-merged)
expr_out <- cfg$paths$tumour_nh_expr_pam50
map_out  <- cfg$paths$tumour_nh_map_rds_pam50

stopifnot(
  "tumour_nh_expr_pam50 not found in config" = !is.null(expr_out),
  "tumour_nh_map_rds_pam50 not found in config" = !is.null(map_out),
  "expr_out must be a single string" = is.character(expr_out) && length(expr_out) == 1,
  "map_out must be a single string" = is.character(map_out) && length(map_out) == 1
)

# Create output directories
dir.create(dirname(expr_out), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(map_out), showWarnings = FALSE, recursive = TRUE)

cell_line_to_sample_id <- cell_line_of[dsmz_original_ids]
cell_line_to_sample_id[is.na(cell_line_to_sample_id)] <- dsmz_original_ids[is.na(cell_line_to_sample_id)]

saveRDS(expr_mat, expr_out)
saveRDS(cell_line_to_sample_id, map_out)

cat(sprintf("Matrix: %d samples × %d genes (%d DSMZ | %d TCGA)\n",
            nrow(expr_mat), ncol(expr_mat),
            sum(dataset_vec == "DSMZ"), sum(dataset_vec == "TCGA")))
cat("Expression matrix saved to:", expr_out, "\n")
cat("Mapping saved to:", map_out, "\n")

