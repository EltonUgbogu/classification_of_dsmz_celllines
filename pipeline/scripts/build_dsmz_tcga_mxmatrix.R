#!/usr/bin/env Rscript

# build_dsmz_tcga_mxmatrix.R
# Build joint DSMZ+TCGA MX500 expression matrix + mapping for tumour neighbourhoods.
# Assumes vst_joint_rds is a true joint matrix with original sample IDs.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(yaml)
  library(optparse)
})

source("scripts/lib_config.R")

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]"),
  make_option("--profile", type = "character", default = NULL,
              help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("=== Building joint DSMZ+TCGA MX500 matrix ===\n")

# Load config using profiled config loader (merges defaults + profile)
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)

vst_path   <- cfg$paths$vst_joint_rds
mx_path    <- cfg$features$mx_final_gene_list
dsmz_meta_path <- cfg$paths$dsmz_meta_csv
out_path   <- cfg$paths$tumour_nh_expr_mx
input_root <- cfg$paths$tumour_nh_input_root
map_out    <- file.path(input_root, "cell_line_to_original_sample_id_MX.rds")

cat("[INFO] Using VST matrix: ", vst_path, "\n", sep = "")
cat("[INFO] Using MX list   : ", mx_path, "\n", sep = "")
cat("[INFO] Using DSMZ meta : ", dsmz_meta_path, "\n", sep = "")
cat("[INFO] Output matrix   : ", out_path, "\n", sep = "")
cat("[INFO] Output mapping  : ", map_out, "\n", sep = "")

stopifnot(file.exists(vst_path))
stopifnot(file.exists(mx_path))
stopifnot(file.exists(dsmz_meta_path))

vst_joint <- readRDS(vst_path)
cat("[INFO] VST dims:", nrow(vst_joint), "genes x", ncol(vst_joint), "samples\n")

# DSMZ/TCGA sanity check (informative only)
n_dsmz <- sum(grepl("^NG-", colnames(vst_joint)))
n_tcga <- sum(grepl("^TCGA-", colnames(vst_joint)))
cat("[INFO] Detected DSMZ samples (NG- prefix):", n_dsmz, "\n")
cat("[INFO] Detected TCGA samples (TCGA- prefix):", n_tcga, "\n")

if (n_dsmz == 0) {
  stop("No DSMZ samples found in VST matrix. ",
       "Please run scripts/rebuild_joint_vst.R first to create a proper joint matrix.")
}

# Read MX gene list
mx_genes <- scan(mx_path, what = character(), quiet = TRUE)
cat("[INFO] MX list length:", length(mx_genes), "\n")

mx_genes <- intersect(mx_genes, rownames(vst_joint))
cat("[INFO] MX genes present in VST:", length(mx_genes), "\n")
if (length(mx_genes) < 100) {
  stop("Too few MX genes found in VST matrix: ", length(mx_genes))
}

# Build samples × genes matrix
expr_mx <- t(vst_joint[mx_genes, , drop = FALSE])

# Extract DSMZ sample IDs (columns in original VST, now rownames in expr_mx)
dsmz_sample_ids <- rownames(expr_mx)[grepl("^NG-", rownames(expr_mx))]

# Read DSMZ metadata and build mapping (same logic as PAM50 script)
meta <- read.csv(dsmz_meta_path, stringsAsFactors = FALSE)
stopifnot("sample_name" %in% colnames(meta))
cell_col <- intersect(c("Cell_line", "Cell_Line", "CellLine"), colnames(meta))[1]
stopifnot(!is.na(cell_col))

cell_line_of <- setNames(meta[[cell_col]], meta$sample_name)

# Map DSMZ sample IDs to cell lines
cell_line_to_sample_id <- cell_line_of[dsmz_sample_ids]
# For samples not found in metadata, use the sample ID itself
cell_line_to_sample_id[is.na(cell_line_to_sample_id)] <- dsmz_sample_ids[is.na(cell_line_to_sample_id)]

# Create output directories
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(map_out), recursive = TRUE, showWarnings = FALSE)

# Save outputs
saveRDS(expr_mx, out_path)
saveRDS(cell_line_to_sample_id, map_out)

cat("\n[OK] Saved MX matrix:\n")
cat("  ", out_path, "\n")
cat("[OK] Saved mapping:\n")
cat("  ", map_out, "\n")
cat("[INFO] Final dims:", nrow(expr_mx), "samples x", ncol(expr_mx), "genes\n")
cat("[INFO] DSMZ:", sum(grepl("^NG-", rownames(expr_mx))),
    "| TCGA:", sum(grepl("^TCGA-", rownames(expr_mx))), "\n")
cat("[INFO] Mapping entries:", length(cell_line_to_sample_id), "DSMZ samples\n")
