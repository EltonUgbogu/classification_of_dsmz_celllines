#!/usr/bin/env Rscript

# build_dsmz_tcga_hvgmatrix.R
# Build joint DSMZ+TCGA HVG500 expression matrix + mapping for tumour neighbourhoods.
# Assumes vst_joint_rds is a true joint matrix with original sample IDs.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(yaml)
})

cat("=== Building joint DSMZ+TCGA HVG500 matrix ===\n")

cfg <- yaml::read_yaml("config/config.yaml")

vst_path   <- cfg$paths$vst_joint_rds
hvg_path   <- cfg$features$hvg_final_gene_list
dsmz_meta_path <- cfg$paths$dsmz_meta_csv
out_path   <- cfg$paths$tumour_nh_expr_hvg
input_root <- cfg$paths$tumour_nh_input_root
map_out    <- file.path(input_root, "cell_line_to_original_sample_id_hvg.rds")

cat("[INFO] Using VST matrix: ", vst_path, "\n", sep = "")
cat("[INFO] Using HVG list  : ", hvg_path, "\n", sep = "")
cat("[INFO] Using DSMZ meta : ", dsmz_meta_path, "\n", sep = "")
cat("[INFO] Output matrix   : ", out_path, "\n", sep = "")
cat("[INFO] Output mapping  : ", map_out, "\n", sep = "")

stopifnot(file.exists(vst_path))
stopifnot(file.exists(hvg_path))
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

# Read HVG gene list
hvg_genes <- scan(hvg_path, what = character(), quiet = TRUE)
cat("[INFO] HVG list length:", length(hvg_genes), "\n")

hvg_genes <- intersect(hvg_genes, rownames(vst_joint))
cat("[INFO] HVG genes present in VST:", length(hvg_genes), "\n")
if (length(hvg_genes) < 100) {
  stop("Too few HVG genes found in VST matrix: ", length(hvg_genes))
}

# Build samples × genes matrix
expr_hvg <- t(vst_joint[hvg_genes, , drop = FALSE])

# Extract DSMZ sample IDs (columns in original VST, now rownames in expr_hvg)
dsmz_sample_ids <- rownames(expr_hvg)[grepl("^NG-", rownames(expr_hvg))]

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
saveRDS(expr_hvg, out_path)
saveRDS(cell_line_to_sample_id, map_out)

cat("\n[OK] Saved HVG matrix:\n")
cat("  ", out_path, "\n")
cat("[OK] Saved mapping:\n")
cat("  ", map_out, "\n")
cat("[INFO] Final dims:", nrow(expr_hvg), "samples x", ncol(expr_hvg), "genes\n")
cat("[INFO] DSMZ:", sum(grepl("^NG-", rownames(expr_hvg))),
    "| TCGA:", sum(grepl("^TCGA-", rownames(expr_hvg))), "\n")
cat("[INFO] Mapping entries:", length(cell_line_to_sample_id), "DSMZ samples\n")

