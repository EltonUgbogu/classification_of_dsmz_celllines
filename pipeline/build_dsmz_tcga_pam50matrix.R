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

# --------------------------------------------------------------------
# 0. CLI + config
# --------------------------------------------------------------------
option_list <- list(
  make_option(
    "--config",
    type = "character",
    default = "config/config.yaml",
    help = "Path to config.yaml [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Using config file:", opt$config, "\n")
cfg <- yaml::read_yaml(opt$config)

`%||%` <- function(x, y) if (!is.null(x)) x else y

unsup_root  <- cfg$paths$unsup_root
input_root  <- cfg$paths$tumour_nh_input_root

tcga_mat <- readRDS(cfg$paths$tcga_brca_pam50_expr %||% cfg$paths$tcga_pam50_expr)
dsmz_mat <- readRDS(cfg$paths$dsmz_bcc_pam50_expr %||% cfg$paths$dsmz_pam50_expr)

# align genes
common_genes <- intersect(rownames(tcga_mat), rownames(dsmz_mat))
stopifnot(length(common_genes) >= 40)
tcga_mat <- tcga_mat[common_genes, , drop = FALSE]
dsmz_mat <- dsmz_mat[common_genes, , drop = FALSE]

# DSMZ metadata
meta <- read.csv(cfg$paths$dsmz_meta_csv, stringsAsFactors = FALSE)
stopifnot("sample_name" %in% colnames(meta))
cell_col <- intersect(c("Cell_line", "Cell_Line", "CellLine"), colnames(meta))[1]
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

dir.create(input_root, showWarnings = FALSE, recursive = TRUE)

cell_line_to_sample_id <- cell_line_of[dsmz_original_ids]
cell_line_to_sample_id[is.na(cell_line_to_sample_id)] <- dsmz_original_ids[is.na(cell_line_to_sample_id)]

expr_out <- file.path(input_root, "BRCA_TCGA-DSMZ_PAM50_samples_x_genes.rds")
map_out  <- file.path(input_root, "cell_line_to_original_sample_id_pam50.rds")

saveRDS(expr_mat, expr_out)
saveRDS(cell_line_to_sample_id, map_out)

cat(sprintf("Matrix: %d samples × %d genes (%d DSMZ | %d TCGA)\n",
            nrow(expr_mat), ncol(expr_mat),
            sum(dataset_vec == "DSMZ"), sum(dataset_vec == "TCGA")))
cat("Expression matrix saved to:", expr_out, "\n")
cat("Mapping saved to:", map_out, "\n")

