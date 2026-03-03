#!/usr/bin/env Rscript
# map_colrownames_to_rds.R
# -------------------------------------------------
# Re-attaches row/colnames from CSVs to the numeric
# matrix RDS created in step 1. No 'arrow' needed.
# -------------------------------------------------

suppressPackageStartupMessages({
  library(tibble)
})

RDS_IN   <- "/home/chu25/data/rbl/rbl_raw_count_matrix_nodim.rds"
ROW_CSV  <- "/home/chu25/data/rbl/rbl_sample_ids.csv"
COL_CSV  <- "/home/chu25/data/rbl/rbl_gene_ids.csv"
RDS_OUT  <- "/home/chu25/data/rbl/rbl_raw_count_matrix.rds"

message("[INFO] Loading numeric-only matrix RDS: ", RDS_IN)
mat <- readRDS(RDS_IN)

if (!is.matrix(mat)) {
  stop("[ERROR] RDS object is not a matrix")
}

# ---------- Load dimnames ----------
message("[INFO] Loading rownames from: ", ROW_CSV)
row_df <- read.csv(ROW_CSV, stringsAsFactors = FALSE)
if (!"sample" %in% colnames(row_df)) {
  stop("[ERROR] ROW_CSV must have a 'sample' column")
}
samples <- as.character(row_df$sample)

message("[INFO] Loading colnames from: ", COL_CSV)
col_df <- read.csv(COL_CSV, stringsAsFactors = FALSE)
if (!"gene" %in% colnames(col_df)) {
  stop("[ERROR] COL_CSV must have a 'gene' column")
}
genes <- as.character(col_df$gene)

# ---------- Sanity checks ----------
if (length(samples) != nrow(mat)) {
  stop(sprintf(
    "[ERROR] Length of samples (%d) != nrow(mat) (%d)",
    length(samples), nrow(mat)
  ))
}
if (length(genes) != ncol(mat)) {
  stop(sprintf(
    "[ERROR] Length of genes (%d) != ncol(mat) (%d)",
    length(genes), ncol(mat)
  ))
}

# ---------- Attach dimnames ----------
rownames(mat) <- samples
colnames(mat) <- genes

message("[INFO] Reattached dimnames: ",
        nrow(mat), " samples × ", ncol(mat), " genes")
message("   First 5 samples: ", paste(head(rownames(mat), 5), collapse = ", "))
message("   First 5 genes: ", paste(head(colnames(mat), 5), collapse = ", "))

# ---------- Save FINAL portable matrix ----------
message("[INFO] Saving fixed matrix to: ", RDS_OUT)
saveRDS(mat, RDS_OUT, compress = "gzip")
message("[SUCCESS] Saved fixed RDS!")

# ---------- Final check ----------
mat_check <- readRDS(RDS_OUT)
cat("[CHECK] Fixed matrix dims: ",
    paste(dim(mat_check), collapse = " × "), "\n")
cat("[CHECK] First 5 samples: ",
    paste(head(rownames(mat_check), 5), collapse = ", "), "\n")
cat("[CHECK] First 5 genes: ",
    paste(head(colnames(mat_check), 5), collapse = ", "), "\n")

if (is.character(rownames(mat_check)) && is.character(colnames(mat_check))) {
  message("[CHECK SUCCESS] Row/Col names are plain character vectors")
} else {
  warning("[CHECK FAIL] Row/Col names are not plain characters")
}
