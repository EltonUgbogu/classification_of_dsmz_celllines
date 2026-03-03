#!/usr/bin/env Rscript
# convert_parquet_to_rds_rbl.R
# -------------------------------------------------
# Converts Parquet → RDS (numeric only) + separate
# CSVs for sample IDs (rows) and gene IDs (cols).
# This is meant to be run in an environment with
# the 'arrow' package available.
# -------------------------------------------------

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

# ---------- CONFIG ----------
PARQUET_PATH <- "/home/chu25/data/rbl/final_rbl_count_matrix.parquet"
RDS_PATH     <- "/home/chu25/data/rbl/rbl_raw_count_matrix_nodim.rds"  # <- no dimnames
ROW_CSV      <- "/home/chu25/data/rbl/rbl_sample_ids.csv"
COL_CSV      <- "/home/chu25/data/rbl/rbl_gene_ids.csv"
META_PATH    <- "/home/chu25/data/rbl/final_rbl_sample_metadata.csv"

message("[INFO] Loading Parquet: ", PARQUET_PATH)
df <- arrow::read_parquet(PARQUET_PATH)
message("   Loaded df: ", nrow(df), " × ", ncol(df))

# ---------- Identify sample ID column ----------
possible_idx <- c("__index_level_0__", "sample", "index", "Row.names")
idx_col <- intersect(possible_idx, colnames(df))[1]
if (is.na(idx_col)) stop("[ERROR] No index column found in Parquet file")

sample_ids <- as.character(df[[idx_col]])   # force plain character
df <- df %>% select(-all_of(idx_col))

# ---------- Build matrix (with dimnames, used only for checks) ----------
mat_full <- as.matrix(df)
rownames(mat_full) <- sample_ids
colnames(mat_full) <- as.character(colnames(df))

message("[INFO] Full matrix (with dimnames): ",
        nrow(mat_full), " samples × ", ncol(mat_full), " genes")
message("   First 5 samples: ", paste(head(rownames(mat_full), 5), collapse=", "))
message("   First 5 genes: ", paste(head(colnames(mat_full), 5), collapse=", "))

# ---------- Save row/col names to CSV ----------
row_df <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
col_df <- data.frame(gene   = as.character(colnames(df)), stringsAsFactors = FALSE)

write.csv(row_df, ROW_CSV, row.names = FALSE)
write.csv(col_df, COL_CSV, row.names = FALSE)
message("[INFO] Saved rownames to: ", ROW_CSV)
message("[INFO] Saved colnames to: ", COL_CSV)

# ---------- Strip dimnames and save numeric-only matrix ----------
mat_nodim <- unname(mat_full)   # removes rownames + colnames completely

message("[INFO] Saving numeric-only RDS (no dimnames): ", RDS_PATH)
saveRDS(mat_nodim, RDS_PATH, compress = "gzip")
message("[SUCCESS] Saved RDS (no dimnames)! Size: ",
        round(file.size(RDS_PATH)/1e6, 2), " MB")

# ---------- Optional: sample order check vs metadata ----------
if (file.exists(META_PATH)) {
  meta <- read.csv(META_PATH, stringsAsFactors = FALSE)
  if (!("sample" %in% colnames(meta))) {
    warning("[WARN] Metadata missing 'sample' column: ", META_PATH)
  } else {
    meta$sample <- as.character(meta$sample)
    if (identical(sample_ids, meta$sample)) {
      message("[SUCCESS] Sample order matches metadata")
    } else {
      message("[WARN] Sample order mismatch between Parquet and metadata")
      message("   Matrix:   ", paste(head(sample_ids), collapse = ", "))
      message("   Metadata: ", paste(head(meta$sample), collapse = ", "))
    }
  }
}

# ---------- Quick check: reload RDS in this env ----------
mat_test <- readRDS(RDS_PATH)
cat("[CHECK] Numeric-only matrix dims: ",
    paste(dim(mat_test), collapse = " × "), "\n")

if (is.null(rownames(mat_test)) && is.null(colnames(mat_test))) {
  message("[CHECK SUCCESS] RDS has *no* dimnames (as intended)")
} else {
  warning("[CHECK FAIL] RDS still has dimnames, this should not happen")
}
