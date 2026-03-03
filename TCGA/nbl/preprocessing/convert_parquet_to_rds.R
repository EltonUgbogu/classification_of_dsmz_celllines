#!/usr/bin/env Rscript
# convert_parquet_to_rds.R
# -------------------------------------------------
# Converts a Parquet count matrix → portable RDS
# Forces all rownames/colnames to plain characters
# -------------------------------------------------

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

# ---------- CONFIG ----------
PARQUET_PATH <- "/home/chu25/data/nbl/final_nbl_count_matrix.parquet"
RDS_PATH     <- "/home/chu25/data/nbl/nbl_raw_count_matrix.rds"
META_PATH    <- "/home/chu25/data/nbl/final_nbl_sample_metadata.csv"

message("[INFO] Loading Parquet: ", PARQUET_PATH)
df <- arrow::read_parquet(PARQUET_PATH)
message("   Loaded df: ", nrow(df), " × ", ncol(df))

# ---------- Identify sample ID column ----------
possible_idx <- c("__index_level_0__", "sample", "index", "Row.names")
idx_col <- intersect(possible_idx, colnames(df))[1]
if (is.na(idx_col)) stop("[ERROR] No index column found")

sample_ids <- as.character(df[[idx_col]])   # force plain character
df <- df %>% select(-all_of(idx_col))

# ---------- Build matrix ----------
mat <- as.matrix(df)
rownames(mat) <- sample_ids
colnames(mat) <- as.character(colnames(df))  # force plain character

# Force dimnames materialization (critical fix)
rownames(mat) <- as.character(rownames(mat))
colnames(mat) <- as.character(colnames(mat))

message("[INFO] Matrix: ", nrow(mat), " samples × ", ncol(mat), " genes")
message("   First 5 samples: ", paste(head(rownames(mat), 5), collapse=", "))
message("   First 5 genes: ", paste(head(colnames(mat), 5), collapse=", "))

# ---------- Save portable RDS ----------
message("[INFO] Saving RDS: ", RDS_PATH)
saveRDS(mat, RDS_PATH, compress="gzip")
message("[SUCCESS] Saved RDS! Size: ", round(file.size(RDS_PATH)/1e6, 2), " MB")

# ---------- Verify sample order ----------
meta <- read.csv(META_PATH, stringsAsFactors = FALSE)
if (!("sample" %in% colnames(meta))) stop("[ERROR] Metadata missing 'sample' column")
meta$sample <- as.character(meta$sample)

if (identical(rownames(mat), meta$sample)) {
  message("[SUCCESS] Sample order matches metadata ✅")
} else {
  message("[WARN] Sample order mismatch ⚠️")
  message("   Matrix: ", paste(head(rownames(mat)), collapse = ", "), "")
  message("   Metadata: ", paste(head(meta$sample), collapse = ", "), "")
}

# ---------- Quick check: reload RDS and verify dimnames ----------
mat_test <- readRDS(RDS_PATH)
cat("[CHECK] Matrix dims: ", paste(dim(mat_test), collapse = " × "), "\n")
cat("[CHECK] First 5 samples: ",
    paste(head(rownames(mat_test), 5), collapse = ", "), "\n")
cat("[CHECK] First 5 genes: ",
    paste(head(colnames(mat_test), 5), collapse = ", "), "\n")

if (is.character(rownames(mat_test)) && is.character(colnames(mat_test))) {
  message("[CHECK SUCCESS] Row/Col names are plain characters")
} else {
  warning("[CHECK FAIL] Row/Col names still not plain characters")
}
