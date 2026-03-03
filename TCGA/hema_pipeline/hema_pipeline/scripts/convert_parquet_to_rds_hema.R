#!/usr/bin/env Rscript
# convert_parquet_to_rds_hema.R
# Parquet → numeric-only RDS (+ row/col ID CSVs) for HEMA dataset.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

# tiny arg parser: --key value
parse_cli <- function(args) {
  res <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) stop("Expected named arg starting with '--', got: ", key)
    if (i == length(args)) stop("No value provided for ", key)
    val <- args[i + 1]
    res[[sub("^--", "", key)]] <- val
    i <- i + 2
  }
  res
}

args <- commandArgs(trailingOnly = TRUE)
opts <- parse_cli(args)

required <- c("parquet", "rds-nodim", "row-csv", "col-csv", "meta")
missing <- setdiff(required, names(opts))
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

PARQUET_PATH <- opts[["parquet"]]
RDS_PATH     <- opts[["rds-nodim"]]
ROW_CSV      <- opts[["row-csv"]]
COL_CSV      <- opts[["col-csv"]]
META_PATH    <- opts[["meta"]]

message("[INFO] Loading Parquet: ", PARQUET_PATH)
df <- arrow::read_parquet(PARQUET_PATH)
message("   Loaded df: ", nrow(df), " × ", ncol(df))

possible_idx <- c("__index_level_0__", "sample", "index", "Row.names")
idx_col <- intersect(possible_idx, colnames(df))[1]
if (is.na(idx_col)) stop("[ERROR] No index column found in Parquet file")

sample_ids <- as.character(df[[idx_col]])
df <- df %>% select(-all_of(idx_col))

mat_full <- as.matrix(df)
rownames(mat_full) <- sample_ids
colnames(mat_full) <- as.character(colnames(df))

message("[INFO] Full matrix (with dimnames): ",
        nrow(mat_full), " samples × ", ncol(mat_full), " genes")

row_df <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
col_df <- data.frame(gene   = as.character(colnames(df)), stringsAsFactors = FALSE)

write.csv(row_df, ROW_CSV, row.names = FALSE)
write.csv(col_df, COL_CSV, row.names = FALSE)
message("[INFO] Saved rownames to: ", ROW_CSV)
message("[INFO] Saved colnames to: ", COL_CSV)

mat_nodim <- unname(mat_full)

message("[INFO] Saving numeric-only RDS (no dimnames): ", RDS_PATH)
saveRDS(mat_nodim, RDS_PATH, compress = "gzip")
message("[SUCCESS] Saved RDS (no dimnames)! Size: ",
        round(file.size(RDS_PATH)/1e6, 2), " MB")

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

mat_test <- readRDS(RDS_PATH)
cat("[CHECK] Numeric-only matrix dims: ",
    paste(dim(mat_test), collapse = " × "), "\n")

if (is.null(rownames(mat_test)) && is.null(colnames(mat_test))) {
  message("[CHECK SUCCESS] RDS has *no* dimnames (as intended)")
} else {
  warning("[CHECK FAIL] RDS still has dimnames, this should not happen")
}
