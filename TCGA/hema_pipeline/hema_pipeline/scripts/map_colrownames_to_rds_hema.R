#!/usr/bin/env Rscript
# map_colrownames_to_rds_hema.R
# Re-attach row/colnames from CSVs to numeric-only matrix RDS for HEMA.

suppressPackageStartupMessages({
  library(tibble)
})

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

required <- c("rds-in", "row-csv", "col-csv", "rds-out")
missing <- setdiff(required, names(opts))
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

RDS_IN  <- opts[["rds-in"]]
ROW_CSV <- opts[["row-csv"]]
COL_CSV <- opts[["col-csv"]]
RDS_OUT <- opts[["rds-out"]]

message("[INFO] Loading numeric-only matrix RDS: ", RDS_IN)
mat <- readRDS(RDS_IN)

if (!is.matrix(mat)) {
  stop("[ERROR] RDS object is not a matrix")
}

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

rownames(mat) <- samples
colnames(mat) <- genes

message("[INFO] Reattached dimnames: ",
        nrow(mat), " samples × ", ncol(mat), " genes")

message("[INFO] Saving fixed matrix to: ", RDS_OUT)
saveRDS(mat, RDS_OUT, compress = "gzip")
message("[SUCCESS] Saved fixed RDS!")

mat_check <- readRDS(RDS_OUT)
cat("[CHECK] Fixed matrix dims: ",
    paste(dim(mat_check), collapse = " × "), "\n")

if (is.character(rownames(mat_check)) && is.character(colnames(mat_check))) {
  message("[CHECK SUCCESS] Row/Col names are plain character vectors")
} else {
  warning("[CHECK FAIL] Row/Col names are not plain characters")
}
