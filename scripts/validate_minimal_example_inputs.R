#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: validate_minimal_example_inputs.R <examples/minimal>")
}

base <- args[[1]]
profiles <- c("brca", "nbl", "rbl")

required_cols <- c(
  "sample_id",
  "profile",
  "cancer_type",
  "sample_type",
  "source"
)

for (profile in profiles) {
  rds_path <- file.path(base, profile, paste0(profile, "_vst_joint.rds"))
  meta_path <- file.path(base, profile, "metadata.tsv")

  if (!file.exists(rds_path)) {
    stop("Missing RDS file: ", rds_path)
  }

  if (!file.exists(meta_path)) {
    stop("Missing metadata file: ", meta_path)
  }

  x <- readRDS(rds_path)
  meta <- read.delim(meta_path, stringsAsFactors = FALSE)

  if (!is.matrix(x)) {
    stop("Expected matrix RDS for ", profile, ", got: ", paste(class(x), collapse = ", "))
  }

  if (!is.numeric(x)) {
    stop("Expected numeric matrix for ", profile)
  }

  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("Matrix must have gene row names and sample column names for ", profile)
  }

  missing_cols <- setdiff(required_cols, colnames(meta))

  if (length(missing_cols) > 0) {
    stop("Missing metadata columns for ", profile, ": ", paste(missing_cols, collapse = ", "))
  }

  if (!all(meta$sample_id %in% colnames(x))) {
    stop("Some metadata sample_id values are not present in matrix columns for ", profile)
  }

  if (!all(colnames(x) %in% meta$sample_id)) {
    stop("Some matrix columns are not present in metadata sample_id for ", profile)
  }

  if (!all(meta$sample_type %in% c("cell_line", "tumour"))) {
    stop("Unexpected sample_type values for ", profile)
  }

  if (sum(meta$sample_type == "cell_line") < 2) {
    stop("Expected at least two synthetic cell-line samples for ", profile)
  }

  if (sum(meta$sample_type == "tumour") < 2) {
    stop("Expected at least two synthetic tumour samples for ", profile)
  }

  cat("Validated minimal example for", profile, "\n")
}

cat("All minimal example inputs validated successfully.\n")
