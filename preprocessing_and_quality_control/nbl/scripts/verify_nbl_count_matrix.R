#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  if (length(args) %% 2 != 0) {
    stop("Expected arguments as --key value pairs.", call. = FALSE)
  }

  keys <- args[seq(1, length(args), by = 2)]
  vals <- args[seq(2, length(args), by = 2)]

  if (!all(grepl("^--", keys))) {
    stop("All argument keys must begin with --.", call. = FALSE)
  }

  setNames(as.list(vals), sub("^--", "", keys))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

required_args <- c(
  "matrix",
  "metadata",
  "output",
  "expected-total",
  "expected-target",
  "expected-srp409177",
  "expected-gse189367"
)
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "), call. = FALSE)
}

matrix_path <- args[["matrix"]]
metadata_path <- args[["metadata"]]
output_path <- args[["output"]]

expected_total <- as.integer(args[["expected-total"]])
expected_counts <- c(
  TARGET_NBL = as.integer(args[["expected-target"]]),
  SRP409177 = as.integer(args[["expected-srp409177"]]),
  GSE189367 = as.integer(args[["expected-gse189367"]])
)

if (file.exists(output_path)) {
  unlink(output_path)
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

normalize_cohort <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

if (!file.exists(matrix_path)) {
  fail("[ERROR] Count matrix does not exist: ", matrix_path)
}
if (!file.exists(metadata_path)) {
  fail("[ERROR] Metadata file does not exist: ", metadata_path)
}

message("[INFO] Reading count object: ", matrix_path)
count_obj <- readRDS(matrix_path)

if (inherits(count_obj, "SummarizedExperiment")) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    fail("[ERROR] Count object is SummarizedExperiment, but the package is unavailable.")
  }
  counts <- SummarizedExperiment::assay(count_obj)
} else if (is.matrix(count_obj)) {
  counts <- count_obj
} else if (is.data.frame(count_obj)) {
  numeric_cols <- vapply(count_obj, is.numeric, logical(1))
  if (!all(numeric_cols)) {
    fail(
      "[ERROR] Count data.frame contains non-numeric columns: ",
      paste(names(count_obj)[!numeric_cols], collapse = ", ")
    )
  }
  counts <- as.matrix(count_obj)
} else {
  fail("[ERROR] Unsupported count object class: ", paste(class(count_obj), collapse = ", "))
}

if (is.null(dim(counts)) || length(dim(counts)) != 2) {
  fail("[ERROR] Count object is not two-dimensional after extraction.")
}
if (is.null(rownames(counts)) || any(!nzchar(rownames(counts)))) {
  fail("[ERROR] Count matrix must have non-empty gene rownames.")
}
if (is.null(colnames(counts)) || any(!nzchar(colnames(counts)))) {
  fail("[ERROR] Count matrix must have non-empty sample column names.")
}

message("[INFO] Count dimensions: ", nrow(counts), " genes x ", ncol(counts), " samples")

if (ncol(counts) != expected_total) {
  fail("[ERROR] Count matrix has ", ncol(counts), " samples; expected ", expected_total, ".")
}

duplicated_samples <- unique(colnames(counts)[duplicated(colnames(counts))])
if (length(duplicated_samples) > 0) {
  fail("[ERROR] Duplicated count matrix sample columns: ", paste(duplicated_samples, collapse = ", "))
}

bad_gene_ids <- rownames(counts)[!grepl("^ENSG[0-9]+(\\.[0-9]+)?$", rownames(counts))]
if (length(bad_gene_ids) > 0) {
  fail(
    "[ERROR] ", length(bad_gene_ids),
    " count matrix rownames are not Ensembl-like ENSG identifiers. Examples: ",
    paste(head(bad_gene_ids, 10), collapse = ", ")
  )
}

if (!is.numeric(counts) && !is.integer(counts)) {
  fail("[ERROR] Count matrix storage mode is not numeric/integer: ", storage.mode(counts))
}
if (anyNA(counts)) {
  fail("[ERROR] Count matrix contains NA values.")
}
if (any(!is.finite(counts))) {
  fail("[ERROR] Count matrix contains non-finite values.")
}
if (any(counts < 0)) {
  fail("[ERROR] Count matrix contains negative counts.")
}

integer_deviation <- max(abs(counts - round(counts)))
if (integer_deviation > 1e-6) {
  fail(
    "[ERROR] Count matrix contains non-integer-like values; max deviation from integer = ",
    signif(integer_deviation, 4), "."
  )
}

message("[INFO] Reading metadata: ", metadata_path)
metadata <- read.csv(metadata_path, check.names = FALSE)

if (nrow(metadata) != expected_total) {
  fail("[ERROR] Metadata has ", nrow(metadata), " rows; expected ", expected_total, ".")
}

cohort_col <- intersect(c("cohort", "source", "project"), names(metadata))
if (length(cohort_col) == 0) {
  fail("[ERROR] Metadata must contain a cohort/source/project column.")
}
cohort_col <- cohort_col[1]

if (!"sample" %in% names(metadata)) {
  if ("aliquot_id" %in% names(metadata)) {
    metadata$sample <- metadata$aliquot_id
  } else {
    fail("[ERROR] Metadata must contain a sample column or aliquot_id column.")
  }
}

metadata[[cohort_col]] <- normalize_cohort(metadata[[cohort_col]])

if ("aliquot_id" %in% names(metadata)) {
  target_rows <- metadata[[cohort_col]] == "TARGET_NBL" &
    !is.na(metadata$aliquot_id) &
    nzchar(metadata$aliquot_id)
  metadata$sample[target_rows] <- metadata$aliquot_id[target_rows]
}

metadata$sample <- as.character(metadata$sample)
if (any(is.na(metadata$sample) | !nzchar(metadata$sample))) {
  fail("[ERROR] Metadata contains missing sample identifiers.")
}

duplicated_meta_samples <- unique(metadata$sample[duplicated(metadata$sample)])
if (length(duplicated_meta_samples) > 0) {
  fail("[ERROR] Duplicated metadata sample identifiers: ", paste(duplicated_meta_samples, collapse = ", "))
}

missing_in_metadata <- setdiff(colnames(counts), metadata$sample)
missing_in_counts <- setdiff(metadata$sample, colnames(counts))
if (length(missing_in_metadata) > 0 || length(missing_in_counts) > 0) {
  fail(
    "[ERROR] Count matrix columns and metadata samples do not match.",
    "\n  Missing in metadata: ", paste(head(missing_in_metadata, 20), collapse = ", "),
    "\n  Missing in counts: ", paste(head(missing_in_counts, 20), collapse = ", ")
  )
}

if (!identical(colnames(counts), metadata$sample)) {
  first_mismatch <- which(colnames(counts) != metadata$sample)[1]
  fail(
    "[ERROR] Count matrix columns and metadata samples have the same set but not the same order.",
    "\n  First mismatch at position ", first_mismatch,
    ": counts=", colnames(counts)[first_mismatch],
    ", metadata=", metadata$sample[first_mismatch]
  )
}

observed_counts <- table(metadata[[cohort_col]])
unexpected_cohorts <- setdiff(names(observed_counts)[observed_counts > 0], names(expected_counts))
if (length(unexpected_cohorts) > 0) {
  fail("[ERROR] Metadata contains unexpected cohorts/sources: ", paste(unexpected_cohorts, collapse = ", "))
}

for (cohort in names(expected_counts)) {
  observed <- unname(observed_counts[cohort])
  if (is.na(observed)) {
    observed <- 0L
  }
  expected <- expected_counts[[cohort]]
  if (observed != expected) {
    fail("[ERROR] Metadata cohort count for ", cohort, " is ", observed, "; expected ", expected, ".")
  }
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
writeLines(
  c(
    "OK",
    paste0("matrix=", matrix_path),
    paste0("metadata=", metadata_path),
    paste0("samples=", expected_total),
    paste0(names(expected_counts), "=", expected_counts)
  ),
  output_path
)

message("[OK] Final NBL count matrix verification passed: ", output_path)
