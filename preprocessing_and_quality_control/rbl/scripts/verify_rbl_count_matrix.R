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
  "expected-gse111168",
  "expected-gse196420",
  "expected-gse268136"
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
  GSE111168_primary_tumours = as.integer(args[["expected-gse111168"]]),
  GSE196420_primary_tumours = as.integer(args[["expected-gse196420"]]),
  GSE268136_primary_tumours = as.integer(args[["expected-gse268136"]])
)

if (file.exists(output_path)) {
  unlink(output_path)
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
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

if (!"cohort" %in% names(metadata)) {
  fail("[ERROR] Metadata must contain a cohort column.")
}

if (!"sample" %in% names(metadata)) {
  alternatives <- intersect(c("sample_id", "run", "Run", "srr", "SRR"), names(metadata))
  if (length(alternatives) == 0) {
    fail("[ERROR] Metadata must contain a sample column or one of: sample_id, run, Run, srr, SRR.")
  }
  metadata$sample <- metadata[[alternatives[1]]]
}

metadata$sample <- as.character(metadata$sample)
metadata$cohort <- trimws(as.character(metadata$cohort))

if (any(is.na(metadata$sample) | !nzchar(metadata$sample))) {
  fail("[ERROR] Metadata contains missing sample identifiers.")
}
if (any(is.na(metadata$cohort) | !nzchar(metadata$cohort))) {
  fail("[ERROR] Metadata contains missing cohort labels.")
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

observed_counts <- table(metadata$cohort)
unexpected_cohorts <- setdiff(names(observed_counts)[observed_counts > 0], names(expected_counts))
if (length(unexpected_cohorts) > 0) {
  fail("[ERROR] Metadata contains unexpected cohorts: ", paste(unexpected_cohorts, collapse = ", "))
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

message("[OK] Final RBL count matrix verification passed: ", output_path)
