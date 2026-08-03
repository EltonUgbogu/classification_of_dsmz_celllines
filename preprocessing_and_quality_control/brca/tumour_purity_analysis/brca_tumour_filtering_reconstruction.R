options(stringsAsFactors = FALSE)

for (path in unname(snakemake@output)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

log_path <- snakemake@output[["analysis_log"]]
log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

read_matrix_like <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("%s not found: %s", label, path))
  obj <- readRDS(path)
  if (methods::is(obj, "SummarizedExperiment")) {
    obj <- SummarizedExperiment::assay(obj)
  } else if (is.list(obj) && !is.null(obj$counts)) {
    obj <- obj$counts
  }
  if (!is.matrix(obj) && !is.data.frame(obj)) {
    stop(sprintf("%s must be matrix/data.frame-like; class=%s", label, paste(class(obj), collapse = ",")))
  }
  mat <- as.matrix(obj)
  if (!is.numeric(mat) || is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(sprintf("%s must be a named numeric matrix", label))
  }
  mat
}

tcga_sample_type <- function(ids) {
  ifelse(grepl("^TCGA-[^-]+-[^-]+-[0-9]{2}", ids), substr(ids, 14L, 15L), NA_character_)
}

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

mode <- as.character(snakemake@params[["mode"]])
expected_raw <- as.integer(snakemake@params[["expected_raw"]])
expected_retained <- as.integer(snakemake@params[["expected_retained"]])
primary_code <- as.character(snakemake@params[["primary_code"]])
original_purity_provenance <- isTRUE(snakemake@params[["original_purity_provenance"]])

if (!identical(mode, "active_sample_set_reconstruction")) {
  stop("Unsupported BRCA tumour-filtering mode: ", mode)
}
if (original_purity_provenance) {
  stop(paste0(
    "This reconstruction script must not be used when genuine original purity-scoring ",
    "provenance is declared available"
  ))
}
if (!file.exists(snakemake@input[["input_validation"]])) {
  stop("Count-input validation marker is missing")
}

cat("[INFO] BRCA tumour-filtered sample-set reconstruction\n")
cat("[INFO] mode=active_sample_set_reconstruction\n")
cat("[INFO] purity_scoring_performed=false\n")
cat("[INFO] original_purity_scoring_provenance_available=false\n")
cat("[INFO] The protected active matrix is read only to recover the established sample IDs.\n")

raw_counts <- read_matrix_like(snakemake@input[["tumour_rds"]], "raw TCGA-BRCA counts")
active_reference <- read_matrix_like(
  snakemake@input[["active_reference_rds"]],
  "active BRCA reference VST"
)

if (ncol(raw_counts) != expected_raw) {
  stop(sprintf("Expected %d raw TCGA samples; observed %d", expected_raw, ncol(raw_counts)))
}
if (anyDuplicated(colnames(raw_counts))) stop("Raw TCGA count matrix has duplicated sample IDs")

retained_ids <- colnames(active_reference)[grepl("^TCGA-", colnames(active_reference))]
if (length(retained_ids) != expected_retained) {
  stop(sprintf(
    "Expected %d retained TCGA sample IDs in active reference; observed %d",
    expected_retained,
    length(retained_ids)
  ))
}
if (anyDuplicated(retained_ids)) stop("Active reference contains duplicated TCGA sample IDs")

missing_from_raw <- setdiff(retained_ids, colnames(raw_counts))
if (length(missing_from_raw)) {
  stop(sprintf(
    "Active retained TCGA IDs missing from raw count matrix: %s",
    paste(missing_from_raw, collapse = ", ")
  ))
}

raw_types <- tcga_sample_type(colnames(raw_counts))
retained_types <- tcga_sample_type(retained_ids)
if (anyNA(retained_types) || any(retained_types != primary_code)) {
  bad <- retained_ids[is.na(retained_types) | retained_types != primary_code]
  stop(sprintf(
    "Retained set contains non-primary or unclassifiable TCGA samples: %s",
    paste(bad, collapse = ", ")
  ))
}

retained_counts <- raw_counts[, retained_ids, drop = FALSE]
if (!identical(colnames(retained_counts), retained_ids)) {
  stop("Retained count matrix order does not match the active reference TCGA order")
}
saveRDS(retained_counts, snakemake@output[["retained_counts"]], compress = TRUE)

retained_manifest <- data.frame(
  retained_order = seq_along(retained_ids),
  sample_id = retained_ids,
  tcga_sample_type = retained_types,
  retained = TRUE,
  selection_mode = mode,
  selection_basis = "TCGA sample IDs present in protected active BRCA matrix",
  purity_score_available = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(retained_manifest, snakemake@output[["retained_manifest"]])

removed_ids <- setdiff(colnames(raw_counts), retained_ids)
removed_types <- tcga_sample_type(removed_ids)
removed_reason <- ifelse(
  is.na(removed_types),
  "unclassifiable_tcga_sample_type",
  ifelse(
    removed_types != primary_code,
    paste0("not_primary_tumour_sample_type_", primary_code),
    "primary_tumour_absent_from_active_reference"
  )
)
removed_manifest <- data.frame(
  raw_order = match(removed_ids, colnames(raw_counts)),
  sample_id = removed_ids,
  tcga_sample_type = removed_types,
  retained = FALSE,
  removal_reason = removed_reason,
  selection_mode = mode,
  purity_score_available = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
removed_manifest <- removed_manifest[order(removed_manifest$raw_order), , drop = FALSE]
write_tsv(removed_manifest, snakemake@output[["removed_manifest"]])

raw_primary <- sum(!is.na(raw_types) & raw_types == primary_code)
removed_primary <- sum(!is.na(removed_types) & removed_types == primary_code)
sample_type_counts <- table(raw_types, useNA = "ifany")
sample_type_summary <- paste(
  paste(names(sample_type_counts), as.integer(sample_type_counts), sep = ":"),
  collapse = ";"
)

report <- data.frame(
  metric = c(
    "workflow_stage",
    "selection_mode",
    "purity_scoring_performed",
    "original_purity_scoring_provenance_available",
    "raw_tcga_samples",
    "raw_primary_tumour_samples",
    "retained_tumour_samples",
    "removed_samples_total",
    "removed_primary_tumour_samples",
    "removed_non_primary_samples",
    "raw_tcga_sample_type_counts",
    "retained_matrix_genes",
    "retained_matrix_samples"
  ),
  value = c(
    "tumour_filtered_sample_set_reconstruction",
    mode,
    "false",
    "false",
    ncol(raw_counts),
    raw_primary,
    ncol(retained_counts),
    length(removed_ids),
    removed_primary,
    sum(is.na(removed_types) | removed_types != primary_code),
    sample_type_summary,
    nrow(retained_counts),
    ncol(retained_counts)
  ),
  interpretation = c(
    "A reproducible lineage bridge from current count inputs to the established active sample set",
    "Sample IDs are reconstructed from the protected active BRCA matrix",
    "No ESTIMATE, tidyestimate, or substitute purity score was run",
    "The repository cannot link original purity scores to the active 818-sample set",
    "All columns in the declared raw TCGA count matrix",
    paste0("TCGA sample type ", primary_code),
    "Exact established active TCGA tumour set",
    "Raw columns not present in the established active set",
    "Primary tumours absent from the active reference",
    "Metastatic, normal, or unclassifiable raw columns",
    "Counts by TCGA sample-type code",
    "Count-level rows retained before later gene harmonisation",
    "Must equal the declared expected retained tumour count"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(report, snakemake@output[["filtering_report"]])

if (ncol(retained_counts) != expected_retained) {
  stop(sprintf(
    "Retained count matrix has %d samples; expected %d",
    ncol(retained_counts),
    expected_retained
  ))
}
if (raw_primary - removed_primary != expected_retained) {
  stop("Primary-tumour accounting does not reconcile to the retained sample count")
}

cat(sprintf("[INFO] Raw TCGA columns: %d\n", ncol(raw_counts)))
cat(sprintf("[INFO] Raw primary tumours: %d\n", raw_primary))
cat(sprintf("[INFO] Retained established tumours: %d\n", ncol(retained_counts)))
cat(sprintf("[INFO] Removed primary tumours: %d\n", removed_primary))
cat(sprintf("[INFO] Removed non-primary samples: %d\n", length(removed_ids) - removed_primary))
cat("[SUCCESS] BRCA tumour-filtered sample-set reconstruction completed\n")
