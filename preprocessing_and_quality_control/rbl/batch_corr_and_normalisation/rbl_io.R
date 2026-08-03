DSMZ_RAW_COUNT_ERROR <- paste(
  "Configured DSMZ input is not a raw integer count matrix.",
  "Do not use transformed VST/log expression for batch-correction input."
)

stop_invalid_dsmz_counts <- function(detail) {
  stop(sprintf("%s %s", DSMZ_RAW_COUNT_ERROR, detail), call. = FALSE)
}

validate_count_matrix_values <- function(counts, path) {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop_invalid_dsmz_counts(sprintf(
      "Expected matrix/data.frame-like input at %s.",
      path
    ))
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop_invalid_dsmz_counts("Row and column names are required.")
  }
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop_invalid_dsmz_counts("Duplicated row or column names were detected.")
  }
  if (anyNA(counts)) {
    stop_invalid_dsmz_counts(sprintf("NA values were detected in %s.", path))
  }
  if (any(!is.finite(counts))) {
    stop_invalid_dsmz_counts(sprintf("Non-finite values were detected in %s.", path))
  }
  if (any(counts < 0)) {
    stop_invalid_dsmz_counts(sprintf("Negative values were detected in %s.", path))
  }
  if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    stop_invalid_dsmz_counts(sprintf("Non-integer values were detected in %s.", path))
  }
  invisible(TRUE)
}

read_dsmz_tsv_counts <- function(counts_path) {
  raw <- utils::read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(raw) < 2) {
    stop(sprintf("[ERROR] DSMZ TSV must contain a gene identifier column and at least one sample column: %s", counts_path))
  }

  gene_col <- colnames(raw)[1]
  gene_ids <- raw[[1]]
  gene_ids_chr <- as.character(gene_ids)
  numeric_gene_ids <- suppressWarnings(as.numeric(gene_ids_chr))
  first_col_is_numeric <- all(!is.na(numeric_gene_ids))
  if (tolower(gene_col) != "gene_id" && first_col_is_numeric) {
    stop(sprintf("[ERROR] First DSMZ TSV column must be gene_id or another non-numeric gene identifier column: %s", gene_col))
  }
  if (any(is.na(gene_ids_chr)) || any(!nzchar(gene_ids_chr))) {
    stop_invalid_dsmz_counts("The gene identifier column contains missing or empty values.")
  }
  count_df <- raw[, -1, drop = FALSE]
  if (anyDuplicated(gene_ids_chr) || anyDuplicated(colnames(count_df))) {
    stop_invalid_dsmz_counts("Duplicated gene identifiers or sample columns were detected.")
  }
  counts <- as.matrix(data.frame(lapply(count_df, function(x) suppressWarnings(as.numeric(x))), check.names = FALSE))
  colnames(counts) <- colnames(count_df)
  rownames(counts) <- gene_ids_chr
  storage.mode(counts) <- "double"
  validate_count_matrix_values(counts, counts_path)
  harmonize_count_matrix(counts)
}

read_dsmz_rds_counts <- function(counts_path) {
  raw <- tryCatch(
    readRDS(counts_path),
    error = function(e) stop_invalid_dsmz_counts(sprintf(
      "readRDS failed for %s: %s.",
      counts_path,
      conditionMessage(e)
    ))
  )
  if (methods::is(raw, "SummarizedExperiment")) {
    raw <- SummarizedExperiment::assay(raw)
  }
  if (is.list(raw) && !is.null(raw$counts)) {
    raw <- raw$counts
  }
  if (is.data.frame(raw) && all(c("Ensembl_ID", "gene_name") %in% colnames(raw))) {
    annotation <- c("Ensembl_ID", "gene_name", "Ensembl_ID_with_version")
    count_cols <- setdiff(colnames(raw), annotation)
    payload <- as.matrix(raw[, count_cols, drop = FALSE])
    storage.mode(payload) <- "numeric"
    rownames(payload) <- rownames(raw)
    validate_count_matrix_values(payload, counts_path)
    return(build_dsmz_matrix(raw))
  }
  if (is.matrix(raw) || is.data.frame(raw)) {
    validate_count_matrix_values(raw, counts_path)
    return(harmonize_count_matrix(raw))
  }
  stop_invalid_dsmz_counts(
    "Expected SummarizedExperiment, matrix/data.frame, list(counts=...), or an annotated DSMZ count table."
  )
}

validate_dsmz_metadata_matches <- function(counts, metadata, sample_col = "sample_name") {
  if (!sample_col %in% colnames(metadata)) {
    stop(sprintf("[ERROR] DSMZ metadata sample column not found: %s", sample_col))
  }
  sample_ids <- metadata[[sample_col]]
  match_counts <- vapply(colnames(counts), function(sample_id) {
    sum(sample_ids == sample_id, na.rm = TRUE)
  }, integer(1))
  if (any(match_counts != 1L)) {
    bad <- names(match_counts)[match_counts != 1L]
    details <- paste(sprintf("%s:%d", bad, match_counts[bad]), collapse = ", ")
    stop(sprintf("[ERROR] DSMZ count columns must match exactly one metadata row in %s. Offending samples: %s", sample_col, details))
  }
  invisible(TRUE)
}

load_dsmz_data <- function(counts_path, meta_path, sample_col = "sample_name") {
  cat("[INFO] Loading DSMZ data...\n")
  if (!file.exists(counts_path)) {
    stop_invalid_dsmz_counts(sprintf("Configured file does not exist: %s.", counts_path))
  }
  if (!file.exists(meta_path)) stop(sprintf("[ERROR] DSMZ metadata CSV not found: %s", meta_path))
  meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)

  ext <- tolower(tools::file_ext(counts_path))
  counts <- switch(
    ext,
    rds = read_dsmz_rds_counts(counts_path),
    tsv = read_dsmz_tsv_counts(counts_path),
    txt = read_dsmz_tsv_counts(counts_path),
    stop(sprintf("[ERROR] Unsupported DSMZ count matrix extension '.%s': %s", ext, counts_path))
  )
  validate_dsmz_metadata_matches(counts, meta, sample_col = sample_col)
  cat(sprintf("[INFO] DSMZ count matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))
  list(counts = counts, meta = meta)
}

harmonize_count_matrix <- function(counts) {
  counts <- as.matrix(counts)
  if (is.null(rownames(counts))) stop("Count matrix must have rownames")
  rownames(counts) <- sub("\\..*$", "", rownames(counts))
  if (any(duplicated(rownames(counts)))) {
    counts <- rowsum(counts, group = rownames(counts), reorder = TRUE)
  }
  storage.mode(counts) <- "double"
  counts
}

load_tumour_rbl_data <- function(file_path) {
  cat("[INFO] Loading tumour data...\n")
  if (!file.exists(file_path)) stop(sprintf("[ERROR] RBL tumour input not found: %s", file_path))
  obj <- readRDS(file_path)

  if (methods::is(obj, "SummarizedExperiment")) {
    counts <- SummarizedExperiment::assay(obj)
  } else if (is.list(obj) && !is.null(obj$counts)) {
    counts <- obj$counts
  } else if (is.matrix(obj) || is.data.frame(obj)) {
    counts <- obj
  } else {
    stop("Unsupported tumour input type. Expected SummarizedExperiment, matrix/data.frame, or list(counts=...).")
  }

  counts <- harmonize_count_matrix(counts)
  cat(sprintf("[INFO] Tumour matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))
  counts
}

load_dsmz_rbl_data <- function(counts, metadata, filter_col = "Disease", filter_values = "Retinoblastoma", sample_col = "sample_name", outdir = NULL) {
  cat("[INFO] Filtering DSMZ metadata...\n")
  if (!filter_col %in% colnames(metadata)) {
    stop(sprintf("Metadata column %s not found", filter_col))
  }
  keep_meta <- metadata[metadata[[filter_col]] %in% filter_values, , drop = FALSE]
  if (nrow(keep_meta) == 0) {
    stop(sprintf("No DSMZ samples matched %s in %s", paste(filter_values, collapse = ", "), filter_col))
  }
  aligned <- align_counts_metadata(keep_meta, counts, col_name = sample_col, outdir = outdir)
  cat(sprintf("[INFO] DSMZ filtered matrix: %d genes x %d samples\n", nrow(aligned$counts), ncol(aligned$counts)))
  list(counts = aligned$counts, metadata = aligned$meta)
}
