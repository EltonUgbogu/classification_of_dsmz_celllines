DSMZ_RAW_COUNT_ERROR <- paste(
  "Configured DSMZ input is not a raw integer count matrix.",
  "Do not use transformed VST/log expression for batch-correction input."
)

stop_invalid_dsmz_counts <- function(detail) {
  stop(sprintf("%s %s", DSMZ_RAW_COUNT_ERROR, detail), call. = FALSE)
}

validate_dsmz_count_object <- function(raw, counts_path) {
  if (!(is.matrix(raw) || is.data.frame(raw))) {
    stop_invalid_dsmz_counts(sprintf(
      "Expected matrix/data.frame-like RDS at %s; class=%s.",
      counts_path,
      paste(class(raw), collapse = ",")
    ))
  }
  if (is.null(rownames(raw)) || is.null(colnames(raw))) {
    stop_invalid_dsmz_counts("Row and column names are required.")
  }
  if (anyDuplicated(rownames(raw)) || anyDuplicated(colnames(raw))) {
    stop_invalid_dsmz_counts("Duplicated row or column names were detected.")
  }

  annotation_cols <- intersect(
    c("Ensembl_ID", "gene_name", "Ensembl_ID_with_version"),
    colnames(raw)
  )
  count_cols <- setdiff(colnames(raw), annotation_cols)
  if (!length(count_cols)) {
    stop_invalid_dsmz_counts("No count columns were found.")
  }
  if (is.data.frame(raw) && !all(vapply(raw[count_cols], is.numeric, logical(1)))) {
    stop_invalid_dsmz_counts("All sample columns must be numeric.")
  }
  counts <- as.matrix(raw[, count_cols, drop = FALSE])
  if (!is.numeric(counts)) {
    stop_invalid_dsmz_counts("All sample columns must be numeric.")
  }
  if (anyNA(counts) || any(!is.finite(counts))) {
    stop_invalid_dsmz_counts("NA or non-finite values were detected.")
  }
  if (any(counts < 0)) {
    stop_invalid_dsmz_counts("Negative values were detected.")
  }
  if (any(abs(counts - round(counts)) > 1e-8)) {
    stop_invalid_dsmz_counts("Non-integer values were detected.")
  }
  invisible(TRUE)
}

load_dsmz_data <- function(counts_path, meta_path) {
  cat("[INFO] Loading DSMZ data...\n")
  if (!file.exists(counts_path)) {
    stop_invalid_dsmz_counts(sprintf("Configured file does not exist: %s.", counts_path))
  }
  if (!file.exists(meta_path)) stop(sprintf("[ERROR] DSMZ metadata CSV not found: %s", meta_path))
  raw <- tryCatch(
    readRDS(counts_path),
    error = function(e) stop_invalid_dsmz_counts(sprintf(
      "readRDS failed for %s: %s.",
      counts_path,
      conditionMessage(e)
    ))
  )
  validate_dsmz_count_object(raw, counts_path)
  meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Ensembl_ID", "gene_name") %in% colnames(raw))) {
    stop_invalid_dsmz_counts(
      "The configured NBL DSMZ table must include Ensembl_ID and gene_name annotations."
    )
  }
  stopifnot("sample_name" %in% colnames(meta))
  list(raw = raw, meta = meta)
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

load_tumour_nbl_data <- function(
    file_path,
    meta_path = NULL,
    expected_total = NULL,
    expected_cohorts = NULL) {
  cat("[INFO] Loading tumour data...\n")
  if (!file.exists(file_path)) stop(sprintf("[ERROR] NBL tumour input not found: %s", file_path))
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
  if (!is.numeric(counts)) {
    stop("[ERROR] NBL tumour matrix must contain numeric counts")
  }
  if (anyNA(counts) || any(!is.finite(counts))) {
    stop("[ERROR] NBL tumour matrix contains missing or non-finite values")
  }
  if (any(counts < 0)) {
    stop("[ERROR] NBL tumour matrix contains negative values")
  }
  integer_like <- all(vapply(
    seq_len(ncol(counts)),
    function(j) all(abs(counts[, j] - round(counts[, j])) <= 1e-6),
    logical(1)
  ))
  if (!integer_like) {
    stop("[ERROR] NBL tumour matrix is not count-like: non-integer values detected")
  }

  if (!is.null(expected_total) && expected_total > 0 && ncol(counts) != expected_total) {
    stop(sprintf(
      "[ERROR] NBL tumour matrix has %d samples; expected %d",
      ncol(counts),
      expected_total
    ))
  }

  if (!is.null(meta_path) && nzchar(meta_path)) {
    if (!file.exists(meta_path)) {
      stop(sprintf("[ERROR] NBL tumour metadata CSV not found: %s", meta_path))
    }
    metadata <- if (grepl("\\.tsv$", meta_path, ignore.case = TRUE)) {
      utils::read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      utils::read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    sample_columns <- intersect(c("sample", "sample_id", "aliquot_id"), colnames(metadata))
    if (length(sample_columns) == 0) {
      stop("[ERROR] NBL tumour metadata lacks a sample, sample_id, or aliquot_id column")
    }
    sample_ids <- as.character(metadata[[sample_columns[1]]])
    if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
      stop("[ERROR] NBL tumour metadata contains missing, empty, or duplicated sample IDs")
    }
    if (nrow(metadata) != ncol(counts)) {
      stop(sprintf(
        "[ERROR] NBL tumour metadata has %d rows but the matrix has %d samples",
        nrow(metadata),
        ncol(counts)
      ))
    }
    if (!identical(colnames(counts), sample_ids)) {
      stop("[ERROR] NBL tumour matrix columns do not exactly match metadata sample IDs in order")
    }

    if (!is.null(expected_cohorts) && length(expected_cohorts) > 0) {
      cohort_columns <- intersect(c("cohort", "source", "project"), colnames(metadata))
      if (length(cohort_columns) == 0) {
        stop("[ERROR] NBL tumour metadata lacks a cohort, source, or project column")
      }
      observed <- table(as.character(metadata[[cohort_columns[1]]]))
      unexpected <- setdiff(names(observed), names(expected_cohorts))
      if (length(unexpected) > 0) {
        stop(sprintf(
          "[ERROR] Unexpected NBL tumour cohorts: %s",
          paste(unexpected, collapse = ", ")
        ))
      }
      for (cohort in names(expected_cohorts)) {
        observed_n <- if (cohort %in% names(observed)) unname(observed[[cohort]]) else 0L
        if (observed_n != expected_cohorts[[cohort]]) {
          stop(sprintf(
            "[ERROR] NBL tumour cohort %s has %d samples; expected %d",
            cohort,
            observed_n,
            expected_cohorts[[cohort]]
          ))
        }
      }
    }
    cat(sprintf(
      "[INFO] Tumour metadata: %d samples; exact matrix order confirmed\n",
      nrow(metadata)
    ))
  }

  cat(sprintf("[INFO] Tumour matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))
  counts
}

load_dsmz_nbl_data <- function(counts, metadata, filter_col = "Disease", filter_values = "Neuroblastoma", sample_col = "sample_name", outdir = NULL) {
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
