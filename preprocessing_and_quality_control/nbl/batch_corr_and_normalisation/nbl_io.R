load_dsmz_data <- function(counts_path, meta_path) {
  cat("[INFO] Loading DSMZ data...\n")
  if (!file.exists(counts_path)) stop(sprintf("[ERROR] DSMZ RDS not found: %s", counts_path))
  if (!file.exists(meta_path)) stop(sprintf("[ERROR] DSMZ metadata CSV not found: %s", meta_path))
  raw <- readRDS(counts_path)
  meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("Ensembl_ID", "gene_name") %in% colnames(raw)))
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

load_tumour_nbl_data <- function(file_path) {
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
