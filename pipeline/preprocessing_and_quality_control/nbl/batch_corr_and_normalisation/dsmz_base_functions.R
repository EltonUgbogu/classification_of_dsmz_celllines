build_dsmz_matrix <- function(dsmz_raw) {
  cat("[INFO] Building DSMZ count matrix...\n")
  annot <- c("Ensembl_ID", "gene_name", "Ensembl_ID_with_version")
  sample_cols <- setdiff(colnames(dsmz_raw), annot)
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(x) as.numeric(as.character(x)))
  M <- as.matrix(dsmz_raw[, sample_cols, drop = FALSE])
  rownames(M) <- sub("\\..*$", "", dsmz_raw$Ensembl_ID)
  if (any(duplicated(rownames(M)))) {
    M <- rowsum(M, group = rownames(M), reorder = TRUE)
  }
  keep <- colSums(!is.na(M)) > 0
  M <- M[, keep, drop = FALSE]
  storage.mode(M) <- "double"
  cat(sprintf("[INFO] Final DSMZ matrix: %d genes x %d samples\n", nrow(M), ncol(M)))
  M
}

align_counts_metadata <- function(meta_data, counts_data, col_name = "sample_name", outdir = NULL) {
  if (!col_name %in% colnames(meta_data)) {
    stop(sprintf("Column %s not found in metadata", col_name))
  }
  meta_data$sample_id <- meta_data[[col_name]]
  matched <- intersect(meta_data$sample_id, colnames(counts_data))
  md_only <- setdiff(meta_data$sample_id, colnames(counts_data))
  ct_only <- setdiff(colnames(counts_data), meta_data$sample_id)
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    if (length(md_only)) {
      utils::write.csv(data.frame(sample_name = md_only), file.path(outdir, "unmatched_metadata_samples.csv"), row.names = FALSE)
    }
    if (length(ct_only)) {
      utils::write.csv(data.frame(sample_name = ct_only), file.path(outdir, "unmatched_count_columns.csv"), row.names = FALSE)
    }
  }
  meta_data <- meta_data[meta_data$sample_id %in% matched, , drop = FALSE]
  meta_data <- meta_data[match(matched, meta_data$sample_id), , drop = FALSE]
  counts_data <- counts_data[, matched, drop = FALSE]
  stopifnot(identical(colnames(counts_data), meta_data$sample_id))
  list(meta = meta_data, counts = counts_data)
}
