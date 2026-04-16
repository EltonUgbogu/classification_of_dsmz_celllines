common_genes <- function(nbl_tumour_counts, dsmz_counts, min_shared = NULL) {
  cat("[INFO] Finding common genes...\n")
  common <- intersect(rownames(nbl_tumour_counts), rownames(dsmz_counts))
  cat(sprintf("[INFO] Shared genes: %d\n", length(common)))
  if (!is.null(min_shared) && length(common) < min_shared) {
    stop(sprintf("Only %d shared genes found; expected at least %d", length(common), min_shared))
  }
  common <- sort(common)
  list(
    nbl_tumour_counts = nbl_tumour_counts[common, , drop = FALSE],
    dsmz_counts = dsmz_counts[common, , drop = FALSE],
    genes = common
  )
}

merge_counts <- function(nbl_tumour_counts, dsmz_counts, dataset_names = c("NBL_TUMOUR", "DSMZ_NBL"), min_shared = NULL) {
  if (length(dataset_names) < 2) {
    stop("dataset_names must have length 2")
  }
  common <- common_genes(nbl_tumour_counts, dsmz_counts, min_shared = min_shared)
  Xc_raw <- cbind(common$nbl_tumour_counts, common$dsmz_counts)
  storage.mode(Xc_raw) <- "double"

  batch <- factor(
    c(rep(dataset_names[1], ncol(common$nbl_tumour_counts)), rep(dataset_names[2], ncol(common$dsmz_counts))),
    levels = dataset_names[1:2]
  )
  names(batch) <- colnames(Xc_raw)

  sample_id <- colnames(Xc_raw)
  names(sample_id) <- sample_id

  list(
    Xc_raw = Xc_raw,
    batch = batch,
    sample_id = sample_id,
    common_genes = common$genes,
    tumour_counts = common$nbl_tumour_counts,
    dsmz_counts = common$dsmz_counts
  )
}
