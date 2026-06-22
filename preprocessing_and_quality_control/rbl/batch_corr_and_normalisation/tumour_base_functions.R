harmonize_gene_ids <- function(tumour_counts) {
  tumour_counts <- as.matrix(tumour_counts)
  rownames(tumour_counts) <- sub("\\..*$", "", rownames(tumour_counts))
  if (any(duplicated(rownames(tumour_counts)))) {
    tumour_counts <- rowsum(tumour_counts, rownames(tumour_counts), reorder = TRUE)
  }
  storage.mode(tumour_counts) <- "double"
  tumour_counts
}

verify_purity_data <- function(tumour_se, purity_file = NULL) {
  p <- NULL
  if (methods::is(tumour_se, "SummarizedExperiment") && "purity" %in% colnames(SummarizedExperiment::colData(tumour_se))) {
    p <- as.numeric(SummarizedExperiment::colData(tumour_se)$purity)
    names(p) <- colnames(SummarizedExperiment::assay(tumour_se))
  } else if (!is.null(purity_file) && file.exists(purity_file)) {
    tab <- utils::read.delim(purity_file, stringsAsFactors = FALSE)
    stopifnot(all(c("sample", "purity") %in% colnames(tab)))
    p <- setNames(tab$purity, tab$sample)
  }
  p
}

purity_adjust_log <- function(tumour_norm, purity = NULL, outdir = NULL) {
  if (is.null(purity)) return(list(V = tumour_norm, dropped = character(0)))
  keep_samp <- intersect(colnames(tumour_norm), names(purity)[!is.na(purity)])
  M <- tumour_norm[, keep_samp, drop = FALSE]
  pu <- purity[keep_samp]
  if (ncol(M) < 3) {
    return(list(V = tumour_norm, dropped = character(0)))
  }
  ct <- apply(M, 1, function(v) suppressWarnings(stats::cor.test(as.numeric(v), pu, method = "spearman")))
  pv <- vapply(ct, `[[`, numeric(1), "p.value")
  rh <- vapply(ct, function(x) unname(x$estimate), numeric(1))
  padj <- stats::p.adjust(pv, "BH")
  drop <- names(which(padj < 0.01 & rh < -0.4))
  out <- tumour_norm
  if (length(drop) > 0) {
    common_keep <- setdiff(rownames(M), drop)
    out <- out[common_keep, , drop = FALSE]
    if (!is.null(outdir)) {
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
      writeLines(drop, file.path(outdir, "purity_negatively_associated_genes.txt"))
    }
  }
  list(V = out, dropped = drop)
}
