# ─────────────────────────────────────────────────────────────────────────────
# Normalization helpers – VST & ComBat-Seq (using Exp_mat consistently)
# Place in R/normalization.R
# ─────────────────────────────────────────────────────────────────────────────

#' @title VST-normalize a gene-expression matrix
#'
#' @description
#' Applies DESeq2's variance-stabilizing transformation (blind by default).
#' Input and output are **matrices** (genes × samples).
#'
#' @param Exp_mat Integer count matrix (genes in rows, samples in columns).
#' @param coldata Data frame with sample metadata (rownames = sample names).
#' @param blind Logical. Use blind VST? (default: `TRUE` – recommended for QC).
#'
#' @return Numeric matrix of VST values (same dimensions as `Exp_mat`).
#' @export
vst_normalize <- function(Exp_mat, coldata, blind = TRUE) {
  cat("[INFO] VST normalizing Exp_mat...\n")

  # ── Input validation ───────────────────────────────────────────────────────
  if (!is.matrix(Exp_mat) && !is.data.frame(Exp_mat))
    stop("`Exp_mat` must be a matrix or data.frame")
  Exp_mat <- as.matrix(round(Exp_mat))
  if (!is.numeric(Exp_mat))
    stop("`Exp_mat` must contain numeric/integer values")

  if (nrow(coldata) != ncol(Exp_mat) ||
      !identical(rownames(coldata), colnames(Exp_mat)))
    stop("`coldata` must have rownames matching colnames(Exp_mat) and same length")

  # ── DESeq2 VST ────────────────────────────────────────────────────────────
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = Exp_mat,
    colData   = coldata,
    design    = ~1
  )
  vsd <- DESeq2::vst(dds, blind = blind)

  SummarizedExperiment::assay(vsd)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Batch-correct counts using ComBat-Seq
#'
#' @description
#' Removes batch effects from **raw count data** using `sva::ComBat_seq()`.
#' Designed for use after `merge_counts()` on the merged `Xc_raw` matrix.
#'
#' @param Exp_mat Integer count matrix (genes × samples).
#' @param batch_factor Factor of batch labels (length = ncol(Exp_mat)).
#' @param seed Random seed for reproducibility (default: 42).
#'
#' @return Batch-corrected count matrix (same dimensions as `Exp_mat`).
#' @export
batch_adjust <- function(Exp_mat, batch_factor, seed = 42) {
  cat("[INFO] Batch adjusting Exp_mat with ComBat-Seq...\n")

  # ── Input validation ───────────────────────────────────────────────────────
  if (!is.matrix(Exp_mat) && !is.data.frame(Exp_mat))
    stop("`Exp_mat` must be a matrix or data.frame")
  Exp_mat <- as.matrix(Exp_mat)
  if (!is.numeric(Exp_mat))
    stop("`Exp_mat` must contain numeric values")

  if (length(batch_factor) != ncol(Exp_mat))
    stop("`batch_factor` length must equal ncol(Exp_mat)")
  batch_factor <- as.factor(batch_factor)

  # ── ComBat-Seq ────────────────────────────────────────────────────────────
  set.seed(seed)
  corrected <- sva::ComBat_seq(
    counts       = Exp_mat,
    batch        = batch_factor,
    group        = NULL,
    covar_mod    = NULL,
    full_mod     = TRUE,
    shrink       = FALSE,
    shrink.disp  = FALSE
  )

  return(corrected)
}