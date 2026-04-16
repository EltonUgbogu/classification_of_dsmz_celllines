vst_normalize <- function(Exp_mat, coldata, blind = TRUE) {
  cat("[INFO] VST normalizing matrix...\n")
  if (!is.matrix(Exp_mat) && !is.data.frame(Exp_mat)) stop("Exp_mat must be a matrix or data.frame")
  Exp_mat <- as.matrix(round(Exp_mat))
  if (!is.numeric(Exp_mat)) stop("Exp_mat must contain numeric values")
  if (any(Exp_mat < 0, na.rm = TRUE)) stop("Exp_mat contains negative values")
  if (nrow(coldata) != ncol(Exp_mat) || !identical(rownames(coldata), colnames(Exp_mat))) {
    stop("coldata rownames must match colnames(Exp_mat)")
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = Exp_mat, colData = coldata, design = ~1)
  vsd <- DESeq2::vst(dds, blind = blind)
  SummarizedExperiment::assay(vsd)
}

batch_adjust <- function(Exp_mat, batch_factor, seed = 42) {
  cat("[INFO] Batch adjusting matrix with ComBat-Seq...\n")
  if (!is.matrix(Exp_mat) && !is.data.frame(Exp_mat)) stop("Exp_mat must be a matrix or data.frame")
  counts <- as.matrix(round(Exp_mat))
  counts[counts < 0] <- 0
  storage.mode(counts) <- "integer"
  if (length(batch_factor) != ncol(counts)) stop("batch_factor length must equal ncol(Exp_mat)")
  batch_factor <- as.factor(batch_factor)
  set.seed(seed)
  sva::ComBat_seq(counts = counts, batch = batch_factor, group = NULL, covar_mod = NULL, full_mod = TRUE, shrink = FALSE, shrink.disp = FALSE)
}
