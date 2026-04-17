choose_minPts <- function(n_obs) {
  val <- as.integer(round(log10(max(10, n_obs)) * 5))
  min(50L, max(5L, val))
}

choose_minCluster <- function(n, frac = 0.03, min_floor = 8, max_cap = 150) {
  p <- max(min_floor, round(frac * n))
  as.integer(min(p, max_cap))
}

top_var_genes <- function(mat, n) {
  iqr <- matrixStats::rowIQRs(mat)
  names(sort(iqr, decreasing = TRUE))[seq_len(min(n, length(iqr)))]
}

make_pcs <- function(M, n_hvg = 3000, max_pc = 30) {
  v <- matrixStats::rowVars(M)
  keep <- which(is.finite(v) & v > 0)
  if (length(keep) < 2) stop("[FATAL] Too few variable genes for PCA.")
  M <- M[keep, , drop = FALSE]
  hvgs_all <- names(sort(matrixStats::rowIQRs(M), decreasing = TRUE))
  hvgs <- hvgs_all[seq_len(min(n_hvg, length(hvgs_all)))]
  X <- t(scale(M[hvgs, , drop = FALSE]))
  X <- X[, colSums(!is.na(t(X))) > 0, drop = FALSE]
  if (nrow(X) < 2 || ncol(X) < 1) stop("[FATAL] PCA input has insufficient dimension.")
  pr <- prcomp(X, center = FALSE, scale. = FALSE)
  if (is.null(pr$x) || ncol(pr$x) < 1) stop("[FATAL] PCA returned 0 components.")
  kPC <- as.integer(min(max_pc, ncol(pr$x)))
  list(PC = pr$x[, 1:kPC, drop = FALSE], kPC = kPC, pr = pr, hvgs = hvgs)
}

make_pcs_matrix <- function(V_adj, n_hvg = 3000, max_pc = 30) {
  if (!is.matrix(V_adj) && !is.data.frame(V_adj)) {
    stop("V_adj must be a matrix or data.frame")
  }
  mat <- as.matrix(V_adj)
  if (ncol(mat) < 2 || nrow(mat) < 2) stop("V_adj must have at least 2 rows and 2 columns")
  pcs_obj <- make_pcs(mat, n_hvg = n_hvg, max_pc = max_pc)
  pcs_obj$PC
}

log_dims <- function(tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  cat(sprintf("[DIM] %-20s: %8s genes x %6s samples\n", tag, format(ng, big.mark = ","), format(ns, big.mark = ",")))
}

write_dims_line <- function(con, tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  writeLines(sprintf("%-24s\tgenes=%d\tsamples=%d", tag, ng, ns), con)
}
