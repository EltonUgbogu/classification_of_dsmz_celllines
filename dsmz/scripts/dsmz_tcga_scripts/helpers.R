# Helpers: low-level utility functions
# Place in R/helpers.R

# Library imports for helper code (not necessary to call library() here; main script should load them)
strip_ensver <- function(x) sub("\\..*$","", x)
as_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

# Simple logging helpers
log_info  <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")
log_warn  <- function(...) cat("[WARN] ", sprintf(...), "\n", sep = "")
log_error <- function(...) cat("[ERROR]", sprintf(...), "\n", sep = "")

#' @brief Calculates minPts for clustering based on dataset size.
#' @param n Integer, total samples.
#' @param frac Numeric, fraction of samples for minPts (default: 0.02).
#' @param min_floor Integer, minimum minPts (default: 5).
#' @param max_cap Integer, maximum minPts (default: 50).
#' @return Integer, calculated minPts.
# Determines the minPts parameter for HDBSCAN based on sample size
choose_minPts <- function(n_obs) {
  # Apply heuristic: minPts is 5 times the log10 of sample size (min 10)
  val <- as.integer(round(log10(max(10, n_obs)) * 5))
  # Clamp the value to ensure minPts is between 5 and 50 (inclusive)
  min(50L, max(5L, val))
}


# Determines minCluster size for dynamic tree cut with flexible parameters
choose_minCluster <- function(n, frac = 0.03, min_floor = 8, max_cap = 150) {
  # Calculate minCluster as fraction of sample size, respecting min_floor
  p <- max(min_floor, round(frac * n))
  # Return minCluster, capped at max_cap
  as.integer(min(p, max_cap))
}



# mat - matrix of data
# n - number of top variable genes to return
top_var_genes <- function(mat, n) {
  iqr <- matrixStats::rowIQRs(mat)
  names(sort(iqr, decreasing = TRUE))[seq_len(min(n, length(iqr)))]
}

make_pcs <- function(M, n_hvg = 3000, max_pc = 30) {
  v <- matrixStats::rowVars(M); keep <- which(is.finite(v) & v > 0)
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
  list(PC = pr$x[, 1:kPC, drop = FALSE],
       kPC = kPC,
       pr = pr,
       hvgs = hvgs)
}


# convenience: return only the PC matrix (used by I/O and clustering code)
make_pcs_matrix <- function(V_adj, n_hvg = 3000, max_pc = 30, center = TRUE, scale. = TRUE) {
  if (!is.matrix(V_adj) && !is.data.frame(V_adj)) {
    stop("V_adj must be a matrix or data.frame")
  }
  mat <- as.matrix(V_adj)
  if (ncol(mat) < 2 || nrow(mat) < 2) stop("V_adj must have at least 2 rows and 2 columns")
  pcs_obj <- helpers_make_pcs(mat, n_hvg = n_hvg, max_pc = max_pc)
  pcs_obj$PC
}



make_umap <- function(PC, seed = 42, n_neighbors = 20, min_dist = 0.3,
                     metric = "cosine") {
  set.seed(seed)
  emb <- uwot::umap(PC, n_neighbors = n_neighbors, min_dist = min_dist,
                    metric = metric, verbose = FALSE)
  emb <- as.data.frame(emb)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$sample   <- rownames(PC)
  emb
}



safe_rds_save <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
}

log_dims <- function(tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  cat(sprintf("[DIM] %-20s: %8s genes x %6s samples\n", tag, format(ng, big.mark=","), format(ns, big.mark=",")))
}

write_dims_line <- function(con, tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  writeLines(sprintf("%-24s\tgenes=%d\tsamples=%d", tag, ng, ns), con)
}


safe_rds_save <- function(object, file) {
  ensure_dir(dirname(file))
  saveRDS(object, file = file)
}

