# Helpers: low-level utility functions
# Place in R/helpers.R

# Library imports for helper code (not necessary to call library() here; main script should load them)
strip_ensver <- function(x) sub("\\..*$","", x)
as_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

#' @brief Calculates minPts for clustering based on dataset size.
#' @param n Integer, total samples.
#' @param frac Numeric, fraction of samples for minPts (default: 0.02).
#' @param min_floor Integer, minimum minPts (default: 5).
#' @param max_cap Integer, maximum minPts (default: 50).
#' @return Integer, calculated minPts.

choose_minPts <- function(n, frac = 0.02, min_floor = 5, max_cap = 50) {
  p <- max(min_floor, round(frac * n)) #' @brief Sets minPts as max of min_floor and rounded fraction of n.
  p <- min(p, max_cap, n - 1L) #' @brief Caps minPts at max_cap and ensures less than n.
  return(as.integer(p)) #' @brief Returns minPts as integer.
}



choose_minCluster <- function(n, frac = 0.03, min_floor = 8, max_cap = 150) {
  p <- max(min_floor, round(frac * n))
  as.integer(min(p, max_cap))
}

make_ens2sym <- function(dsmz_raw) {
  df <- as_df(dsmz_raw)
  ens2sym <- df %>%
    transmute(
      ensembl = strip_ensver(Ensembl_ID),
      symbol  = as.character(gene_name)
    ) %>%
    filter(!is.na(ensembl) & nzchar(ensembl)) %>%
    arrange(ensembl, dplyr::desc(nchar(symbol)), symbol) %>%
    distinct(ensembl, .keep_all = TRUE) %>%
    mutate(symbol_unique = make.unique(ifelse(is.na(symbol) | symbol == "", ensembl, symbol)))
  ens2sym <- as.data.frame(ens2sym, stringsAsFactors = FALSE)
  rownames(ens2sym) <- ens2sym$ensembl
  ens2sym
}

top_var_genes <- function(mat, n = 3000) {
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

# alias to allow wrappers in other files without naming collisions
helpers_make_pcs <- function(M, n_hvg = 3000, max_pc = 30) make_pcs(M, n_hvg = n_hvg, max_pc = max_pc)

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

zscore_by_gene <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- matrixStats::rowSds(mat, na.rm = TRUE); sdv[sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}

safe_pdf <- function(path, expr) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  pdf(path)
  on.exit({
    if (dev.cur() != 1) dev.off()
  }, add = TRUE)
  force(expr)
}

build_dsmz_matrix <- function(dsmz_raw) {
  stopifnot(all(c("Ensembl_ID","gene_name") %in% names(dsmz_raw)))
  annot <- c("Ensembl_ID","gene_name","Ensembl_ID_with_version")
  sample_cols <- setdiff(colnames(dsmz_raw), annot)
  dsmz_raw[sample_cols] <- lapply(dsmz_raw[sample_cols], function(v) as.numeric(as.character(v)))
  M <- as.matrix(dsmz_raw[, sample_cols, drop = FALSE])
  rownames(M) <- strip_ensver(dsmz_raw$Ensembl_ID)
  if (any(duplicated(rownames(M)))) {
    M <- rowsum(M, rownames(M), reorder = TRUE)
  }
  M
}

logCPM <- function(x, prior.count = 1, lib.size = NULL) {
  if (!is.matrix(x)) x <- as.matrix(x)
  edgeR::cpm(x, log = TRUE, prior.count = prior.count, lib.size = lib.size)
}

# Map symbols -> ensembl lookup (expects ens2sym_union data.frame with ensembl,symbol)
make_sym2ens <- function(ens2sym_union) {
  setNames(ens2sym_union$ensembl, toupper(ens2sym_union$symbol))
}

# relabel helper for HDBSCAN outputs
relabel_hdb_numeric <- function(x, keep_unassigned = TRUE) {
  x <- as.character(x)
  unass <- x == "Unassigned"
  labs  <- sort(unique(x[!unass]))
  map   <- setNames(as.character(seq_along(labs)), labs)
  y     <- ifelse(unass & keep_unassigned, NA, map[x])
  factor(y, levels = as.character(seq_len(min(5, length(labs)))))
}

aggregate_median <- function(M, groups) {
  grp <- droplevels(factor(groups))
  lev <- levels(grp)
  out <- sapply(lev, function(g) {
    rowMedians(M[, grp == g, drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  rownames(out) <- rownames(M)
  colnames(out) <- lev
  out
}

# centroid helpers
get_centroid <- function(X, idx) colMeans(X[idx, , drop = FALSE], na.rm = TRUE)
get_k_nearest <- function(X, centroid, k = 20) {
  d <- sqrt(rowSums((X - matrix(centroid, nrow = nrow(X), ncol = ncol(X), byrow = TRUE))^2))
  ord <- order(d, decreasing = FALSE)
  keep <- ord[seq_len(min(k, length(ord)))]
  data.frame(sample = rownames(X)[keep], dist = d[keep], stringsAsFactors = FALSE)
}