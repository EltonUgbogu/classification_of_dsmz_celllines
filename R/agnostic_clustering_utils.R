#!/usr/bin/env Rscript
# agnostic_clustering_utils.R
# Shared helpers for 12-way agnostic clustering
# (PCA/raw × HC/k-means × cell/tumour/integrated)

suppressPackageStartupMessages({
  library(cluster)
  library(Matrix)
})

info <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")

# -------------------------------------------------------------------
# Load expression matrix from RDS (matrix or SummarizedExperiment-like)
# -------------------------------------------------------------------
load_expr_mat <- function(path) {
  if (is.null(path)) return(NULL)
  info("Loading expression from: %s", path)
  obj <- readRDS(path)

  if (is.matrix(obj) || inherits(obj, "Matrix")) {
    mat <- as.matrix(obj)
  } else if (inherits(obj, "SummarizedExperiment") || inherits(obj, "DESeqTransform")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment package is required but not installed.")
    }
    mat <- SummarizedExperiment::assay(obj)
  } else {
    stop("Unsupported object type in RDS: ", paste(class(obj), collapse = ", "))
  }

  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Expression matrix must have rownames (features) and colnames (samples).")
  }

  mat
}

# -------------------------------------------------------------------
# Build sample × gene matrix and dataset labels
# kind: pca_hc_cell, pca_hc_tumour, pca_hc_cell_tumour, etc.
# -------------------------------------------------------------------
build_sample_matrix <- function(kind,
                                cell_rds = NULL,
                                tumour_rds = NULL,
                                cell_label = "CELL",
                                tumour_label = "TUMOUR") {

  is_cell        <- grepl("_cell$", kind)
  is_tumour      <- grepl("_tumour$", kind)
  is_cell_tumour <- grepl("cell_tumour$", kind)

  if (is_cell_tumour) {
    if (is.null(cell_rds) || is.null(tumour_rds)) {
      stop("For *_cell_tumour kinds, both --cell_rds and --tumour_rds are required.")
    }

    cell_gx   <- load_expr_mat(cell_rds)   # genes × samples
    tumour_gx <- load_expr_mat(tumour_rds)

    common_genes <- intersect(rownames(cell_gx), rownames(tumour_gx))
    if (!length(common_genes)) {
      stop("No common genes between cell and tumour matrices.")
    }
    info("Common genes (CELL ∩ TUMOUR): %d", length(common_genes))

    cell_gx   <- cell_gx[common_genes, , drop = FALSE]
    tumour_gx <- tumour_gx[common_genes, , drop = FALSE]

    cell_mat   <- t(cell_gx)   # samples × genes
    tumour_mat <- t(tumour_gx)

    X <- rbind(cell_mat, tumour_mat)
    dataset <- factor(
      c(rep(cell_label,   nrow(cell_mat)),
        rep(tumour_label, nrow(tumour_mat))),
      levels = c(cell_label, tumour_label)
    )

    # Make sample IDs explicit
    rownames(X) <- c(
      paste0(cell_label,   ":", rownames(cell_mat)),
      paste0(tumour_label, ":", rownames(tumour_mat))
    )
  } else if (is_cell) {
    if (is.null(cell_rds)) stop("cell_rds is required for *_cell kinds.")
    cell_gx <- load_expr_mat(cell_rds)       # genes × samples
    X <- t(cell_gx)                          # samples × genes
    dataset <- factor(rep(cell_label, nrow(X)), levels = cell_label)
    rownames(X) <- ifelse(
      is.na(rownames(X)), colnames(cell_gx), rownames(X)
    )
  } else if (is_tumour) {
    if (is.null(tumour_rds)) stop("tumour_rds is required for *_tumour kinds.")
    tumour_gx <- load_expr_mat(tumour_rds)   # genes × samples
    X <- t(tumour_gx)                        # samples × genes
    dataset <- factor(rep(tumour_label, nrow(X)), levels = tumour_label)
    rownames(X) <- ifelse(
      is.na(rownames(X)), colnames(tumour_gx), rownames(X)
    )
  } else {
    stop("Cannot infer input type from kind: ", kind)
  }

  if (any(is.na(rownames(X)))) {
    stop("Sample matrix has NA rownames after processing.")
  }

  list(mat = X, dataset = dataset)
}

# -------------------------------------------------------------------
# Filter zero-variance genes (baseline QC for all methods)
# -------------------------------------------------------------------
filter_zero_var_genes <- function(X, min_sd = 0) {
  sds <- apply(X, 2, sd)
  keep <- which(is.finite(sds) & sds > min_sd)
  if (length(keep) < 5) stop("Too few genes after variance filtering: ", length(keep))
  Xf <- X[, keep, drop = FALSE]
  info("Variance filter: kept %d / %d genes (min_sd=%s)", ncol(Xf), ncol(X), min_sd)
  Xf
}

# -------------------------------------------------------------------
# PCA helper
# -------------------------------------------------------------------
compute_pca_scores <- function(X, n_pcs = 20, scale_features = FALSE) {
  info("Running PCA on matrix: %d samples × %d genes", nrow(X), ncol(X))

  # X is samples x genes
  if (any(!is.finite(X))) {
    stop("Non-finite values in matrix (NA/Inf). Clean before PCA.")
  }

  # Drop zero-variance columns (genes)
  sds <- apply(X, 2, sd)
  keep <- which(sds > 0 & !is.na(sds))
  if (length(keep) < 5) stop("Too few non-constant genes after filtering: ", length(keep))

  Xf <- X[, keep, drop = FALSE]
  info("After filtering zero-variance genes: %d samples × %d genes", nrow(Xf), ncol(Xf))

  # For VST, scaling is usually unnecessary; center is enough
  p <- prcomp(Xf, center = TRUE, scale. = scale_features)
  npc <- min(n_pcs, ncol(p$x))
  info("Using %d PCs (requested %d).", npc, n_pcs)
  list(pca = p, scores = p$x[, seq_len(npc), drop = FALSE])
}

# -------------------------------------------------------------------
# Hierarchical clustering with silhouette-based K selection
# -------------------------------------------------------------------
hc_optimal <- function(X,
                       max_k = 8,
                       dist_method = "euclidean",
                       seed = 42) {
  set.seed(seed)

  if (dist_method == "correlation") {
    info("Distance: 1 - Pearson correlation")
    cm <- cor(t(X), method = "pearson", use = "pairwise.complete.obs")
    cm[!is.finite(cm)] <- 0
    d <- as.dist(1 - cm)
  } else {
    info("Distance: %s", dist_method)
    d <- dist(X, method = dist_method)
  }

  hc <- hclust(d, method = "ward.D2")

  best_k <- NA_integer_
  best_sil <- -Inf
  best_clusters <- NULL

  for (k in 2:max_k) {
    cl <- cutree(hc, k = k)
    sil <- cluster::silhouette(cl, d)
    mean_sil <- mean(sil[, "sil_width"])
    info("HC k=%d -> mean silhouette = %.4f", k, mean_sil)

    if (mean_sil > best_sil ||
        (isTRUE(all.equal(mean_sil, best_sil)) && (is.na(best_k) || k < best_k))) {
      best_sil <- mean_sil
      best_k <- k
      best_clusters <- cl
    }
  }

  info("HC best k=%d (mean silhouette=%.4f)", best_k, best_sil)

  list(
    hc              = hc,
    dist            = d,
    best_k          = best_k,
    best_silhouette = best_sil,
    clusters        = best_clusters
  )
}

# -------------------------------------------------------------------
# k-means clustering with silhouette-based K selection
# -------------------------------------------------------------------
kmeans_optimal <- function(X,
                           k_min = 2,
                           k_max = 8,
                           seed  = 42) {
  set.seed(seed)
  d <- dist(X)  # Euclidean

  best_k <- NA_integer_
  best_sil <- -Inf
  best_clusters <- NULL
  best_km <- NULL

  for (k in k_min:k_max) {
    km <- stats::kmeans(X, centers = k, nstart = 20)
    sil <- cluster::silhouette(km$cluster, d)
    mean_sil <- mean(sil[, "sil_width"])
    info("k-means k=%d -> mean silhouette = %.4f", k, mean_sil)

    if (mean_sil > best_sil ||
        (isTRUE(all.equal(mean_sil, best_sil)) && (is.na(best_k) || k < best_k))) {
      best_sil <- mean_sil
      best_k <- k
      best_clusters <- km$cluster
      best_km <- km
    }
  }

  info("k-means best k=%d (mean silhouette=%.4f)", best_k, best_sil)

  list(
    kmeans          = best_km,
    dist            = d,
    best_k          = best_k,
    best_silhouette = best_sil,
    clusters        = best_clusters
  )
}

# -------------------------------------------------------------------
# Main driver used by all 12 wrappers
# -------------------------------------------------------------------
run_agnostic_clustering <- function(kind,
                                    cell_rds         = NULL,
                                    tumour_rds       = NULL,
                                    outdir,
                                    n_pcs            = 20,
                                    max_k_hc         = 8,
                                    k_min            = 2,
                                    k_max            = 8,
                                    dist_method      = "euclidean",
                                    seed             = 42,
                                    cluster_rds_path,
                                    cell_label       = "CELL",
                                    tumour_label     = "TUMOUR") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  info("==== Running agnostic clustering: %s ====", kind)

  # 1) Build sample matrix
  built <- build_sample_matrix(
    kind        = kind,
    cell_rds    = cell_rds,
    tumour_rds  = tumour_rds,
    cell_label  = cell_label,
    tumour_label = tumour_label
  )
  X       <- built$mat
  dataset <- built$dataset

  # 1b) Minimal QC baseline: drop constant genes (track-3 definition)
  X <- filter_zero_var_genes(X, min_sd = 0)

  # 2) Which method?
  method_type <- if (grepl("^pca_hc", kind)) {
    "pca_hc"
  } else if (grepl("^pca_kmeans", kind)) {
    "pca_kmeans"
  } else if (grepl("^hc_", kind)) {
    "hc"
  } else if (grepl("^kmeans_", kind)) {
    "kmeans"
  } else {
    stop("Unknown 'kind' prefix for method type: ", kind)
  }

  # 3) PCA or raw?
  if (method_type %in% c("pca_hc", "pca_kmeans")) {
    pca_res <- compute_pca_scores(X, n_pcs = n_pcs)
    X_use   <- pca_res$scores
  } else {
    pca_res <- NULL
    X_use   <- X
  }

  # 4) Cluster
  if (method_type %in% c("hc", "pca_hc")) {
    cl_res <- hc_optimal(
      X_use,
      max_k       = max_k_hc,
      dist_method = dist_method,
      seed        = seed
    )
    best_k   <- cl_res$best_k
    clusters <- cl_res$clusters
    hc_res   <- cl_res
    km_res   <- NULL
  } else {
    cl_res <- kmeans_optimal(
      X_use,
      k_min = k_min,
      k_max = k_max,
      seed  = seed
    )
    best_k   <- cl_res$best_k
    clusters <- cl_res$clusters
    km_res   <- cl_res
    hc_res   <- NULL
  }

  clusters <- as.integer(clusters)
  names(clusters) <- rownames(X_use)
  clusters_labeled <- setNames(paste0("C", clusters), names(clusters))

  input_type <- if (grepl("cell_tumour$", kind)) {
    "CELL+TUMOUR"
  } else if (grepl("_cell$", kind)) {
    "CELL"
  } else {
    "TUMOUR"
  }

  res_obj <- list(
    kind             = kind,
    method_type      = method_type,
    input_type       = input_type,
    k                = best_k,
    clusters         = clusters,
    clusters_labeled = clusters_labeled,
    dataset          = dataset,
    pca              = pca_res,
    hc               = hc_res,
    kmeans           = km_res,
    params = list(
      n_pcs      = n_pcs,
      max_k_hc   = max_k_hc,
      k_min      = k_min,
      k_max      = k_max,
      dist       = dist_method,
      seed       = seed
    ),
    timestamp = Sys.time()
  )

  info("Saving cluster object to: %s", cluster_rds_path)
  saveRDS(res_obj, cluster_rds_path)
  invisible(res_obj)
}
