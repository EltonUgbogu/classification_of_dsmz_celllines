#!/usr/bin/env Rscript
# =============================================================================
# hclust_kmeans_utils.R
#
# This module provides the shared statistical infrastructure underpinning the
# 12-way HC/k-means clustering stage. The 12 configurations arise from the
# full factorial combination of three design axes:
#
#   • Feature space  : PCA-reduced scores  ×  raw expression values
#   • Algorithm      : Hierarchical clustering (HC)  ×  k-means
#   • Sample scope   : Cell lines only  ×  Tumours only  ×  Integrated
#
# Each configuration exposes the same statistical pipeline — data ingestion,
# variance filtering, optional dimensionality reduction, distance computation,
# iterative partitioning, and silhouette-based model selection — through a
# single entry-point function, run_hclust_kmeans(), parameterised by
# a 'kind' string that encodes the desired combination.
#
# Statistical rationale:
#   Silhouette width is used as the universal model-selection criterion across
#   both HC and k-means. For a sample i, its silhouette s(i) is defined as:
#
#       s(i) = [ b(i) − a(i) ] / max{ a(i), b(i) }
#
#   where a(i) is the mean intra-cluster dissimilarity and b(i) is the mean
#   dissimilarity to the nearest neighbouring cluster. s(i) ∈ [−1, 1], with
#   values approaching 1 indicating well-separated, compact clusters. The mean
#   silhouette width across all samples provides a scalar summary of partition
#   quality at a given k, enabling a principled sweep over the candidate range.
# =============================================================================

suppressPackageStartupMessages({
  library(cluster)   # Provides silhouette() for partition quality evaluation.
  library(Matrix)    # Supports sparse matrix representations of expression data.
})

# Lightweight logging helper; prefixes all log messages with [INFO].
info <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")


# =============================================================================
# Section 1 — Data Ingestion
# =============================================================================

# -----------------------------------------------------------------------------
# load_expr_mat()
#
# This function deserialises a pre-processed gene expression matrix from an
# RDS file and returns a numeric matrix in genes × samples orientation.
# Three input formats are accepted to accommodate different upstream pipeline
# outputs: (i) a base R matrix, (ii) a sparse Matrix object (appropriate for
# very high-dimensional data), and (iii) a SummarizedExperiment or
# DESeqTransform object, from which the primary assay slot is extracted.
#
# The function enforces the presence of both rownames (features) and colnames
# (samples), as downstream distance and PCA computations require named
# dimensions for unambiguous sample tracking.
# -----------------------------------------------------------------------------
load_expr_mat <- function(path) {
  if (is.null(path)) return(NULL)
  info("Loading expression from: %s", path)
  obj <- readRDS(path)

  if (is.matrix(obj) || inherits(obj, "Matrix")) {
    # Base R matrix or sparse Matrix: coerce directly to dense numeric matrix.
    mat <- as.matrix(obj)

  } else if (inherits(obj, "SummarizedExperiment") || inherits(obj, "DESeqTransform")) {
    # SummarizedExperiment-derived objects: extract the primary assay, which
    # for DESeqTransform objects typically contains VST-stabilised log2 counts.
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment package is required but not installed.")
    }
    mat <- SummarizedExperiment::assay(obj)

  } else {
    stop("Unsupported object type in RDS: ", paste(class(obj), collapse = ", "))
  }

  # Named dimensions are required for all downstream statistical operations.
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Expression matrix must have rownames (features) and colnames (samples).")
  }

  mat
}


# =============================================================================
# Section 2 — Sample Matrix Construction
# =============================================================================

# -----------------------------------------------------------------------------
# build_sample_matrix()
#
# This function constructs the sample × gene input matrix appropriate to the
# requested clustering configuration, as encoded in the 'kind' string.
#
# Statistical considerations:
#   • When operating in integrated (cell_tumour) mode, the feature space is
#     limited to the intersection of genes present in both cohorts. This
#     ensures that all pairwise distances and PC loadings are computed on a
#     consistent, jointly observed feature set, avoiding spurious dissimilarity
#     attributable to missing features rather than true biological divergence.
#   • Sample identifiers are prefixed with their cohort label (e.g., 'CELL:',
#     'TUMOUR:') to prevent namespace collisions in the merged matrix and to
#     enable downstream stratification of clustering results by data source.
#   • The returned dataset factor encodes cohort membership for each sample,
#     permitting downstream scripts to evaluate whether clusters recapitulate
#     or transcend the cell line / tumour boundary.
# -----------------------------------------------------------------------------
build_sample_matrix <- function(kind,
                                cell_rds     = NULL,
                                tumour_rds   = NULL,
                                cell_label   = "CELL",
                                tumour_label = "TUMOUR") {

  # Infer the sample scope from the suffix of the 'kind' string.
  is_cell        <- grepl("_cell$",        kind)
  is_tumour      <- grepl("_tumour$",      kind)
  is_cell_tumour <- grepl("cell_tumour$",  kind)

  if (is_cell_tumour) {
    # Integrated mode: both cohorts are loaded and their feature spaces aligned
    # by gene-set intersection before row-binding into a single sample matrix.
    if (is.null(cell_rds) || is.null(tumour_rds)) {
      stop("For *_cell_tumour kinds, both --cell_rds and --tumour_rds are required.")
    }

    cell_gx   <- load_expr_mat(cell_rds)    # genes × cell line samples
    tumour_gx <- load_expr_mat(tumour_rds)  # genes × tumour samples

    # Feature alignment: limit to genes measured in both cohorts.
    common_genes <- intersect(rownames(cell_gx), rownames(tumour_gx))
    if (!length(common_genes)) {
      stop("No common genes between cell and tumour matrices.")
    }
    info("Common genes (CELL ∩ TUMOUR): %d", length(common_genes))

    cell_gx   <- cell_gx[common_genes, , drop = FALSE]
    tumour_gx <- tumour_gx[common_genes, , drop = FALSE]

    # Transpose from genes × samples to the samples × genes orientation
    # required by distance, PCA, and clustering functions.
    cell_mat   <- t(cell_gx)
    tumour_mat <- t(tumour_gx)

    # Merge cohorts row-wise; the combined matrix spans all samples.
    X <- rbind(cell_mat, tumour_mat)

    # Encode cohort membership as an ordered factor for downstream use.
    dataset <- factor(
      c(rep(cell_label,   nrow(cell_mat)),
        rep(tumour_label, nrow(tumour_mat))),
      levels = c(cell_label, tumour_label)
    )

    # Prefix sample IDs to prevent naming conflicts across cohorts.
    rownames(X) <- c(
      paste0(cell_label,   ":", rownames(cell_mat)),
      paste0(tumour_label, ":", rownames(tumour_mat))
    )

  } else if (is_cell) {
    # Cell-line-only mode: the sample matrix is constructed from the cell
    # line expression data alone, transposed to samples × genes.
    if (is.null(cell_rds)) stop("cell_rds is required for *_cell kinds.")
    cell_gx <- load_expr_mat(cell_rds)
    X       <- t(cell_gx)
    dataset <- factor(rep(cell_label, nrow(X)), levels = cell_label)
    rownames(X) <- ifelse(is.na(rownames(X)), colnames(cell_gx), rownames(X))

  } else if (is_tumour) {
    # Tumour-only mode: analogous to cell-line mode but applied to the tumour
    # cohort, enabling separate investigation of patient sample substructure.
    if (is.null(tumour_rds)) stop("tumour_rds is required for *_tumour kinds.")
    tumour_gx <- load_expr_mat(tumour_rds)
    X         <- t(tumour_gx)
    dataset   <- factor(rep(tumour_label, nrow(X)), levels = tumour_label)
    rownames(X) <- ifelse(is.na(rownames(X)), colnames(tumour_gx), rownames(X))

  } else {
    stop("Cannot infer input type from kind: ", kind)
  }

  # A final integrity check ensures no unnamed samples enter the pipeline,
  # which would compromise the traceability of cluster assignments.
  if (any(is.na(rownames(X)))) {
    stop("Sample matrix has NA rownames after processing.")
  }

  list(mat = X, dataset = dataset)
}


# =============================================================================
# Section 3 — Variance Filtering
# =============================================================================

# -----------------------------------------------------------------------------
# filter_zero_var_genes()
#
# This function removes constant genes (standard deviation ≤ min_sd) from the
# sample matrix prior to dimensionality reduction or distance computation.
#
# Statistical rationale:
#   Constant genes carry no discriminatory information: their contribution to
#   Euclidean distances is zero, and their inclusion in PCA introduces
#   ill-conditioned covariance estimates. By default, min_sd = 0 removes only
#   exactly constant genes; a more stringent threshold may be applied to further
#   attenuate low-information features. A minimum of five genes is enforced
#   as a safeguard against degenerate downstream analyses.
# -----------------------------------------------------------------------------
filter_zero_var_genes <- function(X, min_sd = 0) {
  sds  <- apply(X, 2, sd)
  keep <- which(is.finite(sds) & sds > min_sd)
  if (length(keep) < 5) stop("Too few genes after variance filtering: ", length(keep))
  Xf <- X[, keep, drop = FALSE]
  info("Variance filter: kept %d / %d genes (min_sd=%s)", ncol(Xf), ncol(X), min_sd)
  Xf
}


# =============================================================================
# Section 4 — Dimensionality Reduction
# =============================================================================

# -----------------------------------------------------------------------------
# compute_pca_scores()
#
# This function applies Principal Component Analysis (PCA) to the filtered
# sample × gene matrix and returns the sample scores on the top n_pcs
# principal components.
#
# Statistical rationale:
#   PCA performs an orthogonal linear transformation that maps the p-dimensional
#   gene expression space into a lower-dimensional subspace defined by the
#   eigenvectors of the sample covariance matrix, ordered by descending
#   eigenvalue. Retaining the top n_pcs components (default: 20) concentrates
#   the biologically meaningful variance — which tends to be captured by the
#   leading PCs — while attenuating high-dimensional noise in the tail
#   components. Centering is applied unconditionally; additional feature
#   scaling (scale. = TRUE) is optional and typically unnecessary for
#   variance-stabilised transcript counts, where all features share a
#   comparable dynamic range.
#
#   The resulting n × n_pcs score matrix provides the Euclidean embedding
#   in which subsequent hierarchical agglomeration or k-means partitioning
#   is performed, ensuring that distance computations reflect dominant axes
#   of transcriptomic variation rather than idiosyncratic gene-level noise.
# -----------------------------------------------------------------------------
compute_pca_scores <- function(X, n_pcs = 20, scale_features = FALSE) {
  info("Running PCA on matrix: %d samples × %d genes", nrow(X), ncol(X))

  # Verify numerical integrity before decomposition.
  if (any(!is.finite(X))) {
    stop("Non-finite values in matrix (NA/Inf). Clean before PCA.")
  }

  # A secondary zero-variance check guards against edge cases arising from
  # subsetting in integrated mode that may introduce newly constant features.
  sds  <- apply(X, 2, sd)
  keep <- which(sds > 0 & !is.na(sds))
  if (length(keep) < 5) stop("Too few non-constant genes after filtering: ", length(keep))
  Xf <- X[, keep, drop = FALSE]
  info("After filtering zero-variance genes: %d samples × %d genes", nrow(Xf), ncol(Xf))

  # prcomp() uses a singular value decomposition (SVD) of the centred matrix,
  # which is numerically stable for tall (n >> p) and wide (p >> n) designs.
  p   <- prcomp(Xf, center = TRUE, scale. = scale_features)

  # Cap at the number of computed components to handle low-sample-count cases.
  npc <- min(n_pcs, ncol(p$x))
  info("Using %d PCs (requested %d).", npc, n_pcs)

  list(pca = p, scores = p$x[, seq_len(npc), drop = FALSE])
}


# =============================================================================
# Section 5 — Hierarchical Clustering with Silhouette-Based K Selection
# =============================================================================

# -----------------------------------------------------------------------------
# hc_optimal()
#
# This function performs agglomerative hierarchical clustering on an input
# matrix X (either PC scores or raw features) and selects the optimal number
# of clusters k by maximising the mean silhouette width over k ∈ {2, …, max_k}.
#
# Statistical rationale:
#   1. Distance computation
#      By default, Euclidean distances are computed between sample vectors.
#      When dist_method = "correlation", the dissimilarity is defined as
#      1 − ρ_ij, where ρ_ij is the Pearson correlation between samples i and j.
#      This metric is invariant to linear rescaling of expression profiles,
#      making it appropriate when relative co-expression patterns — rather than
#      absolute magnitude differences — are the primary signal of interest.
#      Non-finite correlation values (arising from zero-variance samples) are
#      set to zero before conversion, preventing propagation of undefined
#      distances into the linkage tree.
#
#   2. Linkage criterion
#      Ward's minimum variance criterion (ward.D2) is used to merge clusters.
#      At each step, the pair of clusters whose union minimises the total
#      within-cluster sum of squared deviations from the centroid is merged.
#      This promotes compact, roughly spherical clusters and is widely
#      regarded as the most appropriate linkage for transcriptomic data where
#      cluster homogeneity — rather than chain-like structure — is expected.
#
#   3. Partition and model selection
#      The dendrogram is cut at each candidate k using cutree(), yielding a
#      hard partition. The mean silhouette width of the resulting partition
#      is computed via cluster::silhouette() and compared across k values.
#      The k that maximises mean silhouette width is retained as the optimal
#      partition, reflecting the k at which clusters are simultaneously most
#      cohesive internally and most separated from their neighbours.
# -----------------------------------------------------------------------------
# -------------------------------------------------------------------
# correlation_chord_dist: Euclidean-compatible correlation dissimilarity
# -------------------------------------------------------------------
# Ward.D2 minimises variance in a Euclidean space and therefore requires a
# Euclidean-compatible dissimilarity. The chord distance
#
#     d(x, y) = sqrt(2 * (1 - r(x, y)))
#
# is the exact Euclidean distance between the centred, unit-norm forms of x
# and y, since ||u_x - u_y||^2 = 2 * (1 - r). It is consequently the
# correlation dissimilarity used for Ward.D2 hierarchical clustering here.
# The raw 1 - r dissimilarity is not a Euclidean metric and is not used with
# Ward. This matches the correlation-geometry transform (centre + unit norm +
# Euclidean distance) applied by the ConsensusClusterPlus branch, so both
# clustering branches operate on the same correlation geometry.
#
# Non-finite correlations (zero-variance profiles) are treated as r = 0.
#
# X: samples x features matrix. Returns a 'dist' object over the samples.
correlation_chord_dist <- function(X) {
  cm <- cor(t(X), method = "pearson", use = "pairwise.complete.obs")
  cm[!is.finite(cm)] <- 0
  as.dist(sqrt(pmax(2 * (1 - cm), 0)))
}

hc_optimal <- function(X,
                       max_k       = 8,
                       dist_method = "euclidean",
                       seed        = 42) {
  set.seed(seed)

  if (dist_method == "correlation") {
    info("Distance: chord sqrt(2 * (1 - Pearson r)) for Ward.D2 compatibility")
    d <- correlation_chord_dist(X)
  } else {
    info("Distance: %s", dist_method)
    d <- dist(X, method = dist_method)
  }

  # Agglomerative clustering under Ward's minimum variance criterion.
  hc <- hclust(d, method = "ward.D2")

  best_k        <- NA_integer_
  best_sil      <- -Inf
  best_clusters <- NULL

  for (k in 2:max_k) {
    cl       <- cutree(hc, k = k)                    # Hard partition at depth k.
    sil      <- cluster::silhouette(cl, d)            # Per-sample silhouette widths.
    mean_sil <- mean(sil[, "sil_width"])              # Scalar summary of partition quality.
    info("HC k=%d -> mean silhouette = %.4f", k, mean_sil)

    # Retain the partition with the highest mean silhouette width.
    if (mean_sil > best_sil) {
      best_sil      <- mean_sil
      best_k        <- k
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


# =============================================================================
# Section 6 — k-Means Clustering with Silhouette-Based K Selection
# =============================================================================

# -----------------------------------------------------------------------------
# kmeans_optimal()
#
# This function applies k-means clustering across a range of k values and
# selects the optimal partition by maximising mean silhouette width, providing
# a consistent model-selection criterion across both clustering algorithms.
#
# Statistical rationale:
#   1. Algorithm
#      k-means partitions n samples into k clusters by iteratively assigning
#      each sample to its nearest centroid (Euclidean distance) and recomputing
#      centroids as within-cluster means, minimising the total within-cluster
#      sum of squares (WCSS). nstart = 20 restarts the algorithm from 20
#      independent random initialisations, mitigating sensitivity to the
#      initial centroid placement and improving the probability of converging
#      to the global WCSS minimum.
#
#   2. Distance for silhouette
#      Although k-means itself uses squared Euclidean distances internally,
#      silhouette widths are computed from the standard (non-squared) Euclidean
#      distance matrix derived from the full feature space, ensuring consistency
#      with the HC evaluation and interpretability of the s(i) ∈ [−1, 1] scale.
#
#   3. Model selection
#      As in hc_optimal(), the k ∈ {k_min, …, k_max} that maximises mean
#      silhouette width is selected. Parallel evaluation of HC and k-means
#      under the same criterion permits direct comparison of the two algorithms'
#      ability to recover stable cluster structure in a given feature space.
# -----------------------------------------------------------------------------
kmeans_optimal <- function(X,
                           k_min = 2,
                           k_max = 8,
                           seed  = 42) {
  set.seed(seed)

  # Euclidean distance matrix for silhouette computation, computed once and
  # reused across all k to avoid redundant O(n²) calculations.
  d <- dist(X)

  best_k        <- NA_integer_
  best_sil      <- -Inf
  best_clusters <- NULL
  best_km       <- NULL

  for (k in k_min:k_max) {
    # nstart = 20 multi-start strategy reduces dependence on initialisation.
    km       <- stats::kmeans(X, centers = k, nstart = 20)
    sil      <- cluster::silhouette(km$cluster, d)
    mean_sil <- mean(sil[, "sil_width"])
    info("k-means k=%d -> mean silhouette = %.4f", k, mean_sil)

    if (mean_sil > best_sil) {
      best_sil      <- mean_sil
      best_k        <- k
      best_clusters <- km$cluster
      best_km       <- km
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


# =============================================================================
# Section 7 — Main Driver
# =============================================================================

# -----------------------------------------------------------------------------
# run_hclust_kmeans()
#
# This is the single entry-point function called by all 12 wrapper scripts.
# It orchestrates the complete statistical pipeline: data loading, QC,
# optional PCA, clustering, model selection, and result serialisation.
#
# The 'kind' string encodes the full configuration as a compound token:
#
#   {method_prefix}_{scope_suffix}
#
#   method_prefix : one of { hc, kmeans, pca_hc, pca_kmeans }
#   scope_suffix  : one of { cell, tumour, cell_tumour }
#
# Pipeline summary:
#   Step 1 — build_sample_matrix()      : data ingestion and feature alignment
#   Step 1b— filter_zero_var_genes()    : baseline QC; removes constant features
#   Step 2 — method type inference      : routes to HC or k-means branch
#   Step 3 — compute_pca_scores() [opt] : dimensionality reduction for pca_* kinds
#   Step 4 — hc_optimal() / kmeans_optimal(): partitioning and model selection
#   Step 5 — result assembly and RDS serialisation
#
# The returned list encapsulates all components needed for downstream
# visualisation (UMAP overlays, silhouette plots) and integration with
# tumour–cell line similarity scoring.
# -----------------------------------------------------------------------------
run_hclust_kmeans <- function(kind,
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
  info("==== Running HC/k-means clustering: %s ====", kind)

  # ------------------------------------------------------------------
  # Step 1: Data ingestion and sample matrix construction.
  # ------------------------------------------------------------------
  built   <- build_sample_matrix(
    kind         = kind,
    cell_rds     = cell_rds,
    tumour_rds   = tumour_rds,
    cell_label   = cell_label,
    tumour_label = tumour_label
  )
  X       <- built$mat
  dataset <- built$dataset

  # Step 1b: Baseline QC — remove exactly constant genes (SD = 0).
  X <- filter_zero_var_genes(X, min_sd = 0)

  # ------------------------------------------------------------------
  # Step 2: Method type inference from 'kind' prefix.
  # This determines whether PCA is applied and which clustering
  # algorithm is invoked.
  # ------------------------------------------------------------------
  method_type <- if (grepl("^pca_hc",     kind)) {
    "pca_hc"
  } else if (grepl("^pca_kmeans", kind)) {
    "pca_kmeans"
  } else if (grepl("^hc_",        kind)) {
    "hc"
  } else if (grepl("^kmeans_",    kind)) {
    "kmeans"
  } else {
    stop("Unknown 'kind' prefix for method type: ", kind)
  }

  # ------------------------------------------------------------------
  # Step 3: Optional PCA dimensionality reduction.
  # For pca_* kinds, the clustering operates on the PC score matrix
  # rather than the full gene expression space.
  # ------------------------------------------------------------------
  if (method_type %in% c("pca_hc", "pca_kmeans")) {
    pca_res <- compute_pca_scores(X, n_pcs = n_pcs)
    X_use   <- pca_res$scores   # n_samples × n_pcs
  } else {
    # Raw mode: distances and cluster assignments are derived directly
    # from the filtered gene expression matrix.
    pca_res <- NULL
    X_use   <- X
  }

  # ------------------------------------------------------------------
  # Step 4: Clustering and silhouette-based model selection.
  # ------------------------------------------------------------------
  if (method_type %in% c("hc", "pca_hc")) {
    cl_res   <- hc_optimal(X_use, max_k = max_k_hc,
                           dist_method = dist_method, seed = seed)
    best_k   <- cl_res$best_k
    clusters <- cl_res$clusters
    hc_res   <- cl_res
    km_res   <- NULL
  } else {
    cl_res   <- kmeans_optimal(X_use, k_min = k_min, k_max = k_max, seed = seed)
    best_k   <- cl_res$best_k
    clusters <- cl_res$clusters
    km_res   <- cl_res
    hc_res   <- NULL
  }

  # ------------------------------------------------------------------
  # Step 5: Result assembly and serialisation.
  # Cluster labels are prefixed with 'C' (e.g., C1, C2) to produce
  # human-readable, factor-compatible identifiers in downstream tables.
  # ------------------------------------------------------------------
  clusters         <- as.integer(clusters)
  names(clusters)  <- rownames(X_use)
  clusters_labeled <- setNames(paste0("C", clusters), names(clusters))

  # Derive a human-readable input-type string for result provenance.
  input_type <- if (grepl("cell_tumour$", kind)) {
    "CELL+TUMOUR"
  } else if (grepl("_cell$", kind)) {
    "CELL"
  } else {
    "TUMOUR"
  }

  # Assemble the result object, preserving all components needed for
  # downstream visualisation, comparison, and reproducibility auditing.
  res_obj <- list(
    kind             = kind,             # Full configuration identifier string.
    method_type      = method_type,      # Inferred algorithm branch.
    input_type       = input_type,       # Sample scope of the analysis.
    k                = best_k,           # Optimal k selected by silhouette criterion.
    clusters         = clusters,         # Integer cluster assignments, named by sample ID.
    clusters_labeled = clusters_labeled, # Character labels (C1, C2, …) for figures.
    dataset          = dataset,          # Cohort factor (CELL / TUMOUR) per sample.
    pca              = pca_res,          # Full PCA object (NULL for raw-feature kinds).
    hc               = hc_res,           # HC result list (NULL for k-means kinds).
    kmeans           = km_res,           # k-means result list (NULL for HC kinds).
    params = list(
      n_pcs    = n_pcs,       # PCs retained in dimensionality reduction.
      max_k_hc = max_k_hc,    # Upper bound of HC k sweep.
      k_min    = k_min,       # Lower bound of k-means k sweep.
      k_max    = k_max,       # Upper bound of k-means k sweep.
      dist     = dist_method, # Distance metric used in HC or silhouette evaluation.
      seed     = seed         # Random seed for reproducibility.
    ),
    timestamp = Sys.time()    # Wall-clock time of execution for audit trail.
  )

  info("Saving cluster object to: %s", cluster_rds_path)
  saveRDS(res_obj, cluster_rds_path)
  invisible(res_obj)
}