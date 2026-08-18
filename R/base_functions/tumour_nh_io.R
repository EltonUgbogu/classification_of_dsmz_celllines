# ==============================================================================
# tumour_nh_io.R
# Path resolution for configured clustering formulations
# ==============================================================================
#
# This module maps a clustering-formulation identifier to the cluster-assignment
# file that formulation writes, and to the tumour-neighbourhood output directory
# it should populate. It is a naming/path mapping only.
#
# It does NOT declare which formulations exist or which are eligible: that is a
# scientific decision owned by configuration
# (patient_referenced_graph.clustering_methods_by_distance) and passed in by the
# caller. Holding a second copy of the eligible set here would make the
# p-consensus denominator depend on whichever list happened to be edited.
#
# Identifier grammar, matching the configured identifiers:
#
#     [CCP_]<ALGORITHM>_<SPACE>_cell_tumour
#
#     CCP_       present for ConsensusClusterPlus formulations, written under
#                consensus/; absent for the HC/k-means formulations, written
#                under hclust_kmeans/
#     ALGORITHM  HCLUST  hierarchical clustering
#                KMEANS  k-means
#     SPACE      expr    expression space
#                pca     PCA-reduced space
#
# Only JOINT cell-line + tumour formulations are representable, so cell-only and
# tumour-only clustering outputs cannot be addressed through this mapping at all.
#
# ==============================================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# nh_method_directory: on-disk directory name for a formulation identifier
# ------------------------------------------------------------------------------
# The clustering rules name their output directories by algorithm and space
# rather than by the identifier used in configuration, so the two vocabularies
# are related here in one place.

nh_method_directory <- function(method_id) {
  is_consensus <- startsWith(method_id, "CCP_")
  remainder <- if (is_consensus) sub("^CCP_", "", method_id) else method_id

  parts <- strsplit(remainder, "_", fixed = TRUE)[[1]]
  if (length(parts) < 4L) {
    stop("Unrecognised clustering formulation identifier: ", method_id)
  }
  algorithm <- parts[[1]]
  space <- parts[[2]]
  scope <- paste(parts[-(1:2)], collapse = "_")

  if (!identical(scope, "cell_tumour")) {
    stop("Only JOINT cell-line + tumour formulations are eligible; got: ", method_id)
  }
  if (!algorithm %in% c("HCLUST", "KMEANS")) {
    stop("Unknown clustering algorithm in identifier '", method_id,
         "'; expected HCLUST or KMEANS.")
  }
  if (!space %in% c("expr", "pca")) {
    stop("Unknown feature space in identifier '", method_id,
         "'; expected expr or pca.")
  }

  # The clustering rules name their output directories in lower case and by
  # algorithm stem, so the configured identifier is translated here.
  stem <- if (identical(algorithm, "HCLUST")) "hc" else "kmeans"

  if (is_consensus) {
    # ConsensusClusterPlus directories: ccp_hc_expr_cell_tumour and siblings.
    list(root = "consensus",
         dirname = paste0("ccp_", stem, "_", space, "_cell_tumour"))
  } else {
    # HC/k-means clustering directories: hc_cell_tumour, pca_hc_cell_tumour,
    # kmeans_cell_tumour, pca_kmeans_cell_tumour.
    dirname <- if (identical(space, "pca")) paste0("pca_", stem) else stem
    list(root = "hclust_kmeans",
         dirname = paste0(dirname, "_cell_tumour"))
  }
}

# ------------------------------------------------------------------------------
# get_nh_methods: resolve paths for the configured formulations
# ------------------------------------------------------------------------------
# unsup_root  profile-scoped unsupervised results root
# direction   feature-distance representation
# method_ids  the exact configured formulation identifiers for this
#             representation, supplied by the caller from configuration
#
# Returns a tibble with columns: method_id, path, outdir, distance, exists.
# The `exists` flag is informational for the caller's validation messages;
# callers must not use it to silently shrink the configured set.

get_nh_methods <- function(unsup_root, direction, method_ids) {

  stopifnot(dir.exists(unsup_root))
  if (!grepl("_(euc|corr)$", direction)) {
    stop("direction must end with _euc or _corr, got: ", direction)
  }
  if (missing(method_ids) || length(method_ids) == 0L) {
    stop("get_nh_methods() requires the configured clustering formulation ",
         "identifiers for '", direction, "'. They are declared in ",
         "patient_referenced_graph.clustering_methods_by_distance and must be ",
         "supplied by the caller rather than reconstructed here.")
  }
  method_ids <- as.character(method_ids)
  if (anyDuplicated(method_ids) > 0L) {
    stop("Duplicate clustering formulation identifier(s): ",
         paste(unique(method_ids[duplicated(method_ids)]), collapse = ", "))
  }

  # The distance component is carried by the representation identifier and is
  # reported for downstream labelling; it does not select the formulation set.
  dist_type <- if (grepl("_corr$", direction)) "correlation" else "euclidean"

  records <- lapply(method_ids, function(method_id) {
    loc <- nh_method_directory(method_id)
    tibble::tibble(
      method_id = method_id,
      path = file.path(unsup_root, loc$root, direction, loc$dirname,
                       paste0(loc$dirname, "_clusters_optimal.rds")),
      outdir = file.path(unsup_root, "tumour_neighbourhoods", direction, method_id),
      distance = dist_type
    )
  })

  dplyr::bind_rows(records) %>%
    dplyr::mutate(exists = file.exists(path))
}

# ------------------------------------------------------------------------------
# make_nh_paths: path configuration object (legacy compatibility)
# ------------------------------------------------------------------------------
# Retained because compute_tumour_neighbourhoods.R checks for this factory.

make_nh_paths <- function(unsup_root) {
  list(
    unsup_root = unsup_root,
    tumour_neighbourhoods = file.path(unsup_root, "tumour_neighbourhoods")
  )
}
