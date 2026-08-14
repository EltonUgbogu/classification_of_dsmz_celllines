# ==============================================================================
# tumour_nh_io.R
# Eligible Clustering-Method Catalogue for Tumour Neighbourhood Computation
# ==============================================================================
#
# This module catalogues the JOINT cell-line + tumour clustering outputs of the
# two upstream clustering stages that are eligible for tumour-neighbourhood
# construction and the p_consensus recurrence fraction:
#
#   AGN_* : HC/k-means clustering — hierarchical clustering and k-means applied
#           to each feature-distance representation in expression space or
#           PCA-reduced space ("AGN" is the retained internal compatibility
#           identifier for this branch; outputs live under agnostic_clustering/).
#   CCP_* : ConsensusClusterPlus resampling-based consensus clustering
#           (outputs live under consensus/).
#
# Only JOINT (*_cell_tumour) outputs are catalogued. Cell-only and tumour-only
# clustering outputs are never eligible. The catalogue maps method identifiers
# to expected file paths; the exact contributing set is declared in
# patient_referenced_graph.clustering_methods_by_distance and validated by the
# callers — eligibility is configuration-driven, never inferred from whichever
# files happen to exist on disk.
#
# ==============================================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# get_nh_methods: Catalogue eligible clustering methods for a given direction
# ------------------------------------------------------------------------------
# Returns a tibble with columns: method_id, path, outdir, distance, exists

get_nh_methods <- function(unsup_root, direction) {

  stopifnot(dir.exists(unsup_root))
  # Accept any direction ending in _euc or _corr (all feature sets + PAM50 + HVG)
  # Remove the hardcoded constraint to allow Variance_euc, MAD_euc, etc.
  if (!grepl("_(euc|corr)$", direction)) {
    stop("direction must end with _euc or _corr, got: ", direction)
  }

  # Infer distance metric from direction suffix
  dist_type <- if (grepl("_corr$", direction)) "correlation" else "euclidean"

  # Root directories for upstream clustering stages
  agn_root  <- file.path(unsup_root, "agnostic_clustering", direction)
  cons_root <- file.path(unsup_root, "consensus", direction)

  # Helper to construct a method record
  mk <- function(method_id, path) {
    tibble::tibble(
      method_id = method_id,
      path      = path,
      outdir    = file.path(unsup_root, "tumour_neighbourhoods", direction, method_id),
      distance  = dist_type
    )
  }

  methods <- list()

  # Step 1: HC/k-means clustering, JOINT outputs (HC for all directions;
  # k-means for Euclidean directions only, since k-means minimises Euclidean
  # inertia and is undefined under correlation distance).
  methods <- c(methods,
    list(
      mk("AGN_HC_expr_cell_tumour",
         file.path(agn_root, "hc_cell_tumour", "hc_cell_tumour_clusters_optimal.rds")),
      mk("AGN_HC_pca_cell_tumour",
         file.path(agn_root, "pca_hc_cell_tumour", "pca_hc_cell_tumour_clusters_optimal.rds"))
    )
  )
  
  if (dist_type == "euclidean") {
    methods <- c(methods,
      list(
        mk("AGN_KM_expr_cell_tumour",
           file.path(agn_root, "kmeans_cell_tumour", "kmeans_cell_tumour_clusters_optimal.rds")),
        mk("AGN_KM_pca_cell_tumour",
           file.path(agn_root, "pca_kmeans_cell_tumour", "pca_kmeans_cell_tumour_clusters_optimal.rds"))
      )
    )
  }
  
  # Step 2: ConsensusClusterPlus consensus clustering, JOINT outputs (same
  # HC-everywhere / k-means-Euclidean-only eligibility as the HC/k-means branch).
  methods <- c(methods,
    list(
      mk("CCP_HC_expr_cell_tumour",
         file.path(cons_root, "ccp_hc_expr_cell_tumour", "ccp_hc_expr_cell_tumour_clusters_optimal.rds")),
      mk("CCP_HC_pca_cell_tumour",
         file.path(cons_root, "ccp_hc_pca_cell_tumour", "ccp_hc_pca_cell_tumour_clusters_optimal.rds"))
    )
  )
  
  if (dist_type == "euclidean") {
    methods <- c(methods,
      list(
        mk("CCP_KM_expr_cell_tumour",
           file.path(cons_root, "ccp_kmeans_expr_cell_tumour", "ccp_kmeans_expr_cell_tumour_clusters_optimal.rds")),
        mk("CCP_KM_pca_cell_tumour",
           file.path(cons_root, "ccp_kmeans_pca_cell_tumour", "ccp_kmeans_pca_cell_tumour_clusters_optimal.rds"))
      )
    )
  }
  
  # Combine and record file existence. The `exists` flag is informational for
  # the callers' validation errors; callers must not use it to shrink the
  # configured method set silently.
  methods <- dplyr::bind_rows(methods) %>%
    dplyr::mutate(exists = file.exists(path))

  methods
}

# ------------------------------------------------------------------------------
# make_nh_paths: Create path configuration object (legacy compatibility)
# ------------------------------------------------------------------------------
# This function exists for backward compatibility with scripts that expect
# a path configuration object. The actual path logic is handled by get_nh_methods().
# Returns a minimal list structure to satisfy the API check.

make_nh_paths <- function(unsup_root) {
  # Return a minimal structure - not actually used by current pipeline
  list(
    unsup_root = unsup_root,
    tumour_neighbourhoods = file.path(unsup_root, "tumour_neighbourhoods")
  )
}
