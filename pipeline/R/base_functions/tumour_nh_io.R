# ==============================================================================
# tumour_nh_io.R
# Pipeline-Native Method Discovery for Tumour Neighbourhood Computation
# ==============================================================================
#
# This module dynamically discovers clustering outputs from upstream pipeline
# stages (agnostic_clustering and consensus) for tumour neighbourhood analysis.
#
# ==============================================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# get_nh_methods: Discover available clustering methods for a given direction
# ------------------------------------------------------------------------------
# Returns a tibble with columns: method_id, path, outdir, distance, exists

get_nh_methods <- function(unsup_root, direction) {
  
  stopifnot(dir.exists(unsup_root))
  stopifnot(direction %in% c("pam50_euc", "pam50_corr", "hvg_euc", "hvg_corr"))
  
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
  
  # Step 1: Agnostic clustering (HC available for all; KM for Euclidean only)
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
  
  # Step 2: Consensus clustering (CCP)
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
  
  # Combine and check file existence
  methods <- dplyr::bind_rows(methods) %>%
    dplyr::mutate(exists = file.exists(path))
  
  methods
}
