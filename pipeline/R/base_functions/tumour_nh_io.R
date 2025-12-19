# tumour_nh_io.R  (pipeline-native)
# -------------------------------------------------------------------
# Pipeline-native method discovery for tumour neighbourhood computation
# Replaces legacy path registry with dynamic discovery based on actual
# pipeline outputs from Step 1 (agnostic_clustering) and Step 2 (consensus)
# -------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# Return a tibble of clustering methods available for a given direction.
# This matches your pipeline outputs:
#  - Step 1 (agnostic_clustering): hc_* and kmeans_* kinds
#  - Step 2 (consensus): ccp_* kinds
get_nh_methods <- function(unsup_root, direction) {
  stopifnot(dir.exists(unsup_root))
  stopifnot(direction %in% c("pam50_euc","pam50_corr","hvg_euc","hvg_corr"))

  # distance type for tumour-neighbourhood computation
  dist_type <- if (grepl("_corr$", direction)) "correlation" else "euclidean"

  # roots
  agn_root  <- file.path(unsup_root, "agnostic_clustering", direction)
  cons_root <- file.path(unsup_root, "consensus", direction)

  # helper
  mk <- function(method_id, path) {
    tibble::tibble(
      method_id = method_id,
      path      = path,
      outdir    = file.path(unsup_root, "tumour_neighbourhoods", direction, method_id),
      distance  = dist_type
    )
  }

  methods <- list()

  # -------------------------
  # Step 1: agnostic (integrated)
  # -------------------------
  methods <- c(methods,
    list(
      mk("AGN_HC_expr_cell_tumour",
         file.path(agn_root, "hc_cell_tumour", "hc_cell_tumour_clusters_optimal.rds")),
      mk("AGN_HC_pca_cell_tumour",
         file.path(agn_root, "pca_hc_cell_tumour", "pca_hc_cell_tumour_clusters_optimal.rds"))
    )
  )

  # kmeans only exists for *_euc directions (your Snakefile enforces this)
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

  # -------------------------
  # Step 2: CCP consensus (integrated)
  # -------------------------
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

  methods <- dplyr::bind_rows(methods) %>%
    dplyr::mutate(exists = file.exists(path))

  methods
}

# Legacy compatibility: provide empty nh_paths for backwards compatibility
# (some scripts may still reference it, but new code should use get_nh_methods())
nh_paths <- list()

make_nh_paths <- function(unsup_root) {
  # Return empty list for legacy compatibility
  # New code should use get_nh_methods(unsup_root, direction) instead
  list()
}
