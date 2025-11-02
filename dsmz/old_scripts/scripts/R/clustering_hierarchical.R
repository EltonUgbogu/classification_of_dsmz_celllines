run_hierarchical_clustering <- function(
  PC,
  V_adj,
  features,
  metadata_ann,
  V_tcga,
  dataset_lab,
  outdir,
  k_fixed = NULL,
  seed = 42
) {
  # ===========================================================================
  # 1. Input Validation & Setup
  # ===========================================================================
  
  ## --- Basic type and dimension checks ---
  stopifnot(
    is.matrix(PC), 
    is.matrix(V_adj), 
    length(features) > 0,
    is.character(features)
  )
  
  n_samples <- nrow(PC)
  if (n_samples < 5) {
    stop("At least 5 samples required for clustering.")
  }
  
  ## --- Sample alignment between PC and V_adj ---
  if (n_samples != ncol(V_adj)) {
    stop(sprintf(
      "PC has %d samples but V_adj has %d. They must match.",
      n_samples, ncol(V_adj)
    ))
  }
  
  ## --- Sample ID handling ---
  if (is.null(rownames(PC))) {
    rownames(PC) <- paste0("sample_", seq_len(n_samples))
    log_warn("[HC] PC missing rownames. Assigned: sample_1, sample_2, ...")
  }
  if (is.null(colnames(V_adj))) {
    colnames(V_adj) <- rownames(PC)
    log_warn("[HC] V_adj missing colnames. Inferred from PC rownames.")
  }
  
  ## --- Enforce matching and ordering ---
  common_samples <- intersect(rownames(PC), colnames(V_adj))
  if (length(common_samples) == 0) {
    stop("No overlapping sample IDs between PC and V_adj.")
  }
  if (length(common_samples) < n_samples) {
    missing_in_V <- setdiff(rownames(PC), colnames(V_adj))
    missing_in_PC <- setdiff(colnames(V_adj), rownames(PC))
    log_warn("[HC] %d samples in PC but not V_adj: %s", 
             length(missing_in_V), paste(head(missing_in_V), collapse = ", "))
    log_warn("[HC] %d samples in V_adj but not PC: %s", 
             length(missing_in_PC), paste(head(missing_in_PC), collapse = ", "))
  }
  
  # Reorder V_adj to match PC
  V_adj <- V_adj[, rownames(PC), drop = FALSE]
  
  ## --- Dataset label ---
  if (length(dataset_lab) == 1) {
    dataset_lab <- rep(dataset_lab, n_samples)
  } else if (length(dataset_lab) != n_samples) {
    stop("dataset_lab must be length 1 or length equal to nrow(PC).")
  }
  
  ## --- Metadata annotation ---
  if (is.character(metadata_ann) && length(metadata_ann) == 1) {
    if (!file.exists(metadata_ann)) {
      stop(sprintf("metadata_ann file not found: %s", metadata_ann))
    }
  } else if (!is.data.frame(metadata_ann)) {
    stop("metadata_ann must be a file path or data.frame.")
  }
  
  ## --- Feature availability ---
  if (is.null(rownames(V_adj))) {
    stop("V_adj must have rownames (gene symbols/IDs).")
  }
  features_present <- intersect(features, rownames(V_adj))
  if (length(features_present) == 0) {
    stop("None of the provided features are present in V_adj.")
  }
  if (length(features_present) < length(features)) {
    log_warn("[HC] %d/%d features missing in V_adj. Proceeding with %d.",
             length(features) - length(features_present), length(features), length(features_present))
    features <- features_present
  }
  
  ## --- Output directory ---
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  log_info("[HC] Output directory: %s", normalizePath(outdir))
  
  ## --- Initialize result container ---
  result_list <- list()
  
  # ===========================================================================
  # 2. Fixed-k Hierarchical Clustering (Optional)
  # ===========================================================================
  
  if (!is.null(k_fixed)) {
    if (!is.numeric(k_fixed) || k_fixed < 2 || k_fixed != as.integer(k_fixed)) {
      stop("k_fixed must be an integer >= 2.")
    }
    k_fixed <- as.integer(k_fixed)
    
    log_info("[HC] === Fixed-k Clustering (k = %d) ===", k_fixed)
    
    # Run HC
    hc_fixed <- run_fixed_hc(PC, k = k_fixed, seed = seed)
    stopifnot(
      !is.null(hc_fixed$hc), 
      inherits(hc_fixed$hc, "hclust"),
      length(hc_fixed$clusters) == n_samples
    )
    
    # Label using TCGA centroids
    log_info("[HC] Mapping %d clusters to TCGA subtypes...", k_fixed)
    clusters_labeled <- label_by_tcga(
      clusters = hc_fixed$clusters,
      V_adj = V_adj,
      features = features,
      metadata_ann = metadata_ann,
      V_tcga = V_tcga
    )
    
    # Save table
    df_fixed <- tibble::tibble(
      sample = rownames(PC),
      dataset = dataset_lab,
      cluster = as.character(clusters_labeled)
    )
    csv_path <- file.path(outdir, sprintf("clusters_hierarchical_fixed_k=%d.csv", k_fixed))
    write.csv(df_fixed, csv_path, row.names = FALSE)
    log_info("[HC] Fixed-k assignments saved: %s", csv_path)
    
    # Plot dendrogram
    pdf_path <- file.path(outdir, sprintf("dendrogram_hierarchical_fixed_k=%d.pdf", k_fixed))
    safe_pdf(pdf_path, {
      plot_dendrogram_fixed(hc_fixed$hc, k = k_fixed, main = sprintf("Fixed-k = %d", k_fixed))
    })
    log_info("[HC] Dendrogram saved: %s", pdf_path)
    
    result_list$fixed <- list(
      hc = hc_fixed$hc,
      clusters_raw = hc_fixed$clusters,
      clusters_labeled = clusters_labeled,
      df = df_fixed
    )
  } else {
    log_info("[HC] k_fixed = NULL → skipping fixed-k clustering.")
  }
  
  # ===========================================================================
  # 3. Dynamic Tree Cut Clustering
  # ===========================================================================
  
  log_info("[HC] === Dynamic Tree Cut Clustering ===")
  
  dyn_res <- run_dynamic_cut(PC)
  stopifnot(
    !is.null(dyn_res$hc), inherits(dyn_res$hc, "hclust"),
    length(dyn_res$clusters) == n_samples,
    inherits(dyn_res$dist, "dist")
  )
  
  k_dyn <- length(unique(dyn_res$clusters))
  log_info("[HC] Dynamic cut discovered %d clusters.", k_dyn)
  
  # Label using TCGA
  log_info("[HC] Mapping dynamic clusters to TCGA subtypes...")
  clusters_dyn_labeled <- label_by_tcga(
    clusters = dyn_res$clusters,
    V_adj = V_adj,
    features = features,
    metadata_ann = metadata_ann,
    V_tcga = V_tcga
  )
  
  # Save assignments
  df_dyn <- tibble::tibble(
    sample = rownames(PC),
    dataset = dataset_lab,
    cluster = as.character(clusters_dyn_labeled)
  )
  csv_path_dyn <- file.path(outdir, sprintf("clusters_hierarchical_dynamic_tree_cut=%d.csv", k_dyn))
  write.csv(df_dyn, csv_path_dyn, row.names = FALSE)
  log_info("[HC] Dynamic assignments saved: %s", csv_path_dyn)
  
  # Plot dendrogram
  pdf_path_dyn <- file.path(outdir, sprintf("dendrogram_hierarchical_dynamic_tree_cut=%d.pdf", k_dyn))
  safe_pdf(pdf_path_dyn, {
    plot_dendrogram_dynamic(dyn_res$hc, clusters_dyn_labeled, main = sprintf("Dynamic Cut (k = %d)", k_dyn))
  })
  log_info("[HC] Dynamic dendrogram saved: %s", pdf_path_dyn)
  
  # ===========================================================================
  # 4. Silhouette Analysis (Dynamic Cut Only)
  # ===========================================================================
  
  sil_obj <- NULL
  if (k_dyn >= 2) {
    dist_ok <- inherits(dyn_res$dist, "dist") &&
               attr(dyn_res$dist, "Size") == n_samples
    
    if (dist_ok) {
      sil_obj <- cluster::silhouette(as.integer(factor(dyn_res$clusters)), dyn_res$dist)
      pdf_path_sil <- file.path(outdir, sprintf("silhouette_dynamic_tree_cut=%d.pdf", k_dyn))
      safe_pdf(pdf_path_sil, {
        plot_silhouette(sil_obj, main = sprintf("Silhouette (k = %d)", k_dyn))
      })
      log_info("[HC] Silhouette plot saved: %s", pdf_path_sil)
      log_info("[HC] Mean silhouette width: %.3f", mean(sil_obj[, 3]))
    } else {
      log_warn("[HC] Skipping silhouette: invalid or missing distance matrix.")
    }
  } else {
    log_warn("[HC] Skipping silhouette: k = %d < 2.", k_dyn)
  }
  
  # ===========================================================================
  # 5. Assemble and Return Results
  # ===========================================================================
  
  result_list$dynamic <- list(
    hc = dyn_res$hc,
    dist = dyn_res$dist,
    clusters_raw = dyn_res$clusters,
    clusters_labeled = clusters_dyn_labeled,
    df = df_dyn,
    k = k_dyn,
    silhouette = sil_obj
  )
  
  log_info("[HC] Hierarchical clustering complete. Results in: %s", normalizePath(outdir))
  invisible(result_list)
}