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
  