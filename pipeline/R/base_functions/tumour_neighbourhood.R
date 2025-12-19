# tumour_neighbourhood.R
# Compute local tumour neighbourhoods for each DSMZ cell line
# within a given clustering (method m)
#
# Supports both Euclidean and correlation distance metrics

compute_tumour_neighbourhoods <- function(
  emb_mat,
  cluster_m,
  dataset,
  top_frac  = 0.10,
  top_n_min = 30,
  top_n_max = 200,
  method_id = "unspecified_method",
  distance  = c("euclidean", "correlation")
) {
  # Match distance argument
  distance <- match.arg(distance)
  
  # --------------------- Input validation ---------------------
  if (!is.matrix(emb_mat) || !is.numeric(emb_mat)) {
    stop("emb_mat must be a numeric matrix (samples × features).")
  }
  if (is.null(rownames(emb_mat))) {
    stop("emb_mat must have rownames corresponding to sample IDs.")
  }
  if (missing(cluster_m) || is.null(cluster_m)) {
    stop("cluster_m is required. Supply a named vector of cluster labels.")
  }
  if (is.null(names(cluster_m))) {
    stop("cluster_m must be named with sample IDs matching rownames(emb_mat).")
  }
  if (missing(dataset) || is.null(dataset)) {
    stop("dataset vector is required (values 'DSMZ' or 'TCGA').")
  }
  if (is.null(names(dataset))) {
    stop("dataset must be named with sample IDs.")
  }

  # Align all inputs to common sample IDs
  common_ids <- intersect(rownames(emb_mat), intersect(names(cluster_m), names(dataset)))
  if (length(common_ids) == 0) {
    stop("No overlapping sample IDs between emb_mat, cluster_m, and dataset.")
  }

  emb_mat   <- emb_mat[common_ids, , drop = FALSE]
  cluster_m <- cluster_m[common_ids]
  dataset   <- dataset[common_ids]

  # Identify DSMZ cell lines and TCGA tumours
  dsmz_ids <- names(dataset)[dataset == "DSMZ"]
  tcga_ids <- names(dataset)[dataset == "TCGA"]

  if (length(dsmz_ids) == 0) stop("No DSMZ samples found.")
  if (length(tcga_ids) == 0) stop("No TCGA tumours found.")

  # --------------------- Distance functions ---------------------
  euclidean_dist <- function(query_vec, target_mat) {
    diffs <- sweep(target_mat, 2, query_vec, "-")
    sqrt(rowSums(diffs^2))
  }
  
  correlation_dist <- function(query_vec, target_mat) {
    # Distance = 1 - Pearson correlation
    # Each row of target_mat is correlated with query_vec
    cors <- cor(t(target_mat), query_vec,
                use = "pairwise.complete.obs",
                method = "pearson")
    as.numeric(1 - cors)
  }
  
  # Select distance function based on argument
  dist_fun <- if (distance == "euclidean") euclidean_dist else correlation_dist
  
  cat("[INFO] Using", distance, "distance for tumour neighbourhood computation\n")

  # --------------------- Main loop over DSMZ cell lines ---------------------
  neighbourhoods <- vector("list", length(dsmz_ids))
  names(neighbourhoods) <- dsmz_ids
  long_records <- vector("list", length(dsmz_ids))

  for (cell_line in dsmz_ids) {
    cluster_of_cell <- cluster_m[[cell_line]]
    if (is.na(cluster_of_cell)) {
      warning("Cell line ", cell_line, " has NA cluster → skipped.")
      next
    }

    # TCGA tumours in the same cluster
    candidate_tumours <- tcga_ids[cluster_m[tcga_ids] == cluster_of_cell]

    # Fallback: if cluster too small, use all TCGA tumours
    if (length(candidate_tumours) < 10) {
      message("[", method_id, "] ", cell_line,
              ": only ", length(candidate_tumours),
              " tumours in cluster → using all TCGA tumours.")
      candidate_tumours <- tcga_ids
    }

    # Compute distances using the selected metric
    query_vec <- emb_mat[cell_line, , drop = TRUE]
    target_mat <- emb_mat[candidate_tumours, , drop = FALSE]
    distances <- dist_fun(query_vec, target_mat)

    # Rank tumours by proximity
    ord <- order(distances)
    ranked_tumours <- candidate_tumours[ord]
    ranked_dist   <- distances[ord]

    # Adaptive neighbourhood size
    n_total <- length(ranked_tumours)
    n_frac  <- ceiling(top_frac * n_total)
    n_keep  <- max(top_n_min, n_frac)
    n_keep  <- min(n_keep, top_n_max, n_total)

    top_tumours <- ranked_tumours[seq_len(n_keep)]

    # Store results
    neighbourhoods[[cell_line]] <- top_tumours

    # Long format for consensus analysis
    long_records[[cell_line]] <- data.frame(
      method      = method_id,
      cell_line   = cell_line,
      tumor_id    = ranked_tumours,
      rank        = seq_along(ranked_tumours),
      distance    = ranked_dist,
      in_top      = ranked_tumours %in% top_tumours,
      cluster     = cluster_of_cell,
      stringsAsFactors = FALSE
    )
  }

  long_df <- bind_rows(long_records)

  # --------------------- Return ---------------------
  list(
    neighbourhoods = neighbourhoods,   # list: cell_line → vector of top TCGA tumour barcodes
    long_df       = long_df,           # tidy table for consensus voting
    method_id     = method_id,
    distance      = distance,          # NEW: which distance metric was used
    params        = list(top_frac = top_frac, top_n_min = top_n_min, top_n_max = top_n_max),
    n_dsmz        = length(dsmz_ids),
    n_tcga        = length(tcga_ids)
  )
}