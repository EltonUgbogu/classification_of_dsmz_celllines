run_iterative_hdbscan <- function(
  PC,                     # matrix: samples x PCs (rownames = sample IDs)
  cluster_labels,         # character: final label names (levels of the output factor)
  seed = 42,
  metric = "cosine",
  umap_n_first = 30L,
  umap_n_next  = 15L,
  umap_min_dist = 0.1,
  min_unassigned_stop = 15L,
  minPts_fun = function(n) max(5L, floor(log2(n) + 3)),  # adaptive minPts rule
  reproject = FALSE,               # if TRUE: re-run PCA on the remaining subset each pass
  original_features = NULL,        # required when reproject = TRUE
  plot_iter = FALSE,               # if TRUE: save a UMAP PNG per pass to `plot_dir`
  plot_dir = "hdbscan_debug",
  logger = log_info                # logging function; expects sprintf-style formatting
) {
  # ------------------------------ 0) Setup & validation -----------------------

  # Validate shape and label requirements early for clearer error messages
  stopifnot(is.matrix(PC), nrow(PC) >= 5,
            is.character(cluster_labels), length(cluster_labels) >= 2)

  # Ensure samples have names; use synthetic names if missing
  if (is.null(rownames(PC))) rownames(PC) <- paste0("s", seq_len(nrow(PC)))

  # If reprojecting onto PCA of remainder, we need the original feature matrix
  if (reproject && is.null(original_features))
    stop("`original_features` must be provided when `reproject = TRUE`.")

  # Make runs deterministic: same seed → same UMAP/HDBSCAN randomness
  set.seed(seed)

  # Normalize PC columns so each dimension contributes comparably to distance
  # (You can disable this if your PC matrix is already scaled.)
  PC <- scale(PC)

  # Number of final labels we aim to map to
  k_desired <- length(cluster_labels)

  # Internal assignment vector over samples:
  # - "Unassigned" initially
  # - "1","2",... for passes/peels
  # - "Noise" only if no clusters were ever found by HDBSCAN
  assignments <- rep("Unassigned", nrow(PC))
  names(assignments) <- rownames(PC)

  # Iteration counter (also becomes the numeric label we assign)
  cur_pass <- 1L

  # Per-pass diagnostics
  silhouette_history <- numeric()  # mean silhouette of the *extracted* cluster each pass
  stability_history  <- numeric()  # relative change vs previous mean silhouette

  # Prepare plot directory if needed
  if (plot_iter) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------- Helper: single HDBSCAN pass ----------------------
  # Given a vector of sample indices (`idx`), run:
  #  - Optional re-PCA (if reproject=TRUE)
  #  - UMAP embedding (with fixed neighbors for the pass)
  #  - HDBSCAN on the embedding
  #  - Select the largest non-noise cluster, and compute its mean silhouette
  one_pass <- function(idx, n_neigh, pass_num) {
    # Abort if too few samples to meaningfully cluster
    if (length(idx) < 5) return(list(ok = FALSE))

    logger("[HDBSCAN] Pass %d | %d samples | n_neighbors=%d",
           pass_num, length(idx), n_neigh)

    # Optionally re-fit PCA on the remainder to sharpen structure for the subset
    # (uses the feature-space provided in `original_features`)
    if (reproject && length(idx) < nrow(PC)) {
      pca_rem <- stats::prcomp(original_features[idx, ], rank. = ncol(PC))
      X_sub <- pca_rem$x
      # Keep original names to allow name-based assignment
      rownames(X_sub) <- rownames(original_features)[idx]
    } else {
      # Otherwise, use the existing (scaled) PC submatrix
      X_sub <- PC[idx, , drop = FALSE]
    }

    # UMAP: set spectral init and single-thread SGD for repeatability
    emb <- uwot::umap(
      X_sub,
      n_neighbors   = n_neigh,
      min_dist      = umap_min_dist,
      metric        = metric,
      init          = "spectral",
      n_sgd_threads = 1,
      verbose       = FALSE
    )

    # Adaptive minPts based on subset size (guard against <2)
    mp <- max(2L, as.integer(minPts_fun(length(idx))))
    logger("[HDBSCAN]   minPts = %d (adaptive)", mp)

    # Density clustering on the embedding
    hdb <- dbscan::hdbscan(emb, minPts = mp)

    # Report how much was labeled as noise (cluster id 0)
    noise_frac <- mean(hdb$cluster == 0)
    logger("[HDBSCAN]   Noise: %.1f%%", noise_frac * 100)

    # Identify the largest non-noise cluster (skip if none)
    cnt <- table(hdb$cluster)
    cnt_no_noise <- cnt[names(cnt) != "0"]
    if (length(cnt_no_noise) == 0) return(list(ok = FALSE))

    largest_id   <- as.integer(names(which.max(cnt_no_noise)))
    members_sub  <- which(hdb$cluster == largest_id)
    # Return sample NAMES (not numeric coercions) to support name-based assignment
    members_orig <- rownames(X_sub)[members_sub]

    # --------- Quality metric: mean silhouette of the extracted cluster --------
    # Compute silhouette on all non-noise points in this pass, then average the
    # silhouette values only for the members of the selected cluster.
    mean_sil <- NA_real_
    if (length(unique(hdb$cluster[hdb$cluster != 0])) >= 2) {
      keep <- which(hdb$cluster != 0)
      if (length(keep) >= 3) {
        cl_full <- as.integer(factor(hdb$cluster[keep]))
        dmat <- stats::dist(emb[keep, , drop = FALSE])
        sil <- cluster::silhouette(cl_full, dmat)

        keep_names <- rownames(emb)[keep]
        members_in_keep <- match(rownames(X_sub)[members_sub], keep_names, nomatch = 0)
        if (any(members_in_keep > 0)) {
          mean_sil <- mean(sil[members_in_keep[members_in_keep > 0], 3])
        }
      }
    }

    # Optional diagnostic plot: colored by HDBSCAN cluster id (0=noise)
    if (plot_iter) {
      df <- data.frame(emb, cluster = factor(hdb$cluster))
      p <- ggplot2::ggplot(df, ggplot2::aes(x = X1, y = X2, color = cluster)) +
        ggplot2::geom_point(size = 1.2, alpha = 0.8) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = sprintf("Pass %d | n=%d", pass_num, length(idx)))
      ggplot2::ggsave(file.path(plot_dir, sprintf("pass_%02d.png", pass_num)), p,
                      width = 6, height = 5)
    }

    list(
      ok = TRUE,
      members = members_orig,  # character vector of sample IDs
      emb = emb,               # 2D UMAP embedding (diagnostic)
      hdb = hdb,               # HDBSCAN result (diagnostic)
      silhouette = mean_sil,   # mean silhouette of extracted cluster
      noise_frac = noise_frac  # noise percentage in this pass
    )
  }

  # -------------------------------- 1) First pass -----------------------------
  # Peel the largest cluster from the full set using a larger neighborhood
  first <- one_pass(seq_len(nrow(PC)), umap_n_first, cur_pass)
  if (!first$ok) {
    logger("[HDBSCAN] First pass failed – returning NA")
    return(factor(rep(NA, nrow(PC)), levels = cluster_labels))
  }
  # Assign pass "1" to its members
  assignments[first$members] <- as.character(cur_pass)
  silhouette_history[cur_pass] <- first$silhouette
  cur_pass <- cur_pass + 1L

  # -------------------- 2) Iteratively peel remaining clusters ----------------
  prev_sil <- first$silhouette
  stable_count <- 0  # counts consecutive passes with Δsilhouette < 5%

  while (cur_pass <= k_desired) {
    # Identify samples not yet assigned
    unass <- which(assignments == "Unassigned")
    # Early stop when remainder is too small to cluster reliably
    if (length(unass) < min_unassigned_stop) {
      logger("[HDBSCAN] Below min_unassigned_stop (%d) – stopping", min_unassigned_stop)
      break
    }

    # Run a pass on the remainder with a smaller neighborhood
    pass <- one_pass(unass, umap_n_next, cur_pass)
    if (!pass$ok) break

    # Assign this pass id to the extracted largest cluster
    assignments[pass$members] <- as.character(cur_pass)
    silhouette_history[cur_pass] <- pass$silhouette

    # Stability criterion: relative change in mean silhouette (guard against NA/0)
    if (!is.na(pass$silhouette) && !is.na(prev_sil) && prev_sil != 0) {
      sil_change <- abs(pass$silhouette - prev_sil) / abs(prev_sil)
      stability_history[cur_pass] <- sil_change
      logger("[HDBSCAN]   Silhouette: %.3f (Δ%.1f%%)", pass$silhouette, sil_change * 100)

      # Count consecutive "stable" improvements; early-stop after 2 passes
      stable_count <- if (sil_change < 0.05) stable_count + 1 else 0
      if (stable_count >= 2) {
        logger("[HDBSCAN] Silhouette stabilized (Δ<5%% for 2 passes) – early stop")
        break
      }
    }

    prev_sil <- pass$silhouette
    cur_pass <- cur_pass + 1L
  }

  # -------------------------- 3) Final cleanup / fallback ---------------------
  # Any remaining Unassigned get forced to the last discovered cluster id (if any),
  # otherwise they are labeled as "Noise".
  leftover <- which(assignments == "Unassigned")
  if (length(leftover)) {
    last_id <- max(suppressWarnings(as.integer(assignments)), na.rm = TRUE)
    if (is.na(last_id) || last_id < 1) {
      assignments[leftover] <- "Noise"
      logger("[HDBSCAN] %d samples → Noise", length(leftover))
    } else {
      assignments[leftover] <- as.character(last_id)
      logger("[HDBSCAN] %d samples → cluster %d", length(leftover), last_id)
    }
  }

  # -------------------------- 4) Map to user labels ---------------------------
  # Internal numeric labels ("1","2",...) are mapped to the user-provided
  # `cluster_labels`. If fewer clusters are discovered than provided, we use
  # a prefix of `cluster_labels`. If more are discovered, we truncate.
  discovered <- sort(setdiff(unique(assignments), c("Unassigned", "Noise")))
  n_disc <- length(discovered)

  if (n_disc == 0) {
    # No valid cluster ever discovered
    return(factor(rep(NA, nrow(PC)), levels = cluster_labels))
  }

  # Build label map: discovered numeric ids → user labels
  map <- setNames(
    cluster_labels[seq_len(min(n_disc, length(cluster_labels)))],
    discovered[seq_len(min(n_disc, length(cluster_labels)))]
  )

  # Translate assignments via the map
  mapped <- map[assignments]

  # Any residual (e.g., "Noise") become the last available label (conservative)
  na_idx <- which(is.na(mapped))
  if (length(na_idx) && length(map)) mapped[na_idx] <- utils::tail(map, 1)

  final_labels <- factor(mapped, levels = cluster_labels)

  # -------------------------- 5) Summary & return -----------------------------
  logger("[HDBSCAN] Final: %d clusters discovered → %d assigned",
         n_disc, sum(!is.na(final_labels)))
  logger("[HDBSCAN] Silhouette history: %s",
         paste(round(silhouette_history, 3), collapse = ", "))

  list(
    labels = final_labels,          # factor (levels = cluster_labels)
    assignments = assignments,      # internal numeric pass ids or "Noise"
    silhouette = silhouette_history,
    stability = stability_history,
    n_clusters = n_disc,
    umap_plots_saved = plot_iter
  )
}
