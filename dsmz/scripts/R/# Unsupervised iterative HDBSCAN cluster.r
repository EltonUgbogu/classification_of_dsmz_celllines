# Unsupervised iterative HDBSCAN clustering utilities
# Requires: dbscan, uwot
# Optional: helpers.R (choose_minPts)

run_iterative_hdbscan <- function(
  PC,
  cluster_labels,
  seed = 42,
  metric = "cosine",
  min_unassigned_stop = 10,   # stop iterating if fewer than this many are left
  umap_neighbors_first = 20,
  umap_neighbors_next  = 15,
  umap_min_dist = 0.3,
  # --- new ---
  outdir = NULL,              # if not NULL, CSVs are written here
  metadata = NULL,            # optional data.frame with rownames = sample IDs
  extra_cols = NULL           # optional character vector of metadata column names to include
) {
  # ---------- Safe logging fallbacks ----------
  if (!exists("log_info", mode = "function"))  log_info <- function(fmt, ...) message(sprintf(fmt, ...))
  if (!exists("log_warn", mode = "function"))  log_warn <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

  # ---------- Validation ----------
  if (!is.matrix(PC)) stop("PC must be a matrix.")
  if (nrow(PC) < 5) stop("PC must have >= 5 samples (rows).")
  if (!is.character(cluster_labels) || length(cluster_labels) < 2) {
    stop("cluster_labels must be a character vector with >= 2 labels.")
  }

  if (is.null(rownames(PC))) rownames(PC) <- paste0("sample_", seq_len(nrow(PC)))
  if (is.null(colnames(PC))) colnames(PC) <- paste0("PC", seq_len(ncol(PC)))

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  # Optional metadata checks
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) stop("metadata must be a data.frame if provided.")
    if (is.null(rownames(metadata))) stop("metadata must have rownames that match PC sample IDs.")
    # Reduce to intersecting samples to avoid accidental reordering
    common_ids <- intersect(rownames(PC), rownames(metadata))
    if (length(common_ids) != nrow(PC)) {
      log_warn("[HDBSCAN] metadata rows do not fully match PC rows; using intersection (%d of %d). Non-matching rows will be NA.", length(common_ids), nrow(PC))
    }
  }

  set.seed(seed)

  # ---------- choose_minPts fallback ----------
  choose_minPts_safe <- if (exists("choose_minPts", mode = "function")) {
    get("choose_minPts", mode = "function")
  } else {
    function(n) max(2L, as.integer(round(log2(max(2, n)) + 3)))
  }

  n <- nrow(PC)
  k_desired <- length(cluster_labels)
  clusters_hdb <- rep("Unassigned", n)
  current_cluster <- 1

  # --- diagnostics to store per-sample ---
  umap1 <- rep(NA_real_, n)
  umap2 <- rep(NA_real_, n)
  hdbscan_id <- rep(NA_integer_, n)
  membership_prob <- rep(NA_real_, n)
  outlier_score <- rep(NA_real_, n)
  iter_assigned <- rep(NA_integer_, n)

  # ---------- Helper to run one HDBSCAN pass over a subset ----------
  run_one_pass <- function(idx, neighb) {
    log_info("[HDBSCAN] Running clustering #%d on %d unassigned samples",
             current_cluster, length(idx))
    emb <- uwot::umap(
      PC[idx, , drop = FALSE],
      n_neighbors = neighb,
      min_dist = umap_min_dist,
      metric = metric,
      verbose = TRUE
    )
    mp <- choose_minPts_safe(nrow(emb))
    mp <- max(2L, mp)
    log_info("[HDBSCAN] Using minPts = %d for n = %d samples", mp, nrow(emb))
    hdb <- dbscan::hdbscan(emb, minPts = mp)

    if (!length(hdb$cluster)) return(list(ok = FALSE))
    cluster_counts <- table(hdb$cluster)
    cluster_counts_no_noise <- cluster_counts[names(cluster_counts) != "0"]
    if (length(cluster_counts_no_noise) == 0) return(list(ok = FALSE))

    largest_id_str <- names(which.max(cluster_counts_no_noise))
    largest_id <- as.integer(largest_id_str)

    # positions within the subset and original indices for the chosen largest cluster
    in_subset_pos <- which(hdb$cluster == largest_id)
    original_idx <- idx[in_subset_pos]

    list(
      ok = TRUE,
      emb = emb,
      hdb = hdb,
      in_subset_pos = in_subset_pos,
      members = original_idx,
      chosen_id = largest_id
    )
  }

  # ---------- First pass ----------
  unassigned_idx <- which(clusters_hdb == "Unassigned")
  first <- run_one_pass(unassigned_idx, umap_neighbors_first)

  if (!first$ok) {
    log_warn("[HDBSCAN] No valid clusters found in first iteration. Returning NA factor.")
    res <- factor(rep(NA_character_, nrow(PC)), levels = cluster_labels)
    # Optionally still write empty diagnostics if outdir set
    if (!is.null(outdir)) {
      write.csv(data.frame(sample_id = rownames(PC), cluster_label = NA_character_),
                file = file.path(outdir, "clusters_all_samples.csv"),
                row.names = FALSE)
    }
    return(res)
  }

  clusters_hdb[first$members] <- as.character(current_cluster)
  # Fill diagnostics for assigned members
  umap1[first$members] <- first$emb[first$in_subset_pos, 1]
  umap2[first$members] <- first$emb[first$in_subset_pos, 2]
  hdbscan_id[first$members] <- first$hdb$cluster[first$in_subset_pos]
  membership_prob[first$members] <- first$hdb$membership_prob[first$in_subset_pos]
  outlier_score[first$members] <- first$hdb$outlier_scores[first$in_subset_pos]
  iter_assigned[first$members] <- current_cluster

  current_cluster <- current_cluster + 1

  # ---------- Iterative passes ----------
  while (current_cluster <= k_desired) {
    unassigned_idx <- which(clusters_hdb == "Unassigned")

    if (length(unassigned_idx) < min_unassigned_stop) {
      log_warn("[HDBSCAN] Too few unassigned samples (%d). Assigning remaining to cluster %d",
               length(unassigned_idx), current_cluster)
      if (length(unassigned_idx) > 0) {
        clusters_hdb[unassigned_idx] <- as.character(current_cluster)
        # No UMAP/HDBSCAN diagnostics for these (left as NA)
        iter_assigned[unassigned_idx] <- current_cluster
      }
      break
    }

    pass <- run_one_pass(unassigned_idx, umap_neighbors_next)

    if (!pass$ok) {
      log_warn("[HDBSCAN] No valid clusters found in this iteration. Assigning remaining to cluster %d",
               current_cluster)
      if (length(unassigned_idx) > 0) {
        clusters_hdb[unassigned_idx] <- as.character(current_cluster)
        iter_assigned[unassigned_idx] <- current_cluster
      }
      break
    }

    clusters_hdb[pass$members] <- as.character(current_cluster)
    umap1[pass$members] <- pass$emb[pass$in_subset_pos, 1]
    umap2[pass$members] <- pass$emb[pass$in_subset_pos, 2]
    hdbscan_id[pass$members] <- pass$hdb$cluster[pass$in_subset_pos]
    membership_prob[pass$members] <- pass$hdb$membership_prob[pass$in_subset_pos]
    outlier_score[pass$members] <- pass$hdb$outlier_scores[pass$in_subset_pos]
    iter_assigned[pass$members] <- current_cluster

    current_cluster <- current_cluster + 1
  }

  # ---------- Post: handle any leftover unassigned ----------
  unassigned_idx <- which(clusters_hdb == "Unassigned")
  if (length(unassigned_idx) > 0) {
    last_assigned_id <- current_cluster - 1
    if (last_assigned_id >= 1) {
      log_warn("[HDBSCAN] %d samples remain unassigned. Assigning to cluster %s",
               length(unassigned_idx), as.character(last_assigned_id))
      clusters_hdb[unassigned_idx] <- as.character(last_assigned_id)
      iter_assigned[unassigned_idx] <- last_assigned_id
    } else {
      log_warn("[HDBSCAN] No clusters assigned at all; marking %d samples as 'Noise'",
               length(unassigned_idx))
      clusters_hdb[unassigned_idx] <- "Noise"
    }
  }

  if (length(clusters_hdb) != nrow(PC)) {
    stop("Internal error: clusters_hdb length must equal nrow(PC).")
  }

  # ---------- Map discovered numeric labels ("1","2",...) to provided cluster_labels ----------
  discovered <- sort(unique(clusters_hdb))
  discovered <- setdiff(discovered, c("Unassigned", "Noise"))
  n_discovered <- length(discovered)

  if (n_discovered == 0) {
    log_warn("[HDBSCAN] No valid clusters discovered. Returning all NA.")
    res <- factor(rep(NA_character_, nrow(PC)), levels = cluster_labels)
    if (!is.null(outdir)) {
      write.csv(data.frame(sample_id = rownames(PC), cluster_label = NA_character_),
                file = file.path(outdir, "clusters_all_samples.csv"),
                row.names = FALSE)
    }
    return(res)
  }

  if (n_discovered > k_desired) {
    log_warn("[HDBSCAN] Discovered %d clusters but only %d labels provided. Using first %d labels.",
             n_discovered, k_desired, k_desired)
    numeric_labels <- as.character(seq_len(n_discovered))
    label_map <- setNames(cluster_labels[seq_len(k_desired)], numeric_labels[seq_len(k_desired)])
  } else if (n_discovered < k_desired) {
    log_warn("[HDBSCAN] Discovered %d clusters but %d labels provided. Using first %d labels.",
             n_discovered, k_desired, n_discovered)
    numeric_labels <- as.character(seq_len(n_discovered))
    label_map <- setNames(cluster_labels[seq_len(n_discovered)], numeric_labels)
  } else {
    numeric_labels <- as.character(seq_len(n_discovered))
    label_map <- setNames(cluster_labels, numeric_labels)
  }

  clusters_mapped <- label_map[clusters_hdb]

  # Any residual "Noise"/"Unassigned" -> last available label (if any)
  residual_idx <- which(is.na(clusters_mapped))
  if (length(residual_idx) > 0 && length(label_map) > 0) {
    log_warn("[HDBSCAN] %d samples were 'Noise'/'Unassigned' after mapping. Assigning to last label.", length(residual_idx))
    clusters_mapped[residual_idx] <- label_map[[length(label_map)]]
  }

  # ---------- Build outputs and write CSVs if requested ----------
  sample_id <- rownames(PC)

  df <- data.frame(
    sample_id = sample_id,
    cluster_numeric = clusters_hdb,
    cluster_label = clusters_mapped,
    umap1 = umap1,
    umap2 = umap2,
    hdbscan_cluster_id = hdbscan_id,
    membership_prob = membership_prob,
    outlier_score = outlier_score,
    iteration = iter_assigned,
    stringsAsFactors = FALSE
  )

  # Attach selected metadata columns (if provided)
  if (!is.null(metadata)) {
    # default to all columns if extra_cols not given
    cols_to_take <- if (is.null(extra_cols)) colnames(metadata) else intersect(extra_cols, colnames(metadata))
    if (length(cols_to_take) > 0) {
      meta_aligned <- metadata[match(sample_id, rownames(metadata)), cols_to_take, drop = FALSE]
      df <- cbind(df, meta_aligned)
    } else if (!is.null(extra_cols)) {
      log_warn("[HDBSCAN] None of extra_cols found in metadata; skipping metadata columns.")
    }
  }

  if (!is.null(outdir)) {
    # 1) master CSV
    write.csv(df, file = file.path(outdir, "clusters_all_samples.csv"), row.names = FALSE)

    # 2) per-cluster CSVs (samples only, plus diagnostics + selected metadata)
    ulabels <- unique(df$cluster_label)
    ulabels <- ulabels[!is.na(ulabels)]
    for (lab in ulabels) {
      sub <- df[df$cluster_label == lab, , drop = FALSE]
      # A concise per-cluster file: sample_id first, keep diagnostics and selected metadata
      fname <- file.path(outdir, sprintf("cluster_%s_samples.csv", make.names(as.character(lab))))
      write.csv(sub, file = fname, row.names = FALSE)
    }

    # 3) small summary
    summary_df <- as.data.frame(sort(table(df$cluster_label), decreasing = TRUE))
    colnames(summary_df) <- c("cluster_label", "n_samples")
    write.csv(summary_df, file = file.path(outdir, "cluster_sizes.csv"), row.names = FALSE)

    log_info("[HDBSCAN] Wrote: clusters_all_samples.csv, cluster_*.csv, and cluster_sizes.csv to %s", outdir)
  }

  factor(clusters_mapped, levels = cluster_labels)
}
