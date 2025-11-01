# Dynamic Tree Cut wrapper (Method for detecting clusters in HC dendrogram)
run_dynamic_cut <- function(
  PC,
  min_cluster_size = NULL,
  deep_split = 2,
  pam_stage = TRUE,
  dist_method = "euclidean",
  linkage = "average",
  # --- new ---
  outdir = NULL,                 # folder to write CSVs (if NULL, no files written)
  barcode_colname = "sample_id", # column name for sample IDs in CSVs
  save_unassigned = TRUE         # if TRUE, write a CSV for label 0 (unassigned)
) {
  # Safe logger fallback
  if (!exists("log_info", mode = "function")) {
    log_info <- function(fmt, ...) message(sprintf(fmt, ...))
  }

  # Ensure sample names exist
  if (is.null(rownames(PC))) {
    rownames(PC) <- paste0("sample_", seq_len(nrow(PC)))
  }

  d  <- stats::dist(PC, method = dist_method)
  hc <- stats::hclust(d, method = linkage)
  distM <- as.matrix(d)

  if (is.null(min_cluster_size)) {
    if (exists("choose_minCluster", mode = "function")) {
      min_cluster_size <- choose_minCluster(nrow(PC))
    } else {
      # simple fallback heuristic
      min_cluster_size <- max(5L, floor(sqrt(nrow(PC))))
    }
  }
  log_info("Running cutreeHybrid (deepSplit=%d, pamStage=%s, minClusterSize=%d)",
           deep_split, pam_stage, min_cluster_size)

  dyn <- dynamicTreeCut::cutreeHybrid(
    dendro = hc, distM = distM,
    deepSplit = deep_split,
    pamStage = pam_stage,
    minClusterSize = min_cluster_size
  )

  # dynamicTreeCut labels are integers; 0 often means unassigned
  clusters_int <- as.integer(dyn$labels)
  clusters <- factor(clusters_int)

  # ===== Write CSVs =====
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    # Master mapping (sample -> cluster)
    mapping <- data.frame(
      !!barcode_colname := rownames(PC),
      cluster = clusters_int,
      stringsAsFactors = FALSE
    )
    utils::write.csv(mapping,
      file = file.path(outdir, "clusters_all_samples.csv"),
      row.names = FALSE
    )

    # Cluster sizes summary
    size_tab <- sort(table(clusters_int))
    sizes <- data.frame(cluster = as.integer(names(size_tab)),
                        n_samples = as.integer(size_tab),
                        row.names = NULL)
    utils::write.csv(sizes,
      file = file.path(outdir, "cluster_sizes.csv"),
      row.names = FALSE
    )

    # Per-cluster CSVs (skip 0 unless save_unassigned = TRUE)
    unique_ids <- sort(unique(clusters_int))
    if (!save_unassigned) unique_ids <- setdiff(unique_ids, 0L)

    for (cid in unique_ids) {
      idx <- which(clusters_int == cid)
      if (length(idx) == 0) next

      df <- data.frame(
        !!barcode_colname := rownames(PC)[idx],
        stringsAsFactors = FALSE
      )
      fname <- if (cid == 0L) "cluster_unassigned_samples.csv"
               else sprintf("cluster_%d_samples.csv", cid)
      utils::write.csv(df, file = file.path(outdir, fname), row.names = FALSE)
    }

    log_info("Wrote per-cluster CSVs, clusters_all_samples.csv and cluster_sizes.csv to: %s", outdir)
  }

  list(hc = hc, dist = d, clusters = clusters, raw = dyn)
}
