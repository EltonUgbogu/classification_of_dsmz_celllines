# Fixed-k hierarchical clustering with per-cluster CSV export
run_fixed_hc <- function(
  PC,
  k = 4,
  dist_method = "euclidean",
  linkage = "average",
  seed = NULL,
  outdir = NULL,            # <- set a folder to write CSVs
  barcode_colname = "sample_id"  # column name used in CSVs
) {
  # Safe logger fallbacks
  if (!exists("log_info", mode = "function")) {
    log_info <- function(fmt, ...) message(sprintf(fmt, ...))
  }

  if (!is.null(seed)) set.seed(seed)
  stopifnot(nrow(PC) >= k)

  # Ensure sample names exist
  if (is.null(rownames(PC))) {
    rownames(PC) <- paste0("sample_", seq_len(nrow(PC)))
  }

  log_info("Calculating distance matrix (%s) on %d PCs for %d samples",
           dist_method, ncol(PC), nrow(PC))
  d <- stats::dist(PC, method = dist_method)

  log_info("Building dendrogram with '%s' linkage", linkage)
  hc <- stats::hclust(d, method = linkage)

  log_info("Cutting dendrogram into k = %d clusters", k)
  clusters <- stats::cutree(hc, k = k)

  # ===== Write CSVs (one per cluster) =====
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    # master mapping for reference
    mapping <- data.frame(
      !!barcode_colname := rownames(PC),
      cluster = as.integer(clusters),
      stringsAsFactors = FALSE
    )
    utils::write.csv(mapping,
      file = file.path(outdir, "clusters_all_samples.csv"),
      row.names = FALSE
    )

    # per-cluster sample lists
    cl_ids <- sort(unique(as.integer(clusters)))
    for (cid in cl_ids) {
      idx <- which(clusters == cid)
      df <- data.frame(
        !!barcode_colname := rownames(PC)[idx],
        stringsAsFactors = FALSE
      )
      utils::write.csv(
        df,
        file = file.path(outdir, sprintf("cluster_%d_samples.csv", cid)),
        row.names = FALSE
      )
    }

    # quick size summary (optional but handy)
    sizes <- as.data.frame(sort(table(as.integer(clusters))))
    colnames(sizes) <- c("cluster", "n_samples")
    utils::write.csv(sizes,
      file = file.path(outdir, "cluster_sizes.csv"),
      row.names = FALSE
    )

    log_info("Wrote %d per-cluster CSVs to: %s", length(cl_ids), outdir)
  }

  list(hc = hc, dist = d, clusters = clusters)
}
