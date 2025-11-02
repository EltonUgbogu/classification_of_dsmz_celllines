# ----------------------------------------------------------------------
# Supervised labeling using TCGA-BRCA PAM50 centroids
# Saves samples per cluster/subtype using colnames(V_adj) automatically
# ----------------------------------------------------------------------
label_by_tcga <- function(
  clusters,          # per-sample cluster assignments (length == ncol(V_adj))
  V_adj,             # genes x samples (columns = sample IDs)
  features,          # feature genes (e.g., PAM50/HVGs)
  metadata_ann,      # CSV path OR data.frame with sample_barcode & subtype column
  V_tcga,            # genes x TCGA-BRCA samples
  subtype = "PAM50_Subtype",
  outdir = NULL,                 # folder to write CSVs (if NULL, no files written)
  write_by_cluster = TRUE,       # write cluster_<id>_samples.csv
  write_by_subtype = TRUE        # write subtype_<name>_samples.csv
) {
  # Safe logger fallbacks
  if (!exists("log_info", mode = "function")) log_info <- function(fmt, ...) message(sprintf(fmt, ...))
  if (!exists("log_warn", mode = "function")) log_warn <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

  log_info("=== Starting supervised label_by_tcga ===")
  stopifnot(length(clusters) == ncol(V_adj), !is.null(colnames(V_tcga)))

  # ---- sample IDs from expression matrix ----
  if (is.null(colnames(V_adj))) {
    colnames(V_adj) <- paste0("sample_", seq_len(ncol(V_adj)))
    log_warn("V_adj had no column names; generated synthetic sample IDs.")
  }
  sample_ids <- colnames(V_adj)

  # ---- load TCGA annotations ----
  tcga_df <- if (is.character(metadata_ann)) {
    read.csv(metadata_ann, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(metadata_ann)) {
    metadata_ann
  } else stop("metadata_ann must be a CSV path or data.frame")

  stopifnot("sample_barcode" %in% colnames(tcga_df), subtype %in% colnames(tcga_df))

  matched <- tcga_df[tcga_df$sample_barcode %in% colnames(V_tcga), ]
  if (nrow(matched) == 0) stop("No matches between TCGA-BRCA samples and metadata_ann$sample_barcode")

  labs <- setNames(as.character(matched[[subtype]]), matched$sample_barcode)
  V_tcga <- V_tcga[, matched$sample_barcode, drop = FALSE]

  # Drop missing subtype labels
  keep <- !is.na(labs)
  if (sum(!keep) > 0) {
    log_warn("Dropping %d TCGA-BRCA samples with missing %s.", sum(!keep), subtype)
    V_tcga <- V_tcga[, keep, drop = FALSE]
    labs <- labs[keep]
  }
  tcga_lvls <- sort(unique(labs))
  log_info("PAM50 subtypes present: %s", paste(tcga_lvls, collapse = ", "))

  # ---- centroids on shared features ----
  features_use <- intersect(features, intersect(rownames(V_tcga), rownames(V_adj)))
  if (length(features_use) < 40L) stop(sprintf("Too few shared features for centroids: %d", length(features_use)))
  log_info("Computing centroids on %d shared features…", length(features_use))

  tcga_means <- calc_centroids(V_tcga, features_use, unname(labs))
  disc_means <- calc_centroids(V_adj,  features_use, clusters)

  # ---- correlate cluster centroids to subtype centroids ----
  common_g   <- intersect(rownames(tcga_means), rownames(disc_means))
  tcga_means <- tcga_means[common_g, , drop = FALSE]
  disc_means <- disc_means[common_g, , drop = FALSE]

  corr_mat <- stats::cor(disc_means, tcga_means, method = "spearman", use = "pairwise.complete.obs")
  log_info("Correlation matrix: %d clusters x %d subtypes", nrow(corr_mat), ncol(corr_mat))
  print(round(corr_mat, 3))

  # ---- cluster -> best subtype mapping ----
  best_map <- setNames(tcga_lvls[apply(corr_mat, 1, which.max)], rownames(corr_mat))
  map_df   <- data.frame(Cluster = names(best_map), PAM50_Subtype = unname(best_map), row.names = NULL)
  log_info("=== Cluster → PAM50 Mapping ===\n%s", paste(capture.output(print(map_df, row.names = FALSE)), collapse = "\n"))

  # ---- per-sample labeled factor ----
  cl_chr <- as.character(clusters)
  if (!all(unique(cl_chr) %in% names(best_map))) {
    log_warn("Some cluster IDs lack a mapping; resulting labels will be NA for those.")
  }
  labeled <- factor(unname(best_map[cl_chr]), levels = sort(unique(unname(best_map))))

  # ---- CSV exports ----
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    # Master mapping
    master <- data.frame(
      sample_id      = sample_ids,
      cluster_id     = cl_chr,
      mapped_subtype = as.character(labeled),
      stringsAsFactors = FALSE
    )
    utils::write.csv(master, file = file.path(outdir, "tcga_labeling_master.csv"), row.names = FALSE)

    # Per-cluster lists
    if (isTRUE(write_by_cluster)) {
      cl_ids <- sort(unique(cl_chr))
      for (cid in cl_ids) {
        idx <- which(cl_chr == cid); if (!length(idx)) next
        df <- data.frame(sample_id = sample_ids[idx], stringsAsFactors = FALSE)
        utils::write.csv(df,
          file = file.path(outdir, sprintf("cluster_%s_samples.csv", make.names(cid))),
          row.names = FALSE
        )
      }
    }

    # Per-subtype lists
    if (isTRUE(write_by_subtype)) {
      stypes <- sort(unique(as.character(labeled)))
      for (st in stypes) {
        idx <- which(as.character(labeled) == st); if (!length(idx)) next
        df <- data.frame(sample_id = sample_ids[idx], stringsAsFactors = FALSE)
        utils::write.csv(df,
          file = file.path(outdir, sprintf("subtype_%s_samples.csv", make.names(st))),
          row.names = FALSE
        )
      }
    }

    # Correlation matrix for audit
    corr_df <- data.frame(cluster = rownames(corr_mat), corr_mat, check.names = FALSE)
    utils::write.csv(corr_df, file = file.path(outdir, "cluster_to_subtype_correlations.csv"), row.names = FALSE)

    log_info("Wrote CSVs to %s (master, per-cluster, per-subtype, correlations).", outdir)
  }

  return(labeled)
}
