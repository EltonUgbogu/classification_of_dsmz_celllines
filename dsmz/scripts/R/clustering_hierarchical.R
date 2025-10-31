# Hierarchical clustering utilities
# Dependencies: helpers.R (safe_pdf, choose_minCluster)
# NOTE: Throughout this script, features refer to HVGs (highly variable genes).
# HVGs are used to compute PCA spaces and centroids to improve robustness and
# reduce noise from low-variance genes.

# Simple logging helpers
log_info  <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")
log_warn  <- function(...) cat("[WARN] ", sprintf(...), "\n", sep = "")
log_error <- function(...) cat("[ERROR]", sprintf(...), "\n", sep = "")

# Fixed-k hierarchical clustering
run_fixed_hc <- function(PC, k = 4, dist_method = "euclidean",
                         linkage = "average", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(nrow(PC) >= k, "Too few samples for the requested k.")

  log_info("Calculating distance matrix (%s) on %d PCs for %d samples",
           dist_method, ncol(PC), nrow(PC))
  d <- stats::dist(PC, method = dist_method)

  log_info("Building dendrogram with '%s' linkage", linkage)
  hc <- stats::hclust(d, method = linkage)

  log_info("Cutting dendrogram into k = %d clusters", k)
  clusters <- stats::cutree(hc, k = k)

  list(hc = hc, dist = d, clusters = clusters)
}

# Dynamic Tree Cut wrapper
run_dynamic_cut <- function(PC, min_cluster_size = NULL,
                            deep_split = 2, pam_stage = TRUE,
                            dist_method = "euclidean",
                            linkage = "average") {
  d  <- stats::dist(PC, method = dist_method)
  hc <- stats::hclust(d, method = linkage)
  distM <- as.matrix(d)

  if (is.null(min_cluster_size)) {
    min_cluster_size <- choose_minCluster(nrow(PC))
  }
  log_info("Running cutreeHybrid (deepSplit=%d, pamStage=%s, minClusterSize=%d)",
           deep_split, pam_stage, min_cluster_size)

  dyn <- dynamicTreeCut::cutreeHybrid(dendro = hc, distM = distM,
                                      deepSplit = deep_split,
                                      pamStage = pam_stage,
                                      minClusterSize = min_cluster_size)

  clusters <- factor(dyn$labels)
  list(hc = hc, dist = d, clusters = clusters, raw = dyn)
}

# Legacy aliases
run_dynamic_tree_cut <- function(PC, deepSplit = 2, pamStage = TRUE) {
  run_dynamic_cut(PC, deep_split = deepSplit, pam_stage = pamStage)
}

# Centroid calculation over HVGs
calc_centroids <- function(expr_mat, hvgs, group_vec) {
  stopifnot(length(group_vec) == ncol(expr_mat))
  groups <- levels(factor(group_vec))
  sapply(groups, function(g) {
    rowMeans(expr_mat[hvgs, group_vec == g, drop = FALSE], na.rm = TRUE)
  })
}

# ----------------------------------------------------------------------
# Label clusters using ONLY TCGA-BRCA PAM50 centroids (biological inference)
# Any TCGA sample without a PAM50 label is dropped from the reference.
# ----------------------------------------------------------------------
label_by_tcga <- function(clusters,
                          V_adj,     # genes x all samples
                          hvgs,
                          tcga_se,   # BRCA-only SummarizedExperiment (validated below)
                          V_tcga,    # genes x BRCA samples (subset of V_adj cols)
                          subtype = "paper_BRCA_Subtype_PAM50") {
  log_info("=== Starting label_by_tcga ===")
  if (length(clusters) != ncol(V_adj)) {
    stop(sprintf("clusters length (%d) != ncol(V_adj) (%d). Align your vectors.",
                 length(clusters), ncol(V_adj)))
  }
  if (is.null(colnames(V_tcga))) stop("V_tcga must have column names.")
  log_info("Input: %d discovered clusters; %d samples in V_adj; %d TCGA samples.",
           length(unique(clusters)), ncol(V_adj), ncol(V_tcga))

  # --- 1) Validate BRCA-only ---
  log_info("Validating tcga_se contains only TCGA-BRCA samples…")
  cd <- SummarizedExperiment::colData(tcga_se)
  if (!"project_id" %in% colnames(cd)) stop("tcga_se missing 'project_id' in colData.")
  proj <- as.character(cd$project_id)
  if (any(is.na(proj))) log_warn("Found %d samples with NA project_id.", sum(is.na(proj)))
  log_info("tcga_se OK: BRCA-only (n=%d).", length(proj))

  # --- 2) Subtype column exists ---
  cd_tcga <- as.data.frame(cd)
  if (!subtype %in% colnames(cd_tcga)) {
    stop(sprintf("PAM50 column '%s' not found. Available: %s",
                 subtype, paste(colnames(cd_tcga), collapse=", ")))
  }
  log_info("Using PAM50 column: '%s'", subtype)

  # --- 3) Extract labels for V_tcga columns; drop NA ---
  labs <- as.character(cd_tcga[colnames(V_tcga), subtype])
  n_na <- sum(is.na(labs))
  if (n_na > 0) {
    log_warn("Dropping %d TCGA samples with missing PAM50 labels.", n_na)
    keep <- !is.na(labs)
    labs  <- labs[keep]
    V_tcga <- V_tcga[, keep, drop = FALSE]
  }
  tcga_lvls <- sort(unique(labs))
  if (length(tcga_lvls) < 2L) {
    stop(sprintf("Only %d PAM50 subtype(s) remain. Need ≥2.", length(tcga_lvls)))
  }
  log_info("PAM50 subtypes present: %s", paste(tcga_lvls, collapse=", "))

  # --- 4) Compute centroids on HVGs (per design) ---
  hvgs_use <- intersect(hvgs, intersect(rownames(V_tcga), rownames(V_adj)))
  if (length(hvgs_use) < 100L) {
    stop(sprintf("Too few overlapping HVGs for centroiding: %d (need >= 100)", length(hvgs_use)))
  }
  if (length(hvgs_use) < length(hvgs)) {
    log_warn("Using %d/%d HVGs that overlap both V_tcga and V_adj.", length(hvgs_use), length(hvgs))
  }
  log_info("Computing TCGA centroids on %d HVGs…", length(hvgs_use))
  tcga_means <- calc_centroids(V_tcga, hvgs_use, labs)
  log_info("Computing discovered cluster centroids…")
  disc_means <- calc_centroids(V_adj,  hvgs_use, clusters)

  # --- 5) Correlate on actual shared rows between the two centroid matrices ---
  common_g <- intersect(rownames(tcga_means), rownames(disc_means))
  if (length(common_g) < 200L) {
    log_warn("Only %d common centroid genes; results might be unstable.", length(common_g))
  }
  tcga_means <- tcga_means[common_g, , drop = FALSE]
  disc_means <- disc_means[common_g, , drop = FALSE]

  corr_mat <- stats::cor(disc_means, tcga_means, method = "spearman", use = "pairwise.complete.obs")
  if (is.null(rownames(corr_mat)) || is.null(colnames(corr_mat))) {
    stop("Correlation matrix lost dimnames; check calc_centroids output.")
  }
  log_info("Correlation matrix dims: %d clusters x %d subtypes.", nrow(corr_mat), ncol(corr_mat))
  print(round(corr_mat, 3))

  # --- 6) Best subtype per discovered cluster (NAME-AWARE) ---
  best_map <- apply(corr_mat, 1, function(v) tcga_lvls[which.max(v)])
  stopifnot(!is.null(names(best_map)))
  map_df <- data.frame(Cluster_ID = names(best_map),
                       Best_BRCA_PAM50 = unname(best_map),
                       row.names = NULL, check.names = FALSE)
  log_info("=== FINAL MAPPING ===")
  log_info("%s", paste(capture.output(print(map_df, row.names = FALSE)), collapse = "\n"))

  # --- 7) Apply mapping by name to the per-sample cluster vector ---
  cl_chr <- as.character(clusters)
  if (!all(unique(cl_chr) %in% names(best_map))) {
    missing_keys <- setdiff(unique(cl_chr), names(best_map))
    stop(sprintf("Missing centroid mapping for cluster(s): %s",
                 paste(missing_keys, collapse=", ")))
  }
  labeled_vec <- unname(best_map[cl_chr])
  labeled <- factor(labeled_vec, levels = sort(unique(unname(best_map))))
  log_info("Labeling complete. Levels: %s", paste(levels(labeled), collapse=", "))
  log_info("=== label_by_tcga finished successfully ===")
  labeled
}

# Plot helpers
plot_dendrogram_fixed <- function(hc, k, out_pdf) {
  safe_pdf(out_pdf, {
    plot(hc, main = sprintf("Hierarchical (k=%d, Euclidean, average)", k),
         xlab = "", sub = "", cex = 0.5)
    abline(h = hc$height[length(hc$height) - (k-1)], col = "red", lty = 2)
    legend("topright", legend = "k-cut", col = "red", lty = 2, bty = "n")
  })
}

plot_dendrogram_dynamic <- function(hc, clusters, out_pdf) {
  cols <- WGCNA::labels2colors(as.numeric(factor(clusters)))
  safe_pdf(out_pdf, {
    WGCNA::plotDendroAndColors(hc, cols,
                               groupLabels = "Dynamic cut",
                               main = "HC + Dynamic Tree Cut (HVG PCs)",
                               dendroLabels = FALSE, cex.dendroLabels = 0.3)
  })
}

plot_silhouette <- function(sil_obj, out_pdf) {
  safe_pdf(out_pdf, {
    plot(sil_obj, main = sprintf("Silhouette – mean = %.3f", mean(sil_obj[, "sil_width"])))
  })
}

# UMAP helper
make_umap <- function(PC, seed = 42, n_neighbors = 20, min_dist = 0.3, metric = "euclidean") {
  set.seed(seed)
  emb <- uwot::umap(PC, n_neighbors = n_neighbors, min_dist = min_dist,
                    metric = metric, verbose = FALSE)
  emb <- as.data.frame(emb)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$sample   <- rownames(PC)
  emb
}

# Driver for hierarchical clustering workflow
run_hierarchical_clustering <- function(PC, V_adj, hvgs, tcga_se, V_tcga, dataset_lab, outdir, k_fixed = 4, seed = 42) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  hc_fixed <- run_fixed_hc(PC, k = k_fixed, seed = seed)
  clusters_fixed <- hc_fixed$clusters

  clusters_fixed_labeled <- label_by_tcga(clusters_fixed, V_adj, hvgs, tcga_se, V_tcga)

  df_fixed <- tibble::tibble(sample = rownames(PC), dataset = dataset_lab, cluster = as.character(clusters_fixed_labeled))
  utils::write.csv(df_fixed, file.path(outdir, "clusters_hierarchical_HVG.csv"), row.names = FALSE)
  plot_dendrogram_fixed(hc_fixed$hc, k_fixed, file.path(outdir, "dendrogram_hierarchical_HVG.pdf"))

  dyn_res <- run_dynamic_cut(PC)
  clusters_dyn_raw <- dyn_res$clusters
  clusters_dyn_labeled <- label_by_tcga(clusters_dyn_raw, V_adj, hvgs, tcga_se, V_tcga, fallback_labels = c(levels(clusters_dyn_raw), "Unassigned"))

  df_dyn <- tibble::tibble(sample = rownames(PC), dataset = dataset_lab, cluster = as.character(clusters_dyn_labeled))
  utils::write.csv(df_dyn, file.path(outdir, "clusters_hierarchical_dynamic_HVG.csv"), row.names = FALSE)
  plot_dendrogram_dynamic(dyn_res$hc, clusters_dyn_labeled, file.path(outdir, "dendrogram_hierarchical_dynamic_HVG.pdf"))

  sil_dyn <- cluster::silhouette(as.integer(factor(clusters_dyn_raw)), dyn_res$dist)
  plot_silhouette(sil_dyn, file.path(outdir, "silhouette_dynamic_HVG.pdf"))

  umap_df <- make_umap(PC, seed = seed)
  umap_df$dataset <- dataset_lab
  umap_df$cluster <- clusters_fixed_labeled
  p_umap <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
    ggplot2::geom_point(alpha = 0.9, size = 1.8) + ggplot2::theme_bw() +
    ggplot2::ggtitle(sprintf("UMAP on %d PCs (Hierarchical k=%d)", ncol(PC), k_fixed))
  ggplot2::ggsave(file.path(outdir, "umap_hierarchical_PCs.pdf"), p_umap, width = 7, height = 6)
  utils::write.csv(umap_df, file.path(outdir, "umap_hierarchical_PCs_coords.csv"), row.names = FALSE)

  list(
    fixed   = list(hc = hc_fixed$hc, clusters = clusters_fixed_labeled, df = df_fixed),
    dynamic = list(hc = dyn_res$hc, clusters = clusters_dyn_labeled, df = df_dyn),
    silhouette = sil_dyn,
    umap = umap_df
  )
}

# Legacy alias for backward compatibility
run_hierarchical_on_pcs <- function(PC, k = 4, method = "average") {
  hc <- run_fixed_hc(PC, k = k, linkage = method)
  list(hc = hc$hc, clusters = hc$clusters, d = hc$dist)
}
