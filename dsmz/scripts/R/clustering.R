# Clustering: PCA, k-means, hierarchical, dynamic tree cut, HDBSCAN flows
# Place in R/clustering.R
# Dependencies: helpers.R (ensure_dir, safe_pdf, choose_minCluster, choose_minPts, get_centroid, get_k_nearest, aggregate_median)

# ============================================================================
# Helper / Utility Functions
# ============================================================================

# ----------------------------------------------------------------------
# 1.1  Logging helpers (keeps the main code tidy)
# ----------------------------------------------------------------------
log_info  <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")
log_warn  <- function(...) cat("[WARN] ", sprintf(...), "\n", sep = "")
log_error <- function(...) cat("[ERROR]", sprintf(...), "\n", sep = "")

# ============================================================================
# K-means Clustering Functions
# ============================================================================

# Performs k-means clustering on principal components with quality control and visualization
run_kmeans_on_pcs <- function(PC, outdir, k_grid = 2:8, seed = 42, 
                              sample_is_tcga = NULL, V_adj = NULL) {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Validate inputs
  stopifnot(is.matrix(PC), all(is.finite(PC)), is.character(outdir))
  if (!is.null(V_adj)) {
    stopifnot(is.matrix(V_adj), all(is.finite(V_adj)), ncol(V_adj) == nrow(PC))
    # Log dimensions and variance check for adjusted matrix
    cat("V_adj dims:", nrow(V_adj), "genes x", ncol(V_adj), "samples\n")
    nzv <- sum(matrixStats::rowVars(V_adj) > 0)
    cat("Genes with non-zero variance:", nzv, "\n")
  }
  
  # Warn for large sample sets
  if (nrow(PC) > 5000) {
    cat("[WARN] Large sample set (%d) for k-means; consider subsampling.\n", nrow(PC))
  }
  
  # Ensure output directory exists
  ensure_dir(outdir)
  
  # Initialize list to store k-means results
  res_list <- list()
  
  # Compute clustering metrics for each k in k_grid
  metrics <- lapply(k_grid, function(k) {
    # Run k-means with specified centers and iterations
    km <- kmeans(PC, centers = k, nstart = 50, iter.max = 100)
    wcss <- sum(km$withinss)  # Within-cluster sum of squares
    # Compute silhouette scores
    sil <- cluster::silhouette(km$cluster, stats::dist(PC))
    sil_mean <- mean(sil[, "sil_width"])
    # Compute Calinski-Harabasz index
    ch <- tryCatch(fpc::calinhara(PC, km$cluster), error = function(e) NA_real_)
    # Compute Davies-Bouldin index
    db <- tryCatch({
      m <- clusterCrit::intCriteria(as.matrix(PC), as.integer(km$cluster), 
                                    "davies_bouldin")
      as.numeric(m$davies_bouldin)
    }, error = function(e) NA_real_)
    # Store results for this k
    res_list[[as.character(k)]] <<- list(km = km, sil = sil)
    data.frame(k = k, WCSS = wcss, Silhouette = sil_mean, CH = ch, DB = db)
  })
  metrics <- do.call(rbind, metrics)
  
  # Save clustering metrics
  utils::write.csv(metrics, file.path(outdir, "kmeans_metrics_HVG.csv"), 
                   row.names = FALSE)
  
  # Generate elbow plot for WCSS
  safe_pdf(file.path(outdir, "kmeans_elbow_HVG.pdf"), {
    plot(metrics$k, metrics$WCSS, type = "b", xlab = "k", 
         ylab = "Within-Cluster SS", main = "Elbow (HVG space)")
  })
  
  # Generate panel of clustering metrics
  safe_pdf(file.path(outdir, "kmeans_metrics_panel_HVG.pdf"), {
    par(mfrow = c(2, 2))
    plot(metrics$k, metrics$Silhouette, type = "b", xlab = "k", 
         ylab = "Avg silhouette", main = "Silhouette (higher is better)")
    plot(metrics$k, metrics$CH, type = "b", xlab = "k", 
         ylab = "Calinski-Harabasz", main = "Calinski-Harabasz (higher)")
    plot(metrics$k, metrics$DB, type = "b", xlab = "k", 
         ylab = "Davies-Bouldin", main = "Davies-Bouldin (lower)")
    plot(metrics$k, metrics$WCSS, type = "b", xlab = "k", 
         ylab = "WCSS", main = "Elbow (lower)")
    par(mfrow = c(1, 1))
  })
  
  # Select best k based on silhouette, then CH if tied
  best_k <- with(metrics, {
    cand <- k[Silhouette == max(Silhouette, na.rm = TRUE)]
    if (length(cand) > 1) cand[which.max(metrics$CH[match(cand, k)])] else cand
  })
  message(sprintf("[kmeans] Best k by silhouette→CH = %d", best_k))
  
  # Extract final k-means and silhouette results
  km_final <- res_list[[as.character(best_k)]]$km
  sil_final <- res_list[[as.character(best_k)]]$sil
  
  # Create cluster labels (C1, C2, ...)
  clusters_kmeans <- factor(km_final$cluster, labels = paste0("C", seq_len(best_k)))
  
  # Create dataset labels (TCGA or DSMZ) if sample_is_tcga is provided
  dataset_lab <- if (!is.null(sample_is_tcga)) {
    ifelse(sample_is_tcga, "TCGA", "DSMZ")
  } else {
    rep("Unknown", nrow(PC))
  }
  
  # Create cluster summary data frame
  cluster_df <- tibble::tibble(
    sample = rownames(PC),
    dataset = dataset_lab,
    cluster = as.character(clusters_kmeans)
  )
  utils::write.csv(cluster_df, file.path(outdir, "clusters_kmeans.csv"), 
                   row.names = FALSE)
  
  # Generate silhouette plot for best k
  safe_pdf(file.path(outdir, sprintf("silhouette_k=%d_HVG.pdf", best_k)), {
    plot(sil_final, main = sprintf("Silhouette (k = %d, HVG space)", best_k), 
         cex.names = 0.8)
  })
  
  # Generate silhouette plots for all k in k_grid
  sil_dir <- file.path(outdir, "silhouettes_per_k")
  dir.create(sil_dir, showWarnings = FALSE, recursive = TRUE)
  for (k in k_grid) {
    silk <- res_list[[as.character(k)]]$sil
    if (!is.null(silk)) {
      safe_pdf(file.path(sil_dir, sprintf("silhouette_k=%d.pdf", k)), {
        plot(silk, main = sprintf("Silhouette (k = %d, HVG space)", k), 
             cex.names = 0.7)
      })
    }
  }
  
  # Log final results
  cat("Optimal k (k-means):", best_k, 
      "with mean silhouette:", round(mean(sil_final[, "sil_width"]), 3), "\n")
  print(table(clusters_kmeans, dataset_lab))
  
  # Return results
  return(list(
    best_k = best_k,
    res_list = res_list,
    metrics = metrics,
    clusters = clusters_kmeans,
    cluster_df = cluster_df
  ))
}

# ============================================================================
# Hierarchical Clustering Functions
# ============================================================================

# ----------------------------------------------------------------------
# 2.1  Hierarchical clustering (fixed k)
# ----------------------------------------------------------------------
run_fixed_hc <- function(PC, k = 5, dist_method = "euclidean",
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

# ----------------------------------------------------------------------
# 2.2  Dynamic tree cut (WGCNA::cutreeHybrid)
# ----------------------------------------------------------------------
run_dynamic_cut <- function(PC, min_cluster_size = NULL,
                            deep_split = 2, pam_stage = TRUE,
                            dist_method = "euclidean",
                            linkage = "average") {
  d  <- stats::dist(PC, method = dist_method)
  hc <- stats::hclust(d, method = linkage)
  distM <- as.matrix(d)
  
  # adaptive minClusterSize if not supplied
  if (is.null(min_cluster_size)) {
    min_cluster_size <- choose_minCluster(nrow(PC))
  }
  log_info("Running cutreeHybrid (deepSplit=%d, pamStage=%s, minClusterSize=%d)",
           deep_split, pam_stage, min_cluster_size)
  
  dyn <- dynamicTreeCut::cutreeHybrid(dendro = hc, distM = distM,
                      deepSplit = deep_split,
                      pamStage = pam_stage,
                      minClusterSize = min_cluster_size)
  
  clusters <- factor(dyn$labels)               # 0 = unassigned
  list(hc = hc, dist = d, clusters = clusters, raw = dyn)
}

# Legacy alias for backward compatibility
run_dynamic_tree_cut <- function(PC, deepSplit = 2, pamStage = TRUE) {
  run_dynamic_cut(PC, deep_split = deepSplit, pam_stage = pamStage)
}

# ----------------------------------------------------------------------
# 2.3  Centroid calculation (common HVGs)
# ----------------------------------------------------------------------
calc_centroids <- function(expr_mat, hvgs, group_vec) {
  # expr_mat: genes x samples, hvgs: character vector of gene names
  # group_vec: factor or character with same length as ncol(expr_mat)
  stopifnot(length(group_vec) == ncol(expr_mat))
  
  groups <- levels(factor(group_vec))
  cents <- sapply(groups, function(g) {
    rowMeans(expr_mat[hvgs, group_vec == g, drop = FALSE], na.rm = TRUE)
  })
  cents
}

# ----------------------------------------------------------------------
# 2.4  Label clusters by TCGA PAM50 correlation
# ----------------------------------------------------------------------
label_by_tcga <- function(clusters, V_adj, hvgs,
                          tcga_se, V_tcga,
                          subtype_col_opts = c("PAM50","Subtype",
                                              "BRCA_Subtype","molecular_subtype"),
                          fallback_labels = c("HER2-high","HER2-low",
                                             "LumA","LumB","Basal")) {
  # 1) locate TCGA annotation column
  tcga_col <- subtype_col_opts[subtype_col_opts %in% colnames(SummarizedExperiment::colData(tcga_se))][1]
  if (is.na(tcga_col)) {
    log_warn("No TCGA subtype column found – using fallback labels")
    cl_levels <- levels(factor(clusters))
    n_levels <- length(cl_levels)
    if (n_levels <= length(fallback_labels)) {
      return(factor(clusters, labels = fallback_labels[seq_len(n_levels)]))
    } else {
      return(factor(clusters))
    }
  }
  log_info("Found TCGA annotation column: '%s'", tcga_col)
  
  # 2) TCGA centroids
  tcga_lab   <- as.character(SummarizedExperiment::colData(tcga_se)[colnames(V_tcga), tcga_col])
  tcga_lvls  <- sort(unique(na.omit(tcga_lab)))
  tcga_means <- calc_centroids(V_tcga, hvgs, tcga_lab)
  
  # 3) discovered centroids
  disc_means <- calc_centroids(V_adj, hvgs, clusters)
  
  # 4) common genes & correlation
  common_g <- intersect(rownames(tcga_means), hvgs)
  if (length(common_g) < 500) {
    log_warn("Only %d common genes for correlation – results may be noisy", length(common_g))
  }
  corr_mat <- stats::cor(disc_means[common_g, , drop = FALSE],
                  tcga_means[common_g, , drop = FALSE],
                  method = "spearman")
  log_info("Correlation matrix (rows=clusters, cols=TCGA):")
  print(corr_mat)
  
  # 5) best-match mapping
  best_tcga   <- apply(corr_mat, 1, function(v) tcga_lvls[which.max(v)])
  final_labs  <- fallback_labels[match(best_tcga, tcga_lvls)]
  final_labs[is.na(final_labs)] <- names(best_tcga)[is.na(final_labs)]
  
  # mapping table
  cl_levels <- levels(factor(clusters))
  map_df <- data.frame(Cluster_ID = cl_levels,
                       Best_TCGA = best_tcga[cl_levels],
                       Applied_Label = final_labs[cl_levels],
                       stringsAsFactors = FALSE)
  log_info("Final mapping:\n%s", paste(capture.output(print(map_df, row.names = FALSE)), collapse = "\n"))
  
  factor(clusters, labels = final_labs[cl_levels])
}

# ----------------------------------------------------------------------
# 2.5  Plotting helpers
# ----------------------------------------------------------------------
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
    plot(sil_obj,
         main = sprintf("Silhouette – mean = %.3f", mean(sil_obj[, "sil_width"])))
  })
}

# ----------------------------------------------------------------------
# 2.6  UMAP visualization
# ----------------------------------------------------------------------
make_umap <- function(PC, seed = 42, n_neighbors = 20, min_dist = 0.3,
                      metric = "euclidean") {
  set.seed(seed)
  emb <- uwot::umap(PC, n_neighbors = n_neighbors, min_dist = min_dist,
                    metric = metric, verbose = FALSE)
  emb <- as.data.frame(emb)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$sample   <- rownames(PC)
  emb
}

# ============================================================================
# Driver Function - Hierarchical Clustering Pipeline
# ============================================================================

# One-line call for the whole hierarchical clustering pipeline
run_hierarchical_clustering <- function(PC,               # samples x PCs (HVG-derived)
                                        V_adj,             # genes x samples (full adjusted)
                                        hvgs,              # character vector of HVGs
                                        tcga_se, V_tcga,   # reference TCGA objects
                                        dataset_lab,       # "TCGA" / "DSMZ" vector
                                        outdir,
                                        k_fixed = 5,
                                        seed = 42) {
  ensure_dir(outdir)                     # create outdir if missing
  
  ## ---------- 1) Fixed-k hierarchical clustering ----------
  hc_fixed <- run_fixed_hc(PC, k = k_fixed, seed = seed)
  clusters_fixed <- hc_fixed$clusters
  
  ## ---------- 2) Label with TCGA PAM50 ----------
  clusters_fixed_labeled <- label_by_tcga(clusters_fixed, V_adj, hvgs,
                                          tcga_se, V_tcga)
  
  ## ---------- 3) Save fixed-k results ----------
  df_fixed <- tibble::tibble(sample = rownames(PC),
                     dataset = dataset_lab,
                     cluster = as.character(clusters_fixed_labeled))
  utils::write.csv(df_fixed,
            file.path(outdir, "clusters_hierarchical_HVG.csv"),
            row.names = FALSE)
  plot_dendrogram_fixed(hc_fixed$hc, k_fixed,
                        file.path(outdir, "dendrogram_hierarchical_HVG.pdf"))
  
  ## ---------- 4) Dynamic tree cut ----------
  dyn_res <- run_dynamic_cut(PC)
  clusters_dyn_raw <- dyn_res$clusters
  clusters_dyn_labeled <- label_by_tcga(clusters_dyn_raw, V_adj, hvgs,
                                        tcga_se, V_tcga,
                                        fallback_labels = c(levels(clusters_dyn_raw),
                                                          "Unassigned"))
  
  ## ---------- 5) Save dynamic results ----------
  df_dyn <- tibble::tibble(sample = rownames(PC),
                   dataset = dataset_lab,
                   cluster = as.character(clusters_dyn_labeled))
  utils::write.csv(df_dyn,
            file.path(outdir, "clusters_hierarchical_dynamic_HVG.csv"),
            row.names = FALSE)
  plot_dendrogram_dynamic(dyn_res$hc, clusters_dyn_labeled,
                          file.path(outdir, "dendrogram_hierarchical_dynamic_HVG.pdf"))
  
  ## ---------- 6) Silhouette QC (dynamic) ----------
  sil_dyn <- cluster::silhouette(as.integer(factor(clusters_dyn_raw)),
                                 dyn_res$dist)
  plot_silhouette(sil_dyn,
                  file.path(outdir, "silhouette_dynamic_HVG.pdf"))
  
  ## ---------- 7) UMAP on the same PC space ----------
  umap_df <- make_umap(PC, seed = seed)
  umap_df$dataset <- dataset_lab
  umap_df$cluster <- clusters_fixed_labeled          # use fixed-k labels
  p_umap <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = cluster, shape = dataset)) +
    ggplot2::geom_point(alpha = 0.9, size = 1.8) + ggplot2::theme_bw() +
    ggplot2::ggtitle(sprintf("UMAP on %d PCs (Hierarchical k=%d)", ncol(PC), k_fixed))
  ggplot2::ggsave(file.path(outdir, "umap_hierarchical_PCs.pdf"), p_umap,
         width = 7, height = 6)
  utils::write.csv(umap_df, file.path(outdir, "umap_hierarchical_PCs_coords.csv"),
            row.names = FALSE)
  
  ## ---------- Return everything ----------
  list(
    fixed   = list(hc = hc_fixed$hc, clusters = clusters_fixed_labeled,
                   df = df_fixed),
    dynamic = list(hc = dyn_res$hc, clusters = clusters_dyn_labeled,
                   df = df_dyn),
    silhouette = sil_dyn,
    umap = umap_df
  )
}

# Legacy alias for backward compatibility
run_hierarchical_on_pcs <- function(PC, k = 5, method = "average") {
  hc <- run_fixed_hc(PC, k = k, linkage = method)
  list(hc = hc$hc, clusters = hc$clusters, d = hc$dist)
}

# ============================================================================
# HDBSCAN Functions
# ============================================================================

# Iterative HDBSCAN clustering to assign samples to five clusters
run_iterative_hdbscan <- function(PC, cluster_labels, seed = 42, 
                                 metric = "euclidean") {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # First UMAP embedding for initial clustering
  emb1 <- uwot::umap(PC, n_neighbors = 20, min_dist = 0.3, metric = metric, 
                     verbose = TRUE)
  # Determine minPts for HDBSCAN based on sample size
  minPts1 <- choose_minPts(nrow(emb1))
  # Log minPts and sample size
  cat(sprintf("[HDBSCAN PCs] Using minPts = %d for n = %d samples\n", minPts1, 
              nrow(emb1)))
  # Run HDBSCAN clustering on first embedding
  hdb1 <- dbscan::hdbscan(emb1, minPts = minPts1)
  # Initialize cluster assignments as "Unassigned"
  clusters_hdb <- rep("Unassigned", nrow(PC))
  # Identify largest non-noise cluster (excluding cluster 0) as cluster "1"
  cluster1_idx <- which(hdb1$cluster == 
                       which.max(table(hdb1$cluster[hdb1$cluster != 0])))
  clusters_hdb[cluster1_idx] <- "1"
  
  # Identify samples not assigned to cluster "1"
  non_cluster1_idx <- which(clusters_hdb == "Unassigned")
  if (length(non_cluster1_idx) < 10) {
    # Warn if too few samples remain for further clustering
    cat("[WARN] Too few non-cluster1 samples (%d); skipping next clustering.\n",
        length(non_cluster1_idx))
  } else {
    # Second UMAP embedding for non-cluster1 samples
    emb2 <- uwot::umap(PC[non_cluster1_idx, , drop = FALSE], n_neighbors = 15, 
                       min_dist = 0.3, metric = "euclidean", verbose = TRUE)
    minPts2 <- choose_minPts(nrow(emb2))
    cat(sprintf("[HDBSCAN PCs] Using minPts = %d for n = %d non-cluster1 samples\n",
                minPts2, nrow(emb2)))
    # Run HDBSCAN on second embedding
    hdb2 <- dbscan::hdbscan(emb2, minPts = minPts2)
    # Assign largest non-noise cluster as cluster "2"
    cluster2_idx <- non_cluster1_idx[which(hdb2$cluster == 
                                         which.max(table(hdb2$cluster[hdb2$cluster != 0])))]
    clusters_hdb[cluster2_idx] <- "2"
  }
  
  # Identify samples not assigned to cluster "1" or "2"
  non_cluster2_idx <- which(clusters_hdb == "Unassigned")
  if (length(non_cluster2_idx) < 10) {
    # Warn if too few samples remain
    cat("[WARN] Too few non-cluster2 samples (%d); skipping next clustering.\n",
        length(non_cluster2_idx))
  } else {
    # Third UMAP embedding for remaining samples
    emb3 <- uwot::umap(PC[non_cluster2_idx, , drop = FALSE], n_neighbors = 15, 
                       min_dist = 0.3, metric = "euclidean", verbose = TRUE)
    minPts3 <- choose_minPts(nrow(emb3))
    cat(sprintf("[HDBSCAN PCs] Using minPts = %d for n = %d non-cluster2 samples\n",
                minPts3, nrow(emb3)))
    # Run HDBSCAN on third embedding
    hdb3 <- dbscan::hdbscan(emb3, minPts = minPts3)
    # Assign largest non-noise cluster as cluster "3"
    cluster3_idx <- non_cluster2_idx[which(hdb3$cluster == 
                                         which.max(table(hdb3$cluster[hdb3$cluster != 0])))]
    clusters_hdb[cluster3_idx] <- "3"
  }
  
  # Identify remaining unassigned samples for final clustering
  non_cluster3_idx <- which(clusters_hdb == "Unassigned")
  if (length(non_cluster3_idx) < 2) {
    # Warn if too few samples remain; assign all as cluster "4"
    cat("[WARN] Too few remaining samples (%d); assigning all as cluster 4.\n",
        length(non_cluster3_idx))
    clusters_hdb[non_cluster3_idx] <- "4"
  } else {
    # Final UMAP embedding for remaining samples
    emb_lum <- uwot::umap(PC[non_cluster3_idx, , drop = FALSE], n_neighbors = 10, 
                          min_dist = 0.3, metric = "euclidean", verbose = TRUE)
    # Run k-means clustering to split into clusters "4" and "5"
    km_lum <- stats::kmeans(emb_lum, centers = 2, nstart = 50)
    lum_clusters <- ifelse(km_lum$cluster == 1, "4", "5")
    clusters_hdb[non_cluster3_idx] <- lum_clusters
  }
  
  # Convert cluster assignments to factor with specified levels
  clusters_hdb <- factor(clusters_hdb, levels = cluster_labels)
  # Return the final cluster assignments
  return(clusters_hdb)
}
