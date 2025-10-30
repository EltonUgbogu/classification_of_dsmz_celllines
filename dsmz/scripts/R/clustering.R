# Clustering: PCA, k-means, hierarchical, dynamic tree cut, HDBSCAN flows
# Place in R/clustering.R

# ---- helpers ---------------------------------------------------------------


# Heuristic: assume rows are observations (cells), columns are features (genes)
# Select top n_hvg features by variance, then run PCA and return PC matrix
 


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

run_hierarchical_on_pcs <- function(PC, k = 5, method = "average") {
  d <- stats::dist(PC, method = "euclidean")
  hc <- stats::hclust(d, method = method)
  clusters <- cutree(hc, k = k)
  list(hc = hc, clusters = clusters, d = d)
}

run_dynamic_tree_cut <- function(PC, deepSplit = 2, pamStage = TRUE) {
  d_pc <- stats::dist(PC, method = "euclidean")
  hc_pc <- stats::hclust(d_pc, method = "average")
  distM_pc <- as.matrix(d_pc)
  minCl_pc <- choose_minCluster(nrow(PC))
  dyn_pc <- dynamicTreeCut::cutreeHybrid(dendro = hc_pc, distM = distM_pc,
                         deepSplit = deepSplit, pamStage = pamStage,
                         minClusterSize = minCl_pc)
  list(dyn = dyn_pc, hc = hc_pc, d = d_pc)
}


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
  hdb1 <- hdbscan(emb1, minPts = minPts1)
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
    hdb2 <- hdbscan(emb2, minPts = minPts2)
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
    hdb3 <- hdbscan(emb3, minPts = minPts3)
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
    km_lum <- kmeans(emb_lum, centers = 2, nstart = 50)
    lum_clusters <- ifelse(km_lum$cluster == 1, "4", "5")
    clusters_hdb[non_cluster3_idx] <- lum_clusters
  }
  
  # Convert cluster assignments to factor with specified levels
  clusters_hdb <- factor(clusters_hdb, levels = cluster_labels)
  # Return the final cluster assignments
  return(clusters_hdb)
}