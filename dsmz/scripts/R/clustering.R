# Clustering: PCA, k-means, hierarchical, dynamic tree cut, HDBSCAN flows
# Place in R/clustering.R

# ---- helpers ---------------------------------------------------------------


# Heuristic: assume rows are observations (cells), columns are features (genes)
# Select top n_hvg features by variance, then run PCA and return PC matrix
make_pcs <- function(V_adj, n_hvg = 3000, max_pc = 30, center = TRUE, scale. = TRUE) {
  if (!is.matrix(V_adj) && !is.data.frame(V_adj)) {
    stop("V_adj must be a matrix or data.frame")
  }
  mat <- as.matrix(V_adj)
  if (ncol(mat) < 2 || nrow(mat) < 2) stop("V_adj must have at least 2 rows and 2 columns")

  # pick highly variable columns (features)
  feature_sd <- apply(mat, 2, stats::sd, na.rm = TRUE)
  feature_rank <- order(feature_sd, decreasing = TRUE, na.last = NA)
  keep <- head(feature_rank, n = min(n_hvg, ncol(mat)))
  mat_hvg <- mat[, keep, drop = FALSE]

  pca <- stats::prcomp(mat_hvg, center = center, scale. = scale.)
  pcs <- pca$x[, seq_len(min(max_pc, ncol(pca$x))), drop = FALSE]
  rownames(pcs) <- rownames(mat)
  pcs
}

choose_minCluster <- function(n_obs) {
  # Reasonable lower bound for dynamic tree cut; clamp to [10, max(10, n/20)]
  max(10L, as.integer(round(n_obs / 20)))
}


run_kmeans_on_pcs <- function(PC, outdir, k_grid = 2:8, seed = 42) {
  set.seed(seed)
  res_list <- list()
  metrics <- lapply(k_grid, function(k) {
    km <- kmeans(PC, centers = k, nstart = 50, iter.max = 100)
    wcss <- sum(km$withinss)
    sil <- cluster::silhouette(km$cluster, stats::dist(PC))
    sil_mean <- mean(sil[, "sil_width"])
    ch <- tryCatch(fpc::calinhara(PC, km$cluster), error = function(e) NA_real_)
    db <- tryCatch({
      m <- clusterCrit::intCriteria(as.matrix(PC), as.integer(km$cluster), "davies_bouldin")
      as.numeric(m$davies_bouldin)
    }, error = function(e) NA_real_)
    res_list[[as.character(k)]] <<- list(km = km, sil = sil)
    data.frame(k = k, WCSS = wcss, Silhouette = sil_mean, CH = ch, DB = db)
  })
  metrics <- do.call(rbind, metrics)
  ensure_dir(outdir)
  utils::write.csv(metrics, file.path(outdir, "kmeans_metrics_HVG.csv"), row.names = FALSE)
  best_k <- with(metrics, {
    cand <- k[Silhouette == max(Silhouette, na.rm = TRUE)]
    if (length(cand) > 1) cand[which.max(metrics$CH[match(cand, k)])] else cand
  })
  list(best_k = best_k, res_list = res_list, metrics = metrics)
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