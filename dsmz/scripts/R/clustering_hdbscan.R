# HDBSCAN clustering utilities
# Dependencies: helpers.R (choose_minPts)

run_iterative_hdbscan <- function(PC, cluster_labels, seed = 42, metric = "euclidean") {
  set.seed(seed)

  emb1 <- uwot::umap(PC, n_neighbors = 20, min_dist = 0.3, metric = metric, verbose = TRUE)
  minPts1 <- choose_minPts(nrow(emb1))
  cat(sprintf("[HDBSCAN PCs] Using minPts = %d for n = %d samples\n", minPts1, nrow(emb1)))
  hdb1 <- dbscan::hdbscan(emb1, minPts = minPts1)
  clusters_hdb <- rep("Unassigned", nrow(PC))
  cluster1_idx <- which(hdb1$cluster == which.max(table(hdb1$cluster[hdb1$cluster != 0])))
  clusters_hdb[cluster1_idx] <- "1"

  non_cluster1_idx <- which(clusters_hdb == "Unassigned")
  if (length(non_cluster1_idx) < 10) {
    cat("[WARN] Too few non-cluster1 samples (", length(non_cluster1_idx), "); skipping next clustering.\n", sep = "")
  } else {
    emb2 <- uwot::umap(PC[non_cluster1_idx, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "euclidean", verbose = TRUE)
    minPts2 <- choose_minPts(nrow(emb2))
    cat(sprintf("[HDBSCAN PCs] Using minPts = %d for n = %d non-cluster1 samples\n", minPts2, nrow(emb2)))
    hdb2 <- dbscan::hdbscan(emb2, minPts = minPts2)
    cluster2_idx <- non_cluster1_idx[which(hdb2$cluster == which.max(table(hdb2$cluster[hdb2$cluster != 0])))]
    clusters_hdb[cluster2_idx] <- "2"
  }

  non_cluster2_idx <- which(clusters_hdb == "Unassigned")
  if (length(non_cluster2_idx) < 10) {
    cat("[WARN] Too few non-cluster2 samples (", length(non_cluster2_idx), "); skipping next clustering.\n", sep = "")
  } else {
    emb3 <- uwot::umap(PC[non_cluster2_idx, , drop = FALSE], n_neighbors = 15, min_dist = 0.3, metric = "euclidean", verbose = TRUE)
    minPts3 <- choose_minPts(nrow(emb3))
    cat(sprintf("[HDBSCAN PCs] Using minPts = %d for n = %d non-cluster2 samples\n", minPts3, nrow(emb3)))
    hdb3 <- dbscan::hdbscan(emb3, minPts = minPts3)
    cluster3_idx <- non_cluster2_idx[which(hdb3$cluster == which.max(table(hdb3$cluster[hdb3$cluster != 0])))]
    clusters_hdb[cluster3_idx] <- "3"
  }

  non_cluster3_idx <- which(clusters_hdb == "Unassigned")
  if (length(non_cluster3_idx) < 2) {
    cat("[WARN] Too few remaining samples (", length(non_cluster3_idx), "); assigning all as cluster 4.\n", sep = "")
    clusters_hdb[non_cluster3_idx] <- "4"
  } else {
    emb_lum <- uwot::umap(PC[non_cluster3_idx, , drop = FALSE], n_neighbors = 10, min_dist = 0.3, metric = "euclidean", verbose = TRUE)
    km_lum <- stats::kmeans(emb_lum, centers = 2, nstart = 50)
    lum_clusters <- ifelse(km_lum$cluster == 1, "4", "5")
    clusters_hdb[non_cluster3_idx] <- lum_clusters
  }

  clusters_hdb <- factor(clusters_hdb, levels = cluster_labels)
  clusters_hdb
}
