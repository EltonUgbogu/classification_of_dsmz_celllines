# K-means clustering utilities
# Dependencies: helpers.R (ensure_dir, safe_pdf)

run_kmeans_on_pcs <- function(PC, outdir, k_grid = 2:8, seed = 42, 
                              sample_is_tcga = NULL, V_adj = NULL) {
  set.seed(seed)
  stopifnot(is.matrix(PC), all(is.finite(PC)), is.character(outdir))
  if (!is.null(V_adj)) {
    stopifnot(is.matrix(V_adj), all(is.finite(V_adj)), ncol(V_adj) == nrow(PC))
    cat("V_adj dims:", nrow(V_adj), "genes x", ncol(V_adj), "samples\n")
    nzv <- sum(matrixStats::rowVars(V_adj) > 0)
    cat("Genes with non-zero variance:", nzv, "\n")
  }
  if (nrow(PC) > 5000) {
    cat("[WARN] Large sample set (", nrow(PC), ") for k-means; consider subsampling.\n", sep = "")
  }
  ensure_dir(outdir)
  res_list <- list()
  metrics <- lapply(k_grid, function(k) {
    km <- stats::kmeans(PC, centers = k, nstart = 50, iter.max = 100)
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
  utils::write.csv(metrics, file.path(outdir, "kmeans_metrics_HVG.csv"), row.names = FALSE)

  safe_pdf(file.path(outdir, "kmeans_elbow_HVG.pdf"), {
    plot(metrics$k, metrics$WCSS, type = "b", xlab = "k", ylab = "Within-Cluster SS", main = "Elbow (HVG space)")
  })
  safe_pdf(file.path(outdir, "kmeans_metrics_panel_HVG.pdf"), {
    op <- par(mfrow = c(2, 2)); on.exit(par(op), add = TRUE)
    plot(metrics$k, metrics$Silhouette, type = "b", xlab = "k", ylab = "Avg silhouette", main = "Silhouette (higher is better)")
    plot(metrics$k, metrics$CH, type = "b", xlab = "k", ylab = "Calinski-Harabasz", main = "Calinski-Harabasz (higher)")
    plot(metrics$k, metrics$DB, type = "b", xlab = "k", ylab = "Davies-Bouldin", main = "Davies-Bouldin (lower)")
    plot(metrics$k, metrics$WCSS, type = "b", xlab = "k", ylab = "WCSS", main = "Elbow (lower)")
  })

  best_k <- with(metrics, {
    cand <- k[Silhouette == max(Silhouette, na.rm = TRUE)]
    if (length(cand) > 1) cand[which.max(metrics$CH[match(cand, k)])] else cand
  })
  message(sprintf("[kmeans] Best k by silhouette→CH = %d", best_k))

  km_final <- res_list[[as.character(best_k)]]$km
  sil_final <- res_list[[as.character(best_k)]]$sil
  clusters_kmeans <- factor(km_final$cluster, labels = paste0("C", seq_len(best_k)))
  dataset_lab <- if (!is.null(sample_is_tcga)) ifelse(sample_is_tcga, "TCGA", "DSMZ") else rep("Unknown", nrow(PC))

  cluster_df <- tibble::tibble(
    sample = rownames(PC),
    dataset = dataset_lab,
    cluster = as.character(clusters_kmeans)
  )
  utils::write.csv(cluster_df, file.path(outdir, "clusters_kmeans.csv"), row.names = FALSE)

  safe_pdf(file.path(outdir, sprintf("silhouette_k=%d_HVG.pdf", best_k)), {
    plot(sil_final, main = sprintf("Silhouette (k = %d, HVG space)", best_k), cex.names = 0.8)
  })
  sil_dir <- file.path(outdir, "silhouettes_per_k")
  dir.create(sil_dir, showWarnings = FALSE, recursive = TRUE)
  for (k in k_grid) {
    silk <- res_list[[as.character(k)]]$sil
    if (!is.null(silk)) {
      safe_pdf(file.path(sil_dir, sprintf("silhouette_k=%d.pdf", k)), {
        plot(silk, main = sprintf("Silhouette (k = %d, HVG space)", k), cex.names = 0.7)
      })
    }
  }

  cat("Optimal k (k-means):", best_k, "with mean silhouette:", round(mean(sil_final[, "sil_width"]), 3), "\n")
  print(table(clusters_kmeans, dataset_lab))

  list(best_k = best_k, res_list = res_list, metrics = metrics, clusters = clusters_kmeans, cluster_df = cluster_df)
}
