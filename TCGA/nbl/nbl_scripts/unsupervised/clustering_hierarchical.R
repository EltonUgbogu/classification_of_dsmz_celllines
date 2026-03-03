# hierarchical_helpers.R
run_fixed_hc <- function(Exp_mat, k, seed = 42, dist_method = "euclidean", 
                hclust_method = "ward.D2") {
  stopifnot(is.matrix(Exp_mat) || is.data.frame(Exp_mat))
  Exp_mat <- as.matrix(Exp_mat)
  n_samples <- nrow(Exp_mat)

  if (is.null(k) || !is.numeric(k) || k < 2 || k != as.integer(k))
    stop("k must be an integer >= 2")
  if (n_samples < k)
    stop(sprintf("Too few samples: %d < k=%d", n_samples, k))

  set.seed(seed)
  d <- dist(Exp_mat, method = dist_method)
  hc <- hclust(d, method = hclust_method)
  cl <- cutree(hc, k = k)

  list(hc = hc, dist = d, clusters = cl)
}

run_dynamic_cut <- function(Exp_mat,
                            min_cluster_size = NULL,
                            deep_split = 2,
                            pam_stage = TRUE,
                            dist_method = "euclidean",
                            linkage = "average") {
  d <- dist(Exp_mat, method = dist_method)
  hc <- hclust(d, method = linkage)
  distM <- as.matrix(d)

  if (is.null(min_cluster_size)) {
    min_cluster_size <- max(5, ceiling(nrow(Exp_mat) * 0.05))
  }

  log_info("Running cutreeHybrid (deepSplit=%d, pamStage=%s, minClusterSize=%d)",
           deep_split, pam_stage, min_cluster_size)

  dyn <- dynamicTreeCut::cutreeHybrid(
    dendro = hc, distM = distM,
    deepSplit = deep_split, pamStage = pam_stage,
    minClusterSize = min_cluster_size
  )
  clusters <- factor(dyn$labels)

  list(hc = hc, dist = d, clusters = clusters, raw = dyn)
}


find_optimal_k <- function(Exp_mat, max_k = 10, seed = 42, dist_method = "euclidean", 
        linkage = "ward.D2") {
  if (max_k < 2) stop("max_k must be >= 2")
  if (nrow(Exp_mat) < max_k) max_k <- nrow(Exp_mat) - 1

  set.seed(seed)
  d <- dist(Exp_mat, method = dist_method)
  hc <- hclust(d, method = linkage)

  sil_scores <- numeric(max_k - 1)
  cluster_list <- vector("list", max_k - 1)
  ks <- 2:max_k

  log_info("[OPTIMAL K] Testing k = 2 to %d", max_k)

  for (i in seq_along(ks)) {
    k <- ks[i]
    clusters <- cutree(hc, k = k)
    sil <- cluster::silhouette(as.integer(clusters), d)
    sil_scores[i] <- mean(sil[, 3])
    cluster_list[[i]] <- clusters
    log_info("  k = %2d → mean silhouette = %.3f | sizes: %s",
             k, sil_scores[i], paste(table(clusters), collapse = ", "))
  }

  best_idx <- which.max(sil_scores)
  best_k <- ks[best_idx]
  best_clusters <- cluster_list[[best_idx]]

  log_info("[OPTIMAL K] Best k = %d (silhouette = %.3f)", best_k, sil_scores[best_idx])

  list(
    k = best_k,
    silhouette_mean = sil_scores[best_idx],
    silhouette_obj = cluster::silhouette(as.integer(best_clusters), d),
    clusters = best_clusters,
    hc = hc,
    dist = d,
    all_ks = ks,                    # ← ADDED
    all_silhouettes = sil_scores    # ← ADDED
  )
}


run_hierarchical_clustering <- function(Exp_mat,
                                        title,
                                        max_k = 10,
                                        seed = 42,
                                        outdir = getwd(),
                                        use_optimal_k = TRUE) {
  prefix <- gsub("[^A-Za-z0-9_-]", "_", title)
  prefix <- paste0(prefix, "_")
  result_list <- list()

  # === Optimal k via silhouette ===
  if (use_optimal_k) {
    opt <- find_optimal_k(Exp_mat, max_k = max_k, seed = seed)
    k_opt <- opt$k
    log_info("[HC] Using optimal k = %d", k_opt)

    # Save optimal clusters
    df_opt <- tibble::tibble(sample = rownames(Exp_mat), cluster = as.character(opt$clusters))
    csv_opt <- file.path(outdir, sprintf("%sclusters_optimal_k=%d.csv", prefix, k_opt))
    write.csv(df_opt, csv_opt, row.names = FALSE)
    log_info("[HC] Optimal clusters saved: %s", csv_opt)

    # Dendrogram
    pdf_dend <- file.path(outdir, sprintf("%sdendrogram_optimal_k=%d.pdf", prefix, k_opt))
    plot_dendrogram_fixed(opt$hc, k = k_opt, out_pdf = pdf_dend,
                          main = sprintf("%s - Optimal k=%d (silhouette=%.3f)", title, k_opt, 
                          opt$silhouette_mean))

    # Silhouette plot
    pdf_sil <- file.path(outdir, sprintf("%ssilhouette_optimal_k=%d.pdf", prefix, k_opt))
    plot_silhouette(opt$silhouette_obj, out_pdf = pdf_sil,
                    main = sprintf("%s - Silhouette (k=%d, mean=%.3f)", title, k_opt, 
                    opt$silhouette_mean))

    result_list$optimal <- list(
      k = k_opt,
      silhouette_mean = opt$silhouette_mean,
      clusters_labeled = opt$clusters,
      df = df_opt,
      hc = opt$hc,
      dist = opt$dist,
      silhouette_obj = opt$silhouette_obj
    )
  }

  # === Dynamic Tree Cut (always) ===
  log_info("[HC] Running dynamic tree cut...")
  dyn_res <- run_dynamic_cut(Exp_mat)
  k_dyn <- length(unique(dyn_res$clusters))
  log_info("[HC] Dynamic cut: k = %d | sizes: %s", k_dyn, paste(table(dyn_res$clusters), 
          collapse = ", "))

  df_dyn <- tibble::tibble(sample = rownames(Exp_mat), cluster = as.character(dyn_res$clusters))
  csv_dyn <- file.path(outdir, sprintf("%sclusters_dynamic_k=%d.csv", prefix, k_dyn))
  write.csv(df_dyn, csv_dyn, row.names = FALSE)

  result_list$dynamic <- list(
    k = k_dyn,
    clusters_labeled = dyn_res$clusters,
    df = df_dyn,
    hc = dyn_res$hc
  )

  invisible(result_list)
}