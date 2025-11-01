run_kmeans_on_pcs <- function(
  PC,
  outdir,
  k_grid = 2:8,
  seed = 42,
  sample_is_tcga = NULL,
  V_adj = NULL
) {
  # ===========================================================================
  # 1. INITIALIZATION & ENVIRONMENT SETUP
  # ===========================================================================
  
  ## Set random seed for full reproducibility (kmeans, dist)
  set.seed(seed)
  
  ## --- Fallback definitions for helper functions (if not in helpers.R) ---
  if (!exists("ensure_dir", mode = "function")) {
    ensure_dir <- function(p) {
      if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  if (!exists("safe_pdf", mode = "function")) {
    safe_pdf <- function(path, expr, width = 7, height = 5) {
      oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
      pdf(file = path, width = width, height = height)
      on.exit(try(dev.off(), silent = TRUE), add = TRUE)
      force(expr)
    }
  }
  
  # ===========================================================================
  # 2. INPUT VALIDATION & SANITIZATION
  # ===========================================================================
  
  ## --- Core matrix checks ---
  if (!is.matrix(PC)) stop("PC must be a numeric matrix (samples x PCs).")
  if (!is.numeric(PC)) storage.mode(PC) <- "double"
  if (!all(is.finite(PC))) {
    stop("PC contains non-finite values (NA/NaN/Inf). Clean or impute before clustering.")
  }
  if (ncol(PC) < 2) {
    stop("At least 2 principal components required for meaningful clustering and silhouette.")
  }
  
  ## --- Sample ID assignment (critical for output tracking) ---
  if (is.null(rownames(PC))) {
    rownames(PC) <- sprintf("S%05d", seq_len(nrow(PC)))
    message("[INFO] Assigned synthetic sample IDs: S00001, S00002, ...")
  }
  
  ## --- Optional: V_adj (expression matrix) diagnostics ---
  if (!is.null(V_adj)) {
    if (!is.matrix(V_adj)) stop("V_adj must be a numeric matrix (genes x samples).")
    if (!is.numeric(V_adj)) storage.mode(V_adj) <- "double"
    if (!all(is.finite(V_adj))) stop("V_adj contains non-finite values.")
    
    ## Sample alignment check
    if (!is.null(colnames(V_adj))) {
      common <- intersect(colnames(V_adj), rownames(PC))
      if (length(common) != nrow(PC)) {
        warning(sprintf(
          "[WARN] Only %d/%d samples in V_adj match PC rownames. Proceeding without reordering.",
          length(common), nrow(PC)
        ))
      }
    } else {
      if (ncol(V_adj) != nrow(PC)) {
        stop("When V_adj has no colnames, ncol(V_adj) must equal nrow(PC).")
      }
    }
    
    ## Variance diagnostic (optional, requires matrixStats)
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      nzv <- sum(matrixStats::rowVars(V_adj) > 0)
      message(sprintf("Genes with non-zero variance in V_adj: %d", nzv))
    } else {
      message("[INFO] matrixStats not available; skipping variance check.")
    }
  }
  
  ## --- TCGA flag validation ---
  if (!is.null(sample_is_tcga)) {
    if (!is.logical(sample_is_tcga)) stop("sample_is_tcga must be logical.")
    if (length(sample_is_tcga) != nrow(PC)) {
      stop("sample_is_tcga length must equal nrow(PC).")
    }
  }
  
  ## --- Large dataset warning ---
  if (nrow(PC) > 5000) {
    message(sprintf("[PERF] Large dataset (n = %d). k-means may be slow; consider subsampling.", nrow(PC)))
  }
  
  ## --- Output directory setup ---
  ensure_dir(outdir)
  sil_dir <- file.path(outdir, "silhouettes_per_k")
  ensure_dir(sil_dir)
  
  # ===========================================================================
  # 3. PRECOMPUTE DISTANCE MATRIX (for silhouette across all k)
  # ===========================================================================
  
  ## Attempt to compute Euclidean distance matrix once
  ## For very large n, this may fail (memory) → fall back to NA silhouettes
  D <- try(stats::dist(PC), silent = TRUE)
  if (inherits(D, "try-error")) {
    warning("[WARN] Failed to compute dist() matrix. Silhouette scores will be NA.")
    D <- NULL
  }
  
  # ===========================================================================
  # 4. K-MEANS LOOP OVER k_grid
  # ===========================================================================
  
  res_list <- list()        # Stores kmeans object + silhouette per k
  metrics_rows <- list()    # Aggregated metrics table
  
  for (k in k_grid) {
    ## Skip invalid k (e.g., k > n_samples)
    if (k >= nrow(PC)) {
      warning(sprintf("[SKIP] k = %d >= n_samples (%d). Skipping.", k, nrow(PC)))
      next
    }
    
    message(sprintf("[k-means] Running k = %d (nstart = 50)...", k))
    
    ## Run k-means with multiple restarts
    km <- kmeans(PC, centers = k, nstart = 50, iter.max = 100)
    wcss <- sum(km$withinss)
    
    ## --- Silhouette (requires D) ---
    sil_mean <- NA_real_
    sil_obj <- NULL
    if (!is.null(D)) {
      sil_try <- try(cluster::silhouette(km$cluster, D), silent = TRUE)
      if (!inherits(sil_try, "try-error")) {
        sil_obj <- sil_try
        sil_mean <- mean(sil_try[, "sil_width"])
      }
    }
    
    ## --- Calinski-Harabasz (fpc package) ---
    ch <- tryCatch(
      fpc::calinhara(PC, km$cluster),
      error = function(e) NA_real_
    )
    
    ## --- Davies-Bouldin (clusterCrit package) ---
    db <- tryCatch({
      if (!requireNamespace("clusterCrit", quietly = TRUE)) {
        NA_real_
      } else {
        crit <- clusterCrit::intCriteria(as.matrix(PC), as.integer(km$cluster), "davies_bouldin")
        as.numeric(crit$davies_bouldin)
      }
    }, error = function(e) NA_real_)
    
    ## Store results
    res_list[[as.character(k)]] <- list(km = km, sil = sil_obj)
    metrics_rows[[as.character(k)]] <- data.frame(
      k = k,
      WCSS = wcss,
      Silhouette = sil_mean,
      CH = ch,
      DB = db,
      stringsAsFactors = FALSE
    )
    
    ## Per-k silhouette plot
    if (!is.null(sil_obj)) {
      safe_pdf(file.path(sil_dir, sprintf("silhouette_k=%d.pdf", k)), {
        plot(sil_obj,
             main = sprintf("Silhouette Plot (k = %d)", k),
             cex.names = 0.7,
             col = rainbow(k))
      })
    }
  }
  
  # ===========================================================================
  # 5. AGGREGATE METRICS & SAVE
  # ===========================================================================
  
  if (length(metrics_rows) == 0) {
    stop("No valid k values were processed. Check k_grid and input data.")
  }
  
  metrics <- do.call(rbind, metrics_rows)
  rownames(metrics) <- NULL
  write.csv(metrics, file.path(outdir, "kmeans_metrics.csv"), row.names = FALSE)
  message("[OUTPUT] Metrics table saved: kmeans_metrics.csv")
  
  # ===========================================================================
  # 6. DIAGNOSTIC PLOTS
  # ===========================================================================
  
  ## Elbow plot (WCSS)
  safe_pdf(file.path(outdir, "kmeans_elbow.pdf"), {
    plot(metrics$k, metrics$WCSS,
         type = "b", pch = 19, col = "blue",
         xlab = "Number of clusters (k)",
         ylab = "Within-cluster sum of squares",
         main = "Elbow Method (lower is better)")
    grid()
  })
  
  ## Multi-panel metrics
  safe_pdf(file.path(outdir, "kmeans_metrics_panel.pdf"), width = 10, height = 8, {
    op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
    
    plot(metrics$k, metrics$Silhouette, type = "b", pch = 19, col = "darkgreen",
         xlab = "k", ylab = "Avg Silhouette", main = "Silhouette (higher better)")
    grid()
    
    plot(metrics$k, metrics$CH, type = "b", pch = 19, col = "purple",
         xlab = "k", ylab = "Calinski-Harabasz", main = "CH Index (higher better)")
    grid()
    
    plot(metrics$k, metrics$DB, type = "b", pch = 19, col = "red",
         xlab = "k", ylab = "Davies-Bouldin", main = "DB Index (lower better)")
    grid()
    
    plot(metrics$k, metrics$WCSS, type = "b", pch = 19, col = "orange",
         xlab = "k", ylab = "WCSS", main = "Elbow (lower better)")
    grid()
  })
  
  # ===========================================================================
  # 7. SELECT OPTIMAL k (ROBUST HIERARCHY)
  # ===========================================================================
  
  best_k <- NA_integer_
  
  ## Primary: Maximize silhouette
  if (any(is.finite(metrics$Silhouette))) {
    max_sil <- max(metrics$Silhouette, na.rm = TRUE)
    candidates <- metrics$k[metrics$Silhouette == max_sil]
    
    ## Tie-break: Maximize CH among silhouette winners
    if (length(candidates) > 1 && any(is.finite(metrics$CH))) {
      best_k <- candidates[which.max(metrics$CH[match(candidates, metrics$k)])]
    } else {
      best_k <- candidates[1]
    }
  }
  ## Fallback 1: Maximize CH
  else if (any(is.finite(metrics$CH))) {
    best_k <- metrics$k[which.max(metrics$CH)]
  }
  ## Fallback 2: Largest k (conservative elbow)
  else {
    best_k <- max(metrics$k)
    warning(sprintf(
      "[FALLBACK] No valid silhouette or CH → using largest k = %d", best_k
    ))
  }
  
  message(sprintf(
    "[RESULT] Optimal k = %d | Silhouette = %.3f | CH = %.1f | DB = %.3f",
    best_k,
    metrics$Silhouette[metrics$k == best_k],
    metrics$CH[metrics$k == best_k],
    metrics$DB[metrics$k == best_k]
  ))
  
  # ===========================================================================
  # 8. FINAL CLUSTER ASSIGNMENTS & OUTPUT
  # ===========================================================================
  
  km_final <- res_list[[as.character(best_k)]]$km
  sil_final <- res_list[[as.character(best_k)]]$sil
  
  ## Final cluster labels
  clusters_kmeans <- factor(km_final$cluster, labels = paste0("C", seq_len(best_k)))
  
  ## Dataset label (TCGA vs DSMZ)
  dataset_lab <- if (is.null(sample_is_tcga)) {
    rep("Unknown", nrow(PC))
  } else {
    ifelse(sample_is_tcga, "TCGA", "DSMZ")
  }
  
  ## Final assignment table
  cluster_df <- tibble::tibble(
    sample = rownames(PC),
    dataset = dataset_lab,
    cluster = as.character(clusters_kmeans)
  )
  write.csv(cluster_df, file.path(outdir, "clusters_kmeans.csv"), row.names = FALSE)
  message("[OUTPUT] Final clusters saved: clusters_kmeans.csv")
  
  ## Final silhouette plot
  if (!is.null(sil_final)) {
    safe_pdf(file.path(outdir, sprintf("silhouette_k=%d_final.pdf", best_k)), {
      plot(sil_final,
           main = sprintf("Final Silhouette (k = %d)", best_k),
           cex.names = 0.8,
           col = rainbow(best_k))
    })
    message(sprintf(
      "[QUALITY] Mean silhouette width = %.3f", mean(sil_final[, "sil_width"])
    ))
  }
  
  ## Print contingency table
  print(table(clusters_kmeans, dataset_lab))
  
  # ===========================================================================
  # 9. PER-K REPRODUCIBILITY 
  # ===========================================================================
  
  assign_dir <- file.path(outdir, "cluster_assignments_per_k")
  ensure_dir(assign_dir)
  
  for (k in names(res_list)) {
    kmk <- res_list[[k]]$km
    dfk <- data.frame(
      sample = rownames(PC),
      cluster = paste0("C", kmk$cluster),
      stringsAsFactors = FALSE
    )
    write.csv(dfk, file.path(assign_dir, sprintf("clusters_k=%s.csv", k)), row.names = FALSE)
  }
  message("[OUTPUT] Per-k assignments saved in: ", assign_dir)
  
  # ===========================================================================
  # 10. RETURN STRUCTURED RESULT
  # ===========================================================================
  
  invisible(list(
    best_k = best_k,
    res_list = res_list,
    metrics = metrics,
    clusters = clusters_kmeans,
    cluster_df = cluster_df
  ))
}