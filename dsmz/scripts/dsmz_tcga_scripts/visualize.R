# Visualization helpers
# Place in R/visualize.R

# M - VST normalized matrix of expression data
# batch_labels - vector of batch labels
# title - title of the plot
# path_pdf - path to save the plot

make_pca_plot <- function(M, 
                          batch_labels,
                          title = "PCA Plot", 
                          path_pdf,
                          center = TRUE, 
                          scale = FALSE) {  # ← scale = FALSE by default for VST
  cat("[INFO] Making PCA plot...\n")
  # Input validation
  stopifnot(is.matrix(M) || is.data.frame(M))
  stopifnot(length(batch_labels) == ncol(M))

  pr <- prcomp(t(M), center = center, scale. = scale)
  cat("[INFO] Calculating variance explained...\n")
  var_expl <- (pr$sdev^2) / sum(pr$sdev^2) * 100
  cat("[INFO] Creating data frame...\n")
  df <- data.frame(
    PC1 = pr$x[, 1],
    PC2 = pr$x[, 2],
    batch = factor(batch_labels),
    sample = rownames(pr$x),
    stringsAsFactors = FALSE
  )
  cat("[INFO] Saving PCA plot to PDF...\n")
  # Save to PDF safely
  safe_pdf(path_pdf, {
    layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))

    # PCA scatter
    plot(df$PC1, df$PC2,
         pch = 19, cex = 1.2,
         col = as.numeric(df$batch),
         xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
         ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
         main = title,
         frame.plot = FALSE)
    legend("topright",
           legend = levels(df$batch),
           col = seq_along(levels(df$batch)),
           pch = 19, cex = 0.9, bty = "n")

    # Variance explained barplot
    barplot(var_expl[1:min(20, length(var_expl))],
            names.arg = paste0("PC", 1:min(20, length(var_expl))),
            las = 2, cex.names = 0.7,
            ylab = "% Variance Explained",
            main = "Top 20 PCs")

    layout(1)
  })

  message(sprintf("[INFO] PCA plot saved: %s", path_pdf))
  message(sprintf("   PC1: %.1f%%, PC2: %.1f%%", var_expl[1], var_expl[2]))

  return(df)
}

# M - matrix of data (merged counts matrix of TCGA and DSMZ)
# title - title of the plot
# path_pdf - path to save the plot
# batch_labels - vector of batch labels e,g c("TCGA", "DSMZ")
# seed - random seed (default: 42)
# n_neighbors - number of neighbors (default: 20)
# min_dist - minimum distance (default: 0.3)
# metric - metric to use for UMAP (default: cosine)
make_umap <- function(M, title, path_pdf, batch_labels, seed = 42, n_neighbors = 20, min_dist = 0.3,
                      metric = "cosine", scale = FALSE) {
  cat("[INFO] Making UMAP...\n")
  set.seed(seed)
  
  if (length(batch_labels) != ncol(M)) {
    stop(sprintf("length(batch_labels) must equal number of samples (ncol(M)): %d != %d", length(batch_labels), ncol(M)))
  }
  cat(sprintf("[INFO] Running UMAP with %d samples and %d features\n", ncol(M), nrow(M)))
  cat(sprintf("[INFO] Running UMAP with scale = %s\n", scale))
  emb <- uwot::umap(t(M),
                    n_neighbors = n_neighbors,
                    min_dist = min_dist,
                    metric = metric,
                    scale = scale,          # Use uwot's internal scaling
                    verbose = FALSE)
  
  emb <- as.data.frame(emb)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$sample <- colnames(M)
  emb$batch <- factor(batch_labels)
  
  cat(sprintf("[INFO] Saving UMAP as pdf: %s\n", path_pdf))
  safe_pdf(path_pdf, {
    ggplot(emb, aes(UMAP1, UMAP2, color = batch)) + 
      geom_point(alpha = 0.85, size = 1.5) + 
      theme_bw() + 
      ggtitle(title)
  })
  
  cat(sprintf("[INFO] UMAP saved to PDF: %s\n", path_pdf))
  invisible(emb)
}

plot_umap_on_PCs <- function(PCs, batch_labels, title, n_neighbors = 20, min_dist = 0.3, metric = "cosine", seed = 42) {
  set.seed(seed)
  emb <- uwot::umap(PC, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric)
  emb <- as.data.frame(emb); colnames(emb) <- c("UMAP1","UMAP2")
  emb$sample <- rownames(PCs)
  emb$batch <- factor(batch_labels)
  p <- ggplot(emb, aes(UMAP1, UMAP2, color = batch)) + geom_point(alpha = 0.85, size = 1.5) + theme_bw() + ggtitle(title)
  ggsave(path_pdf, p, width = 7, height = 6)
  emb
}

plot_umap_with_annotation <- function(PC_or_X, ann, out_file, n_neighbors = 20, min_dist = 0.3, metric = "cosine", title = NULL) {
  set.seed(42)
  emb <- uwot::umap(PC_or_X, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric)
  emb <- as.data.frame(emb); colnames(emb) <- c("UMAP1","UMAP2")
  emb$sample <- rownames(PC_or_X)
  if (!is.null(ann)) emb <- cbind(emb, ann[emb$sample, , drop = FALSE])
  p <- ggplot(emb, aes(UMAP1, UMAP2, color = ann[,1])) + geom_point(alpha = 0.85, size = 1.5) + theme_bw() + ggtitle(title)
  ggsave(out_file, p, width = 7, height = 6)
  emb
}

save_heatmap <- function(mat, ann_col = NULL, out_pdf, main = "Heatmap", cluster_cols = TRUE, cluster_rows = TRUE) {
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  pdf(out_pdf, width = 10, height = 8)
  pheatmap::pheatmap(mat, annotation_col = ann_col, main = main, cluster_cols = cluster_cols, cluster_rows = cluster_rows, show_rownames = TRUE, show_colnames = FALSE)
  dev.off()
}

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