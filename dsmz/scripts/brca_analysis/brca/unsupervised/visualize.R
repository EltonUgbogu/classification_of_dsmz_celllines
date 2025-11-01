# Visualization helpers
# Place in R/visualize.R

make_pca_plot <- function(M, lab_batch, title, path_pdf, n_top = 5000) {
  sel <- top_var_genes(M, n = min(n_top, nrow(M)))
  X   <- t(scale(M[sel, , drop = FALSE]))
  pr  <- prcomp(X, center = FALSE, scale. = FALSE)
  var_expl <- (pr$sdev^2) / sum(pr$sdev^2) * 100
  df <- data.frame(PC1 = pr$x[,1], PC2 = pr$x[,2],
                   batch = lab_batch,
                   sample = rownames(pr$x),
                   stringsAsFactors = FALSE)
  safe_pdf(path_pdf, {
    layout(matrix(c(1,2), nrow = 1))
    plot(df$PC1, df$PC2, pch = 19, cex = 0.7, col = as.factor(df$batch),
         xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
         ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
         main = title)
    legend("topright", legend = levels(as.factor(df$batch)),
           col = seq_along(levels(as.factor(df$batch))), pch = 19, cex = 0.8)
    barplot(var_expl[1:20], las = 2, main = "Explained variance (top 20 PCs)",
            ylab = "% variance", xlab = "PC")
    layout(1)
  })
  invisible(list(scores = df, var_expl = var_expl))
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