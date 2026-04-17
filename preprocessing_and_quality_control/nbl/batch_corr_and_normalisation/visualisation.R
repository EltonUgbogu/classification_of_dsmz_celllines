safe_pdf <- function(file, expr, width = 7, height = 5) {
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

make_pca_plot <- function(Exp_mat, batch_labels, title, path_pdf, center = TRUE, scale = FALSE, top_pcs = 20) {
  Exp_mat <- as.matrix(Exp_mat)
  row_vars <- matrixStats::rowSds(Exp_mat, na.rm = TRUE)
  keep_genes <- is.finite(row_vars) & row_vars > 1e-8
  Exp_mat <- Exp_mat[keep_genes, , drop = FALSE]
  if (nrow(Exp_mat) == 0) stop("No variable genes left for PCA")
  pr <- stats::prcomp(t(Exp_mat), center = center, scale. = scale)
  var_expl <- (pr$sdev^2) / sum(pr$sdev^2) * 100
  top_n <- min(top_pcs, length(var_expl))
  df <- data.frame(PC1 = pr$x[, 1], PC2 = pr$x[, 2], batch = as.factor(batch_labels), sample = colnames(Exp_mat))
  safe_pdf(path_pdf, {
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
    graphics::plot(df$PC1, df$PC2,
      pch = 19,
      col = as.integer(df$batch),
      xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
      ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
      main = title,
      frame.plot = FALSE
    )
    graphics::legend("topright", legend = levels(df$batch), col = seq_along(levels(df$batch)), pch = 19, bty = "n")
    graphics::barplot(var_expl[seq_len(top_n)], names.arg = paste0("PC", seq_len(top_n)), las = 2, cex.names = 0.7, ylab = "% Variance Explained")
    graphics::layout(1)
  })
  invisible(df)
}

make_umap <- function(Exp_mat, title, path_pdf, batch_labels, seed = 42, n_neighbors = 20, min_dist = 0.3, metric = "cosine", scale = FALSE, center = FALSE, transpose = TRUE) {
  Exp_mat <- as.matrix(Exp_mat)
  if (transpose) Exp_mat <- t(Exp_mat)
  if (center) Exp_mat <- scale(Exp_mat, center = TRUE, scale = FALSE)
  set.seed(seed)
  emb <- uwot::umap(Exp_mat, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric, scale = scale, verbose = FALSE)
  emb <- as.data.frame(emb)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$sample <- rownames(Exp_mat)
  emb$batch <- as.factor(batch_labels)
  p <- ggplot2::ggplot(emb, ggplot2::aes(x = UMAP1, y = UMAP2, color = batch)) +
    ggplot2::geom_point(alpha = 0.85, size = 1.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title) +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  safe_pdf(path_pdf, print(p), width = 7, height = 6)
  invisible(emb)
}

pick_genes_for_dispersion <- function(counts, n = 10000, seed = 42) {
  set.seed(seed)
  g <- nrow(counts)
  if (n >= g) return(seq_len(g))
  vars <- matrixStats::rowVars(counts)
  means <- rowMeans(counts)
  hi <- head(order(vars, decreasing = TRUE), n %/% 3)
  lo <- head(order(means, decreasing = FALSE), n %/% 3)
  mid <- sample(setdiff(seq_len(g), c(hi, lo)), n %/% 3)
  unique(c(hi, lo, mid))
}

plot_dispersion <- function(counts, coldata, title, out_pdf, max_genes = 10000L, seed = 42L) {
  dds <- DESeq2::DESeqDataSetFromMatrix(round(counts), coldata, design = ~1)
  g <- nrow(dds)
  if (!is.null(max_genes) && max_genes < g) {
    pick <- pick_genes_for_dispersion(counts, n = max_genes, seed = seed)
    dds <- dds[pick, ]
  }
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds, fitType = "parametric")
  safe_pdf(out_pdf, DESeq2::plotDispEsts(dds, main = title), width = 7, height = 5)
  invisible(out_pdf)
}

plot_mean_sd <- function(V_mat, title, out_pdf) {
  df <- data.frame(mean = rowMeans(V_mat, na.rm = TRUE), sd = matrixStats::rowSds(V_mat, na.rm = TRUE))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean, y = sd)) +
    ggplot2::geom_point(alpha = 0.7, size = 0.4, colour = "black") +
    ggplot2::geom_smooth(method = "loess", span = 0.3, colour = "red", se = FALSE, linewidth = 0.8) +
    ggplot2::labs(x = "Mean of VST values", y = "Row SD", title = paste(title, "VST: Mean vs SD")) +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(out_pdf, p, width = 7, height = 5, dpi = 300)
  invisible(out_pdf)
}

qc_vst_diagnostics <- function(counts, coldata, title, outdir, full_disp = TRUE) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  V <- vst_normalize(counts, coldata)
  plot_mean_sd(V, title, file.path(outdir, paste0(title, "_Mean_vs_SD.pdf")))
  plot_dispersion(counts, coldata, title, file.path(outdir, paste0(title, "_Dispersion.pdf")), max_genes = if (full_disp) NULL else 10000L)
  invisible(V)
}
