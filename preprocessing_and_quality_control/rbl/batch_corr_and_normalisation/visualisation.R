## Cohort preprocessing figures, shared verbatim by NBL and RBL.
##
## Every computation in this file is unchanged. PCA is still `prcomp` over all
## rows with non-zero standard deviation, UMAP is still `uwot::umap` with the
## same neighbourhood parameters and seed, the dispersion panel is still the
## DESeq2 parametric fit over the same gene selection, and the mean-SD panel is
## still a loess trend with span 0.3 over every finite row. Sample inclusion,
## normalisation and batch correction are untouched.
##
## Only the rendering changed. Colours, point styling, legend wording,
## typography and the gridless panel frame now come from the shared figure
## module (`preprocessing_and_quality_control/scripts/R/batch_effect_and_embedding_figures.R`),
## so a patient tumour is the same blue and a DSMZ cell line the same orange in
## every cohort, and no panel carries a background grid.
##
## The shared module must be loaded before any function here is called; the
## preprocessing scripts source it from the rule's `shared_figure_module` input.

make_pca_plot <- function(Exp_mat, batch_labels, title, path_pdf, center = TRUE, scale = FALSE, top_pcs = 20, subtitle = NULL) {
  Exp_mat <- as.matrix(Exp_mat)
  row_vars <- matrixStats::rowSds(Exp_mat, na.rm = TRUE)
  keep_genes <- is.finite(row_vars) & row_vars > 1e-8
  Exp_mat <- Exp_mat[keep_genes, , drop = FALSE]
  if (nrow(Exp_mat) == 0) stop("No variable genes left for PCA")
  pr <- stats::prcomp(t(Exp_mat), center = center, scale. = scale)
  var_expl <- (pr$sdev^2) / sum(pr$sdev^2) * 100
  top_n <- min(top_pcs, length(var_expl))
  df <- data.frame(
    PC1 = pr$x[, 1],
    PC2 = pr$x[, 2],
    batch = as.factor(batch_labels),
    sample = colnames(Exp_mat)
  )

  ## Scores panel on the left, scree on the right, one shared legend beneath
  ## both, so the legend can never sit on top of a point.
  style <- figure_series_style(batch_labels)
  figure_open(path_pdf, width = 9.8, height = 6.6)
  on.exit(figure_close(), add = TRUE)
  graphics::layout(
    matrix(c(1, 2, 3, 3), nrow = 2L, byrow = TRUE),
    widths = c(3, 1.15),
    heights = c(1, 0.14)
  )
  graphics::par(oma = c(0.4, 0, 3.2, 0))

  figure_reset_par()
  graphics::par(mar = c(4.2, 4.8, 2.0, 1.2), pty = "s")
  limits <- figure_square_limits(df$PC1, df$PC2)
  figure_panel_frame(
    df$PC1,
    df$PC2,
    xlab = sprintf("PC1 (%.1f%% variance)", var_expl[1]),
    ylab = sprintf("PC2 (%.1f%% variance)", var_expl[2]),
    xlim = limits$xlim,
    ylim = limits$ylim,
    scale_style = "numeric"
  )
  figure_draw_series(df$PC1, df$PC2, style)

  ## Scree panel: the same variance-explained values as before, drawn as bars
  ## on a gridless frame instead of a default barplot.
  figure_reset_par()
  graphics::par(mar = c(4.2, 4.6, 2.0, 1.0), pty = "m")
  scree <- var_expl[seq_len(top_n)]
  positions <- seq_len(top_n)
  figure_panel_frame(
    xlim = c(0.4, top_n + 0.6),
    ylim = c(0, max(scree) * 1.12),
    xlab = "Principal component",
    ylab = "% variance explained",
    x_axis = FALSE
  )
  graphics::rect(
    positions - 0.38, 0, positions + 0.38, scree,
    col = FIGURE_INK_MUTED, border = "#ffffff", lwd = 0.7
  )
  graphics::mtext(
    paste0("PC", positions),
    side = 1, line = 0.6, at = positions,
    cex = 0.58, col = FIGURE_INK_MUTED, las = 2
  )

  figure_series_legend(style)
  graphics::mtext(
    title, side = 3, line = 1.6, outer = TRUE, adj = 0.02,
    font = 2, cex = 1.35, col = FIGURE_INK
  )
  if (!is.null(subtitle)) {
    graphics::mtext(
      subtitle, side = 3, line = 0.3, outer = TRUE, adj = 0.02,
      cex = 0.95, col = FIGURE_INK_MUTED
    )
  }
  invisible(df)
}

make_umap <- function(Exp_mat, title, path_pdf, batch_labels, seed = 42, n_neighbors = 20, min_dist = 0.3, metric = "cosine", scale = FALSE, center = FALSE, transpose = TRUE, subtitle = NULL) {
  Exp_mat <- as.matrix(Exp_mat)
  if (transpose) Exp_mat <- t(Exp_mat)
  if (center) Exp_mat <- scale(Exp_mat, center = TRUE, scale = FALSE)
  set.seed(seed)
  emb <- uwot::umap(
    Exp_mat,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    scale = scale,
    verbose = FALSE
  )
  emb <- as.data.frame(emb)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$sample <- rownames(Exp_mat)
  emb$batch <- as.factor(batch_labels)

  ## "bare" axes: UMAP coordinates carry no direct numerical meaning, so the
  ## panel keeps the two axis lines and drops ticks and tick labels.
  plot_embedding_figure(
    list(list(embedding = list(
      coords = as.matrix(emb[, c("UMAP1", "UMAP2")]),
      xlab = "UMAP 1",
      ylab = "UMAP 2"
    ))),
    batch_labels,
    path_pdf,
    overall_title = title,
    overall_subtitle = subtitle,
    scale_style = "bare"
  )
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

  ## DESeq2 draws the panel; we supply the palette and our own legend so the
  ## three dispersion quantities never borrow a series colour's meaning.
  figure_open(out_pdf, width = 7.8, height = 5.9)
  on.exit(figure_close(), add = TRUE)
  figure_reset_par()
  graphics::par(mar = c(4.4, 4.8, 4.4, 1.6))
  DESeq2::plotDispEsts(
    dds,
    main = "",
    genecol = grDevices::adjustcolor(FIGURE_INK_MUTED, alpha.f = 0.40),
    fitcol = FIGURE_TREND,
    finalcol = grDevices::adjustcolor(FIGURE_SERIES_COLOURS[[1L]], alpha.f = 0.60),
    legend = FALSE,
    xlab = "Mean of normalised counts",
    ylab = "Dispersion"
  )
  graphics::mtext(
    title, side = 3, line = 2.2, adj = 0, font = 2,
    cex = 1.15 * graphics::par("cex"), col = FIGURE_INK
  )
  graphics::mtext(
    sprintf("DESeq2 parametric dispersion fit  |  %s genes", figure_count(nrow(dds))),
    side = 3, line = 1.1, adj = 0,
    cex = 0.85 * graphics::par("cex"), col = FIGURE_INK_MUTED
  )
  graphics::legend(
    "bottomleft",
    legend = c("Gene-wise estimate", "Fitted trend", "Final (shrunk) estimate"),
    pch = c(16, NA, 16),
    lty = c(NA, 1, NA),
    lwd = c(NA, 2.4, NA),
    col = c(FIGURE_INK_MUTED, FIGURE_TREND, FIGURE_SERIES_COLOURS[[1L]]),
    bty = "n",
    cex = 0.9,
    text.col = FIGURE_INK
  )
  invisible(out_pdf)
}

plot_mean_sd <- function(V_mat, title, out_pdf) {
  means <- rowMeans(V_mat, na.rm = TRUE)
  sds <- matrixStats::rowSds(V_mat, na.rm = TRUE)
  keep <- is.finite(means) & is.finite(sds)
  means <- means[keep]
  sds <- sds[keep]

  figure_open(out_pdf, width = 7.8, height = 5.9)
  on.exit(figure_close(), add = TRUE)
  figure_reset_par()
  graphics::par(mar = c(4.4, 4.8, 4.4, 1.6))
  figure_panel_frame(
    means,
    sds,
    xlab = "Mean of VST values",
    ylab = "Row standard deviation",
    main = paste(title, "VST: mean versus standard deviation"),
    subtitle = sprintf(
      "%s genes  |  colour = local point density  |  orange line = loess trend (span 0.3)",
      figure_count(length(means))
    )
  )
  graphics::points(
    means, sds,
    pch = 16, cex = 0.5,
    col = grDevices::adjustcolor(figure_density_colours(means, sds), alpha.f = 0.55)
  )
  ## Same loess trend as before (span 0.3); drawn directly instead of through
  ## ggplot2::geom_smooth so the panel shares the base-R figure style.
  trend <- stats::loess(sds ~ means, span = 0.3)
  ordered <- order(means)
  graphics::lines(
    means[ordered], stats::fitted(trend)[ordered],
    col = FIGURE_TREND, lwd = 2.8
  )
  invisible(out_pdf)
}

qc_vst_diagnostics <- function(counts, coldata, title, outdir, full_disp = TRUE) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  V <- vst_normalize(counts, coldata)
  plot_mean_sd(V, title, file.path(outdir, paste0(title, "_Mean_vs_SD.pdf")))
  plot_dispersion(counts, coldata, title, file.path(outdir, paste0(title, "_Dispersion.pdf")), max_genes = if (full_disp) NULL else 10000L)
  invisible(V)
}
