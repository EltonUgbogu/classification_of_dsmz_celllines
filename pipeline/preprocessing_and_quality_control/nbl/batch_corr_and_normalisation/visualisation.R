# ─────────────────────────────────────────────────────────────────────────────
# Visualization helpers
# Place in R/visualize.R
# ─────────────────────────────────────────────────────────────────────────────

#' @title Safe PDF writer
#'
#' @description Opens a PDF device, evaluates expression, and closes on exit.
#'
#' @param file Path to PDF.
#' @param expr Expression to evaluate (e.g., plot code).
#' @param width,height Dimensions in inches.
#'
#' @return NULL
#' @nexport
safe_pdf <- function(file, expr, width = 7, height = 5) {
  pdf(file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

log_info <- function(...) cat(sprintf("[INFO] %s\n", sprintf(...)))
log_warn <- function(...) cat(sprintf("[WARN] %s\n", sprintf(...)))

# ─────────────────────────────────────────────────────────────────────────────
#' @title Robust PCA + Scree Plot (PDF)
#'
#' @description
#' Creates a two-panel PCA plot (left = PC1 vs PC2 coloured by batch,
#' right = % variance explained for top 20 PCs).
#' Constant/zero-variance genes are **silently removed** before `prcomp()`.
#'
#' @param Exp_mat Gene-by-sample matrix (numeric). Rownames = genes, colnames = samples.
#' @param batch_labels Factor/character vector of length `ncol(Exp_mat)` (e.g. "DSMZ","TCGA").
#' @param title Main title for the PCA panel.
#' @param path_pdf Full path to the output PDF.
#' @param center,scale Passed to `prcomp()` (default: `center=TRUE, scale=FALSE` for VST).
#' @param top_pcs How many PCs to show in the barplot (default 20).
#' @param quiet Suppress progress messages? Default `FALSE`.
#'
#' @return Invisibly a data-frame with PC1, PC2, batch, sample, and % variance per PC.
#' @export
make_pca_plot <- function(
  Exp_mat,
  batch_labels,
  title       = "PCA Plot",
  path_pdf,
  center      = TRUE,
  scale       = FALSE,
  top_pcs     = 20,
  quiet       = FALSE
) {
  ## ------------------------------------------------------------------ ##
  ## 0. INPUT VALIDATION
  ## ------------------------------------------------------------------ ##
  if (!is.matrix(Exp_mat) && !is.data.frame(Exp_mat))
    stop("Exp_mat must be a matrix or data.frame")
  Exp_mat <- as.matrix(Exp_mat)
  if (!is.numeric(Exp_mat))
    stop("Exp_mat must contain numeric values")
  if (length(batch_labels) != ncol(Exp_mat))
    stop("batch_labels length != ncol(Exp_mat)")
  batch_labels <- as.factor(batch_labels)

  ## ------------------------------------------------------------------ ##
  ## 1. REMOVE CONSTANT / ZERO-VARIANCE GENES
  ## ------------------------------------------------------------------ ##
  row_vars <- matrixStats::rowSds(Exp_mat, na.rm = TRUE)
  row_vars[is.na(row_vars)] <- 0
  keep_genes <- row_vars > 1e-8
  n_removed <- sum(!keep_genes)
  if (n_removed > 0) {
    if (!quiet) message("[INFO] Removing ", n_removed, " constant/zero-variance genes")
    Exp_mat <- Exp_mat[keep_genes, , drop = FALSE]
  }
  if (nrow(Exp_mat) == 0)
    stop("No variable genes left after filtering")

  ## ------------------------------------------------------------------ ##
  ## 2. PCA
  ## ------------------------------------------------------------------ ##
  if (!quiet) message("[INFO] Running prcomp() (center=", center, ", scale=", scale, ")")
  pr <- prcomp(t(Exp_mat), center = center, scale. = scale)
  var_expl <- (pr$sdev^2) / sum(pr$sdev^2) * 100
  top_n <- min(top_pcs, length(var_expl))

  ## ------------------------------------------------------------------ ##
  ## 3. BUILD PLOT DATA FRAME
  ## ------------------------------------------------------------------ ##
  df <- data.frame(
    PC1    = pr$x[, 1],
    PC2    = pr$x[, 2],
    batch  = batch_labels,
    sample = colnames(Exp_mat),
    stringsAsFactors = FALSE
  )
  df$pct_var <- var_expl[seq_along(var_expl)]

  ## ------------------------------------------------------------------ ##
  ## 4. PDF OUTPUT (two-panel layout)
  ## ------------------------------------------------------------------ ##
  if (!quiet) message("[INFO] Writing PDF: ", path_pdf)
  safe_pdf(path_pdf, {
    layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))

    ## ---- Left panel: PCA scatter ----
    plot(df$PC1, df$PC2,
         pch = 19, cex = 1.1,
         col = as.numeric(df$batch),
         xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
         ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
         main = title,
         frame.plot = FALSE)
    legend("topright",
           legend = levels(df$batch),
           col = seq_along(levels(df$batch)),
           pch = 19, cex = 0.9, bty = "n")

    ## ---- Right panel: Top-N PCs barplot ----
    barplot(var_expl[1:top_n],
            names.arg = paste0("PC", 1:top_n),
            las = 2, cex.names = 0.7,
            ylab = "% Variance Explained",
            main = paste("Top", top_n, "PCs"))
    layout(1)
  })

  ## ------------------------------------------------------------------ ##
  ## 5. FINAL MESSAGES
  ## ------------------------------------------------------------------ ##
  message(sprintf("[SUCCESS] PCA plot saved: %s", path_pdf))
  message(sprintf(" PC1: %.1f%%, PC2: %.1f%%", var_expl[1], var_expl[2]))
  if (n_removed > 0) {
    message(sprintf(" Removed %d constant genes before PCA", n_removed))
  }

  invisible(df)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Robust UMAP Plot (PDF)
#'
#' @description
#' Creates a UMAP plot colored by batch and saves to PDF.
#' Optionally transposes, centers, and scales the input.
#'
#' @param Exp_mat Gene-by-sample matrix (numeric). Rownames = genes, colnames = samples.
#' @param title Main title for the plot.
#' @param path_pdf Full path to the output PDF.
#' @param batch_labels Factor/character vector of length `ncol(Exp_mat)`.
#' @param seed Random seed (default: 42).
#' @param n_neighbors UMAP neighbors (default: 20).
#' @param min_dist UMAP min_dist (default: 0.3).
#' @param metric Distance metric (default: "cosine").
#' @param scale If TRUE, uwot scales internally (default: FALSE).
#' @param center If TRUE, externally center columns (default: FALSE).
#' @param transpose If TRUE (default), transposes so samples = rows.
#'
#' @return Invisibly a data-frame with UMAP1, UMAP2, sample, batch.
#' @export
make_umap <- function(
  Exp_mat,
  title,
  path_pdf,
  batch_labels,
  seed         = 42,
  n_neighbors  = 20,
  min_dist     = 0.3,
  metric       = "cosine",
  scale        = FALSE,
  center       = FALSE,
  transpose    = TRUE
) {
  cat("[INFO] Making UMAP...\n")
  set.seed(seed)

  ## ------------------------------------------------------------------
  ## Input validation
  ## ------------------------------------------------------------------
  if (!is.matrix(Exp_mat) && !is.data.frame(Exp_mat))
    stop("Exp_mat must be a matrix or data.frame")
  Exp_mat <- as.matrix(Exp_mat)

  if (transpose) Exp_mat <- t(Exp_mat)

  if (nrow(Exp_mat) != length(batch_labels))
    stop(sprintf("Number of samples mismatch: %d (rows) != %d (batch_labels)",
                 nrow(Exp_mat), length(batch_labels)))

  ## ------------------------------------------------------------------
  ## Optional external centering
  ## ------------------------------------------------------------------
  if (center) {
    cat("[INFO] Centering matrix (external)...\n")
    Exp_mat <- scale(Exp_mat, center = TRUE, scale = FALSE)
  }

  cat(sprintf("[INFO] Running UMAP: %d samples, %d features, scale = %s\n",
              nrow(Exp_mat), ncol(Exp_mat), scale))

  ## ------------------------------------------------------------------
  ## UMAP
  ## ------------------------------------------------------------------
  emb <- uwot::umap(
    X           = Exp_mat,
    n_neighbors = n_neighbors,
    min_dist    = min_dist,
    metric      = metric,
    scale       = scale,
    verbose     = FALSE,
    ret_model   = FALSE
  )
  emb <- as.data.frame(emb)
  colnames(emb) <- paste0("UMAP", seq_len(ncol(emb)))
  emb$sample <- rownames(Exp_mat)
  emb$batch  <- factor(batch_labels)

  ## ------------------------------------------------------------------
  ## Plot & save
  ## ------------------------------------------------------------------
  cat(sprintf("[INFO] Saving UMAP to %s\n", path_pdf))
  safe_pdf(path_pdf, {
    ggplot2::ggplot(emb, ggplot2::aes(x = UMAP1, y = UMAP2, color = batch)) +
      ggplot2::geom_point(alpha = 0.85, size = 1.5) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = title) +
      ggplot2::theme(legend.title = ggplot2::element_blank())
  })
  cat(sprintf("[INFO] Done – PDF written to %s\n", path_pdf))

  invisible(emb)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title VST Diagnostics (Mean-SD + Dispersion)
#'
#' @description Runs VST and saves mean-SD and dispersion plots.
#'
#' @param counts Raw count matrix.
#' @param coldata Sample metadata.
#' @param title Dataset name (e.g. "TCGA-BRCA").
#' @param outdir Output directory.
#' @param full_disp Use all genes for dispersion? (default TRUE).
#'
#' @return NULL
#' @export
qc_vst_diagnostics <- function(
  counts,
  coldata,
  title,
  outdir,
  full_disp = TRUE
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # 1. VST
  V <- vst_normalize(counts, coldata)

  # 2. Mean-SD
  plot_mean_sd(
    V_mat   = V,
    title   = title,
    out_pdf = file.path(outdir, paste0(title, "_Mean_vs_SD.pdf"))
  )

  # 3. Dispersion
  n_genes <- if (full_disp) NULL else 10000L
  plot_dispersion(
    counts     = counts,
    coldata    = coldata,
    title      = title,
    out_pdf    = file.path(outdir, paste0(title, "_Dispersion.pdf")),
    max_genes  = n_genes
  )
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Smart stratified subsampling for dispersion
#'
#' @description Selects high/medium/low variance genes for dispersion plot.
#'
#' @param counts Count matrix.
#' @param n Number of genes to pick.
#' @param seed Random seed.
#'
#' @return Integer vector of row indices.
#' @noRd
pick_genes_for_dispersion <- function(counts, n = 10000, seed = 42) {
  set.seed(seed)
  g <- nrow(counts)
  if (n >= g) return(seq_len(g))

  vars  <- matrixStats::rowVars(counts)
  means <- rowMeans(counts)

  hi  <- head(order(vars,  decreasing = TRUE),  n %/% 3)
  lo  <- head(order(means, decreasing = FALSE), n %/% 3)
  mid <- sample(setdiff(seq_len(g), c(hi, lo)), n %/% 3)

  unique(c(hi, lo, mid))
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Dispersion plot (full or subsampled)
#'
#' @description Estimates dispersion with DESeq2 and plots.
#'
#' @param counts Raw counts.
#' @param coldata Metadata.
#' @param title Plot title.
#' @param out_pdf Output PDF path.
#' @param max_genes Max genes to use (NULL = all).
#' @param seed Random seed for subsampling.
#'
#' @return NULL
#' @export
plot_dispersion <- function(
  counts,
  coldata,
  title,
  out_pdf,
  max_genes = 10000L,
  seed      = 42L
) {
  dds <- DESeq2::DESeqDataSetFromMatrix(round(counts), coldata, design = ~1)
  g <- nrow(dds)

  if (is.null(max_genes) || max_genes >= g) {
    message("[INFO] Using ALL ", g, " genes for dispersion")
    dds_use <- dds
  } else {
    pick <- pick_genes_for_dispersion(counts, n = max_genes, seed = seed)
    dds_use <- dds[pick, ]
    message(sprintf("[INFO] Subsampled %d genes (hi/med/lo)", length(pick)))
  }

  dds_use <- DESeq2::estimateSizeFactors(dds_use)
  dds_use <- DESeq2::estimateDispersions(dds_use, fitType = "parametric")

  safe_pdf(out_pdf, {
    DESeq2::plotDispEsts(
      dds_use,
      main = sprintf("%s: Dispersion (%s genes)", title,
                     if (is.null(max_genes)) "all" else format(nrow(dds_use), big.mark = ",")),
      cex = 0.6, cex.axis = 0.9, cex.lab = 1.1, cex.main = 1.2
    )
  })
  message("[OUTPUT] ", out_pdf)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title VST Mean-vs-SD Plot
#'
#' @description Scatter of mean vs SD after VST.
#'
#' @param V_mat VST-transformed matrix.
#' @param title Plot title.
#' @param out_pdf Output PDF path.
#'
#' @return NULL
#' @export
plot_mean_sd <- function(V_mat, title, out_pdf) {
  mean_v <- rowMeans(V_mat, na.rm = TRUE)
  sd_v   <- matrixStats::rowSds(V_mat, na.rm = TRUE)
  df <- data.frame(mean = mean_v, sd = sd_v)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean, y = sd)) +
    ggplot2::geom_point(alpha = 0.7, size = 0.4, colour = "black") +
    ggplot2::geom_smooth(method = "loess", span = 0.3, colour = "red", se = FALSE, size = 1.2) +
    ggplot2::labs(
      x = "Mean of counts",
      y = "Row SD",
      title = paste(title, "VST: Mean vs SD")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold", size = 13),
      axis.title   = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(out_pdf, p, width = 7, height = 5, dpi = 600)
  message("[OUTPUT] ", out_pdf)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title UMAP on PCA space
#'
#' @description Projects PCs into 2D UMAP.
#'
#' @param PCs PC matrix (samples x PCs).
#' @param batch_labels Batch vector.
#' @param title Plot title.
#' @param path_pdf Output PDF.
#' @param n_neighbors,min_dist,metric,seed UMAP params.
#'
#' @return UMAP embedding data.frame.
#' @export
plot_umap_on_PCs <- function(
  PCs,
  batch_labels,
  title,
  path_pdf,
  n_neighbors = 20,
  min_dist    = 0.3,
  metric      = "cosine",
  seed        = 42
) {
  set.seed(seed)
  emb <- uwot::umap(PCs, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric)
  emb <- as.data.frame(emb)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$sample <- rownames(PCs)
  emb$batch  <- factor(batch_labels)

  p <- ggplot2::ggplot(emb, ggplot2::aes(UMAP1, UMAP2, color = batch)) +
    ggplot2::geom_point(alpha = 0.85, size = 1.5) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title)

  ggplot2::ggsave(path_pdf, p, width = 7, height = 6)
  invisible(emb)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title UMAP with annotation
#'
#' @description UMAP colored by first column of annotation.
#'
#' @param PC_or_X Input matrix (samples x features).
#' @param ann Annotation data.frame (samples in rows).
#' @param out_file Output PDF.
#' @param n_neighbors,min_dist,metric,title,seed UMAP params.
#'
#' @return UMAP embedding with annotation.
#' @export
plot_umap_with_annotation <- function(
  PC_or_X,
  ann,
  out_file,
  n_neighbors = 20,
  min_dist    = 0.3,
  metric      = "cosine",
  title       = NULL,
  seed        = 42
) {
  set.seed(seed)
  emb <- uwot::umap(PC_or_X, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric)
  emb <- as.data.frame(emb)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$sample <- rownames(PC_or_X)
  if (!is.null(ann)) {
    emb <- cbind(emb, ann[emb$sample, , drop = FALSE])
  }

  p <- ggplot2::ggplot(emb, ggplot2::aes(UMAP1, UMAP2, color = ann[,1])) +
    ggplot2::geom_point(alpha = 0.85, size = 1.5) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title)

  ggplot2::ggsave(out_file, p, width = 7, height = 6)
  invisible(emb)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Save heatmap with pheatmap
#'
#' @description Saves clustered heatmap to PDF.
#'
#' @param mat Matrix to plot.
#' @param ann_col Column annotation.
#' @param out_pdf Output path.
#' @param main Title.
#' @param cluster_cols,cluster_rows Logical.
#'
#' @return NULL
#' @export
# visualize.R
library(pheatmap)
library(RColorBrewer)

save_heatmap_pub <- function(mat,
                             annotation_row = NULL,
                             annotation_col = NULL,
                             out_pdf,
                             main = "Heatmap",
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             show_rownames = FALSE,
                             show_colnames = FALSE,
                             width = 10,
                             height = 8) {
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  pdf(out_pdf, width = width, height = height, pointsize = 10)
  on.exit(dev.off(), add = TRUE)

  pheatmap::pheatmap(
    mat = mat,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    main = main,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    fontsize_row = 7,
    fontsize_col = 7,
    fontsize = 10,
    color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
    border_color = NA,
    cellwidth = 3,
    cellheight = 8,
    treeheight_row = 50,
    treeheight_col = 50,
    legend = TRUE,
    annotation_legend = TRUE,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    silent = TRUE
  )
}

log_info <- function(...) cat(sprintf("[INFO] %s\n", sprintf(...)))
log_warn <- function(...) cat(sprintf("[WARN] %s\n", sprintf(...)))

# ─────────────────────────────────────────────────────────────────────────────
#' @title Fixed-k dendrogram
#'
#' @description Plots dendrogram with k-cut line.
#'
#' @param hc hclust object.
#' @param k Number of clusters.
#' @param out_pdf Output PDF.
#'
#' @return NULL
#' @export
plot_dendrogram_fixed <- function(hc, k, out_pdf, main = NULL) {
  safe_pdf(out_pdf, {
    plot(hc,
         main = if (is.null(main)) sprintf("Hierarchical (k=%d)", k) else main,
         xlab = "", sub = "", cex = 0.5, hang = -1)
    abline(h = hc$height[length(hc$height) - (k - 1)], col = "red", lty = 2)
    legend("topright", legend = "k-cut", col = "red", lty = 2, bty = "n")
  })
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Dynamic tree cut dendrogram
#'
#' @description Plots dendrogram with WGCNA colors.
#'
#' @param hc hclust object.
#' @param clusters Cluster assignments.
#' @param out_pdf Output PDF.
#'
#' @return NULL
#' @export
plot_dendrogram_dynamic <- function(hc, clusters, out_pdf, main = NULL) {
  cols <- WGCNA::labels2colors(as.numeric(factor(clusters)))
  safe_pdf(out_pdf, {
    WGCNA::plotDendroAndColors(
      dendro = hc,
      colors = cols,
      groupLabels = "Dynamic cut",
      main = if (is.null(main)) "Dynamic Tree Cut" else main,
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05
    )
  })
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Silhouette plot
#'
#' @description Plots silhouette from cluster::silhouette.
#'
#' @param sil_obj silhouette object.
#' @param out_pdf Output PDF.
#'
#' @return NULL
#' @export


# visualize.R
plot_silhouette_curve <- function(opt_result, out_pdf, width = 6, height = 4) {
  # --- VALIDATE INPUT ---
  if (is.null(opt_result$all_ks) || length(opt_result$all_ks) == 0) {
    log_warn("No silhouette data to plot. Skipping.")
    return(invisible(NULL))
  }

  df <- data.frame(
    k = opt_result$all_ks,
    silhouette = opt_result$all_silhouettes
  )

  # --- PRE-COMPUTE LABEL (NO k IN aes) ---
  label_text <- sprintf("k=%d\nsil=%.3f", opt_result$k, round(opt_result$silhouette_mean, 3))

  # --- SAFE Y POSITION ---
  y_max <- max(df$silhouette, na.rm = TRUE)
  y_pos <- y_max + 0.05 * (y_max - min(df$silhouette, na.rm = TRUE))
  if (is.infinite(y_pos)) y_pos <- 0.5  # fallback

  p <- ggplot(df, aes(x = k, y = silhouette)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(size = 3, color = "steelblue") +
    geom_vline(xintercept = opt_result$k, linetype = "dashed", color = "red", linewidth = 1) +
    annotate(
      "text",
      x = opt_result$k,
      y = y_pos,
      label = label_text,
      color = "red", size = 4, vjust = 0, hjust = 0.5
    ) +
    labs(
      title = "Silhouette Score vs Number of Clusters",
      x = "Number of clusters (k)",
      y = "Mean silhouette width"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(out_pdf, p, width = width, height = height, dpi = 300, device = "pdf")
  log_info("Silhouette curve saved: %s", out_pdf)
}

plot_silhouette <- function(sil_obj, out_pdf, main = NULL) {
  safe_pdf(out_pdf, {
    plot(sil_obj,
         main = if (is.null(main))
           sprintf("Silhouette - mean = %.3f", mean(sil_obj[, "sil_width"]))
         else main,
         col = "steelblue")
  })
}