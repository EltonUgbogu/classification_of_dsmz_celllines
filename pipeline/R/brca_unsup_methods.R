# brca_unsup_methods.R
# Helper functions for method-aware configuration and filename generation
# + Shared plotting helpers for PCA/UMAP (thesis-figure-ready)

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(ggplot2)
  library(uwot)
  library(ggrepel)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ----------------------------------------------------------------------
# Get method ID from command line arguments
# ----------------------------------------------------------------------
get_method_id <- function(default = "pam50_euc") {
  args <- commandArgs(trailingOnly = TRUE)
  m_arg <- args[grepl("^--method", args)]
  if (length(m_arg) == 0) return(default)
  m <- sub("^--method(=)?", "", m_arg[1])
  if (m == "") default else m
}

# ----------------------------------------------------------------------
# Map feature/distance to tags for filenames
# ----------------------------------------------------------------------
make_method_tags <- function(feature, distance) {
  feature_tag <- toupper(feature)
  dist_tag <- switch(
    tolower(distance),
    "correlation" = "CORR", "corr" = "CORR",
    "euclidean" = "EUC", "euc" = "EUC",
    toupper(distance)
  )
  list(feature_tag = feature_tag, dist_tag = dist_tag)
}

# ----------------------------------------------------------------------
# Get full method context (feature, distance, prefixes)
# ----------------------------------------------------------------------
get_method_context <- function(
  cfg, method_id, 
  tumour_prefix = "BRCA_TCGA", 
  overlay_prefix = "BRCA_TCGA-DSMZ"
) {
  if (is.null(cfg$methods[[method_id]])) {
    stop(sprintf("Unknown method_id '%s' – check config.yaml$methods", method_id))
  }
  
  m_cfg <- cfg$methods[[method_id]]
  feature <- m_cfg$feature
  distance <- m_cfg$distance
  
  tags <- make_method_tags(feature, distance)
  
  cluster_prefix <- sprintf("%s_%s_%s", tumour_prefix, tags$feature_tag, tags$dist_tag)
  overlay_prefix_full <- sprintf("%s_%s_%s", overlay_prefix, tags$feature_tag, tags$dist_tag)
  
  list(
    method_id = method_id,
    feature = feature,
    distance = distance,
    feature_tag = tags$feature_tag,
    dist_tag = tags$dist_tag,
    cluster_prefix = cluster_prefix,
    overlay_prefix = overlay_prefix_full
  )
}

## ======================================================================
## Shared dimred helpers (PCA + UMAP) for TCGA / DSMZ / joint plots
## ======================================================================

shorten_cell_line_name <- function(x) {
  x %>%
    gsub("^NG-[0-9]+_", "", .) %>%   # remove NG-XXXXX_ prefix
    gsub("_lib[0-9_]+$", "", .)      # remove _lib... suffix
}

pam50_colours <- c(
  "Basal"  = "#4DAF4A",
  "Her2"   = "#FF7F00",
  "LumA"   = "#377EB8",
  "LumB"   = "#E41A1C",
  "Normal" = "#984EA3"
)

cluster_colours <- c(
  "1" = "#E41A1C",
  "2" = "#377EB8",
  "3" = "#4DAF4A",
  "4" = "#FF7F00",
  "5" = "#984EA3",
  "6" = "#FFFF33",
  "7" = "#A65628",
  "8" = "#F781BF"
)

theme_dimred <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      plot.title = element_text(
        size = base_size + 3,
        face = "bold",
        colour = "grey10",
        hjust = 0,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        size = base_size,
        colour = "grey40",
        hjust = 0,
        margin = margin(b = 10)
      ),
      axis.line = element_line(linewidth = 0.4, colour = "grey30"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
      axis.text = element_text(size = base_size - 1, colour = "grey20"),
      axis.title = element_text(size = base_size, colour = "grey10"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.position = "right",
      legend.justification = "top",
      legend.background = element_rect(
        fill = "white",
        colour = "grey80",
        linewidth = 0.3
      ),
      legend.margin = margin(6, 8, 6, 8),
      legend.title = element_text(
        size = base_size,
        face = "bold",
        colour = "grey20"
      ),
      legend.text = element_text(
        size = base_size - 1,
        colour = "grey30"
      ),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.size = unit(0.5, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      plot.margin = margin(15, 15, 10, 10)
    )
}

# ----------------------------------------------------------------------
# Generic PCA scatter helper (TCGA ± DSMZ)
# mat: genes × samples matrix
# meta_df: data.frame with sample, dataset_type, cluster, PAM50_Subtype (any subset)
# ----------------------------------------------------------------------
save_pca_by_cluster <- function(mat,
                                meta_df,
                                out_pdf,
                                title = "PCA: samples",
                                subtitle = NULL,
                                label_dsmz = FALSE) {
  stopifnot(all(colnames(mat) %in% meta_df$sample))
  
  mat <- mat[, meta_df$sample, drop = FALSE]
  pca <- prcomp(t(mat), scale. = FALSE)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  
  pc_df <- as.data.frame(pca$x[, 1:2])
  colnames(pc_df) <- c("PC1", "PC2")
  pc_df$sample <- rownames(pca$x)
  
  plot_df <- pc_df %>%
    left_join(meta_df, by = "sample") %>%
    mutate(
      dataset_type = dataset_type %||% "Tumour"
    )
  
  if (is.null(subtitle)) {
    subtitle <- sprintf(
      "PC1: %.1f%%, PC2: %.1f%% variance explained",
      100 * var_expl[1], 100 * var_expl[2]
    )
  }
  
  tumour_df <- plot_df %>% filter(dataset_type == "Tumour")
  dsmz_df   <- plot_df %>% filter(dataset_type == "Cell line")
  
  p <- ggplot() +
    geom_point(
      data = tumour_df,
      aes(x = PC1, y = PC2, colour = cluster, shape = dataset_type),
      size = 1.8, alpha = 0.8
    )
  
  if (nrow(dsmz_df) > 0) {
    p <- p +
      geom_point(
        data = dsmz_df,
        aes(x = PC1, y = PC2, shape = dataset_type),
        size = 2.0, stroke = 0.7, colour = "black"
      )
    
    if (label_dsmz) {
      dsmz_df <- dsmz_df %>%
        mutate(cell_line_short = shorten_cell_line_name(sample))
      
      p <- p +
        ggrepel::geom_label_repel(
          data = dsmz_df,
          aes(x = PC1, y = PC2, label = cell_line_short),
          size = 2.4,
          fill = "white",
          label.size = 0.1,
          label.r = unit(0.1, "lines"),
          colour = "black",
          max.overlaps = Inf,
          segment.size = 0.2,
          segment.alpha = 0.6
        )
    }
  }
  
  present_clusters <- as.character(sort(unique(as.integer(as.character(tumour_df$cluster)))))
  plot_colours <- cluster_colours[present_clusters]
  if (any(is.na(plot_colours))) {
    extra <- scales::hue_pal()(sum(is.na(plot_colours)))
    plot_colours[is.na(plot_colours)] <- extra
  }
  
  p <- p +
    scale_colour_manual(values = plot_colours, name = "Cluster") +
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell line" = 4),
      name   = "Sample type"
    ) +
    theme_dimred(base_size = 11) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "PC1",
      y        = "PC2"
    )
  
  ggsave(out_pdf, p, width = 7, height = 5.5, units = "in", dpi = 300)
  ggsave(sub("\\.pdf$", ".png", out_pdf), p,
         width = 7, height = 5.5, units = "in", dpi = 300)
  
  invisible(p)
}

# ----------------------------------------------------------------------
# Generic UMAP helper (single dataset OR TCGA+DSMZ overlay)
# mat: genes × samples; meta_df: sample, dataset_type, cluster, PAM50_Subtype
# ----------------------------------------------------------------------
save_umap_by_cluster <- function(mat,
                                 meta_df,
                                 out_pdf,
                                 title = "UMAP: samples",
                                 subtitle = "Coloured by cluster",
                                 n_neighbors = 20,
                                 min_dist = 0.3,
                                 metric = "cosine",
                                 label_dsmz = FALSE) {
  stopifnot(all(colnames(mat) %in% meta_df$sample))
  
  mat <- mat[, meta_df$sample, drop = FALSE]
  
  set.seed(42)
  um <- uwot::umap(
    t(mat),
    n_neighbors = n_neighbors,
    min_dist    = min_dist,
    metric      = metric
  )
  
  colnames(um) <- c("UMAP1", "UMAP2")
  rownames(um) <- colnames(mat)
  
  um_df <- as.data.frame(um) %>%
    mutate(sample = rownames(um)) %>%
    left_join(meta_df, by = "sample") %>%
    mutate(
      dataset_type = dataset_type %||% "Tumour"
    )
  
  tumour_df <- um_df %>% filter(dataset_type == "Tumour")
  dsmz_df   <- um_df %>% filter(dataset_type == "Cell line")
  
  present_clusters <- as.character(sort(unique(as.integer(as.character(tumour_df$cluster)))))
  plot_colours <- cluster_colours[present_clusters]
  if (any(is.na(plot_colours))) {
    extra <- scales::hue_pal()(sum(is.na(plot_colours)))
    plot_colours[is.na(plot_colours)] <- extra
  }
  
  p <- ggplot() +
    geom_point(
      data = tumour_df,
      aes(x = UMAP1, y = UMAP2, colour = cluster, shape = dataset_type),
      size = 1.8,
      alpha = 0.75
    )
  
  if (nrow(dsmz_df) > 0) {
    p <- p +
      geom_point(
        data = dsmz_df,
        aes(x = UMAP1, y = UMAP2, shape = dataset_type),
        size = 2,
        stroke = 0.7,
        colour = "black"
      )
    
    if (label_dsmz) {
      dsmz_df <- dsmz_df %>%
        mutate(cell_line_short = shorten_cell_line_name(sample))
      
      p <- p +
        ggrepel::geom_label_repel(
          data = dsmz_df,
          aes(x = UMAP1, y = UMAP2, label = cell_line_short),
          size = 2.4,
          fill = "white",
          label.size = 0.1,
          label.r = unit(0.1, "lines"),
          colour = "black",
          max.overlaps = Inf,
          segment.size = 0.2,
          segment.alpha = 0.6
        )
    }
  }
  
  p <- p +
    scale_colour_manual(values = plot_colours, name = "Cluster") +
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell line" = 4),
      name   = "Sample type"
    ) +
    coord_fixed(ratio = 1) +
    theme_dimred(base_size = 11) +
    guides(
      colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
      shape  = guide_legend(override.aes = list(size = 3))
    ) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "UMAP 1",
      y        = "UMAP 2"
    )
  
  ggsave(out_pdf, p, width = 7, height = 5.5, units = "in", dpi = 300)
  ggsave(sub("\\.pdf$", ".png", out_pdf), p,
         width = 7, height = 5.5, units = "in", dpi = 300)
  
  invisible(p)
}
