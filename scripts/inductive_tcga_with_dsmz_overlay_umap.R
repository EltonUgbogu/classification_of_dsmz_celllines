#!/usr/bin/env Rscript

# inductive_tcga_with_dsmz_overlay_umap.R
#
# Project DSMZ BCC cell lines onto the inductive TCGA BRCA UMAP space
# and overlay them on cluster-coloured and PAM50-coloured UMAPs.

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(uwot)
  library(scales)
  library(patchwork)
  library(grid)
  library(ggrepel)
  library(yaml)
  library(optparse)
})

# shorten_cell_line_name is now in R/brca_unsup_methods.R

source("scripts/lib_config.R")

# ----------------------------------------------------------------------
# 0) CLI: get method parameter
# ----------------------------------------------------------------------
option_list <- list(
  make_option("--method", type = "character", default = "pam50_euc",
              help = "Method ID: pam50_corr, pam50_euc, MX_corr, or MX_euc [default: %default]"),
  make_option(c("-c", "--config"), type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default: %default]"),
  make_option("--profile", type = "character", default = NULL,
              help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')")
)

opt <- parse_args(OptionParser(option_list = option_list))
method_id <- opt$method

# Load method-aware helper
source("R/brca_unsup_methods.R")

# ----------------------------------------------------------------------
# 0a) Load config.yaml
# ----------------------------------------------------------------------
# Load config using profiled config loader (merges defaults + profile)
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)
unsup_root <- cfg$paths$unsup_root
data_root <- cfg$paths$data_root

# Get method context
m_ctx <- get_method_context(cfg, method_id)

cat("[INFO] Method:", method_id,
    "| feature:", m_ctx$feature,
    "| distance:", m_ctx$distance, "\n")
cat("[INFO] Overlay prefix:", m_ctx$overlay_prefix, "\n")

# ----------------------------------------------------------------------
# 0b) Paths and parameters
# ----------------------------------------------------------------------
results_dir <- file.path(unsup_root, "inductive_tcga")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Load expression matrices based on method feature
if (m_ctx$feature == "PAM50") {
  tcga_vst_rds <- cfg$paths$tcga_pam50_expr
  dsmz_vst_rds <- cfg$paths$dsmz_pam50_expr
} else if (m_ctx$feature == "MX") {
  tcga_vst_rds <- cfg$paths$tcga_mx_expr
  dsmz_vst_rds <- cfg$paths$dsmz_mx_expr
} else {
  stop("Unknown feature type: ", m_ctx$feature)
}

cat("[INFO] Using TCGA expression matrix:", tcga_vst_rds, "\n")
cat("[INFO] Using DSMZ expression matrix:", dsmz_vst_rds, "\n")

# Read canonical consensus clusters file (bestk, not k3-specific)
# Now produced by tcga_inductive_clustering_methods.R under inductive_clustering_methods/
consensus_csv <- file.path(
  results_dir,
  "inductive_clustering_methods",
  "BRCA_TCGA_PAM50_CONSENSUS_clusters_bestk.csv"
)
subtype_csv <- cfg$paths$tcga_subtype_csv

# Output paths (method-aware, K-agnostic filenames)
out_umap_cluster <- file.path(
  results_dir,
  sprintf("%s_UMAP_OVERLAY_INDUCTIVE_SAMPLES.pdf", m_ctx$overlay_prefix)
)
out_umap_pam50 <- file.path(
  results_dir,
  sprintf("%s_UMAP_OVERLAY_PAM50_SAMPLES.pdf", m_ctx$overlay_prefix)
)
out_umap_combined <- file.path(
  results_dir,
  sprintf("%s_UMAP_OVERLAY_COMBINED_SAMPLES.pdf", m_ctx$overlay_prefix)
)

# ----------------------------------------------------------------------
# 1) Load TCGA VST + consensus clusters
# ----------------------------------------------------------------------
cat("[INFO] Loading TCGA VST matrix...\n")

if (!file.exists(tcga_vst_rds))
  stop("TCGA VST matrix not found: ", tcga_vst_rds)

tcga_vst <- readRDS(tcga_vst_rds)
cat("  TCGA VST dims:", dim(tcga_vst), "\n")

# Load consensus clusters
if (!file.exists(consensus_csv))
  stop("Consensus CSV not found: ", consensus_csv)

clusters <- read.csv(consensus_csv, header = TRUE)
if (ncol(clusters) == 2) {
  colnames(clusters) <- c("sample", "cluster")
} else {
  stop("Consensus CSV must have exactly two columns: sample, cluster")
}
clusters$sample <- as.character(clusters$sample)

# Align VST with clusters
common <- intersect(colnames(tcga_vst), clusters$sample)
tcga_vst <- tcga_vst[, common, drop = FALSE]
clusters <- clusters[match(common, clusters$sample), ]
stopifnot(identical(colnames(tcga_vst), clusters$sample))

cat("  Samples with cluster assignment:", length(common), "\n")

# ----------------------------------------------------------------------
# 1a) Load DSMZ BCC VST and align genes
# ----------------------------------------------------------------------
cat("[INFO] Loading DSMZ BCC VST matrix...\n")

if (!file.exists(dsmz_vst_rds))
  stop("DSMZ BCC VST matrix not found: ", dsmz_vst_rds)

dsmz_vst <- readRDS(dsmz_vst_rds)
cat("  DSMZ BCC VST dims:", dim(dsmz_vst), "\n")

# Harmonise by intersecting genes
common_genes <- intersect(rownames(tcga_vst), rownames(dsmz_vst))
cat("  Common genes (TCGA x DSMZ):", length(common_genes), "\n")

if (length(common_genes) < 50) {
  stop("Too few common genes; check rownames/annotation.")
}

tcga_vst <- tcga_vst[common_genes, , drop = FALSE]
dsmz_vst <- dsmz_vst[common_genes, , drop = FALSE]

# ----------------------------------------------------------------------
# 1b) Load PAM50 subtypes
# ----------------------------------------------------------------------
cat("[INFO] Loading PAM50 subtypes...\n")

if (!file.exists(subtype_csv))
  stop("Subtype CSV not found: ", subtype_csv)

pam50_df <- read.csv(subtype_csv, stringsAsFactors = FALSE) %>%
  transmute(
    sample = sample_barcode,
    PAM50_Subtype = PAM50_Subtype
  )

pam50_df <- pam50_df[pam50_df$sample %in% colnames(tcga_vst), , drop = FALSE]
cat("  Samples with PAM50 annotation:", nrow(pam50_df), "\n")

# ----------------------------------------------------------------------
# 2) Inductive UMAP: fit on TCGA, transform DSMZ
# ----------------------------------------------------------------------
cat("[INFO] Computing inductive UMAP...\n")

set.seed(42)
um_model <- umap(
  t(tcga_vst),
  n_neighbors = 20,
  min_dist = 0.3,
  metric = "cosine",
  ret_model = TRUE
)

# TCGA embedding
um_tcga <- um_model$embedding
colnames(um_tcga) <- c("UMAP1", "UMAP2")
rownames(um_tcga) <- colnames(tcga_vst)

um_df_tcga <- as.data.frame(um_tcga) %>%
  mutate(
    sample = rownames(um_tcga),
    cluster = factor(
      clusters$cluster[match(rownames(um_tcga), clusters$sample)]
    )
  ) %>%
  left_join(pam50_df, by = "sample") %>%
  mutate(
    dataset = "TCGA BRCA tumour",
    dataset_type = "Tumour"
  )

cat("  TCGA samples embedded:", nrow(um_df_tcga), "\n")

# Project DSMZ into the same UMAP space
dsmz_um <- umap_transform(t(dsmz_vst), um_model)
colnames(dsmz_um) <- c("UMAP1", "UMAP2")
rownames(dsmz_um) <- colnames(dsmz_vst)

um_df_dsmz <- as.data.frame(dsmz_um) %>%
  mutate(
    sample = rownames(dsmz_um),
    cell_line_short = shorten_cell_line_name(sample),
    cluster = NA,
    PAM50_Subtype = NA_character_,
    dataset = "DSMZ BCC cell line",
    dataset_type = "Cell line"
  )

cat("  DSMZ samples projected:", nrow(um_df_dsmz), "\n")

# Combined UMAP dataframe
um_df <- bind_rows(um_df_tcga, um_df_dsmz)

# ----------------------------------------------------------------------
# 3) Colour palettes and theme (now in R/brca_unsup_methods.R)
# ----------------------------------------------------------------------
# pam50_colours, cluster_colours, and theme_dimred are loaded from shared helpers
theme_umap <- theme_dimred  # alias for backward compatibility

# ----------------------------------------------------------------------
# 4) Plotting functions with DSMZ overlay
# ----------------------------------------------------------------------

plot_umap_cluster <- function(um_df,
                              out_pdf = NULL,
                              title = "TCGA BRCA + DSMZ BCC overlay",
                              subtitle = "Coloured by inductive cluster",
                              k = NULL) {

  plot_df <- um_df %>%
    filter(!is.na(UMAP1), !is.na(UMAP2))

  tumour_df <- plot_df %>%
    filter(dataset_type == "Tumour", !is.na(cluster))
  dsmz_df <- plot_df %>%
    filter(dataset_type == "Cell line")

  if (nrow(tumour_df) == 0) {
    stop("No tumour samples with valid cluster assignment found")
  }

  present_clusters <- as.character(
    sort(unique(as.integer(as.character(tumour_df$cluster))))
  )
  plot_colours <- cluster_colours[present_clusters]

  if (any(is.na(plot_colours))) {
    extra <- scales::hue_pal()(sum(is.na(plot_colours)))
    plot_colours[is.na(plot_colours)] <- extra
  }

  if (!is.null(k)) {
    subtitle <- sprintf("Coloured by inductive cluster (k = %d)", k)
  }

  p <- ggplot() +
    geom_point(
      data = tumour_df,
      aes(x = UMAP1, y = UMAP2, colour = cluster, shape = dataset_type),
      size = 1.8,
      alpha = 0.75
    ) +
    geom_point(
      data = dsmz_df,
      aes(x = UMAP1, y = UMAP2, shape = dataset_type),
      size = 2,
      stroke = 0.7,
      colour = "black"
    ) +
    geom_label_repel(
      data = dsmz_df,
      aes(x = UMAP1, y = UMAP2, label = cell_line_short),
      size = 2.4,                # smaller text
      fill = "white",            # white background box
      label.size = 0.1,          # thin border around label
      label.r = unit(0.1, "lines"),
      colour = "black",
      max.overlaps = Inf,
      segment.size = 0.2,
      segment.alpha = 0.6
    ) +
    scale_colour_manual(values = plot_colours, name = "Cluster") +
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell line" = 4),
      name = "Sample type"
    ) +
    coord_fixed(ratio = 1) +
    theme_umap(base_size = 11) +
    guides(
      colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
      shape = guide_legend(override.aes = list(size = 3))
    ) +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2")

  if (!is.null(out_pdf)) {
    ggsave(out_pdf, p, width = 7, height = 5.5, units = "in", dpi = 300)
    out_png <- sub("\\.pdf$", ".png", out_pdf)
    ggsave(out_png, p, width = 7, height = 5.5, units = "in", dpi = 300)
    cat("[SUCCESS] Saved cluster UMAP to:\n  ", out_pdf, "\n")
  }

  return(p)
}

plot_umap_pam50 <- function(um_df,
                            out_pdf = NULL,
                            title = "TCGA BRCA + DSMZ BCC overlay",
                            subtitle = "Coloured by PAM50 subtype") {

  plot_df <- um_df %>%
    filter(!is.na(UMAP1), !is.na(UMAP2))

  tcga_df <- plot_df %>%
    filter(dataset_type == "Tumour", !is.na(PAM50_Subtype)) %>%
    mutate(subtype = factor(PAM50_Subtype, levels = names(pam50_colours)))

  dsmz_df <- plot_df %>%
    filter(dataset_type == "Cell line")

  if (nrow(tcga_df) == 0) {
    stop("No tumour samples with valid PAM50_Subtype found")
  }

  present_subtypes <- levels(droplevels(tcga_df$subtype))
  
  # Warn if unknown PAM50 subtypes are present
  unknown <- setdiff(present_subtypes, names(pam50_colours))
  if (length(unknown) > 0) {
    warning("Unknown PAM50 labels in data: ", paste(unknown, collapse = ", "))
  }
  
  plot_colours <- pam50_colours[present_subtypes]

  p <- ggplot() +
    geom_point(
      data = tcga_df,
      aes(x = UMAP1, y = UMAP2, colour = subtype, shape = dataset_type),
      size = 1.8,
      alpha = 0.75
    ) +
    geom_point(
      data = dsmz_df,
      aes(x = UMAP1, y = UMAP2, shape = dataset_type),
      size = 2,
      stroke = 0.7,
      colour = "black"
    ) +
    geom_label_repel(
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
    ) +
    scale_colour_manual(
      values = plot_colours,
      name = "PAM50 subtype",
      na.value = "grey70"
    ) +
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell line" = 4),
      name = "Sample type"
    ) +
    coord_fixed(ratio = 1) +
    theme_umap(base_size = 11) +
    guides(
      colour = guide_legend(
        title = "PAM50 subtype",
        override.aes = list(size = 3, alpha = 1)
      ),
      shape = guide_legend(override.aes = list(size = 3))
    ) +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2")

  if (!is.null(out_pdf)) {
    ggsave(out_pdf, p, width = 7, height = 5.5, units = "in", dpi = 300)
    out_png <- sub("\\.pdf$", ".png", out_pdf)
    ggsave(out_png, p, width = 7, height = 5.5, units = "in", dpi = 300)
    cat("[SUCCESS] Saved PAM50 UMAP to:\n  ", out_pdf, "\n")
  }

  return(p)
}

plot_umap_combined <- function(um_df,
                               out_pdf = NULL,
                               k = NULL) {

  plot_df <- um_df %>%
    filter(!is.na(UMAP1), !is.na(UMAP2))

  tumour_df <- plot_df %>%
    filter(dataset_type == "Tumour", !is.na(cluster), !is.na(PAM50_Subtype))
  dsmz_df <- plot_df %>%
    filter(dataset_type == "Cell line")

  if (nrow(tumour_df) == 0) {
    stop("No tumour samples with valid cluster and PAM50 found")
  }

  tumour_df <- tumour_df %>%
    mutate(
      subtype = factor(PAM50_Subtype, levels = names(pam50_colours)),
      cluster = factor(cluster)
    )

  present_subtypes <- levels(droplevels(tumour_df$subtype))
  pam50_plot_colours <- pam50_colours[present_subtypes]

  present_clusters <- as.character(
    sort(unique(as.integer(as.character(tumour_df$cluster))))
  )
  cluster_plot_colours <- cluster_colours[present_clusters]

  # Cluster panel
  p1 <- ggplot() +
    geom_point(
      data = tumour_df,
      aes(x = UMAP1, y = UMAP2, colour = cluster, shape = dataset_type),
      size = 1.5,
      alpha = 0.75
    ) +
    geom_point(
      data = dsmz_df,
      aes(x = UMAP1, y = UMAP2, shape = dataset_type),
      size = 1.8,
      stroke = 0.6,
      colour = "black"
    ) +
    geom_label_repel(
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
    ) +
    scale_colour_manual(values = cluster_plot_colours, name = "Cluster") +
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell line" = 4),
      name = "Sample type"
    ) +
    coord_fixed(ratio = 1) +
    theme_umap(base_size = 10) +
    guides(
      colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)),
      shape = guide_legend(override.aes = list(size = 2.5))
    ) +
    labs(title = "Inductive cluster", x = "UMAP 1", y = "UMAP 2")

  # PAM50 panel
  p2 <- ggplot() +
    geom_point(
      data = tumour_df,
      aes(x = UMAP1, y = UMAP2, colour = subtype, shape = dataset_type),
      size = 1.5,
      alpha = 0.75
    ) +
    geom_point(
      data = dsmz_df,
      aes(x = UMAP1, y = UMAP2, shape = dataset_type),
      size = 1.8,
      stroke = 0.6,
      colour = "black"
    ) +
    geom_label_repel(
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
    ) +
    scale_colour_manual(
      values = pam50_plot_colours,
      name = "PAM50 subtype"
    ) +
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell line" = 4),
      name = "Sample type"
    ) +
    coord_fixed(ratio = 1) +
    theme_umap(base_size = 10) +
    guides(
      colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)),
      shape = guide_legend(override.aes = list(size = 2.5))
    ) +
    labs(title = "PAM50 subtype", x = "UMAP 1", y = "UMAP 2")

  combined <- p1 + p2 +
    plot_annotation(
      title = "TCGA BRCA tumours + DSMZ BCC cell lines",
      subtitle = "Cell lines projected into tumour UMAP space",
      theme = theme(
        plot.title = element_text(
          size = 14,
          face = "bold",
          colour = "grey10"
        ),
        plot.subtitle = element_text(size = 11, colour = "grey40")
      )
    )

  if (!is.null(out_pdf)) {
    ggsave(
      out_pdf, combined,
      width = 12, height = 5.5, units = "in", dpi = 300
    )
    out_png <- sub("\\.pdf$", ".png", out_pdf)
    ggsave(
      out_png, combined,
      width = 12, height = 5.5, units = "in", dpi = 300
    )
    cat("[SUCCESS] Saved combined UMAP to:\n  ", out_pdf, "\n")
  }

  return(combined)
}

# ----------------------------------------------------------------------
# 5) Generate plots
# ----------------------------------------------------------------------
cat("\n[INFO] Generating overlay UMAPs...\n")

# Determine k dynamically from clusters
n_clusters <- length(unique(clusters$cluster))
cat("[INFO] Number of clusters detected:", n_clusters, "\n")

plot_umap_cluster(um_df, out_pdf = out_umap_cluster, k = n_clusters)
plot_umap_pam50(um_df, out_pdf = out_umap_pam50)
plot_umap_combined(um_df, out_pdf = out_umap_combined, k = n_clusters)

# ----------------------------------------------------------------------
# 6) Assign DSMZ cell lines to nearest cluster (kNN in UMAP space)
# ----------------------------------------------------------------------
cat("\n[INFO] Assigning DSMZ cell lines to nearest clusters...\n")

# Get TCGA tumour coordinates and clusters
tcga_umap <- um_df %>%
  filter(dataset_type == "Tumour", !is.na(cluster)) %>%
  select(sample, UMAP1, UMAP2, cluster, PAM50_Subtype)

# Get DSMZ cell line coordinates
dsmz_umap <- um_df %>%
  filter(dataset_type == "Cell line") %>%
  select(sample, UMAP1, UMAP2)

# For each cell line, find k nearest tumours and assign by majority vote
k_neighbors <- 10

assign_nearest_cluster <- function(dsmz_row, tcga_df, k = 10) {
  # Compute Euclidean distance in UMAP space
  distances <- sqrt(
    (tcga_df$UMAP1 - dsmz_row$UMAP1)^2 +
    (tcga_df$UMAP2 - dsmz_row$UMAP2)^2
  )

  # Get k nearest neighbours (guard against k > number of tumours)
  k_eff <- min(k, nrow(tcga_df))
  nearest_idx <- order(distances)[1:k_eff]
  nearest_tumours <- tcga_df[nearest_idx, ]

  # Majority vote for cluster
  cluster_votes <- table(nearest_tumours$cluster)
  assigned_cluster <- names(which.max(cluster_votes))

  # Majority vote for PAM50 (among neighbours with known PAM50)
  pam50_known <- nearest_tumours %>%
    filter(!is.na(PAM50_Subtype))
  if (nrow(pam50_known) > 0) {
    pam50_votes <- table(pam50_known$PAM50_Subtype)
    assigned_pam50 <- names(which.max(pam50_votes))
  } else {
    assigned_pam50 <- NA_character_
  }

  # Mean distance to k neighbours
  mean_dist <- mean(distances[nearest_idx])

  data.frame(
    sample_id = dsmz_row$sample,                       # full NG-... ID
    cell_line = shorten_cell_line_name(dsmz_row$sample),  # short name
    UMAP1 = dsmz_row$UMAP1,
    UMAP2 = dsmz_row$UMAP2,
    assigned_cluster = assigned_cluster,
    assigned_PAM50 = assigned_pam50,
    mean_kNN_distance = round(mean_dist, 4),
    stringsAsFactors = FALSE
  )
}

# Apply to all cell lines
dsmz_assignments <- do.call(rbind, lapply(
  seq_len(nrow(dsmz_umap)),
  function(i) assign_nearest_cluster(dsmz_umap[i, ], tcga_umap, k_neighbors)
))

# Reorder columns for better readability and sort
dsmz_assignments <- dsmz_assignments %>%
  select(
    sample_id,
    cell_line,
    assigned_cluster,
    assigned_PAM50,
    mean_kNN_distance,
    UMAP1, UMAP2
  ) %>%
  arrange(assigned_cluster, assigned_PAM50, cell_line)

# Print summary
cat("\nDSMZ cell line cluster assignments (k =", k_neighbors, "nearest):\n")
print(dsmz_assignments, row.names = FALSE)

# Summary by cluster
cat("\nCell lines per cluster:\n")
cluster_summary <- dsmz_assignments %>%
  group_by(assigned_cluster) %>%
  summarise(
    n_cell_lines = n(),
    cell_lines = paste(cell_line, collapse = ", "),
    .groups = "drop"
  )
print(cluster_summary, width = Inf)

# Summary by PAM50
cat("\nCell lines per PAM50 subtype:\n")
pam50_summary <- dsmz_assignments %>%
  group_by(assigned_PAM50) %>%
  summarise(
    n_cell_lines = n(),
    cell_lines = paste(cell_line, collapse = ", "),
    .groups = "drop"
  )
print(pam50_summary, width = Inf)

# Save to CSV
assignments_csv <- file.path(
  results_dir,
  "DSMZ_cell_line_cluster_assignments.csv"
)
write.csv(dsmz_assignments, assignments_csv, row.names = FALSE)
cat("\n[OUTPUT] Cell line assignments saved to:", assignments_csv, "\n")

# ----------------------------------------------------------------------
# 7) Overlay summary (CLI-friendly report)
# ----------------------------------------------------------------------
cat("\n================ OVERLAY SUMMARY ================\n")
cat("Total TCGA tumours in UMAP:  ", nrow(um_df_tcga), "\n")
cat("Total DSMZ BCC cell lines:   ", nrow(um_df_dsmz), "\n")
cat("Consensus clusters (K):      ", n_clusters, "\n\n")

cat("Cell lines per cluster:\n")
print(cluster_summary, width = Inf)

cat("\nCell lines per PAM50 subtype:\n")
print(pam50_summary, width = Inf)
cat("=================================================\n")

cat("\n[SUCCESS] DSMZ BCC overlay on inductive TCGA UMAP complete.\n")
