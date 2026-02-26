#!/usr/bin/env Rscript
# pan_cancer_umap.R
# Convert raw counts → VST → PCA → UMAP and generate summary plots.

suppressPackageStartupMessages({
  library(tidyverse)
  library(uwot)
  library(DESeq2)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(scales)
})

# -------------------------------------------------------------
# 1) CONFIG: paths to raw count matrices
# -------------------------------------------------------------
config_Datasets <- c(
  # neuroblastoma (NBL)
  nbl_counts        = "/home/chu25/dsmz/dsmz_nbl_analysis/preprocessing/results/nbl_tumour_counts.rds",
  dsmz_nbl_counts   = "/home/chu25/dsmz/dsmz_nbl_analysis/preprocessing/results/nbl_dsmz_counts.rds",
  # breast cancer (BRCA)
  dsmz_bcc_counts   = "/home/chu25/dsmz/bcc_analysis/results/dsmz_bcc_counts.rds",
  tcga_brca_counts  = "/home/chu25/dsmz/bcc_analysis/results/tcga_brca_counts.rds",
  # leu/lym (AML + DSMZ LL-100-like)
  dsmz_leu_lym_counts = "/home/chu25/dsmz/leu_lym_dsmz_analysis/results/dsmz_leu_lym_counts.rds",
  tcga_laml_counts    = "/home/chu25/dsmz/leu_lym_dsmz_analysis/results/tcga_laml_counts.rds",
  # retinoblastoma (RBL)
  rbl_dsmz_counts   = "/home/chu25/dsmz/dsmz_rbl_analysis/preprocessing/results/RBL_dsmz_counts.rds",
  rbl_tumour_counts = "/home/chu25/dsmz/dsmz_rbl_analysis/preprocessing/results/RBL_tumour_counts.rds"
)

outdir <- "/home/chu25/colDataAnnData/results/pan_cancer_umap"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Mapping dataset_name → cancer_type / sample_type
dataset_meta <- tibble(
  id   = names(config_Datasets),
  path = unname(config_Datasets)
) %>%
  mutate(
    cancer_type = case_when(
      grepl("nbl", id)          ~ "NBL",
      grepl("brca|bcc", id)     ~ "BRCA",
      grepl("leu_lym|laml", id) ~ "LEU/LYM",
      grepl("rbl", id)          ~ "RBL",
      TRUE ~ "OTHER"
    ),
    sample_type = ifelse(grepl("dsmz", id), "Cell Line", "Tumour")
  )

# -------------------------------------------------------------
# 2) LOAD COUNTS AND CONVERT TO VST (per dataset)
# -------------------------------------------------------------
load_counts_as_vst <- function(path, dataset_id) {
  counts <- readRDS(path)      # genes x samples
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"

  coldata <- data.frame(
    sample     = colnames(counts),
    dataset_id = dataset_id,
    row.names  = colnames(counts)
  )

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = coldata,
    design    = ~ 1
  )

  assay(vst(dds, blind = TRUE))
}

vst_list <- purrr::map2(
  dataset_meta$path,
  dataset_meta$id,
  load_counts_as_vst
)
names(vst_list) <- dataset_meta$id

cat("[INFO] Loaded and VST-normalized", length(vst_list), "datasets.\n")

# -------------------------------------------------------------
# 3) HARMONISE GENES (intersect Ensembl IDs)
# -------------------------------------------------------------
clean_ensg <- function(mat) {
  rn <- rownames(mat)
  rn2 <- sub("\\.\\d+$", "", rn)  # strip Ensembl version
  rownames(mat) <- rn2
  mat
}
vst_list <- map(vst_list, clean_ensg)

common_genes <- Reduce(intersect, map(vst_list, rownames))
cat("[INFO] Common genes across all datasets:", length(common_genes), "\n")

vst_list <- map(vst_list, ~ .x[common_genes, , drop = FALSE])

# -------------------------------------------------------------
# 4) MERGE ALL INTO ONE EXPRESSION MATRIX
# -------------------------------------------------------------
merged_vst <- do.call(cbind, vst_list)

annot <- purrr::imap_dfr(vst_list, function(mat, nm) {
  tibble(
    sample     = colnames(mat),
    dataset_id = nm
  )
}) %>%
  left_join(dataset_meta, by = c("dataset_id" = "id"))

annot <- annot[match(colnames(merged_vst), annot$sample), ]
stopifnot(all(annot$sample == colnames(merged_vst)))

cat("[INFO] Total merged samples:", ncol(merged_vst), "\n")
print(table(annot$cancer_type, annot$sample_type))

# Optional HVG selection (currently disabled)
# library(matrixStats)
# var_genes <- head(order(matrixStats::rowVars(merged_vst), decreasing = TRUE), 2000)
# merged_vst <- merged_vst[var_genes, , drop = FALSE]

# -------------------------------------------------------------
# 5) PCA
# -------------------------------------------------------------
set.seed(42)
cat("[INFO] Running PCA on merged VST matrix...\n")

pca_n_pcs <- 50

pca_res <- prcomp(
  t(merged_vst),          # samples x genes
  center = TRUE,
  scale. = TRUE
)

pcs <- pca_res$x[, 1:pca_n_pcs, drop = FALSE]
cat("[INFO] PCA done. Using", ncol(pcs), "PCs for UMAP.\n")

# -------------------------------------------------------------
# 6) UMAP on PCs (cosine distance)
# -------------------------------------------------------------
set.seed(42)
emb <- uwot::umap(
  X           = pcs,
  n_neighbors = 30,
  min_dist    = 0.3,
  metric      = "cosine",
  scale       = FALSE,
  verbose     = TRUE
)

emb_df <- as.data.frame(emb)
colnames(emb_df) <- c("UMAP1", "UMAP2")
emb_df$sample <- colnames(merged_vst)
emb_df <- left_join(emb_df, annot, by = "sample")

saveRDS(emb_df, file.path(outdir, "pan_cancer_umap_vst_coords.rds"))

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------
prepare_plot_data <- function(emb_df) {
  emb_df %>%
    mutate(
      type = ifelse(sample_type == "Cell Line", "Cell Line", "Tumour"),
      type = factor(type, levels = c("Tumour", "Cell Line")),
      disease = case_when(
        cancer_type == "BRCA"    ~ "Breast Cancer",
        cancer_type == "NBL"     ~ "Neuroblastoma",
        cancer_type == "LEU/LYM" ~ "Leukemia/Lymphoma",
        cancer_type == "RBL"     ~ "Retinoblastoma",
        TRUE                     ~ cancer_type
      ),
      disease = factor(
        disease,
        levels = c("Breast Cancer", "Neuroblastoma",
                   "Leukemia/Lymphoma", "Retinoblastoma")
      )
    )
}

plot_df <- prepare_plot_data(emb_df)

# -----------------------------------------------------------------------------
# COLOUR PALETTE
# -----------------------------------------------------------------------------
disease_palette <- c(
  "Breast Cancer"       = "#E41A1C",
  "Neuroblastoma"       = "#377EB8",
  "Leukemia/Lymphoma"   = "#4DAF4A",
  "Retinoblastoma"      = "#FF7F00"
)

# -----------------------------------------------------------------------------
# PLOT 1: MAIN OVERVIEW
# -----------------------------------------------------------------------------
plot_main_overview <- function(df, palette) {

  tumours    <- df %>% dplyr::filter(type == "Tumour")
  cell_lines <- df %>% dplyr::filter(type == "Cell Line")

  centroids <- tumours %>%
    dplyr::group_by(disease) %>%
    dplyr::summarise(
      UMAP1 = mean(UMAP1, na.rm = TRUE),
      UMAP2 = mean(UMAP2, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot() +
    geom_point(
      data   = tumours,
      aes(UMAP1, UMAP2, shape = type),
      colour = "grey85",
      fill   = "grey85",
      size   = 0.6,
      alpha  = 0.35,
      stroke = 0.1
    ) +
    geom_point(
      data   = cell_lines,
      aes(UMAP1, UMAP2, colour = disease, shape = type),
      size   = 0.6,
      alpha  = 0.95,
      stroke = 0.7
    ) +
    geom_text(
      data     = centroids,
      aes(UMAP1, UMAP2, label = disease),
      size     = 1.8,
      fontface = "bold",
      colour   = "grey20"
    ) +
    scale_colour_manual(
      values = palette,
      name   = "Cancer type"
    ) +
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell Line" = 21),
      name   = "Sample type"
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid       = element_blank(),
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      axis.text        = element_text(size = 9, colour = "grey20"),
      axis.title       = element_text(size = 10, colour = "grey10"),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      legend.key.size  = unit(0.8, "lines"),
      legend.background= element_rect(fill = scales::alpha("white", 0.8),
                                      colour = NA),
      plot.title       = element_text(size = 12, face = "bold", colour = "grey10"),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      plot.caption     = element_text(size = 7, colour = "grey50", hjust = 0),
      plot.margin      = margin(15, 15, 10, 10),
      panel.background = element_rect(fill = "white", colour = NA)
    ) +
    guides(
      colour = guide_legend(order = 1, override.aes = list(size = 3)),
      shape  = guide_legend(order = 2, override.aes = list(size = 3))
    ) +
    labs(
      title    = "Pan-Cancer Tumour–Cell Line Alignment",
      subtitle = "UMAP projection of tumours and DSMZ cell lines (VST normalised)",
      caption  = "Grey clouds: tumours; coloured points: cell lines (by cancer type)",
      x = "UMAP 1",
      y = "UMAP 2"
    )
}

# -----------------------------------------------------------------------------
# PLOT 2: FACETED BY CANCER TYPE
# -----------------------------------------------------------------------------
plot_faceted_by_cancer <- function(df, palette) {
  
  tumours    <- df %>% filter(type == "Tumour")
  cell_lines <- df %>% filter(type == "Cell Line")
  
  xlim <- range(df$UMAP1) + c(-0.5, 0.5) * diff(range(df$UMAP1))
  ylim <- range(df$UMAP2) + c(-0.5, 0.5) * diff(range(df$UMAP2))
  
  ggplot(mapping = aes(UMAP1, UMAP2)) +
    geom_point(
      data  = tumours,
      shape = 16,
      size  = 0.3,
      alpha = 0.15,
      colour = "grey60"
    ) +
    geom_point(
      data = tumours,
      aes(fill = disease),
      shape = 16,
      size  = 1.2,
      alpha = 0.5
    ) +
    geom_point(
      data = cell_lines,
      aes(fill = disease),
      shape  = 21,
      size   = 3,
      alpha  = 0.95,
      colour = "black",
      stroke = 0.5
    ) +
    facet_wrap(~ disease, ncol = 2, scales = "fixed") +
    scale_fill_manual(values = palette, guide = "none") +
    coord_fixed(ratio = 1, xlim = xlim, ylim = ylim) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = "grey70"),
      strip.text       = element_text(face = "bold", size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = 8),
      axis.title       = element_text(size = 9),
      plot.title       = element_text(size = 12, face = "bold"),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.spacing    = unit(1, "lines")
    ) +
    labs(
      title    = "Cancer-Type Specific Views",
      subtitle = "Each panel highlights one cancer type; grey background shows all tumours",
      x = "UMAP 1",
      y = "UMAP 2"
    )
}

# -----------------------------------------------------------------------------
# PLOT 3: CELL LINES WITH LABELS
# -----------------------------------------------------------------------------
plot_cell_lines_labeled <- function(df, palette) {
  
  tumours    <- df %>% filter(type == "Tumour")
  cell_lines <- df %>% filter(type == "Cell Line")
  
  ggplot(mapping = aes(UMAP1, UMAP2)) +
    geom_point(
      data = tumours,
      aes(colour = disease),
      shape = 16,
      size  = 0.8,
      alpha = 0.25
    ) +
    geom_point(
      data = cell_lines,
      aes(fill = disease),
      shape  = 21,
      size   = 3,
      alpha  = 0.9,
      colour = "black",
      stroke = 0.5
    ) +
    scale_fill_manual(values = palette, name = "Cancer type") +
    scale_colour_manual(values = palette, guide = "none") +
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 10) +
    theme(
      legend.position   = "right",
      legend.title      = element_text(face = "bold", size = 9),
      legend.text       = element_text(size = 8),
      axis.line         = element_line(linewidth = 0.4),
      axis.text         = element_text(size = 8),
      axis.title        = element_text(size = 9),
      plot.title        = element_text(size = 11, face = "bold"),
      plot.subtitle     = element_text(size = 8, colour = "grey40"),
      panel.background  = element_rect(fill = "grey98")
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(size = 3, shape = 21, stroke = 0.3)
      )
    ) +
    labs(
      title    = "Cell Line Positions",
      subtitle = "DSMZ cell lines highlighted; tumours shown as background",
      x = "UMAP 1",
      y = "UMAP 2"
    )
}

# -----------------------------------------------------------------------------
# PLOT 4: SAMPLE TYPE SPLIT
# -----------------------------------------------------------------------------
plot_sample_type_split <- function(df, palette) {
  
  ggplot(df, aes(UMAP1, UMAP2, fill = disease)) +
    geom_point(
      aes(size = type, alpha = type),
      shape = 21,
      colour = "black",
      stroke = 0.3
    ) +
    facet_wrap(~ type, ncol = 2) +
    scale_fill_manual(values = palette, name = "Cancer type") +
    scale_size_manual(values = c("Tumour" = 1.2, "Cell Line" = 2.5), guide = "none") +
    scale_alpha_manual(values = c("Tumour" = 0.5, "Cell Line" = 0.9), guide = "none") +
    coord_fixed(ratio = 1) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey15", colour = "grey15"),
      strip.text       = element_text(face = "bold", size = 11, colour = "white"),
      panel.grid       = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 9),
      axis.text        = element_text(size = 8),
      plot.title       = element_text(size = 12, face = "bold"),
      panel.spacing    = unit(1.5, "lines")
    ) +
    guides(
      fill = guide_legend(
        title = "Cancer type",
        override.aes = list(size = 3, alpha = 1),
        nrow = 1
      )
    ) +
    labs(
      title = "Tumours vs Cell Lines",
      x = "UMAP 1",
      y = "UMAP 2"
    )
}

# -----------------------------------------------------------------------------
# PLOT 5: DENSITY HEATMAP
# -----------------------------------------------------------------------------
plot_density_heatmap <- function(df, palette) {
  
  tumours    <- df %>% filter(type == "Tumour")
  cell_lines <- df %>% filter(type == "Cell Line")
  
  ggplot(mapping = aes(UMAP1, UMAP2)) +
    stat_density_2d(
      data    = tumours,
      aes(fill = after_stat(density)),
      geom    = "raster",
      contour = FALSE,
      alpha   = 0.8
    ) +
    scale_fill_viridis_c(
      option = "magma",
      name   = "Tumour\nDensity",
      guide  = guide_colorbar(barwidth = 1, barheight = 6)
    ) +
    geom_point(
      data  = cell_lines,
      aes(colour = disease),
      shape = 16,
      size  = 2.5,
      alpha = 0.95
    ) +
    scale_colour_manual(values = palette, name = "Cell line\ncancer type") +
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 8),
      legend.text     = element_text(size = 7),
      axis.line       = element_line(linewidth = 0.4),
      axis.text       = element_text(size = 8),
      axis.title      = element_text(size = 9),
      plot.title      = element_text(size = 11, face = "bold"),
      plot.subtitle   = element_text(size = 8, colour = "grey40"),
      panel.background= element_rect(fill = "grey10")
    ) +
    labs(
      title    = "Tumour Density Landscape",
      subtitle = "Heatmap shows tumour density; points show cell line positions",
      x = "UMAP 1",
      y = "UMAP 2"
    )
}

# -----------------------------------------------------------------------------
# PLOT 6: SUMMARY STATISTICS PANEL
# -----------------------------------------------------------------------------
plot_summary_stats <- function(df, palette) {
  
  df2 <- df %>%
    dplyr::as_tibble() %>%
    dplyr::filter(!is.na(disease), !is.na(type))
  
  counts <- df2 %>%
    dplyr::group_by(disease, type) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = type,
      values_from = n,
      values_fill = list(n = 0)
    ) %>%
    dplyr::mutate(Total = Tumour + `Cell Line`)
  
  p <- counts %>%
    tidyr::pivot_longer(
      cols      = c(Tumour, `Cell Line`),
      names_to  = "type",
      values_to = "n"
    ) %>%
    dplyr::mutate(
      type = factor(type, levels = c("Tumour", "Cell Line"))
    ) %>%
    ggplot(aes(x = disease, y = n, fill = type)) +
    geom_col(
      position = "dodge",
      width    = 0.7,
      colour   = "grey30",
      linewidth= 0.3
    ) +
    geom_text(
      aes(label = n),
      position = position_dodge(width = 0.7),
      vjust    = -0.5,
      size     = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c("Tumour" = "grey70", "Cell Line" = "grey30"),
      name   = "Sample type"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "top",
      axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
      axis.title      = element_text(size = 10),
      plot.title      = element_text(size = 11, face = "bold")
    ) +
    labs(
      title = "Sample Counts by Cancer Type",
      x     = NULL,
      y     = "Number of Samples"
    )
  
  p
}

# -----------------------------------------------------------------------------
# GENERATE ALL PLOTS
# -----------------------------------------------------------------------------
cat("[INFO] Generating plots...\n")
p1 <- plot_main_overview(plot_df, disease_palette)
p2 <- plot_faceted_by_cancer(plot_df, disease_palette)
p3 <- plot_cell_lines_labeled(plot_df, disease_palette)
p4 <- plot_sample_type_split(plot_df, disease_palette)
p5 <- plot_density_heatmap(plot_df, disease_palette)
p6 <- plot_summary_stats(plot_df, disease_palette)

# -----------------------------------------------------------------------------
# SAVE INDIVIDUAL PLOTS
# -----------------------------------------------------------------------------
save_plot <- function(plot, name, width, height) {
  ggsave(
    file.path(outdir, paste0(name, ".pdf")),
    plot, width = width, height = height, units = "in"
  )
  ggsave(
    file.path(outdir, paste0(name, ".png")),
    plot, width = width, height = height, units = "in", dpi = 300
  )
  cat("[SAVED]", name, "\n")
}

save_plot(p1, "01_main_overview",      15, 6)
save_plot(p2, "02_faceted_by_cancer",  8, 7)
save_plot(p3, "03_cell_lines_labeled", 8, 7)
save_plot(p4, "04_sample_type_split",  9, 5)
save_plot(p5, "05_density_heatmap",    7, 6)
save_plot(p6, "06_summary_counts",     6, 4)

# -----------------------------------------------------------------------------
# COMPOSITE FIGURE
# -----------------------------------------------------------------------------
composite <- (p1 | p6) / (p2) +
  plot_annotation(
    title    = "Pan-Cancer UMAP Analysis: TCGA Tumours × DSMZ Cell Lines",
    subtitle = "VST-normalized expression data projected via UMAP (cosine distance)",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  ) +
  plot_layout(heights = c(1, 1.2))

ggsave(
  file.path(outdir, "00_composite_figure.pdf"),
  composite, width = 12, height = 10, units = "in"
)
ggsave(
  file.path(outdir, "00_composite_figure.png"),
  composite, width = 12, height = 10, units = "in", dpi = 300
)

cat("\n[SUCCESS] All plots saved to:", outdir, "\n")

cat("\n=== DATASET SUMMARY ===\n")
print(table(plot_df$disease, plot_df$type))
cat("\n=== GENERATED FILES ===\n")
cat("  00_composite_figure.pdf/png   - Multi-panel summary\n")
cat("  01_main_overview.pdf/png      - Main overview UMAP\n")
cat("  02_faceted_by_cancer.pdf/png  - Per–cancer-type panels\n")
cat("  03_cell_lines_labeled.pdf/png - Cell line positions\n")
cat("  04_sample_type_split.pdf/png  - Tumour vs cell line panels\n")
cat("  05_density_heatmap.pdf/png    - Tumour density heatmap\n")
cat("  06_summary_counts.pdf/png     - Sample count bar chart\n")

