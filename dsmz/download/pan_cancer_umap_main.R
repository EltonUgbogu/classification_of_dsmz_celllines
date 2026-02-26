#!/usr/bin/env Rscript

# pan_cancer_umap_main.R
#
# Unified pan-cancer UMAP pipeline
# Methods:
#   --method expr  : UMAP on VST expression
#   --method pca   : PCA -> UMAP
#   --method cpca  : VST -> cPCA; ComBat_seq+VST -> project -> UMAP

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(uwot)
  library(ggplot2)
  library(sva)       # ComBat_seq
  library(optparse)
  library(scales)
})

# ─────────────────────────────────────────────────────────────
# 0) CLI
# ─────────────────────────────────────────────────────────────

option_list <- list(
  make_option(
    "--method",
    type    = "character",
    default = "expr",
    help    = "Embedding space: expr, pca, or cpca [default %default]"
  )
)

opt    <- parse_args(OptionParser(option_list = option_list))
method <- match.arg(opt$method, c("expr", "pca", "cpca"))
cat("[INFO] Method:", method, "\n")

# ─────────────────────────────────────────────────────────────
# 1) CONFIG
# ─────────────────────────────────────────────────────────────

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

base_outdir <- "/home/chu25/colDataAnnData/results"
outdir      <- file.path(base_outdir, "pan_cancer_umap_main")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

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

# ─────────────────────────────────────────────────────────────
# 2) HELPERS
# ─────────────────────────────────────────────────────────────

clean_ensg <- function(mat) {
  rn <- rownames(mat)
  rn2 <- sub("\\.\\d+$", "", rn)
  rownames(mat) <- rn2
  mat
}

load_counts <- function(path) {
  counts <- readRDS(path)       # genes x samples
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"
  counts
}

vst_per_dataset <- function(counts, dataset_id) {
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
        levels = c("Breast Cancer",
                   "Neuroblastoma",
                   "Leukemia/Lymphoma",
                   "Retinoblastoma")
      )
    )
}

theme_umap_main <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(linewidth = 0.4, colour = "grey20"),
      axis.ticks        = element_line(linewidth = 0.3, colour = "grey30"),
      axis.text         = element_text(size = 9, colour = "grey20"),
      axis.title        = element_text(size = 10, colour = "grey10"),
      legend.position   = "right",
      legend.title      = element_text(face = "bold", size = 10),
      legend.text       = element_text(size = 9),
      legend.key.size   = unit(0.8, "lines"),
      legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
      plot.title        = element_text(size = 13, face = "bold", colour = "grey10"),
      plot.subtitle     = element_text(size = 10, colour = "grey40"),
      plot.caption      = element_text(size = 7, colour = "grey50", hjust = 0),
      plot.margin       = margin(12, 12, 8, 8),
      panel.background  = element_rect(fill = "white", colour = NA)
    )
}

plot_main_overview <- function(df, palette, subtitle_text) {

  ggplot(df, aes(UMAP1, UMAP2)) +
    # Tumours: small, faint circles
    geom_point(
      data   = df %>% filter(type == "Tumour"),
      aes(colour = disease, shape = type),
      size   = 2.5,
      alpha  = 0.95,
      stroke = 0.5
    ) +
    # Cell lines: larger crosses, more opaque
    geom_point(
      data   = df %>% filter(type == "Cell Line"),
      aes(colour = disease, shape = type),
      size   = 2.0,
      alpha  = 0.95,
      stroke = 0.5
    ) +
    scale_colour_manual(
      values = palette,
      name   = "Cancer type"
    ) +
    # Circles vs crosses
    scale_shape_manual(
      values = c("Tumour" = 16, "Cell Line" = 4),
      name   = "Sample type"
    ) +
    coord_fixed(ratio = 1) +
    theme_umap_main() +
    guides(
      colour = guide_legend(
        order = 1,
        override.aes = list(size = 3, alpha = 1)
      ),
      shape  = guide_legend(
        order = 2,
        override.aes = list(size = 3, alpha = 1)
      )
    ) +
    labs(
      title    = "Pan-Cancer Tumour–Cell Line Alignment",
      subtitle = subtitle_text,
      caption  = "Colour: cancer type; shape: circle = tumour, cross = cell line",
      x = "UMAP 1",
      y = "UMAP 2"
    )
}

save_plot <- function(plot, path_base, width = 18, height = 10) {
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height, units = "in")
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, units = "in", dpi = 300)
  cat("[SAVED]", path_base, ".pdf/png\n")
}

# ─────────────────────────────────────────────────────────────
# 3) LOAD RAW COUNTS, HARMONISE GENES
# ─────────────────────────────────────────────────────────────

counts_list <- purrr::map(dataset_meta$path, load_counts)
names(counts_list) <- dataset_meta$id

counts_list <- map(counts_list, clean_ensg)

common_genes <- Reduce(intersect, map(counts_list, rownames))
cat("[INFO] Common genes across all datasets:", length(common_genes), "\n")

counts_list <- map(counts_list, ~ .x[common_genes, , drop = FALSE])

merged_counts <- do.call(cbind, counts_list)

annot <- purrr::imap_dfr(counts_list, function(mat, nm) {
  tibble(
    sample     = colnames(mat),
    dataset_id = nm
  )
}) %>%
  left_join(dataset_meta, by = c("dataset_id" = "id"))

annot <- annot[match(colnames(merged_counts), annot$sample), ]
stopifnot(all(annot$sample == colnames(merged_counts)))

cat("[INFO] Total merged samples:", ncol(merged_counts), "\n")
print(table(annot$cancer_type, annot$sample_type))

# ─────────────────────────────────────────────────────────────
# 4) PER-DATASET VST (for expr + pca + cpca basis)
# ─────────────────────────────────────────────────────────────

cat("[INFO] Computing per-dataset VST...\n")
vst_list <- map2(counts_list, names(counts_list), vst_per_dataset)
names(vst_list) <- names(counts_list)

vst_list <- map(vst_list, clean_ensg)
vst_list <- map(vst_list, ~ .x[common_genes, , drop = FALSE])

vst_merged <- do.call(cbind, vst_list)
stopifnot(all(colnames(vst_merged) == annot$sample))

# ─────────────────────────────────────────────────────────────
# 5) BUILD EMBEDDING MATRIX X FOR UMAP
# ─────────────────────────────────────────────────────────────

disease_palette <- c(
  "Breast Cancer"       = "#E41A1C",
  "Neuroblastoma"       = "#377EB8",
  "Leukemia/Lymphoma"   = "#4DAF4A",
  "Retinoblastoma"      = "#FF7F00"
)

X_umap      <- NULL
subtitle    <- NULL
embed_name  <- NULL

if (method == "expr") {
  # ── UMAP on VST expression ────────────────────────────────
  X_umap     <- t(vst_merged)
  subtitle   <- "UMAP on per-dataset VST expression (cosine distance)"
  embed_name <- "expr"

} else if (method == "pca") {
  # ── PCA -> UMAP ───────────────────────────────────────────
  cat("[INFO] Running PCA on merged VST...\n")
  pca_n_pcs <- 50
  pca_res <- prcomp(
    t(vst_merged),
    center = TRUE,
    scale. = TRUE
  )
  pcs <- pca_res$x[, 1:pca_n_pcs, drop = FALSE]
  cat("[INFO] PCA done. Using", ncol(pcs), "PCs.\n")
  X_umap     <- pcs
  subtitle   <- "UMAP on top PCs (per-dataset VST → PCA → UMAP)"
  embed_name <- "pca"

} else if (method == "cpca") {
  # ── cPCA (VST) + ComBat_seq (raw) + VST + projection → UMAP
  cat("[INFO] Running cPCA + ComBat_seq pipeline...\n")

  # cPCA basis from per-dataset VST
  data_sxg_pre <- t(vst_merged)   # samples x genes
  global_means <- colMeans(data_sxg_pre)

  center_with_global <- function(mat, means) {
    sweep(mat, 2, means, FUN = "-")
  }

  target_idx     <- annot$sample_type == "Tumour"
  background_idx <- annot$sample_type == "Cell Line"

  ALPHA        <- 0.8
  N_COMPONENTS <- 50

  target_centered <- center_with_global(data_sxg_pre[target_idx, , drop = FALSE], global_means)
  background_centered <- center_with_global(data_sxg_pre[background_idx, , drop = FALSE], global_means)

  cov_target <- (t(target_centered) %*% target_centered) /
                (nrow(target_centered) - 1)
  cov_background <- (t(background_centered) %*% background_centered) /
                    (nrow(background_centered) - 1)

  cov_contrastive <- cov_target - ALPHA * cov_background

  cat("[INFO] Eigen-decomposition of contrastive covariance...\n")
  cpc_eig <- eigen(cov_contrastive, symmetric = TRUE)
  k <- min(N_COMPONENTS, ncol(data_sxg_pre))
  cpc_vectors <- cpc_eig$vectors[, 1:k, drop = FALSE]   # genes x k
  colnames(cpc_vectors) <- paste0("cPC", seq_len(k))

  # ComBat_seq on raw counts
  annot$dataset_id <- factor(annot$dataset_id)
  batch <- annot$dataset_id

  cat("[INFO] Running ComBat_seq on raw counts...\n")
  combat_counts <- ComBat_seq(
    count = merged_counts,
    batch = batch
  )

  # VST on ComBat_seq-corrected counts (global)
  coldata_post <- data.frame(
    sample      = colnames(combat_counts),
    dataset_id  = annot$dataset_id,
    sample_type = annot$sample_type,
    cancer_type = annot$cancer_type,
    row.names   = colnames(combat_counts)
  )

  dds_post <- DESeqDataSetFromMatrix(
    countData = round(pmax(combat_counts, 0)),
    colData   = coldata_post,
    design    = ~ 1
  )
  vst_post <- assay(vst(dds_post, blind = TRUE))

  # Project ComBat_seq+VST into cPCA basis
  data_sxg_post        <- t(vst_post)
  data_sxg_post_center <- center_with_global(data_sxg_post, global_means)

  cpc_scores <- data_sxg_post_center %*% cpc_vectors   # samples x k
  colnames(cpc_scores) <- colnames(cpc_vectors)
  rownames(cpc_scores) <- colnames(combat_counts)

  X_umap     <- cpc_scores
  subtitle   <- "UMAP on contrastive PCs (VST cPCA basis → ComBat_seq + VST projection)"
  embed_name <- "cpca"
}

# ─────────────────────────────────────────────────────────────
# 6) RUN UMAP
# ─────────────────────────────────────────────────────────────

set.seed(42)
cat("[INFO] Running UMAP (method =", method, ")...\n")
emb <- uwot::umap(  
  X           = X_umap,
  n_neighbors = 30,
  min_dist    = 0.3,
  metric      = "cosine",
  scale       = FALSE,
  verbose     = TRUE
)

emb_df <- as.data.frame(emb)
colnames(emb_df) <- c("UMAP1", "UMAP2")
emb_df$sample    <- rownames(X_umap)
emb_df <- left_join(emb_df, annot, by = "sample")

saveRDS(
  emb_df,
  file.path(outdir, paste0("pan_cancer_umap_", embed_name, "_coords.rds"))
)

cat("[INFO] Embedding saved for method:", method, "\n")

plot_df <- prepare_plot_data(emb_df)

# ─────────────────────────────────────────────────────────────
# 7) MAIN VIEW PLOT ONLY
# ─────────────────────────────────────────────────────────────

p_main <- plot_main_overview(plot_df, disease_palette, subtitle)

fname_base <- file.path(outdir, paste0("01_main_overview_", embed_name))
save_plot(p_main, fname_base, width, height)

cat("\n[SUCCESS] Main plot (method =", method, ") written to:\n  ", fname_base, ".pdf/png\n", sep = "")
cat("\n=== DATASET SUMMARY ===\n")
print(table(plot_df$disease, plot_df$type))
