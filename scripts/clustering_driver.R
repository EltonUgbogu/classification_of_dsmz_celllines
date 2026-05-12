#!/usr/bin/env Rscript

# clustering_driver.R
# Generic driver for combos:
#   feature_set: pam50 / MX
#   distance:    euclidean / correlation
#
# Uses:
#   - config$paths$vst_joint_rds
#   - config$features$pam50_gene_list
#   - config$features$mx_final_gene_list
#   - config$clustering$k_grid
#
# Produces (inside --outdir):
#   BRCA_TCGA-DSMZ_<FEAT>_<DIST>_CLUSTERS_kmeans_kK_SAMPLES.csv
#   BRCA_TCGA-DSMZ_<FEAT>_<DIST>_PCA_coordinates_kK_SAMPLES.tsv
#   BRCA_TCGA-DSMZ_<FEAT>_<DIST>_UMAP_coordinates_kK_SAMPLES.tsv
#   BRCA_TCGA-DSMZ_<FEAT>_<DIST>_CLUSTERING_SUMMARY_kmeans_kK.txt
#   DONE is touched by Snakemake (rule output)

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(readr)
  library(tibble)
  library(uwot)
  library(cluster)
})

source("scripts/lib_config.R")

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ------------------------------------------------------------------------------
# 1) CLI
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-c", "--config"), type = "character", help = "Path to config.yaml"),
  make_option("--profile", type = "character", default = NULL,
              help = "Config profile name (default: SNAKEMAKE_PROFILE or 'default')"),
  make_option("--feature_set", type = "character", help = "pam50 or MX"),
  make_option("--distance",    type = "character", help = "euclidean or correlation"),
  make_option("--outdir",      type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$config) || is.null(opt$feature_set) ||
    is.null(opt$distance) || is.null(opt$outdir)) {
  stop("Missing required arguments. Use --config, --feature_set, --distance, --outdir.", call. = FALSE)
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cat("=== clustering_driver.R ===\n")
cat("config:      ", opt$config,      "\n")
cat("feature_set: ", opt$feature_set, "\n")
cat("distance:    ", opt$distance,    "\n")
cat("outdir:      ", opt$outdir,      "\n\n")

# ------------------------------------------------------------------------------
# 2) Load config + data
# ------------------------------------------------------------------------------
# Load config using profiled config loader (merges defaults + profile)
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)

vst_path <- cfg$paths$vst_joint_rds %||% stop("paths$vst_joint_rds missing")
vst_mat <- readRDS(vst_path)  # genes x samples (ENSG with possible version suffixes)

# Strip version suffix from vst rownames if present (ENSG... .12)
rn <- rownames(vst_mat)
rn <- sub("\\.\\d+$", "", rn)
rownames(vst_mat) <- rn

feature_set <- match.arg(opt$feature_set, c("pam50", "MX"))
distance    <- match.arg(opt$distance, c("euclidean", "correlation"))

# PAM50: prefer Ensembl list if provided (matches vst_joint_rds rownames)
pam50_ens_file <- cfg$features$pam50_ensembl_gene_list
if (!is.null(pam50_ens_file) && file.exists(pam50_ens_file)) {
  pam50_genes <- readr::read_lines(pam50_ens_file)
} else {
  pam50_sym_file <- cfg$features$pam50_gene_list %||% stop("features$pam50_gene_list missing")
  pam50_genes <- readr::read_lines(pam50_sym_file)
}

# --- MX list only if needed ---
mx_genes <- NULL
if (feature_set == "MX") {
  mx_file <- cfg$features$mx_final_gene_list %||% stop("features$mx_final_gene_list missing")
  if (!file.exists(mx_file)) {
    stop("MX file missing for feature_set = 'MX': ", mx_file)
  }
  mx_genes <- readr::read_lines(mx_file)
}

# choose current gene set
if (feature_set == "pam50") {
  current_genes <- pam50_genes
  feat_label <- "PAM50"
} else if (feature_set == "MX") {
  current_genes <- mx_genes
  feat_label <- "MX500"
} else {
  stop("Unknown feature_set: ", feature_set)
}

sel_genes <- intersect(current_genes, rownames(vst_mat))

if (length(sel_genes) < 5) {
  stop("Too few genes after intersection for feature_set = ", feature_set,
       ". Got ", length(sel_genes), " genes.")
}

cat("[INFO] Selected genes for ", feature_set, ": ", length(sel_genes), "\n")

vst_sel <- vst_mat[sel_genes, , drop = FALSE]  # genes x samples
expr    <- t(vst_sel)                          # samples x genes

# k-grid from config
k_grid <- cfg$clustering$k_grid %||% c(2L, 3L, 4L, 5L, 6L)
k_grid <- sort(unique(as.integer(k_grid)))
cat("[INFO] k_grid: ", paste(k_grid, collapse = ", "), "\n")

# UMAP params from config (with defaults)
umap_cfg <- cfg$umap %||% list()
n_neighbors <- umap_cfg$n_neighbors %||% 20L
min_dist    <- umap_cfg$min_dist    %||% 0.3
metric_euc  <- umap_cfg$metric_euclidean %||% "euclidean"
metric_cor  <- umap_cfg$metric_correlation %||% "cosine"

# ------------------------------------------------------------------------------
# 3) PCA
# ------------------------------------------------------------------------------
cat("[INFO] Running PCA...\n")
pca <- prcomp(expr, center = TRUE, scale. = TRUE)
pc_scores <- pca$x  # samples x PCs

# Use first 20 PCs or fewer if not available
n_pcs <- min(20L, ncol(pc_scores))
pc_use <- pc_scores[, seq_len(n_pcs), drop = FALSE]

# ------------------------------------------------------------------------------
# 4) Distance matrix for silhouette (EUC or COR)
# ------------------------------------------------------------------------------
cat("[INFO] Computing distance matrix using ", distance, "...\n")

if (distance == "euclidean") {
  dmat <- dist(pc_use, method = "euclidean")
  dist_label <- "EUC"
} else {
  # correlation distance = 1 - cor
  cmat <- stats::cor(t(pc_use), method = "pearson", use = "pairwise.complete.obs")
  dmat <- as.dist(1 - cmat)
  dist_label <- "COR"
}

# ------------------------------------------------------------------------------
# 5) k-means over k_grid with silhouette selection
# ------------------------------------------------------------------------------
cat("[INFO] Running k-means grid search over k = {", paste(k_grid, collapse = ", "),
    "} with nstart=20...\n")

results <- lapply(k_grid, function(k) {
  set.seed(42)
  km <- kmeans(pc_use, centers = k, nstart = 20)
  sil <- cluster::silhouette(km$cluster, dmat)
  mean_sil <- mean(sil[, "sil_width"])
  data.frame(k = k, mean_sil = mean_sil)
})

sil_df <- bind_rows(results)
best_row <- sil_df[order(-sil_df$mean_sil, sil_df$k), ][1, ]
best_k <- best_row$k
best_sil <- best_row$mean_sil

cat("[INFO] Best k by silhouette: k =", best_k,
    " (mean silhouette =", round(best_sil, 3),
    "; ties resolved to smaller k)\n")

# Run final kmeans at best_k
set.seed(42)
km_final <- kmeans(pc_use, centers = best_k, nstart = 20)
clusters <- km_final$cluster
cluster_factor <- factor(clusters)

# ------------------------------------------------------------------------------
# 6) UMAP on PCs
# ------------------------------------------------------------------------------
cat("[INFO] Running UMAP...\n")
umap_metric <- if (distance == "euclidean") metric_euc else metric_cor

set.seed(42)
umap_coords <- uwot::umap(
  pc_use,
  n_neighbors = n_neighbors,
  min_dist    = min_dist,
  metric      = umap_metric
)

# ------------------------------------------------------------------------------
# 7) Output naming + writing
# ------------------------------------------------------------------------------
prefix <- sprintf("BRCA_TCGA-DSMZ_%s_%s", feat_label, dist_label)

clusters_path <- file.path(
  opt$outdir,
  sprintf("%s_CLUSTERS_kmeans_k%d_SAMPLES.csv", prefix, best_k)
)

pca_coord_path <- file.path(
  opt$outdir,
  sprintf("%s_PCA_coordinates_k%d_SAMPLES.tsv", prefix, best_k)
)

umap_coord_path <- file.path(
  opt$outdir,
  sprintf("%s_UMAP_coordinates_k%d_SAMPLES.tsv", prefix, best_k)
)

summary_path <- file.path(
  opt$outdir,
  sprintf("%s_CLUSTERING_SUMMARY_kmeans_k%d.txt", prefix, best_k)
)

# Clusters CSV
cat("[INFO] Writing clusters to: ", clusters_path, "\n")
cluster_df <- tibble(
  sample_id = rownames(expr),
  cluster   = as.integer(cluster_factor)
)
readr::write_csv(cluster_df, clusters_path)

# PCA coordinates
cat("[INFO] Writing PCA coordinates to: ", pca_coord_path, "\n")
pca_df <- as_tibble(pc_use, rownames = "sample_id") %>%
  mutate(cluster = as.integer(cluster_factor))
readr::write_tsv(pca_df, pca_coord_path)

# UMAP coordinates
cat("[INFO] Writing UMAP coordinates to: ", umap_coord_path, "\n")
umap_df <- tibble(
  sample_id = rownames(expr),
  UMAP1     = umap_coords[, 1],
  UMAP2     = umap_coords[, 2],
  cluster   = as.integer(cluster_factor)
)
readr::write_tsv(umap_df, umap_coord_path)

# Summary
cat("[INFO] Writing summary to: ", summary_path, "\n")
summary_lines <- c(
  sprintf("feature_set:   %s", feature_set),
  sprintf("distance:      %s", distance),
  sprintf("feat_label:    %s", feat_label),
  sprintf("dist_label:    %s", dist_label),
  sprintf("n_genes:       %d", length(sel_genes)),
  sprintf("n_samples:     %d", nrow(expr)),
  sprintf("k_grid:        %s", paste(k_grid, collapse = ", ")),
  sprintf("best_k:        %d", best_k),
  sprintf("mean_sil_best: %.4f", best_sil)
)
writeLines(summary_lines, summary_path)

cat("[INFO] Done.\n")
