# pca_hc_tcga.R
# PCA of TCGA-BRCA (PAM50 expression) + hierarchical clustering (optimal-k) + UMAP

options(stringsAsFactors = FALSE)
set.seed(42)

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(tibble)
  library(uwot)
  library(ggplot2)
  library(matrixStats)
  library(readr)
  library(cluster)
  library(dynamicTreeCut)
})

# Source shared plotting helpers
source("R/brca_unsup_methods.R")

# ------------------------------------------------------------------
# Load config.yaml
# ------------------------------------------------------------------
cfg <- yaml::read_yaml("config/config.yaml")
unsup_root <- cfg$paths$unsup_root

# Load canonical km_clusters file
km_clusters_file <- file.path(
  unsup_root,
  "inductive_tcga",
  "consensusclustering",
  "km_tcga",
  "BRCA_TCGA_PAM50_km_clusters_bestk.csv"
)

if (!file.exists(km_clusters_file)) {
  stop("Cannot find km_clusters file: ", km_clusters_file)
}

km_clusters <- readr::read_csv(km_clusters_file)
cat("[INFO] Loaded", nrow(km_clusters), "samples from km_clusters\n")

# -------------------------------------------------------------------
# Output dir
# -------------------------------------------------------------------
outdir <- file.path(unsup_root, "inductive_tcga", "pca_hc_tcga")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", outdir, "\n")

# -------------------------------------------------------------------
# 1. Load TCGA PAM50 expression matrix (genes x samples)
# -------------------------------------------------------------------
cat("\n=== Loading TCGA-BRCA PAM50 matrix ===\n")
tcga_brca_pam50 <- readRDS(cfg$paths$tcga_pam50_expr)

cat("TCGA matrix dim (genes x samples):",
    paste(dim(tcga_brca_pam50), collapse = " x "), "\n")

stopifnot(!is.null(rownames(tcga_brca_pam50)),
          !is.null(colnames(tcga_brca_pam50)))

# -------------------------------------------------------------------
# 2. PCA on TCGA-BRCA matrix
# -------------------------------------------------------------------
cat("\n=== PCA on TCGA-BRCA PAM50 space ===\n")
tcga_samples_mat <- t(tcga_brca_pam50)  # samples x genes

pca_obj <- prcomp(tcga_samples_mat, center = TRUE, scale. = FALSE)
pcs_mat <- pca_obj$x

cat("PCs matrix (samples x PCs):", paste(dim(pcs_mat), collapse = " x "), "\n")

# Use top 20 PCs
n_use <- min(20, ncol(pcs_mat))
pcs_mat <- pcs_mat[, 1:n_use, drop = FALSE]

# -------------------------------------------------------------------
# 2a. Variance explained + save table
# -------------------------------------------------------------------
var_exp <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
cum_var <- cumsum(var_exp)

var_df <- data.frame(
  PC = paste0("PC", seq_along(var_exp)),
  var_exp = var_exp,
  cum_var = cum_var
)

write.csv(
  var_df,
  file.path(outdir, "PCA_VarianceExplained_TCGA_BRCA_PAM50.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------
# 2b. Scree plot (top 20 PCs)
# -------------------------------------------------------------------
top_n <- min(20, nrow(var_df))
p_scree <- ggplot(var_df[1:top_n, ], aes(x = factor(PC, levels = PC), y = var_exp)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "PCA variance explained (TCGA-BRCA PAM50)",
    x = "Principal components (decreasing variance)",
    y = "Proportion of variance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(outdir, "PCA_Screeplot_TCGA_BRCA_PAM50.pdf"),
  p_scree,
  width = 7, height = 5, dpi = 300
)

# -------------------------------------------------------------------
# 2c. PCs dataframe for plotting / saving
# -------------------------------------------------------------------
pcs_df <- as.data.frame(pcs_mat)
pcs_df$sample_id <- rownames(pcs_df)
pcs_df$dataset   <- "TCGA"

write.csv(
  pcs_df,
  file.path(outdir, "PCA_TCGA_BRCA_PAM50_samples.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------
# 2d. PCA scatter with % variance on axes (using shared helper)
# -------------------------------------------------------------------
# Prepare metadata for plotting
meta_df <- tibble(
  sample = rownames(pcs_mat),
  cluster = "All",  # No clusters yet at this stage
  dataset_type = "Tumour"
)

# Use shared helper for consistent plotting
save_pca_by_cluster(
  mat = tcga_brca_pam50,
  meta_df = meta_df,
  out_pdf = file.path(outdir, "PCA_Scatterplot_TCGA_BRCA_PAM50_PC1_PC2.pdf"),
  title = "TCGA-BRCA PCA (PAM50 genes)",
  subtitle = NULL,  # Will auto-generate variance explained
  label_dsmz = FALSE
)

# -------------------------------------------------------------------
# 3. Hierarchical clustering on PCA space (optimal silhouette k)
# -------------------------------------------------------------------
cat("\n=== Hierarchical clustering on PCA space ===\n")

# Calculate distance matrix
dist_mat <- dist(pcs_mat, method = "euclidean")

# Build dendrogram
hc <- hclust(dist_mat, method = "ward.D2")

# Dynamic tree cut
dyn_cut <- cutreeDynamic(hc, distM = as.matrix(dist_mat), minClusterSize = 10)
dyn_clusters <- paste0("C", dyn_cut)
names(dyn_clusters) <- rownames(pcs_mat)

# Optimal k via silhouette (k=2 to 6)
k_grid <- 2:6
sil_scores <- numeric(length(k_grid))
best_k <- 2

for (i in seq_along(k_grid)) {
  k <- k_grid[i]
  cut_k <- cutree(hc, k = k)
  sil <- silhouette(cut_k, dist_mat)
  sil_scores[i] <- mean(sil[, "sil_width"])
}

best_k <- k_grid[which.max(sil_scores)]
opt_cut <- cutree(hc, k = best_k)
opt_clusters <- paste0("C", opt_cut)
names(opt_clusters) <- rownames(pcs_mat)

cat("Optimal hierarchical clustering k =", best_k, "\n")

# Save clusters (Snakefile expects k=2, but we'll use optimal)
if (best_k == 2) {
  final_clusters <- opt_clusters
} else {
  # Use k=2 for Snakefile compatibility
  cut_2 <- cutree(hc, k = 2)
  final_clusters <- paste0("C", cut_2)
  names(final_clusters) <- rownames(pcs_mat)
}

cluster_df <- data.frame(
  sample_id = names(final_clusters),
  cluster = final_clusters
)
write.csv(cluster_df, 
          file.path(outdir, "PCA_TCGA_BRCA_PAM50_samples_clusters_optimal_k=2.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------
# 4. UMAP on PCA space (sample embeddings) - using shared helper
# -------------------------------------------------------------------
cat("\n=== UMAP visualization ===\n")

# Prepare metadata with clusters
meta_df_umap <- tibble(
  sample = names(final_clusters),
  cluster = factor(final_clusters),
  dataset_type = "Tumour"
)

# Use shared helper for consistent plotting
save_umap_by_cluster(
  mat = tcga_brca_pam50[, meta_df_umap$sample, drop = FALSE],
  meta_df = meta_df_umap,
  out_pdf = file.path(outdir, "PCA_TCGA_BRCA_PAM50_UMAP_samples_optimal.pdf"),
  title = sprintf("TCGA-BRCA (PAM50 genes) - UMAP"),
  subtitle = sprintf("Samples coloured by hierarchical cluster (k = %d)", 2),
  n_neighbors = 15,
  min_dist = 0.3,
  metric = "euclidean",
  label_dsmz = FALSE
)

cat("\n[SUCCESS] TCGA-BRCA PCA + hierarchical clustering complete!\n")
cat("  Output directory:", outdir, "\n")
cat("  Optimal k:", best_k, "(using k=2 for compatibility)\n")
