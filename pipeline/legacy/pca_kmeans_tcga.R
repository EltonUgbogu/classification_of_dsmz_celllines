# pca_kmeans_tcga.R
# PCA of TCGA-BRCA (PAM50 expression) + k-means clustering (optimal-k) + UMAP

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
outdir <- file.path(unsup_root, "inductive_tcga", "pca_kmeans_tcga")
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
cat("\n=== PCA on TCGA PAM50 space ===\n")
tcga_samples_mat <- t(tcga_brca_pam50)  # samples x genes

pca_obj <- prcomp(tcga_samples_mat, center = TRUE, scale. = FALSE)
pcs_mat_all <- pca_obj$x

cat("PCs matrix (samples x PCs):",
    paste(dim(pcs_mat_all), collapse = " x "), "\n")

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
pcs_df <- as.data.frame(pcs_mat_all)
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
  sample = rownames(pcs_mat_all),
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
# 2e. Use top 20 PCs for clustering/UMAP
# -------------------------------------------------------------------
cat("\n=== Using top 20 PCs for clustering ===\n")
n_use <- min(20, ncol(pcs_mat_all))
pcs_mat <- pcs_mat_all[, 1:n_use, drop = FALSE]
cat("Using", n_use, "PCs for clustering/UMAP\n")

# -------------------------------------------------------------------
# 3. K-means clustering on PCA space (optimal-k search)
# -------------------------------------------------------------------
cat("\n=== K-means on PCA space ===\n")
kmeans_dir <- file.path(outdir, "patient_clustering_kmeans_tcga")
dir.create(kmeans_dir, recursive = TRUE, showWarnings = FALSE)

k_grid <- 2:6
sil_scores <- numeric(length(k_grid))
km_results <- list()

for (i in seq_along(k_grid)) {
  k <- k_grid[i]
  km <- kmeans(pcs_mat, centers = k, nstart = 50, iter.max = 100)
  km_results[[as.character(k)]] <- km
  
  # Calculate silhouette
  d <- dist(pcs_mat, method = "euclidean")
  sil <- cluster::silhouette(km$cluster, d)
  sil_scores[i] <- mean(sil[, "sil_width"])
}

best_k <- k_grid[which.max(sil_scores)]
km_best <- km_results[[as.character(best_k)]]
clusters_labeled <- paste0("C", km_best$cluster)
names(clusters_labeled) <- rownames(pcs_mat)

cat("Optimal k-means k =", best_k, "\n")

# Save final clusters
cluster_df <- data.frame(
  sample = names(clusters_labeled),
  cluster = clusters_labeled
)
write.csv(cluster_df, file.path(kmeans_dir, "final_clusters.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 4. UMAP on PCA space (sample embeddings) - using shared helper
# -------------------------------------------------------------------
cat("\n=== UMAP visualization ===\n")

# Prepare metadata with clusters
meta_df_umap <- tibble(
  sample = names(clusters_labeled),
  cluster = factor(clusters_labeled),
  dataset_type = "Tumour"
)

# Use shared helper for consistent plotting
save_umap_by_cluster(
  mat = tcga_brca_pam50[, meta_df_umap$sample, drop = FALSE],
  meta_df = meta_df_umap,
  out_pdf = file.path(outdir, sprintf("PCA_TCGA_BRCA_PAM50_UMAP_samples_kmeans_optimal_k%d.pdf", best_k)),
  title = "TCGA-BRCA (PAM50 genes) - UMAP",
  subtitle = sprintf("Samples coloured by k-means cluster (k = %d)", best_k),
  n_neighbors = 15,
  min_dist = 0.3,
  metric = "euclidean",
  label_dsmz = FALSE
)

cat("\n[SUCCESS] TCGA-only PCA + k-means clustering complete!\n")
cat("  Output directory:", outdir, "\n")
cat("  Optimal k:", best_k, "\n")
