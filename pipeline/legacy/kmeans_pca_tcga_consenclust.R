# kmeans_pca_tcga_consenclust.R
# Consensus k-means clustering on TCGA-BRCA PAM50 expression (PCA space) – TCGA-only

set.seed(42)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(ConsensusClusterPlus)
  library(uwot)
  library(ggplot2)
  library(readr)
  library(matrixStats)
  library(tibble)
})

# Source shared plotting helpers
source("R/brca_unsup_methods.R")

# ------------------------------------------------------------------
# Load config.yaml and base_functions helpers
# ------------------------------------------------------------------
cfg <- yaml::read_yaml("config/config.yaml")

if (is.null(cfg$paths$tumour_nh_base_functions)) {
  stop("Config must define paths$tumour_nh_base_functions")
}

base_fun_dir <- cfg$paths$tumour_nh_base_functions

if (!dir.exists(base_fun_dir)) {
  stop("Base-functions directory does not exist: ", base_fun_dir)
}

helper_files <- list.files(
  base_fun_dir,
  pattern = "\\.R$",
  full.names = TRUE
)

if (length(helper_files) == 0L) {
  stop("No .R helper files found in base_functions dir: ", base_fun_dir)
}

cat("[INFO] Loading", length(helper_files), "helper files from:", base_fun_dir, "\n")
invisible(lapply(helper_files, source))

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
outdir <- file.path(unsup_root, "inductive_tcga", "consensusclustering", "kmeans_pca_tcga")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", outdir, "\n")

# -------------------------------------------------------------------
# 1. Load TCGA PAM50 matrix (genes × samples)
# -------------------------------------------------------------------
cat("\n=== Loading TCGA-BRCA PAM50 matrix ===\n")

# Small helper for config fallback
`%||%` <- function(x, y) if (!is.null(x)) x else y

tcga_pam50_path <- cfg$paths$tcga_pam50_expr %||% cfg$paths$tcga_brca_pam50_expr

if (is.null(tcga_pam50_path)) {
  stop("Config must define either paths$tcga_pam50_expr or paths$tcga_brca_pam50_expr")
}

if (!file.exists(tcga_pam50_path)) {
  stop("TCGA PAM50 expr file not found at: ", tcga_pam50_path)
}

expr_mat <- readRDS(tcga_pam50_path)
cat(sprintf("TCGA: %d genes × %d samples\n", nrow(expr_mat), ncol(expr_mat)))

# Remove constant genes
expr_mat <- expr_mat[rowVars(expr_mat) > 0, , drop = FALSE]
cat("Retained", nrow(expr_mat), "variable PAM50 genes\n")

# Create alias for downstream code that expects tcga_brca_pam50
tcga_brca_pam50 <- expr_mat

# -------------------------------------------------------------------
# 2. PCA on TCGA matrix
# -------------------------------------------------------------------
cat("\n=== PCA on TCGA-BRCA PAM50 ===\n")
samples_mat <- t(expr_mat)  # samples × genes
pca_obj <- prcomp(samples_mat, center = TRUE, scale. = FALSE)
pcs_mat <- pca_obj$x

cat(sprintf("PCs matrix (samples × PCs): %s\n", paste(dim(pcs_mat), collapse = " × ")))

# -------------------------------------------------------------------
# 2a. Variance explained + save table
# -------------------------------------------------------------------
var_exp <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
cum_var <- cumsum(var_exp)
var_df <- data.frame(
  PC      = paste0("PC", seq_along(var_exp)),
  var_exp = var_exp,
  cum_var = cum_var
)
write.csv(
  var_df,
  file.path(outdir, "PCA_VarianceExplained_TCGA_PAM50.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------
# 2b. Scree plot
# -------------------------------------------------------------------
top_n <- min(20, nrow(var_df))
p_scree <- ggplot(var_df[1:top_n, ],
                  aes(x = factor(PC, levels = PC), y = var_exp)) +
  geom_col(fill = "#2c7bb6", width = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA variance explained (TCGA-BRCA PAM50)",
    x     = "Principal components",
    y     = "Proportion of variance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(outdir, "PCA_Screeplot_TCGA_PAM50.pdf"),
  p_scree,
  width = 7, height = 5, dpi = 300
)

# -------------------------------------------------------------------
# 2c. PCs dataframe
# -------------------------------------------------------------------
pcs_df <- as.data.frame(pcs_mat)
pcs_df$sample_id <- rownames(pcs_df)
pcs_df$dataset   <- "TCGA"
write.csv(pcs_df, file.path(outdir, "PCA_TCGA_PAM50_samples.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 2d. PCA scatter (using shared helper)
# -------------------------------------------------------------------
# Prepare metadata for plotting
meta_df_pca <- tibble(
  sample = rownames(pcs_mat),
  cluster = "All",  # No clusters yet at this stage
  dataset_type = "Tumour"
)

# Use shared helper for consistent plotting
save_pca_by_cluster(
  mat = tcga_brca_pam50,
  meta_df = meta_df_pca,
  out_pdf = file.path(outdir, "PCA_Scatterplot_TCGA_PAM50_PC1_PC2.pdf"),
  title = "TCGA-BRCA PCA (PAM50 genes)",
  subtitle = NULL,  # Will auto-generate variance explained
  label_dsmz = FALSE
)

# -------------------------------------------------------------------
# 3. ConsensusClusterPlus on PCA space (k-means)
# -------------------------------------------------------------------
cat("\n=== Consensus clustering (k-means) on PCA space ===\n")
data_for_ccp <- t(pcs_mat)  # PCs × samples

ccp_outdir <- file.path(outdir, "consensus")
dir.create(ccp_outdir, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd()
setwd(ccp_outdir)
on.exit(setwd(old_wd), add = TRUE)

ccp_results <- ConsensusClusterPlus(
  data_for_ccp,
  maxK         = 8,
  reps         = 1000,
  pItem        = 0.8,
  pFeature     = 1.0,
  clusterAlg   = "km",
  distance     = "euclidean",
  seed         = 42,
  title        = "TCGA_PCA_kmeans",
  plot         = "png",
  writeTable   = TRUE,
  verbose      = FALSE
)

# -------------------------------------------------------------------
# 3b. Select optimal K via ΔAUC
# -------------------------------------------------------------------
cat("\n=== Selecting optimal K ===\n")
auc_vals <- sapply(2:8, function(k) {
  cm <- ccp_results[[k]]$consensusMatrix
  vals <- cm[upper.tri(cm)]
  e <- ecdf(vals)
  x <- seq(0, 1, length.out = 1000)
  sum(diff(x) * (head(e(x), -1) + tail(e(x), -1)) / 2)
})
delta_auc <- diff(auc_vals)
best_k    <- which.max(delta_auc) + 2
cat("Selected optimal K =", best_k, "(max ΔAUC)\n")

delta_df <- data.frame(
  K         = 3:8,
  delta_AUC = delta_auc
)
p_delta <- ggplot(delta_df, aes(x = factor(K), y = delta_AUC)) +
  geom_col(fill = "#2c7bb6", width = 0.7) +
  geom_text(aes(label = sprintf("%.4f", delta_AUC)),
            vjust = -0.6, size = 3.5) +
  labs(
    title    = "Consensus Clustering (k-means): ΔAUC (TCGA PCA)",
    subtitle = sprintf("Optimal K = %d (maximum ΔAUC)", best_k),
    x        = "Number of clusters K",
    y        = "ΔAUC"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "#d7191c", face = "bold")
  )
ggsave(
  file.path(ccp_outdir, sprintf("DeltaAUC_optimalK%d_TCGA_PCA_kmeans.pdf", best_k)),
  p_delta,
  width = 7, height = 5.5, dpi = 300
)

# Final clusters
final_clusters <- setNames(
  paste0("C", ccp_results[[best_k]]$consensusClass),
  names(ccp_results[[best_k]]$consensusClass)
)

# Save results
results_df <- tibble(
  sample_id = names(final_clusters),
  cluster   = final_clusters,
  dataset   = "TCGA"
)
write_csv(results_df, file.path(ccp_outdir, sprintf("consensus_clusters_K%d.csv", best_k)))

# Also write canonical file (for consistency, though this script picks its own best_k)
canonical_file <- file.path(ccp_outdir, "consensus_clusters_K3.csv")
if (best_k == 3) {
  write_csv(results_df, canonical_file)
} else {
  # If best_k != 3, still write K3 file for Snakefile compatibility
  k3_clusters <- setNames(
    paste0("C", ccp_results[[3]]$consensusClass),
    names(ccp_results[[3]]$consensusClass)
  )
  k3_df <- tibble(
    sample_id = names(k3_clusters),
    cluster   = k3_clusters,
    dataset   = "TCGA"
  )
  write_csv(k3_df, canonical_file)
}

saveRDS(
  list(
    ccp_results = ccp_results,
    best_k      = best_k,
    clusters    = final_clusters,
    pcs_mat     = pcs_mat,
    auc_values  = auc_vals,
    delta_auc   = delta_auc,
    timestamp   = Sys.time()
  ),
  file.path(ccp_outdir, sprintf("full_consensus_clustering_K%d.rds", best_k))
)

# -------------------------------------------------------------------
# 4. UMAP of PCA space coloured by clusters (using shared helper)
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
  out_pdf = file.path(outdir, sprintf("UMAP_consensus_K%d.pdf", best_k)),
  title = "TCGA-BRCA (PAM50 genes) - UMAP",
  subtitle = sprintf("Samples coloured by consensus k-means cluster (K = %d)", best_k),
  n_neighbors = 15,
  min_dist = 0.3,
  metric = "euclidean",
  label_dsmz = FALSE
)

cat("\nConsensus k-means clustering on PCA space completed (TCGA-only).\n")
cat("Optimal K :", best_k, "\n")
cat("Results saved in:", outdir, "\n")
