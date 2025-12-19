# km_tcga_consensusclustering.R
# Consensus clustering on TCGA-BRCA PAM50 only (k-means)

set.seed(42)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ConsensusClusterPlus)
  library(matrixStats)
  library(uwot)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# Source shared plotting helpers
source("R/brca_unsup_methods.R")

stop_with <- function(...) stop(paste0(...), call. = FALSE)
info      <- function(...) cat("[INFO] ", sprintf(...), "\n", sep = "")
stage     <- function(...) cat("\n[", sprintf(...), "]\n", sep = "")

# Load config
suppressPackageStartupMessages(library(yaml))
cfg <- yaml::read_yaml("config/config.yaml")

# Set output directory
outdir <- file.path(cfg$paths$unsup_root, "inductive_tcga", "consensusclustering", "km_tcga")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

stage("LOAD TCGA PAM50 expression")

tcga_mat <- readRDS(cfg$paths$tcga_pam50_expr)
info("TCGA dim (genes x samples): %s", paste(dim(tcga_mat), collapse=" x "))

stopifnot(!is.null(rownames(tcga_mat)), !is.null(colnames(tcga_mat)))

stage("PREPROCESS")

v <- rowVars(tcga_mat, na.rm = TRUE)
expr_mat <- tcga_mat[v > 0 & is.finite(v), , drop = FALSE]
info("After removing constant genes: %s", paste(dim(expr_mat), collapse=" x "))

if (anyDuplicated(colnames(expr_mat))) {
  stop_with("Duplicate TCGA sample IDs found. Fix colnames before CCP.")
}

stage("CCP k-means")

old_wd <- getwd()
setwd(outdir)
on.exit(setwd(old_wd), add = TRUE)

ccp_res <- ConsensusClusterPlus(
  expr_mat,
  maxK       = 8,
  reps       = 1000,
  pItem      = 0.8,
  pFeature   = 1.0,
  clusterAlg = "km",
  distance   = "euclidean",
  seed       = 42,
  title      = "TCGA_BRCA_PAM50_KM",
  plot       = "png",
  writeTable = TRUE,
  verbose    = TRUE
)

saveRDS(ccp_res, file.path(outdir, "tcga_km_ccp_res.rds"))
info("Saved CCP object: tcga_km_ccp_res.rds")

stage("ICL")

icl <- calcICL(ccp_res, title="TCGA_BRCA_PAM50_KM", plot="png")
write.csv(icl$clusterConsensus,
          file.path(outdir, "clusterConsensus.csv"), row.names=FALSE)
write.csv(icl$itemConsensus,
          file.path(outdir, "itemConsensus.csv"), row.names=FALSE)
info("ICL tables written.")

stage("EXPORT labels all K")

labels_df <- bind_rows(lapply(2:8, function(k){
  cl <- ccp_res[[k]]$consensusClass
  data.frame(sample_id = names(cl),
             cluster   = paste0("C", cl),
             k         = k,
             stringsAsFactors = FALSE)
}))
write.csv(labels_df, file.path(outdir, "consensus_labels_all_k.csv"), row.names=FALSE)

stage("SELECT best K")

auc_values <- sapply(2:8, function(k) {
  cm <- ccp_res[[k]]$consensusMatrix
  vals <- cm[upper.tri(cm)]
  ec <- ecdf(vals)
  x <- seq(0, 1, length.out=500)
  sum(diff(x) * (ec(head(x,-1)) + ec(tail(x,-1))) / 2)
})
delta_auc <- diff(auc_values)
best_k <- which.max(delta_auc) + 2

info("AUC (K=2..8): %s", paste(round(auc_values,4), collapse=" "))
info("Delta AUC: %s", paste(round(delta_auc,4), collapse=" "))
info("Selected best K: %d", best_k)

final_clusters <- ccp_res[[best_k]]$consensusClass
final_clusters <- setNames(paste0("C", final_clusters), names(final_clusters))

# Create cluster data.frame
cluster_df <- tibble(
  sample_id = names(final_clusters),
  cluster   = final_clusters
)

# Write canonical "best-k" file (independent of actual K value)
out_clusters <- file.path(outdir, "BRCA_TCGA_PAM50_km_clusters_bestk.csv")
write_csv(cluster_df, out_clusters)
info("Saved best-k clusters to: %s", out_clusters)

# Write metadata file with chosen K
meta_df <- tibble(
  method = "kmeans_tcga",
  best_k = best_k
)
out_meta <- file.path(outdir, "BRCA_TCGA_PAM50_km_best_k.tsv")
write_tsv(meta_df, out_meta)
info("Saved best_k metadata to: %s", out_meta)

# Also keep the K-specific file for debugging/reference
write_csv(
  cluster_df,
  file.path(outdir, sprintf("tcga_km_clusters_K%d.csv", best_k))
)

stage("UMAP")

pca_obj <- prcomp(t(expr_mat), center=TRUE, scale.=FALSE)
pcs <- pca_obj$x[, 1:min(20, ncol(pca_obj$x)), drop=FALSE]

um <- uwot::umap(pcs, n_neighbors=15, min_dist=0.3)

umap_df <- data.frame(
  UMAP1 = um[,1],
  UMAP2 = um[,2],
  sample_id = rownames(pcs),
  cluster = final_clusters[rownames(pcs)]
)

# Use shared helper for consistent plotting
meta_df_umap <- tibble(
  sample = names(final_clusters),
  cluster = factor(final_clusters),
  dataset_type = "Tumour"
)

save_umap_by_cluster(
  mat = expr_mat[, meta_df_umap$sample, drop = FALSE],
  meta_df = meta_df_umap,
  out_pdf = file.path(outdir, sprintf("UMAP_tcga_km_K%d.pdf", best_k)),
  title = "TCGA-BRCA (PAM50 genes) - UMAP",
  subtitle = sprintf("Samples coloured by consensus k-means cluster (K = %d)", best_k),
  n_neighbors = 15,
  min_dist = 0.3,
  metric = "euclidean",
  label_dsmz = FALSE
)

info("Done. Results in: %s", outdir)
