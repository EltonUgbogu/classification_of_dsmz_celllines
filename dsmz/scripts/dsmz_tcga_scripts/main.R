# main.R
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(uwot)
  library(ggplot2)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(matrixStats)
  library(sva)
  library(pheatmap)
})

# --- Source your functions ---
src <- function(...) lapply(file.path("R", c(...)), source, chdir = TRUE)
src("io.R",
    "tcga_base_functions.R",
    "dsmz_base_functions.R",
    "tcga_dsmz_base_functions.R",
    "normalization.R",
    "visualize.R",
    "helpers.R",
    "clustering_fixed.R",
    "clustering_dynamic.R",
    "clustering_hdbscan.R",
    "labeling_tcga.R")


tcga <- load_tcga_data(config$tcga_se_rds)
dsmz <- load_dsmz_data(config$dsmz_rds, config$dsmz_meta_csv)
outdir <- "/home/chu25/dsmz/results/unsupervised"
purity <- verify_purity_data(tcga$se, config$purity_tsv)
tcga_counts <- harmonize_gene_ids(tcga$counts)
dsmz_counts <- build_dsmz_matrix(dsmz$raw)
common <- common_genes(tcga_counts, dsmz_counts, min_shared=NULL)
merged_exp_data <- merge_counts(tcga_counts, dsmz_counts, dataset_names=c("TCGA", "DSMZ"))
tcga_vst <- vst_per_dataset(common$tcga_counts, coldata = data.frame(row.names = colnames(common$tcga_counts)), design = ~ 1)
dsmz_vst <- vst_per_dataset(common$dsmz_counts, coldata = data.frame(row.names = colnames(common$dsmz_counts)), design = ~ 1)
batch_corrected_counts <- batch_adjust(merged_exp_data$Xc_raw, merged_exp_data$batch, seed=config$seed)
coldata <- tibble(sample = colnames(merged_exp_data$Xc_raw), 
ifelse(colnames(merged_exp_data$Xc_raw) %in% colnames(tcga_vst), "TCGA", "DSMZ")) %>% 
 column_to_rownames("sample")
vst_batch_corrected_counts <- vst_per_dataset(batch_corrected_counts, coldata = coldata, design = ~ 1)
tcga_vst_batch_corrected <- vst_batch_corrected_counts[, colnames(tcga_vst), drop = FALSE]
dsmz_vst_batch_corrected <- vst_batch_corrected_counts[, colnames(dsmz_vst), drop = FALSE]
pca_pre_BC <- make_pca_plot(merged_exp_data$Xc_raw, merged_exp_data$batch, "PCA Pre-Batch Correction", file.path(outdir, "pca_pre_BC.pdf"), scale = TRUE)
pca_post_BC <- make_pca_plot(vst_batch_corrected_counts, coldata$dataset, "PCA Post-Batch Correction", file.path(outdir, "pca_post_BC.pdf"), scale = FALSE)
umap_pre_BC <- make_umap(merged_exp_data$Xc_raw, "UMAP Pre-Batch Correction", file.path(outdir, "umap_pre_BC.pdf"), batch_labels = merged_exp_data$batch, scale = TRUE)
umap_post_BC <- make_umap(vst_batch_corrected_counts, "UMAP Post-Batch Correction", file.path(outdir, "umap_post_BC.pdf"), batch_labels = coldata$dataset, scale = FALSE)
saveRDS(dsmz_vst_batch_corrected, file.path(outdir, "dsmz_vst_batch_corrected.rds"))
saveRDS(tcga_vst_batch_corrected, file.path(outdir, "tcga_vst_batch_corrected.rds"))


