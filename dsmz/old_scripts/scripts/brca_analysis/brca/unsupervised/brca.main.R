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
  library(SummarizedExperiment)
  library(tidyverse)
  library(DESeq2)
})

# --- Source your functions ---
files <- c(
  "io.R",
  "brca_io.R",
  "brca_base_functions.R",
  "tcga_base_functions.R",
  "dsmz_base_functions.R",
  "tcga_dsmz_base_functions.R",
  "normalization.R",
  "visualize.R",
  "helpers.R"
)
invisible(lapply(files, source, chdir = TRUE))

# === Load raw data ===
tcga_brca_counts_raw <- load_tcga_brca_data(config$tcga_se_rds, config$tcga_project_filter)
dsmz_brca_counts_raw <- load_dsmz_breast_data(config$dsmz_rds, config$dsmz_meta_csv, config$dsmz_organ_filter)

# === Subtype filtering (returns list) ===
tcga_brca_subtype_res <- load_tcga_brca_data_by_pam50(tcga_brca_counts_raw, config$tcga_brca_subtype_metadata)
tcga_brca_subtype_counts <- tcga_brca_subtype_res$tcga_subtype_counts
tcga_subtype_labels      <- tcga_brca_subtype_res$labels

# === Save intermediates ===
saveRDS(tcga_brca_counts_raw, file.path(config$outdir, "tcga_brca_counts.rds"))
saveRDS(dsmz_brca_counts_raw, file.path(config$outdir, "dsmz_brca_counts.rds"))
saveRDS(tcga_brca_subtype_counts, file.path(config$outdir, "tcga_brca_subtype_counts.rds"))
saveRDS(tcga_subtype_labels, file.path(config$outdir, "tcga_brca_subtype_labels.rds"))

# === Purity (uncomment only if verify_purity_data exists and config$purity_tsv is defined) ===
# purity <- verify_purity_data(tcga_brca_counts_raw, config$purity_tsv)

# === Harmonize gene IDs ===
tcga_brca_counts <- harmonize_gene_ids(tcga_brca_subtype_counts)
dsmz_brca_counts <- build_dsmz_matrix(dsmz_brca_counts_raw)  # dsmz_brca_counts_raw = filtered matrix with Ensembl_ID

# === Common genes ===
common <- common_genes(tcga_brca_counts, dsmz_brca_counts, min_shared = NULL)

# === Merge counts ===
merged_exp_data <- merge_counts(
  common$tcga_brca_subtype_counts,
  common$dsmz_brca_counts,
  dataset_names = c("TCGA-BRCA", "DSMZ-Breast")
)

# === VST per dataset (keep your function) ===
tcga_vst <- vst_per_dataset(
  countData = common$tcga_brca_subtype_counts,
  coldata = data.frame(row.names = colnames(common$tcga_brca_subtype_counts)),
  design = ~ 1
)

dsmz_vst <- vst_per_dataset(
  countData = common$dsmz_brca_counts,
  coldata = data.frame(row.names = colnames(common$dsmz_brca_counts)),
  design = ~ 1
)

# === Batch correction ===
seed <- 42  # <-- define seed
batch_corrected_counts <- batch_adjust(merged_exp_data$Xc_raw, merged_exp_data$batch, seed = seed)

# === coldata for merged data ===
coldata <- tibble(
  sample = colnames(merged_exp_data$Xc_raw),
  dataset = ifelse(
    colnames(merged_exp_data$Xc_raw) %in% colnames(common$tcga_brca_subtype_counts),
    "TCGA-BRCA",
    "DSMZ-Breast"
  )
) %>% column_to_rownames("sample")

# === VST on batch-corrected counts ===
vst_batch_corrected_counts <- vst_per_dataset(
  countData = batch_corrected_counts,
  coldata = coldata,
  design = ~ 1
)

# === Extract subset VSTs ===
tcga_vst_brca_subtype_batch_corrected <- vst_batch_corrected_counts[, colnames(tcga_vst), drop = FALSE]
dsmz_vst_batch_corrected <- vst_batch_corrected_counts[, colnames(dsmz_vst), drop = FALSE]

# === Visualization ===
pca_pre_BC <- make_pca_plot(
  merged_exp_data$Xc_raw, merged_exp_data$batch,
  "PCA Pre-Batch Correction",
  file.path(config$outdir, "pca_pre_BC.pdf"),
  scale = TRUE
)

pca_post_BC <- make_pca_plot(
  vst_batch_corrected_counts, coldata$dataset,
  "PCA Post-Batch Correction",
  file.path(config$outdir, "pca_post_BC.pdf"),
  scale = FALSE
)

umap_pre_BC <- make_umap(
  merged_exp_data$Xc_raw,
  "UMAP Pre-Batch Correction",
  file.path(config$outdir, "umap_pre_BC.pdf"),
  batch_labels = merged_exp_data$batch,
  scale = TRUE
)

umap_post_BC <- make_umap(
  vst_batch_corrected_counts,
  "UMAP Post-Batch Correction",
  file.path(config$outdir, "umap_post_BC.pdf"),
  batch_labels = coldata$dataset,
  scale = FALSE
)

# === Save final VSTs ===
saveRDS(dsmz_vst_batch_corrected, file.path(config$outdir, "dsmz_vst_brca_batch_corrected.rds"))
saveRDS(tcga_vst_brca_subtype_batch_corrected, file.path(config$outdir, "tcga_vst_brca_subtype_batch_corrected.rds"))