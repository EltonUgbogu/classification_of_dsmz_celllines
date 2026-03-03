# ==============================================================================
# nbl_main.R – TCGA-BRCA + DSMZ Breast Cell Lines (with PAM50 Subtypes)
# ==============================================================================
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

# --- Source all helper functions ------------------------------------------------
files <- c(
  "nbl_io.R",
  "tcga_base_functions.R",
  "dsmz_base_functions.R",
  "tcga_dsmz_base_functions.R",
  "normalization.R",
  "visualize.R",
  "helpers.R"
)
invisible(lapply(files, source, chdir = TRUE))

# ==============================================================================
# 1. Load raw data
# ==============================================================================
tcga_nbl_counts_raw <- load_tcga_nbl_data(
  config$tcga_se_rds,
  config$tcga_project_filter
)

dsmz_counts_raw <- load_dsmz_data(
  config$dsmz_rds,
  config$dsmz_meta_csv
)

# ==============================================================================
# 2. Subtype filtering (TCGA) & DSMZ breast-cell line selection
# ==============================================================================
tcga_nbl_subtype_res <- load_tcga_nbl_data_by_pam50(
  tcga_nbl_counts_raw,
  config$tcga_nbl_subtype_metadata
)

tcga_nbl_subtype_counts <- tcga_nbl_subtype_res$tcga_subtype_counts
tcga_subtype_labels      <- tcga_nbl_subtype_res$labels

tcga_nbl_counts <- harmonize_gene_ids(tcga_nbl_subtype_counts)

dsmz_counts <- build_dsmz_matrix(dsmz_counts_raw$raw)

dsmz_nbl_counts <- load_dsmz_breast_data(
  dsmz_counts,
  config$dsmz_meta_csv,
  config$dsmz_organ_filter
)

# ==============================================================================
# 3. Save intermediate count matrices
# ==============================================================================
saveRDS(tcga_nbl_counts_raw,
        file.path(config$outdir, "tcga_nbl_counts.rds"))
saveRDS(tcga_nbl_subtype_counts,
        file.path(config$outdir, "tcga_nbl_subtype_counts.rds"))
saveRDS(tcga_subtype_labels,
        file.path(config$outdir, "tcga_nbl_subtype_labels.rds"))
saveRDS(dsmz_nbl_counts,
        file.path(config$outdir, "dsmz_nbl_counts.rds"))

# ==============================================================================
# 4. Identify common genes (for downstream individual VST)
# ==============================================================================
common <- common_genes(
  tcga_counts = tcga_nbl_counts,
  dsmz_counts = dsmz_nbl_counts,
  min_shared  = 5000
)

# ==============================================================================
# 5. Merge the two datasets (full gene union)
# ==============================================================================
merged_exp_data <- merge_counts(
  tcga_counts   = common$tcga_counts,
  dsmz_counts   = common$dsmz_counts,
  dataset_names = c("TCGA-BRCA-subtype", "DSMZ-BCC")
)

# ==============================================================================
# 6. Extract merged data
# ==============================================================================
counts_raw   <- merged_exp_data$Xc_raw
batch_labels <- merged_exp_data$batch
sample_ids   <- merged_exp_data$sample_id

# --- Debug: Print first few entries ---
cat("\n=== DEBUG: Merged Data ===\n")
cat("First 5 rows of counts_raw (gene x sample):\n")
print(counts_raw[1:5, 1:5])
cat("\nFirst 5 batch labels:\n")
print(head(batch_labels, 5))
cat("\nFirst 5 sample_ids:\n")
print(head(sample_ids, 5))

# ==============================================================================
# 7. Build coldata + add DSMZ cell line names + TCGA subtypes
# ==============================================================================
sample_names <- names(sample_ids)

# --- Build tibble ---
coldata <- tibble(
  sample  = sample_ids,
  dataset = batch_labels[sample_names]  # <-- values + names preserved
) %>%
  column_to_rownames("sample")  # rownames = sample names

# --- CRITICAL: Ensure dataset is a NAMED factor ---
coldata$dataset <- batch_labels[rownames(coldata)]  # <-- RE-ASSIGN WITH NAMES

# --- Add DSMZ cell line names ---
dsmz_meta <- read.csv(config$dsmz_meta_csv, stringsAsFactors = FALSE)
dsmz_samples <- rownames(coldata)[coldata$dataset == "DSMZ-BCC"]

if (length(dsmz_samples) > 0) {
  if (!all(c("Cell_Line", "sample_name") %in% colnames(dsmz_meta))) {
    stop("Missing 'Cell_Line' or 'sample_name' in DSMZ metadata")
  }
  cell_line_map <- setNames(dsmz_meta$Cell_Line, dsmz_meta$sample_name)
  coldata$cell_line <- NA_character_
  coldata[dsmz_samples, "cell_line"] <- cell_line_map[dsmz_samples]
}

# --- Add TCGA PAM50 subtypes ---
tcga_samples <- rownames(coldata)[coldata$dataset == "TCGA-BRCA-subtype"]
coldata$subtype <- NA_character_
coldata[tcga_samples, "subtype"] <- tcga_subtype_labels[tcga_samples]

# --- DEBUG ---
cat("\n=== First 5 rows of coldata (FINAL) ===\n")
print(head(coldata, 5))
cat("\nNames of coldata$dataset (first 5):\n")
print(head(names(coldata$dataset), 5))
cat("Class: ", class(coldata$dataset), "\n")


# ==============================================================================
# 9. Batch correction (ComBat-Seq)
# ==============================================================================
batch_corrected_counts <- batch_adjust(
  Exp_mat      = counts_raw,
  batch_factor = batch_labels
)

# ==============================================================================
# 10. colData for individual datasets
# ==============================================================================
tcga_nbl_subtype_coldata <- data.frame(row.names = colnames(common$tcga_counts))
dsmz_nbl_coldata         <- data.frame(row.names = colnames(common$dsmz_counts))

# ==============================================================================
# 11. VST on original (pre-merge) matrices
# ==============================================================================
tcga_nbl_subtype_vst <- vst_normalize(
  Exp_mat = common$tcga_counts,
  coldata = tcga_nbl_subtype_coldata
)

dsmz_nbl_vst <- vst_normalize(
  Exp_mat = common$dsmz_counts,
  coldata = dsmz_nbl_coldata
)
# merge vst normalized dataset
vst_non_batch_corrected_counts <- cbind(tcga_nbl_subtype_vst, dsmz_nbl_vst)
# ==============================================================================
# 12. VST on batch-corrected merged matrix
# ==============================================================================
vst_batch_corrected_counts <- vst_normalize(
  Exp_mat = batch_corrected_counts,
  coldata = coldata
)

# ==============================================================================
# 13. Extract batch-corrected VSTs per dataset
# ==============================================================================
tcga_vst_nbl_subtype_batch_corrected <- vst_batch_corrected_counts[
  , colnames(tcga_nbl_subtype_vst), drop = FALSE
]

dsmz_vst_nbl_batch_corrected <- vst_batch_corrected_counts[
  , colnames(dsmz_nbl_vst), drop = FALSE
]

# ==============================================================================
# 15. Save batch-corrected VST matrices
# ==============================================================================
saveRDS(dsmz_vst_nbl_batch_corrected,
        file.path(config$outdir, "dsmz_vst_nbl_batch_corrected.rds"))
saveRDS(tcga_vst_nbl_subtype_batch_corrected,
        file.path(config$outdir, "tcga_vst_nbl_subtype_batch_corrected.rds"))

# ==============================================================================
# 16. Visualisation
# ==============================================================================

pca_pre_BC <- make_pca_plot(
  Exp_mat      = vst_non_batch_corrected_counts,
  batch_labels = batch_labels,
  title        = "PCA Pre-Batch Correction",
  path_pdf     = file.path(config$outdir, "pca_pre_BC.pdf"),
  center       = TRUE,
  scale        = TRUE
)

pca_post_BC <- make_pca_plot(
  Exp_mat      = vst_batch_corrected_counts,
  batch_labels = batch_labels,
  title        = "PCA Post-Batch Correction",
  path_pdf     = file.path(config$outdir, "pca_post_BC.pdf"),
  center       = TRUE,
  scale        = FALSE
)

umap_pre_BC <- make_umap(
  Exp_mat      = counts_raw,
  title        = "UMAP Pre-Batch Correction",
  path_pdf     = file.path(config$outdir, "umap_pre_BC.pdf"),
  batch_labels = batch_labels,
  center       = TRUE,
  scale        = TRUE
)

umap_post_BC <- make_umap(
  Exp_mat      = vst_batch_corrected_counts,
  title        = "UMAP Post-Batch Correction",
  path_pdf     = file.path(config$outdir, "umap_post_BC.pdf"),
  batch_labels = batch_labels,
  center       = TRUE,
  scale        = FALSE
)

# ==============================================================================
# 17. QC diagnostics
# ==============================================================================
qc_vst_diagnostics(
  counts    = common$tcga_counts,
  coldata   = tcga_nbl_subtype_coldata,
  title     = "TCGA-BRCA_ST_RAW",
  outdir    = config$outdir,
  full_disp = TRUE
)

qc_vst_diagnostics(
  counts    = dsmz_nbl_counts,
  coldata   = dsmz_nbl_coldata,
  title     = "DSMZ-BCC_RAW",
  outdir    = config$outdir,
  full_disp = TRUE
)

# ==============================================================================
# Done
# ==============================================================================
message("[SUCCESS] nbl_main.R completed – all outputs in: ", config$outdir)