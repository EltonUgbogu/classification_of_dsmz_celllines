root_dir <- normalizePath(file.path(snakemake@scriptdir, "..", ".."), mustWork = TRUE)
source(file.path(root_dir, "helpers.R"))
source(file.path(root_dir, "dsmz_base_functions.R"))
source(file.path(root_dir, "rbl_io.R"))
source(file.path(root_dir, "rbl_tumour_dsmz_base_functions.R"))
source(file.path(root_dir, "normalisation.R"))
source(file.path(root_dir, "visualisation.R"))
source(file.path(root_dir, "tumour_base_functions.R"))

if (length(snakemake@log) > 0) {
  log_file <- snakemake@log[[1]]
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)
}

for (path in unname(snakemake@output)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
if (length(snakemake@log) > 0) {
  dir.create(dirname(snakemake@log[[1]]), recursive = TRUE, showWarnings = FALSE)
}

params <- snakemake@params
inputs <- snakemake@input
outputs <- snakemake@output

cat("[INFO] Starting RBL preprocess workflow\n")
cat(sprintf("[INFO] Tumour input: %s\n", inputs[["tumour_rds"]]))
cat(sprintf("[INFO] DSMZ input: %s\n", inputs[["dsmz_count_matrix"]]))
cat(sprintf("[INFO] DSMZ metadata: %s\n", inputs[["dsmz_sample_metadata"]]))

# Load and filter inputs.
tumour_counts <- load_tumour_rbl_data(inputs[["tumour_rds"]])
dsmz_data <- load_dsmz_data(
  inputs[["dsmz_count_matrix"]],
  inputs[["dsmz_sample_metadata"]],
  sample_col = params[["dsmz_sample_col"]]
)
dsmz_filtered <- load_dsmz_rbl_data(
  counts = dsmz_data$counts,
  metadata = dsmz_data$meta,
  filter_col = params[["dsmz_filter_col"]],
  filter_values = unlist(params[["dsmz_filter_values"]]),
  sample_col = params[["dsmz_sample_col"]],
  outdir = file.path(dirname(outputs[["dsmz_counts_rds"]]), "audit")
)

saveRDS(tumour_counts, outputs[["tumour_counts_rds"]])
saveRDS(dsmz_filtered$counts, outputs[["dsmz_counts_rds"]])

# Merge on common genes.
merged <- merge_counts(
  tumour_counts,
  dsmz_filtered$counts,
  dataset_names = c(params[["tumour_label"]], params[["dsmz_label"]]),
  min_shared = params[["min_shared"]]
)
saveRDS(merged$Xc_raw, outputs[["merged_counts_rds"]])

# Build coldata.
coldata <- data.frame(
  dataset = unname(merged$batch),
  row.names = names(merged$batch),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
coldata$cell_line <- NA_character_
if (params[["dsmz_cellline_col"]] %in% colnames(dsmz_filtered$metadata)) {
  cell_line_map <- setNames(dsmz_filtered$metadata[[params[["dsmz_cellline_col"]]]], dsmz_filtered$metadata$sample_id)
  dsmz_samples <- rownames(coldata)[coldata$dataset == params[["dsmz_label"]]]
  coldata[dsmz_samples, "cell_line"] <- cell_line_map[dsmz_samples]
}
coldata$subtype <- NA_character_
subtype_file <- params[["subtype_metadata"]]
if (!is.null(subtype_file) && nzchar(subtype_file) && file.exists(subtype_file)) {
  subtype_tab <- utils::read.csv(subtype_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (all(c(params[["subtype_sample_col"]], params[["subtype_col"]]) %in% colnames(subtype_tab))) {
    subtype_map <- setNames(subtype_tab[[params[["subtype_col"]]]], subtype_tab[[params[["subtype_sample_col"]]]])
    tumour_samples <- rownames(coldata)[coldata$dataset == params[["tumour_label"]]]
    coldata[tumour_samples, "subtype"] <- subtype_map[tumour_samples]
  }
}
utils::write.table(coldata, file = outputs[["coldata_tsv"]], sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Batch correction.
batch_corrected_counts <- batch_adjust(merged$Xc_raw, merged$batch)
saveRDS(batch_corrected_counts, outputs[["batch_corrected_counts_rds"]])

# VST before and after batch correction.
tumour_coldata <- coldata[colnames(merged$tumour_counts), , drop = FALSE]
dsmz_coldata <- coldata[colnames(merged$dsmz_counts), , drop = FALSE]

tumour_vst <- vst_normalize(merged$tumour_counts, tumour_coldata)
dsmz_vst <- vst_normalize(merged$dsmz_counts, dsmz_coldata)
joint_vst_pre <- cbind(tumour_vst, dsmz_vst)
joint_vst_post <- vst_normalize(batch_corrected_counts, coldata)

tumour_vst_post <- joint_vst_post[, colnames(merged$tumour_counts), drop = FALSE]
dsmz_vst_post <- joint_vst_post[, colnames(merged$dsmz_counts), drop = FALSE]

saveRDS(tumour_vst, outputs[["tumour_vst_rds"]])
saveRDS(dsmz_vst, outputs[["dsmz_vst_rds"]])
saveRDS(joint_vst_pre, outputs[["joint_vst_pre_bc_rds"]])
saveRDS(joint_vst_post, outputs[["joint_vst_post_bc_rds"]])
saveRDS(tumour_vst_post, outputs[["tumour_vst_post_bc_rds"]])
saveRDS(dsmz_vst_post, outputs[["dsmz_vst_post_bc_rds"]])

collapse_col <- "DSMZ_Cell_line_norm"
if (!collapse_col %in% colnames(dsmz_filtered$metadata)) {
  stop(sprintf("[ERROR] DSMZ metadata column required for cell-line collapse not found: %s", collapse_col))
}
if (!"sample_id" %in% colnames(dsmz_filtered$metadata)) {
  stop("[ERROR] DSMZ filtered metadata is missing sample_id")
}
cell_line_map <- stats::setNames(dsmz_filtered$metadata[[collapse_col]], dsmz_filtered$metadata$sample_id)
cell_line_key <- unname(cell_line_map[colnames(dsmz_vst_post)])
if (anyNA(cell_line_key) || any(!nzchar(cell_line_key))) {
  missing_samples <- colnames(dsmz_vst_post)[is.na(cell_line_key) | !nzchar(cell_line_key)]
  stop(sprintf("[ERROR] Missing DSMZ cell-line collapse key for samples: %s", paste(missing_samples, collapse = ", ")))
}
collapsed_celllines <- unique(cell_line_key)
dsmz_vst_post_collapsed <- vapply(
  collapsed_celllines,
  function(cell_line) {
    rowMeans(dsmz_vst_post[, cell_line_key == cell_line, drop = FALSE])
  },
  numeric(nrow(dsmz_vst_post))
)
rownames(dsmz_vst_post_collapsed) <- rownames(dsmz_vst_post)
colnames(dsmz_vst_post_collapsed) <- collapsed_celllines
saveRDS(dsmz_vst_post_collapsed, outputs[["dsmz_vst_post_bc_collapsed_cellline_rds"]])
cat(sprintf("[INFO] DSMZ post-BC libraries: %d; collapsed cell lines: %d\n", ncol(dsmz_vst_post), ncol(dsmz_vst_post_collapsed)))

# Plots.
make_pca_plot(joint_vst_pre, merged$batch[colnames(joint_vst_pre)], "PCA Pre-Batch Correction", outputs[["pca_pre_pdf"]], center = TRUE, scale = FALSE)
make_pca_plot(joint_vst_post, merged$batch[colnames(joint_vst_post)], "PCA Post-Batch Correction", outputs[["pca_post_pdf"]], center = TRUE, scale = FALSE)
make_umap(joint_vst_pre, "UMAP Pre-Batch Correction", outputs[["umap_pre_pdf"]], merged$batch[colnames(joint_vst_pre)], center = TRUE, scale = FALSE)
make_umap(joint_vst_post, "UMAP Post-Batch Correction", outputs[["umap_post_pdf"]], merged$batch[colnames(joint_vst_post)], center = TRUE, scale = FALSE)

# QC outputs tracked by Snakemake.
plot_mean_sd(tumour_vst, params[["tumour_label"]], outputs[["tumour_qc_pdf"]])
plot_mean_sd(dsmz_vst, params[["dsmz_label"]], outputs[["dsmz_qc_pdf"]])

# Additional dispersion diagnostics.
plot_dispersion(merged$tumour_counts, tumour_coldata, paste(params[["tumour_label"]], "Dispersion"), file.path(dirname(outputs[["tumour_qc_pdf"]]), "qc_tumour_dispersion.pdf"), max_genes = 10000L)
plot_dispersion(merged$dsmz_counts, dsmz_coldata, paste(params[["dsmz_label"]], "Dispersion"), file.path(dirname(outputs[["dsmz_qc_pdf"]]), "qc_dsmz_dispersion.pdf"), max_genes = 10000L)

cat("[SUCCESS] RBL preprocess workflow completed\n")
