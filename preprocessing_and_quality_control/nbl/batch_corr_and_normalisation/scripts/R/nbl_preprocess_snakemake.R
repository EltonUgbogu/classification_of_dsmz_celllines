root_dir <- normalizePath(file.path(snakemake@scriptdir, "..", ".."), mustWork = TRUE)
source(file.path(root_dir, "helpers.R"))
source(file.path(root_dir, "dsmz_base_functions.R"))
source(file.path(root_dir, "nbl_io.R"))
source(file.path(root_dir, "nbl_tumour_dsmz_base_functions.R"))
source(file.path(root_dir, "normalisation.R"))
source(file.path(root_dir, "visualisation.R"))
source(file.path(root_dir, "tcga_base_functions.R"))
## Palette, typography, gridless panel frame and atomic PDF device shared with
## BRCA and RBL. visualisation.R draws through these, so it must be loaded
## before any figure is produced.
source(snakemake@input[["shared_figure_module"]])

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

read_manifest_ids <- function(path) {
  manifest <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  if (!"sample_id" %in% colnames(manifest)) {
    stop("[ERROR] Sample manifest must contain sample_id: ", path)
  }
  ids <- as.character(manifest$sample_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("[ERROR] Sample manifest contains invalid sample IDs: ", path)
  }
  ids
}

cat("[INFO] Starting NBL purity-filtered preprocess workflow\n")
cat("[INFO] Authoritative order: retained raw tumour counts + raw DSMZ counts -> joint ComBat-seq -> joint DESeq2 VST\n")
cat(sprintf("[INFO] Tumour input: %s\n", inputs[["tumour_rds"]]))
tumour_meta_input <- inputs[["tumour_meta_csv"]]
tumour_meta_path <- if (length(tumour_meta_input) > 0) {
  as.character(tumour_meta_input[[1]])
} else {
  ""
}
cat(sprintf(
  "[INFO] Tumour metadata: %s\n",
  if (nzchar(tumour_meta_path)) tumour_meta_path else "<not configured>"
))
cat(sprintf("[INFO] DSMZ input: %s\n", inputs[["dsmz_rds"]]))
cat(sprintf("[INFO] DSMZ metadata: %s\n", inputs[["dsmz_meta_csv"]]))

# Load and validate the raw inputs selected by the upstream tumour-only purity
# rule. The retained manifest is the explicit sample contract for ComBat-seq.
tumour_counts <- load_tumour_nbl_data(
  inputs[["tumour_rds"]],
  meta_path = NULL,
  expected_total = params[["expected_tumour_samples"]],
  expected_cohorts = NULL
)
retained_ids <- read_manifest_ids(inputs[["retained_samples"]])
excluded_ids <- read_manifest_ids(inputs[["excluded_samples"]])
if (!identical(colnames(tumour_counts), retained_ids)) {
  stop("[ERROR] Retained raw NBL tumour columns do not exactly equal the retained manifest")
}
if (length(intersect(colnames(tumour_counts), excluded_ids))) {
  stop("[ERROR] An excluded NBL tumour occurs in the ComBat-seq tumour input")
}
tumour_metadata <- if (grepl("\\.tsv$", tumour_meta_path, ignore.case = TRUE)) {
  utils::read.delim(tumour_meta_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  utils::read.csv(tumour_meta_path, stringsAsFactors = FALSE, check.names = FALSE)
}
tumour_sample_col <- intersect(
  c("sample", "sample_id", "aliquot_id"),
  colnames(tumour_metadata)
)[1]
tumour_cohort_col <- intersect(
  c("cohort", "source", "project"),
  colnames(tumour_metadata)
)[1]
if (is.na(tumour_sample_col) || is.na(tumour_cohort_col)) {
  stop("[ERROR] NBL tumour metadata requires sample and cohort columns")
}
tumour_metadata <- tumour_metadata[
  match(colnames(tumour_counts), tumour_metadata[[tumour_sample_col]]),
  ,
  drop = FALSE
]
if (anyNA(tumour_metadata[[tumour_sample_col]]) ||
    !identical(as.character(tumour_metadata[[tumour_sample_col]]), colnames(tumour_counts))) {
  stop("[ERROR] NBL tumour metadata could not be aligned exactly to raw count columns")
}
dsmz_data <- load_dsmz_data(inputs[["dsmz_rds"]], inputs[["dsmz_meta_csv"]])
dsmz_counts_all <- build_dsmz_matrix(dsmz_data$raw)
dsmz_filtered <- load_dsmz_nbl_data(
  counts = dsmz_counts_all,
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
  specimen_type = c(
    rep("tumour", ncol(merged$tumour_counts)),
    rep("cell_line", ncol(merged$dsmz_counts))
  ),
  cohort = c(
    as.character(tumour_metadata[[tumour_cohort_col]]),
    rep("DSMZ", ncol(merged$dsmz_counts))
  ),
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
coldata_out <- data.frame(
  sample_id = rownames(coldata),
  coldata,
  row.names = NULL,
  check.names = FALSE
)
utils::write.table(
  coldata_out,
  file = outputs[["coldata_tsv"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Batch correction.
batch_corrected_counts <- batch_adjust(merged$Xc_raw, merged$batch)
saveRDS(batch_corrected_counts, outputs[["batch_corrected_counts_rds"]])

# VST is fitted jointly before and jointly after batch correction so the
# diagnostic matrices share one sample order and one fitted transformation.
joint_vst_pre <- vst_normalize(merged$Xc_raw, coldata)
joint_vst_post <- vst_normalize(batch_corrected_counts, coldata)

tumour_vst <- joint_vst_pre[, colnames(merged$tumour_counts), drop = FALSE]
dsmz_vst <- joint_vst_pre[, colnames(merged$dsmz_counts), drop = FALSE]
tumour_vst_post <- joint_vst_post[, colnames(merged$tumour_counts), drop = FALSE]
dsmz_vst_post <- joint_vst_post[, colnames(merged$dsmz_counts), drop = FALSE]

saveRDS(tumour_vst, outputs[["tumour_vst_rds"]])
saveRDS(dsmz_vst, outputs[["dsmz_vst_rds"]])
saveRDS(joint_vst_pre, outputs[["joint_vst_pre_bc_rds"]])
saveRDS(joint_vst_post, outputs[["joint_vst_post_bc_rds"]])
saveRDS(tumour_vst_post, outputs[["tumour_vst_post_bc_rds"]])
saveRDS(dsmz_vst_post, outputs[["dsmz_vst_post_bc_rds"]])

## ---------------------------------------------------------------------------
## Batch-effect quantification now lives in the dedicated
## `nbl_batch_effect_quantification` rule, which reads the joint VST matrices
## written above. Keeping it separate means the derived statistics, feature
## manifests and paired PCA figures can be deleted and rebuilt on their own,
## without regenerating a count or VST matrix.
## ---------------------------------------------------------------------------

# Plots. Titles and subtitles follow the same pattern as BRCA and RBL so the
# three cohorts read identically. The standalone qualitative figures reuse the
# same ordered pre-correction top3000 manifest before and after correction.
plot_labels_pre <- merged$batch[colnames(joint_vst_pre)]
plot_labels_post <- merged$batch[colnames(joint_vst_post)]
embedding_subtitle <- sprintf(
  "%s patient tumours and %s DSMZ cell-line profiles",
  figure_count(sum(plot_labels_pre == params[["tumour_label"]])),
  figure_count(sum(plot_labels_pre == params[["dsmz_label"]]))
)
qc_feature_spaces <- build_feature_spaces(joint_vst_pre, joint_vst_post, params[["qc_top_genes"]])
qc_ids <- qc_feature_spaces$spaces$top3000$ids
qc_stage_pca <- list(
  before = compute_stage_pca(joint_vst_pre, qc_ids),
  after = compute_stage_pca(joint_vst_post, qc_ids)
)
embedding_footer <- paste(
  qc_feature_spaces$spaces$top3000$label,
  " | PCA center=TRUE scale.=FALSE",
  sprintf(" | UMAP seed=%d metric=cosine n_neighbors=%d min_dist=%.1f",
          params[["seed"]], min(20L, ncol(joint_vst_pre) - 1L), 0.3)
)
umap_pre_embedding <- compute_feature_space_umap(
  joint_vst_pre, qc_ids, params[["seed"]], snakemake@threads,
  context = "pca_umap_qualitative_pre"
)
umap_post_embedding <- compute_feature_space_umap(
  joint_vst_post, qc_ids, params[["seed"]], snakemake@threads,
  context = "pca_umap_qualitative_post"
)

plot_embedding_figure(
  list(list(embedding = qc_stage_pca$before)),
  plot_labels_pre,
  outputs[["pca_pre_pdf"]],
  overall_title = "NBL profiles before ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "numeric"
)
plot_embedding_figure(
  list(list(embedding = qc_stage_pca$after)),
  plot_labels_post,
  outputs[["pca_post_pdf"]],
  overall_title = "NBL profiles after ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "numeric"
)
plot_embedding_figure(
  list(list(embedding = umap_pre_embedding)),
  plot_labels_pre,
  outputs[["umap_pre_pdf"]],
  overall_title = "NBL profiles before ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "bare"
)
plot_embedding_figure(
  list(list(embedding = umap_post_embedding)),
  plot_labels_post,
  outputs[["umap_post_pdf"]],
  overall_title = "NBL profiles after ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "bare"
)

# QC outputs tracked by Snakemake.
plot_mean_sd(tumour_vst_post, "NBL tumours: post-correction", outputs[["tumour_qc_pdf"]])
plot_mean_sd(dsmz_vst_post, "DSMZ NBL cell lines: post-correction", outputs[["dsmz_qc_pdf"]])

# Additional dispersion diagnostics.
plot_dispersion(
  merged$tumour_counts,
  coldata[colnames(merged$tumour_counts), , drop = FALSE],
  "NBL tumours: raw-count dispersion after purity filtering",
  outputs[["tumour_dispersion_pdf"]],
  max_genes = 10000L
)
plot_dispersion(
  merged$dsmz_counts,
  coldata[colnames(merged$dsmz_counts), , drop = FALSE],
  "DSMZ NBL cell lines: raw-count dispersion",
  outputs[["dsmz_dispersion_pdf"]],
  max_genes = 10000L
)

if ("provenance_tsv" %in% names(outputs)) {
  sha256_file <- function(path) {
    strsplit(system2("sha256sum", path, stdout = TRUE), "[[:space:]]+")[[1]][1]
  }
  provenance <- data.frame(
    field = c(
      "generated_at",
      "cohort",
      "processing_stage",
      "tumour_raw_count_input",
      "tumour_raw_count_sha256",
      "retained_sample_manifest",
      "retained_sample_manifest_sha256",
      "excluded_sample_manifest",
      "excluded_sample_manifest_sha256",
      "purity_validation",
      "purity_validation_sha256",
      "tumour_metadata_input",
      "tumour_metadata_sha256",
      "dsmz_raw_count_input",
      "dsmz_raw_count_sha256",
      "batch_correction_method",
      "normalisation_method",
      "pca_center_scale_policy",
      "paired_feature_space_qc",
      "umap_seed",
      "umap_metric",
      "umap_n_neighbors",
      "umap_min_dist",
      "tumour_samples",
      "dsmz_samples",
      "shared_genes",
      "joint_vst_post_bc_sha256"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      "NBL",
      "raw tumour counts -> tumour-only purity -> retained raw tumour counts + raw DSMZ -> ComBat-seq -> joint DESeq2 VST",
      inputs[["tumour_rds"]],
      sha256_file(inputs[["tumour_rds"]]),
      inputs[["retained_samples"]],
      sha256_file(inputs[["retained_samples"]]),
      inputs[["excluded_samples"]],
      sha256_file(inputs[["excluded_samples"]]),
      inputs[["purity_validation"]],
      sha256_file(inputs[["purity_validation"]]),
      tumour_meta_path,
      sha256_file(tumour_meta_path),
      inputs[["dsmz_rds"]],
      sha256_file(inputs[["dsmz_rds"]]),
      "sva::ComBat_seq on raw joint counts",
      "DESeq2::vst",
      "prcomp(t(vst_matrix[feature_ids, , drop = FALSE]), center = TRUE, scale. = FALSE)",
      "top3000 selected once from the pre-ComBat-seq joint VST matrix and reused unchanged before and after correction for standalone PCA/UMAP QC",
      params[["seed"]],
      "cosine",
      min(20L, ncol(joint_vst_pre) - 1L),
      0.3,
      ncol(merged$tumour_counts),
      ncol(merged$dsmz_counts),
      nrow(merged$Xc_raw),
      sha256_file(outputs[["joint_vst_post_bc_rds"]])
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    provenance,
    file = outputs[["provenance_tsv"]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

cat("[SUCCESS] NBL preprocess workflow completed\n")
