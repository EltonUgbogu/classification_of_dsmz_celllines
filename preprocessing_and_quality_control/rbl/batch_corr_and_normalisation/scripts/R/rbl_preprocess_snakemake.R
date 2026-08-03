root_dir <- normalizePath(file.path(snakemake@scriptdir, "..", ".."), mustWork = TRUE)
source(file.path(root_dir, "helpers.R"))
source(file.path(root_dir, "dsmz_base_functions.R"))
source(file.path(root_dir, "rbl_io.R"))
source(file.path(root_dir, "rbl_tumour_dsmz_base_functions.R"))
source(file.path(root_dir, "normalisation.R"))
source(file.path(root_dir, "visualisation.R"))
source(file.path(root_dir, "tumour_base_functions.R"))
## Palette, typography, gridless panel frame and atomic PDF device shared with
## BRCA and NBL. visualisation.R draws through these, so it must be loaded
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

validate_retained_counts <- function(counts) {
  if (!is.numeric(counts) || is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("[ERROR] Retained RBL tumour object must be a named numeric matrix")
  }
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("[ERROR] Retained RBL tumour counts contain duplicated identifiers")
  }
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      any(abs(counts - round(counts)) > 1e-8)) {
    stop("[ERROR] Retained RBL tumour counts are not valid non-negative integer counts")
  }
  invisible(TRUE)
}

cat("[INFO] Starting RBL purity-filtered preprocess workflow\n")
cat("[INFO] Authoritative order: retained raw tumour counts + raw DSMZ counts -> joint ComBat-seq -> joint DESeq2 VST\n")
cat(sprintf("[INFO] Tumour input: %s\n", inputs[["tumour_rds"]]))
cat(sprintf("[INFO] Tumour metadata: %s\n", inputs[["tumour_meta_csv"]]))
cat(sprintf("[INFO] DSMZ input: %s\n", inputs[["dsmz_count_matrix"]]))
cat(sprintf("[INFO] DSMZ metadata: %s\n", inputs[["dsmz_sample_metadata"]]))

# Load the retained raw tumour counts produced by the upstream tumour-only
# purity rule and enforce its manifest before any joint processing.
tumour_counts <- load_tumour_rbl_data(inputs[["tumour_rds"]])
validate_retained_counts(tumour_counts)
retained_ids <- read_manifest_ids(inputs[["retained_samples"]])
excluded_ids <- read_manifest_ids(inputs[["excluded_samples"]])
if (!identical(colnames(tumour_counts), retained_ids)) {
  stop("[ERROR] Retained raw RBL tumour columns do not exactly equal the retained manifest")
}
if (length(intersect(colnames(tumour_counts), excluded_ids))) {
  stop("[ERROR] An excluded RBL tumour occurs in the ComBat-seq tumour input")
}
tumour_metadata <- utils::read.csv(
  inputs[["tumour_meta_csv"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
tumour_sample_col <- intersect(c("sample", "sample_id", "run"), colnames(tumour_metadata))[1]
tumour_cohort_col <- intersect(c("cohort", "source", "project"), colnames(tumour_metadata))[1]
if (is.na(tumour_sample_col) || is.na(tumour_cohort_col)) {
  stop("[ERROR] RBL tumour metadata requires sample and cohort columns")
}
if (as.integer(params[["expected_tumour_samples"]]) > 0L &&
    ncol(tumour_counts) != as.integer(params[["expected_tumour_samples"]])) {
  stop(sprintf(
    "[ERROR] RBL raw tumour matrix has %d samples; expected %d",
    ncol(tumour_counts),
    as.integer(params[["expected_tumour_samples"]])
  ))
}
tumour_metadata <- tumour_metadata[
  match(colnames(tumour_counts), tumour_metadata[[tumour_sample_col]]),
  ,
  drop = FALSE
]
if (anyNA(tumour_metadata[[tumour_sample_col]]) ||
    !identical(as.character(tumour_metadata[[tumour_sample_col]]), colnames(tumour_counts))) {
  stop("[ERROR] RBL tumour metadata could not be aligned exactly to raw count columns")
}
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

# VST is fitted jointly before and jointly after batch correction.
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
## `rbl_batch_effect_quantification` rule, which reads the joint VST matrices
## written above. Keeping it separate means the derived statistics, feature
## manifests and paired PCA figures can be deleted and rebuilt on their own,
## without regenerating a count or VST matrix.
## ---------------------------------------------------------------------------

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
if (as.integer(params[["expected_collapsed_dsmz"]]) > 0L &&
    ncol(dsmz_vst_post_collapsed) != as.integer(params[["expected_collapsed_dsmz"]])) {
  stop(sprintf(
    "[ERROR] Collapsed DSMZ object has %d profiles; expected %d",
    ncol(dsmz_vst_post_collapsed),
    as.integer(params[["expected_collapsed_dsmz"]])
  ))
}

# Plots. Titles and subtitles follow the same pattern as BRCA and NBL so the
# three cohorts read identically; the underlying computations are unchanged.
plot_labels_pre <- merged$batch[colnames(joint_vst_pre)]
plot_labels_post <- merged$batch[colnames(joint_vst_post)]
embedding_subtitle <- sprintf(
  "%s patient tumours and %s DSMZ cell-line profiles",
  figure_count(sum(plot_labels_pre == params[["tumour_label"]])),
  figure_count(sum(plot_labels_pre == params[["dsmz_label"]]))
)

make_pca_plot(joint_vst_pre, plot_labels_pre, "RBL profiles before ComBat-seq", outputs[["pca_pre_pdf"]], center = TRUE, scale = FALSE, subtitle = embedding_subtitle)
make_pca_plot(joint_vst_post, plot_labels_post, "RBL profiles after ComBat-seq", outputs[["pca_post_pdf"]], center = TRUE, scale = FALSE, subtitle = embedding_subtitle)
make_umap(joint_vst_pre, "RBL profiles before ComBat-seq", outputs[["umap_pre_pdf"]], plot_labels_pre, center = TRUE, scale = FALSE, subtitle = embedding_subtitle)
make_umap(joint_vst_post, "RBL profiles after ComBat-seq", outputs[["umap_post_pdf"]], plot_labels_post, center = TRUE, scale = FALSE, subtitle = embedding_subtitle)

# QC outputs tracked by Snakemake.
plot_mean_sd(tumour_vst_post, "RBL tumours: post-correction", outputs[["tumour_qc_pdf"]])
plot_mean_sd(dsmz_vst_post, "DSMZ RBL cell lines: post-correction", outputs[["dsmz_qc_pdf"]])

# Additional dispersion diagnostics.
plot_dispersion(
  merged$tumour_counts,
  coldata[colnames(merged$tumour_counts), , drop = FALSE],
  "RBL tumours: raw-count dispersion after purity filtering",
  outputs[["tumour_dispersion_pdf"]],
  max_genes = 10000L
)
plot_dispersion(
  merged$dsmz_counts,
  coldata[colnames(merged$dsmz_counts), , drop = FALSE],
  "DSMZ RBL cell lines: raw-count dispersion",
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
      "tumour_samples",
      "dsmz_samples",
      "shared_genes",
      "joint_vst_post_bc_sha256"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      "RBL",
      "raw tumour counts -> tumour-only purity -> retained raw tumour counts + raw DSMZ -> ComBat-seq -> joint DESeq2 VST",
      inputs[["tumour_rds"]],
      sha256_file(inputs[["tumour_rds"]]),
      inputs[["retained_samples"]],
      sha256_file(inputs[["retained_samples"]]),
      inputs[["excluded_samples"]],
      sha256_file(inputs[["excluded_samples"]]),
      inputs[["purity_validation"]],
      sha256_file(inputs[["purity_validation"]]),
      inputs[["tumour_meta_csv"]],
      sha256_file(inputs[["tumour_meta_csv"]]),
      inputs[["dsmz_count_matrix"]],
      sha256_file(inputs[["dsmz_count_matrix"]]),
      "sva::ComBat_seq on raw joint counts",
      "DESeq2::vst",
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

cat("[SUCCESS] RBL preprocess workflow completed\n")
