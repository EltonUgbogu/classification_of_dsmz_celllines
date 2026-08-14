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

## RBL sequences some cell lines more than once, so the DSMZ block holds more
## libraries than distinct cell lines. Averaging the replicate libraries of a
## cell line in VST space gives one profile per biological cell line; the same
## key and the same column order are used pre- and post-correction so the two
## embedding figures stay comparable.
collapse_dsmz_by_cellline <- function(mat) {
  collapsed <- vapply(
    collapsed_celllines,
    function(cell_line) rowMeans(mat[, cell_line_key == cell_line, drop = FALSE]),
    numeric(nrow(mat))
  )
  rownames(collapsed) <- rownames(mat)
  colnames(collapsed) <- collapsed_celllines
  collapsed
}

collapse_metadata_by_cellline <- function(metadata) {
  sample_groups <- split(colnames(dsmz_vst_post), cell_line_key)
  collapsed_rows <- lapply(
    collapsed_celllines,
    function(cell_line) {
      rows <- metadata[match(sample_groups[[cell_line]], metadata$sample_id), , drop = FALSE]
      template <- rows[1L, , drop = FALSE]
      template$sample_id <- cell_line
      template[[collapse_col]] <- cell_line
      if ("Cell_Line" %in% colnames(template)) template$Cell_Line <- cell_line
      template
    }
  )
  collapsed <- do.call(rbind, collapsed_rows)
  rownames(collapsed) <- NULL
  collapsed
}

collapse_counts_by_cellline <- function(counts) {
  collapsed <- vapply(
    collapsed_celllines,
    function(cell_line) rowMeans(counts[, cell_line_key == cell_line, drop = FALSE]),
    numeric(nrow(counts))
  )
  rownames(collapsed) <- rownames(counts)
  colnames(collapsed) <- collapsed_celllines
  collapsed
}

assert_expected_replicate_groups <- function() {
  observed_sizes <- stats::setNames(as.integer(table(cell_line_key)), names(table(cell_line_key)))
  expected_sizes <- c(RBL_15 = 2L, RBL_20 = 2L)
  if (!identical(observed_sizes[names(expected_sizes)], expected_sizes)) {
    stop(sprintf(
      "[ERROR] DSMZ replicate collapse did not match the expected biological groups: observed %s",
      paste(sprintf("%s=%d", names(observed_sizes), observed_sizes), collapse = ", ")
    ))
  }
  if (sum(observed_sizes == 1L) != 7L || length(observed_sizes) != 9L) {
    stop(sprintf(
      "[ERROR] DSMZ replicate collapse produced %d biological groups (%d singleton groups); expected 9 total groups with 7 singleton groups",
      length(observed_sizes), sum(observed_sizes == 1L)
    ))
  }
  invisible(observed_sizes)
}

dsmz_vst_post_collapsed <- collapse_dsmz_by_cellline(dsmz_vst_post)
dsmz_vst_pre_collapsed <- collapse_dsmz_by_cellline(dsmz_vst)
assert_expected_replicate_groups()
saveRDS(dsmz_vst_pre_collapsed, outputs[["dsmz_vst_pre_bc_collapsed_cellline_rds"]])
saveRDS(dsmz_vst_post_collapsed, outputs[["dsmz_vst_post_bc_collapsed_cellline_rds"]])
if (!identical(colnames(dsmz_vst_pre_collapsed), colnames(dsmz_vst_post_collapsed))) {
  stop("[ERROR] Pre- and post-ComBat-seq collapsed DSMZ VST objects differ in group ordering")
}
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
# three cohorts read identically. The standalone qualitative figures reuse the
# same ordered pre-correction top3000 manifest before and after correction.
#
# The embeddings are drawn over the replicate-collapsed DSMZ block, so one point
# is one cell line rather than one sequencing library. Every matrix written above
# is untouched: this collapse exists only for the four embedding figures, and the
# saved joint VST objects still carry all libraries.
embedding_pre <- cbind(tumour_vst, dsmz_vst_pre_collapsed)
embedding_post <- cbind(tumour_vst_post, dsmz_vst_post_collapsed)
plot_labels_pre <- c(
  rep(params[["tumour_label"]], ncol(tumour_vst)),
  rep(params[["dsmz_label"]], ncol(dsmz_vst_pre_collapsed))
)
plot_labels_post <- c(
  rep(params[["tumour_label"]], ncol(tumour_vst_post)),
  rep(params[["dsmz_label"]], ncol(dsmz_vst_post_collapsed))
)
embedding_subtitle <- sprintf(
  "%s patient tumours and %s DSMZ cell-line groups",
  figure_count(ncol(tumour_vst)),
  figure_count(ncol(dsmz_vst_pre_collapsed))
)
if (ncol(tumour_vst) != 68L) {
  stop(sprintf("[ERROR] RBL plotting path has %d tumour profiles; expected 68 purity-retained tumours", ncol(tumour_vst)))
}
if (ncol(dsmz_vst_pre_collapsed) != 9L || ncol(dsmz_vst_post_collapsed) != 9L) {
  stop(sprintf(
    "[ERROR] RBL plotting path has %d pre-BC and %d post-BC DSMZ biological groups; expected 9 at both stages",
    ncol(dsmz_vst_pre_collapsed), ncol(dsmz_vst_post_collapsed)
  ))
}
qc_feature_spaces <- build_feature_spaces(embedding_pre, embedding_post, params[["qc_top_genes"]])
qc_ids <- qc_feature_spaces$spaces$top3000$ids
qc_stage_pca <- list(
  before = compute_stage_pca(embedding_pre, qc_ids),
  after = compute_stage_pca(embedding_post, qc_ids)
)
embedding_footer <- paste(
  qc_feature_spaces$spaces$top3000$label,
  " | PCA center=TRUE scale.=FALSE",
  sprintf(" | UMAP seed=%d metric=cosine n_neighbors=%d min_dist=%.1f",
          params[["seed"]], min(20L, ncol(embedding_pre) - 1L), 0.3)
)
umap_pre_embedding <- compute_feature_space_umap(
  embedding_pre, qc_ids, params[["seed"]], snakemake@threads,
  context = "pca_umap_qualitative_pre"
)
umap_post_embedding <- compute_feature_space_umap(
  embedding_post, qc_ids, params[["seed"]], snakemake@threads,
  context = "pca_umap_qualitative_post"
)

plot_embedding_figure(
  list(list(embedding = qc_stage_pca$before)),
  plot_labels_pre,
  outputs[["pca_pre_pdf"]],
  overall_title = "RBL profiles before ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "numeric"
)
plot_embedding_figure(
  list(list(embedding = qc_stage_pca$after)),
  plot_labels_post,
  outputs[["pca_post_pdf"]],
  overall_title = "RBL profiles after ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "numeric"
)
plot_embedding_figure(
  list(list(embedding = umap_pre_embedding)),
  plot_labels_pre,
  outputs[["umap_pre_pdf"]],
  overall_title = "RBL profiles before ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "bare"
)
plot_embedding_figure(
  list(list(embedding = umap_post_embedding)),
  plot_labels_post,
  outputs[["umap_post_pdf"]],
  overall_title = "RBL profiles after ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "bare"
)

# QC outputs tracked by Snakemake.
plot_mean_sd(tumour_vst_post, "RBL tumours: post-correction", outputs[["tumour_qc_pdf"]])
plot_mean_sd(dsmz_vst_post_collapsed, "DSMZ RBL cell-line groups: post-correction", outputs[["dsmz_qc_pdf"]])

# Additional dispersion diagnostics.
plot_dispersion(
  merged$tumour_counts,
  coldata[colnames(merged$tumour_counts), , drop = FALSE],
  "RBL tumours: raw-count dispersion after purity filtering",
  outputs[["tumour_dispersion_pdf"]],
  max_genes = 10000L
)
dsmz_counts_collapsed <- collapse_counts_by_cellline(merged$dsmz_counts)
dsmz_metadata_collapsed <- collapse_metadata_by_cellline(dsmz_filtered$metadata)
dsmz_coldata_collapsed <- data.frame(
  dataset = rep(params[["dsmz_label"]], length(collapsed_celllines)),
  specimen_type = rep("cell_line", length(collapsed_celllines)),
  cohort = rep("DSMZ", length(collapsed_celllines)),
  cell_line = collapsed_celllines,
  subtype = rep(NA_character_, length(collapsed_celllines)),
  row.names = collapsed_celllines,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(
  identical(colnames(dsmz_counts_collapsed), rownames(dsmz_coldata_collapsed)),
  identical(colnames(dsmz_counts_collapsed), dsmz_metadata_collapsed$sample_id)
)
plot_dispersion(
  dsmz_counts_collapsed,
  dsmz_coldata_collapsed,
  "DSMZ RBL cell-line groups: raw-count dispersion",
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
      "plotting_population",
      "pca_center_scale_policy",
      "paired_feature_space_qc",
      "umap_seed",
      "umap_metric",
      "umap_n_neighbors",
      "umap_min_dist",
      "tumour_samples",
      "dsmz_samples",
      "collapsed_dsmz_groups",
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
      "68 purity-retained tumours + 9 DSMZ biological groups for standalone PCA/UMAP figures; 68 tumours + 11 DSMZ libraries preserved in the joint matrices",
      "prcomp(t(vst_matrix[feature_ids, , drop = FALSE]), center = TRUE, scale. = FALSE)",
      "top3000 selected once from the pre-ComBat-seq collapsed plotting VST matrix and reused unchanged before and after correction for standalone PCA/UMAP QC",
      params[["seed"]],
      "cosine",
      min(20L, ncol(embedding_pre) - 1L),
      0.3,
      ncol(merged$tumour_counts),
      ncol(merged$dsmz_counts),
      ncol(dsmz_vst_post_collapsed),
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
