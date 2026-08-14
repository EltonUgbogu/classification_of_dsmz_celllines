options(stringsAsFactors = FALSE)

for (path in unname(snakemake@output)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

log_path <- snakemake@output[["workflow_log"]]
log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
on.exit({
  sink()
  close(log_con)
}, add = TRUE)

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

save_rds_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = TRUE)
  if (!file.rename(temporary, path)) {
    if (!file.copy(temporary, path, overwrite = TRUE)) {
      stop("Failed to promote temporary RDS to output: ", path)
    }
    unlink(temporary)
  }
  invisible(path)
}

read_count_matrix <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("%s not found: %s", label, path))
  obj <- readRDS(path)
  if (methods::is(obj, "SummarizedExperiment")) {
    obj <- SummarizedExperiment::assay(obj)
  } else if (is.list(obj) && !is.null(obj$counts)) {
    obj <- obj$counts
  }
  if (!is.matrix(obj) && !is.data.frame(obj)) {
    stop(sprintf("%s must be matrix/data.frame-like; class=%s", label, paste(class(obj), collapse = ",")))
  }
  mat <- as.matrix(obj)
  if (!is.numeric(mat) || is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(sprintf("%s must be a named numeric matrix", label))
  }
  mat
}

validate_counts <- function(counts, label, chunk_rows = 1000L) {
  if (anyDuplicated(rownames(counts))) stop(label, " has duplicated raw row names")
  if (anyDuplicated(colnames(counts))) stop(label, " has duplicated column names")
  for (start in seq.int(1L, nrow(counts), by = chunk_rows)) {
    idx <- start:min(nrow(counts), start + chunk_rows - 1L)
    block <- counts[idx, , drop = FALSE]
    if (anyNA(block)) stop(label, " contains NA values")
    if (any(!is.finite(block))) stop(label, " contains non-finite values")
    if (any(block < 0)) stop(label, " contains negative values")
    if (any(abs(block - round(block)) > 1e-8)) stop(label, " contains non-integer-like values")
  }
  invisible(TRUE)
}

harmonise_ensembl_counts <- function(counts, label) {
  original_ids <- rownames(counts)
  harmonised_ids <- sub("\\..*$", "", original_ids)
  duplicate_occurrences <- sum(duplicated(harmonised_ids))
  if (duplicate_occurrences) {
    cat(sprintf(
      "[INFO] %s: aggregating %d duplicate occurrence(s) after Ensembl-version harmonisation\n",
      label,
      duplicate_occurrences
    ))
    counts <- rowsum(counts, group = harmonised_ids, reorder = TRUE)
  } else {
    rownames(counts) <- harmonised_ids
  }
  storage.mode(counts) <- "double"
  list(counts = counts, duplicate_occurrences = duplicate_occurrences)
}

sha256_file <- function(path) {
  digest::digest(object = path, algo = "sha256", file = TRUE, serialize = FALSE)
}

vst_normalise <- function(counts, coldata, label) {
  cat(sprintf("[INFO] DESeq2 VST: %s (%d genes x %d samples)\n", label, nrow(counts), ncol(counts)))
  integer_counts <- round(as.matrix(counts))
  # Keep the rounded matrix in double storage. ComBat-seq can produce valid
  # integer-valued counts above R's 32-bit integer limit; coercing those values
  # to integer introduces NA values before DESeq2 sees the matrix.
  if (any(integer_counts > .Machine$integer.max)) {
    cat(sprintf(
      "[INFO] %s contains rounded counts above R's 32-bit integer limit; preserving double storage\n",
      label
    ))
  }
  if (!identical(colnames(integer_counts), rownames(coldata))) {
    stop("coldata row order does not match count columns for ", label)
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = integer_counts,
    colData = coldata,
    design = ~1
  )
  transformed <- tryCatch(
    DESeq2::vst(dds, blind = TRUE, fitType = "parametric"),
    error = function(condition) {
      cat(sprintf(
        "[WARN] Parametric DESeq2 VST failed for %s: %s\n",
        label,
        conditionMessage(condition)
      ))
      cat(sprintf("[INFO] Retrying DESeq2 VST for %s with fitType=local\n", label))
      DESeq2::vst(dds, blind = TRUE, fitType = "local")
    }
  )
  out <- SummarizedExperiment::assay(transformed)
  if (!identical(dim(out), dim(counts))) stop("VST changed matrix dimensions for ", label)
  if (anyNA(out) || any(!is.finite(out))) stop("VST returned invalid values for ", label)
  out
}

combat_seq_adjust <- function(counts, batch, seed) {
  cat(sprintf(
    "[INFO] sva::ComBat_seq on %d genes x %d samples; batch sizes=%s\n",
    nrow(counts),
    ncol(counts),
    paste(names(table(batch)), as.integer(table(batch)), sep = ":", collapse = ",")
  ))
  integer_counts <- round(as.matrix(counts))
  set.seed(seed)
  adjusted <- sva::ComBat_seq(
    counts = integer_counts,
    batch = as.factor(batch),
    group = NULL,
    covar_mod = NULL,
    full_mod = TRUE,
    shrink = FALSE,
    shrink.disp = FALSE
  )
  adjusted <- as.matrix(adjusted)
  rownames(adjusted) <- rownames(counts)
  colnames(adjusted) <- colnames(counts)
  if (!identical(dim(adjusted), dim(counts))) stop("ComBat-seq changed matrix dimensions")
  if (anyNA(adjusted) || any(!is.finite(adjusted)) || any(adjusted < 0)) {
    stop("ComBat-seq returned NA, non-finite, or negative values")
  }
  adjusted
}

## Shared figure and batch-effect library (single implementation for all cohorts).
source(snakemake@input[["shared_figure_module"]])


read_manifest_ids <- function(path) {
  manifest <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  if (!"sample_id" %in% colnames(manifest)) {
    stop("Retained-sample manifest must contain sample_id: ", path)
  }
  ids <- as.character(manifest$sample_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Retained-sample manifest has invalid sample IDs")
  }
  ids
}

cat("[INFO] Starting BRCA purity-filtered raw-count batch-correction and normalisation workflow\n")
cat("[INFO] Authoritative order: retained raw tumour counts + raw DSMZ counts -> joint ComBat-seq -> joint DESeq2 VST\n")
cat("[WARN] Tumour-versus-cell-line source is confounded with the two batch labels; interpret correction accordingly.\n")

expected_tumours <- as.integer(snakemake@params[["expected_tumours"]])
expected_dsmz <- as.integer(snakemake@params[["expected_dsmz"]])
expected_shared <- as.integer(snakemake@params[["expected_shared_genes"]])
min_shared <- as.integer(snakemake@params[["min_shared"]])
tumour_label <- as.character(snakemake@params[["tumour_label"]])
dsmz_label <- as.character(snakemake@params[["dsmz_label"]])
dsmz_sample_col <- as.character(snakemake@params[["dsmz_sample_col"]])
dsmz_code_col <- as.character(snakemake@params[["dsmz_code_col"]])
dsmz_code <- as.character(snakemake@params[["dsmz_code"]])
dsmz_cell_line_col <- as.character(snakemake@params[["dsmz_cell_line_col"]])
seed <- as.integer(snakemake@params[["seed"]])
qc_top_genes <- as.integer(snakemake@params[["qc_top_genes"]])
permanova_permutations <- as.integer(snakemake@params[["permanova_permutations"]])
threads <- as.integer(snakemake@threads)

tumour_raw <- read_count_matrix(
  snakemake@input[["tumour_rds"]],
  "purity-retained raw TCGA-BRCA counts"
)
dsmz_raw <- read_count_matrix(snakemake@input[["dsmz_rds"]], "DSMZ BRCA counts")
validate_counts(tumour_raw, "purity-retained raw TCGA-BRCA counts")
validate_counts(dsmz_raw, "DSMZ BRCA counts")

retained_ids <- read_manifest_ids(snakemake@input[["retained_samples"]])
if (!identical(colnames(tumour_raw), retained_ids)) {
  stop("Purity-retained raw tumour columns do not exactly equal the retained manifest")
}
if (length(intersect(
  colnames(tumour_raw),
  read_manifest_ids(snakemake@input[["excluded_samples"]])
))) {
  stop("An excluded tumour appears in the retained raw tumour matrix")
}
if (expected_tumours > 0L && ncol(tumour_raw) != expected_tumours) {
  stop(sprintf(
    "Expected %d purity-retained TCGA primary tumours; observed %d",
    expected_tumours,
    ncol(tumour_raw)
  ))
}
cat(sprintf(
  "[INFO] Purity-retained raw TCGA primary tumours: %d\n",
  ncol(tumour_raw)
))
if (ncol(dsmz_raw) != expected_dsmz) {
  stop(sprintf("Expected %d DSMZ profiles; observed %d", expected_dsmz, ncol(dsmz_raw)))
}
if (length(intersect(colnames(tumour_raw), colnames(dsmz_raw)))) {
  stop("Tumour and DSMZ matrices have colliding sample IDs")
}

dsmz_meta <- utils::read.csv(
  snakemake@input[["dsmz_meta_csv"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_meta <- c(dsmz_sample_col, dsmz_code_col, dsmz_cell_line_col)
missing_meta <- setdiff(required_meta, colnames(dsmz_meta))
if (length(missing_meta)) {
  stop("DSMZ metadata missing columns: ", paste(missing_meta, collapse = ", "))
}
brca_meta <- dsmz_meta[dsmz_meta[[dsmz_code_col]] %in% dsmz_code, , drop = FALSE]
match_count <- vapply(
  colnames(dsmz_raw),
  function(id) sum(brca_meta[[dsmz_sample_col]] == id, na.rm = TRUE),
  integer(1)
)
if (any(match_count != 1L)) {
  stop(
    "Every DSMZ BRCA count column must match exactly one BRCA metadata row: ",
    paste(names(match_count)[match_count != 1L], collapse = ", ")
  )
}
brca_meta <- brca_meta[match(colnames(dsmz_raw), brca_meta[[dsmz_sample_col]]), , drop = FALSE]
if (!identical(brca_meta[[dsmz_sample_col]], colnames(dsmz_raw))) {
  stop("DSMZ metadata order could not be aligned exactly to count columns")
}
cell_line_names <- brca_meta[[dsmz_cell_line_col]]
if (anyNA(cell_line_names) || any(!nzchar(trimws(cell_line_names)))) {
  stop("DSMZ BRCA metadata has missing cell-line labels")
}
if (length(unique(cell_line_names)) != expected_dsmz) {
  stop("DSMZ BRCA input is expected to contain one profile per cell line")
}

tumour_harmonised <- harmonise_ensembl_counts(tumour_raw, "TCGA-BRCA")
dsmz_harmonised <- harmonise_ensembl_counts(dsmz_raw, "DSMZ BRCA")
shared_genes <- sort(intersect(
  rownames(tumour_harmonised$counts),
  rownames(dsmz_harmonised$counts)
))
if (length(shared_genes) < min_shared) {
  stop(sprintf("Only %d shared genes; configured minimum is %d", length(shared_genes), min_shared))
}
if (length(shared_genes) != expected_shared) {
  stop(sprintf("Expected %d shared genes; observed %d", expected_shared, length(shared_genes)))
}

tumour_counts <- tumour_harmonised$counts[shared_genes, , drop = FALSE]
dsmz_counts <- dsmz_harmonised$counts[shared_genes, , drop = FALSE]
merged_counts <- cbind(tumour_counts, dsmz_counts)
batch <- factor(
  c(rep(tumour_label, ncol(tumour_counts)), rep(dsmz_label, ncol(dsmz_counts))),
  levels = c(tumour_label, dsmz_label)
)
names(batch) <- colnames(merged_counts)

coldata <- data.frame(
  dataset = as.character(batch),
  specimen_type = c(
    rep("tumour", ncol(tumour_counts)),
    rep("cell_line", ncol(dsmz_counts))
  ),
  cohort = c(
    rep("TCGA-BRCA", ncol(tumour_counts)),
    rep("DSMZ", ncol(dsmz_counts))
  ),
  cell_line = c(rep(NA_character_, ncol(tumour_counts)), cell_line_names),
  row.names = colnames(merged_counts),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

save_rds_atomic(tumour_counts, snakemake@output[["tumour_counts_rds"]])
save_rds_atomic(dsmz_counts, snakemake@output[["dsmz_counts_rds"]])
save_rds_atomic(merged_counts, snakemake@output[["merged_counts_rds"]])
write_tsv(
  data.frame(sample_id = rownames(coldata), coldata, row.names = NULL, check.names = FALSE),
  snakemake@output[["coldata_tsv"]]
)

final_order <- colnames(merged_counts)
joint_vst_pre <- vst_normalise(
  merged_counts,
  coldata,
  "joint purity-filtered pre-ComBat-seq"
)
tumour_vst_pre <- joint_vst_pre[, colnames(tumour_counts), drop = FALSE]
dsmz_vst_pre <- joint_vst_pre[, colnames(dsmz_counts), drop = FALSE]

batch_corrected_counts <- combat_seq_adjust(merged_counts, batch, seed = seed)
cat(sprintf(
  "[INFO] ComBat-seq adjusted-count range: %.0f..%.0f; values above R integer limit: %d\n",
  min(batch_corrected_counts),
  max(batch_corrected_counts),
  sum(batch_corrected_counts > .Machine$integer.max)
))
joint_vst_post <- vst_normalise(
  batch_corrected_counts,
  coldata,
  "joint purity-filtered post-ComBat-seq"
)
tumour_vst_post <- joint_vst_post[, colnames(tumour_counts), drop = FALSE]
dsmz_vst_post <- joint_vst_post[, colnames(dsmz_counts), drop = FALSE]

save_rds_atomic(batch_corrected_counts, snakemake@output[["batch_corrected_counts_rds"]])
save_rds_atomic(tumour_vst_pre, snakemake@output[["tumour_vst_rds"]])
save_rds_atomic(dsmz_vst_pre, snakemake@output[["dsmz_vst_rds"]])
save_rds_atomic(joint_vst_pre, snakemake@output[["joint_vst_pre_bc_rds"]])
save_rds_atomic(joint_vst_post, snakemake@output[["joint_vst_post_bc_rds"]])
save_rds_atomic(tumour_vst_post, snakemake@output[["tumour_vst_post_bc_rds"]])
save_rds_atomic(dsmz_vst_post, snakemake@output[["dsmz_vst_post_bc_rds"]])

plot_labels <- coldata[final_order, "dataset"]
embedding_subtitle <- sprintf(
  "%s patient tumours and %s DSMZ cell-line profiles",
  figure_count(sum(plot_labels == tumour_label)),
  figure_count(sum(plot_labels == dsmz_label))
)

## Batch-effect quantification, its feature manifests and the paired
## before/after PCA figures now live in the dedicated
## `brca_batch_effect_quantification` rule, which consumes the joint VST
## matrices written above. Keeping them in a separate rule means the derived
## statistics and figures can be deleted and rebuilt on their own, without
## regenerating a single count or VST matrix.
##
## The per-stage PCA figures below still use the top-variable-gene feature
## space, rebuilt here by the same shared code path so the embedding is
## identical to the one the quantification rule tests.
qc_feature_spaces <- build_feature_spaces(joint_vst_pre, joint_vst_post, qc_top_genes)
qc_ids <- qc_feature_spaces$spaces$top3000$ids
qc_stage_pca <- list(
  before = compute_stage_pca(joint_vst_pre, qc_ids),
  after = compute_stage_pca(joint_vst_post, qc_ids)
)

## Per-stage PCA and UMAP figures reuse the same ordered pre-correction top3000
## manifest, so the standalone qualitative figures stay on the audited paired
## feature space rather than re-ranking genes independently per stage.
embedding_footer <- paste(
  qc_feature_spaces$spaces$top3000$label,
  " | PCA center=TRUE scale.=FALSE",
  sprintf(" | UMAP seed=%d metric=cosine n_neighbors=%d min_dist=%.1f",
          seed, min(20L, ncol(joint_vst_pre) - 1L), 0.3)
)
umap_pre_embedding <- compute_feature_space_umap(
  joint_vst_pre, qc_ids, seed, threads,
  context = "pca_umap_qualitative_pre"
)
umap_post_embedding <- compute_feature_space_umap(
  joint_vst_post, qc_ids, seed, threads,
  context = "pca_umap_qualitative_post"
)

plot_embedding_figure(
  list(list(embedding = qc_stage_pca$before)),
  plot_labels,
  snakemake@output[["pca_pre_pdf"]],
  overall_title = "BRCA profiles before ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "numeric"
)
plot_embedding_figure(
  list(list(embedding = qc_stage_pca$after)),
  plot_labels,
  snakemake@output[["pca_post_pdf"]],
  overall_title = "BRCA profiles after ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "numeric"
)
plot_embedding_figure(
  list(list(embedding = umap_pre_embedding)),
  plot_labels,
  snakemake@output[["umap_pre_pdf"]],
  overall_title = "BRCA profiles before ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "bare"
)
plot_embedding_figure(
  list(list(embedding = umap_post_embedding)),
  plot_labels,
  snakemake@output[["umap_post_pdf"]],
  overall_title = "BRCA profiles after ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "bare"
)
plot_embedding_figure(
  list(
    list(embedding = umap_pre_embedding, main = "A  Before correction"),
    list(embedding = umap_post_embedding, main = "B  After correction")
  ),
  plot_labels,
  snakemake@output[["umap_before_after_pdf"]],
  overall_title = "BRCA profiles before and after ComBat-seq",
  overall_subtitle = embedding_subtitle,
  footer = embedding_footer,
  scale_style = "bare"
)
plot_vst_mean_sd(
  tumour_vst_post,
  "BRCA tumours: post-correction VST mean versus standard deviation",
  snakemake@output[["tumour_qc_pdf"]],
  seed
)
plot_vst_mean_sd(
  dsmz_vst_post,
  "DSMZ BRCA cell lines: post-correction VST mean versus standard deviation",
  snakemake@output[["dsmz_qc_pdf"]],
  seed
)
plot_count_dispersion(
  tumour_counts,
  "BRCA tumours: raw-count dispersion after purity filtering",
  snakemake@output[["tumour_dispersion_pdf"]],
  seed
)
plot_count_dispersion(
  dsmz_counts,
  "DSMZ BRCA cell lines: raw-count dispersion",
  snakemake@output[["dsmz_dispersion_pdf"]],
  seed
)

package_names <- c("DESeq2", "sva", "SummarizedExperiment", "matrixStats", "uwot", "digest")
package_versions <- vapply(
  package_names,
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)
key_output_paths <- c(
  merged_counts_raw = snakemake@output[["merged_counts_rds"]],
  merged_counts_batch_corrected = snakemake@output[["batch_corrected_counts_rds"]],
  joint_vst_post_bc = snakemake@output[["joint_vst_post_bc_rds"]]
)
provenance <- data.frame(
  field = c(
    "generated_at",
    "workflow",
    "reference_implementation",
    "raw_tumour_count_input",
    "retained_sample_manifest",
    "retained_sample_manifest_sha256",
    "purity_validation",
    "purity_validation_sha256",
    "dsmz_count_input",
    "dsmz_metadata_input",
    "tumour_input_sha256",
    "dsmz_input_sha256",
    "dsmz_metadata_sha256",
    "tumour_dimensions",
    "dsmz_dimensions",
    "shared_genes",
    "joint_dimensions",
    "joint_column_order",
    "batch_method",
    "batch_design_warning",
    "vst_method",
    "random_seed",
    "qc_top_variable_genes",
    "pca_center_scale_policy",
    "paired_feature_space_qc",
    "umap_seed",
    "umap_metric",
    "umap_n_neighbors",
    "umap_min_dist",
    "R_version",
    paste0("package_", package_names),
    paste0("output_sha256_", names(key_output_paths))
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    "BRCA raw-count batch correction and normalisation",
    "RBL structural logic adapted to BRCA inputs and outputs",
    snakemake@input[["tumour_rds"]],
    snakemake@input[["retained_samples"]],
    sha256_file(snakemake@input[["retained_samples"]]),
    snakemake@input[["purity_validation"]],
    sha256_file(snakemake@input[["purity_validation"]]),
    snakemake@input[["dsmz_rds"]],
    snakemake@input[["dsmz_meta_csv"]],
    sha256_file(snakemake@input[["tumour_rds"]]),
    sha256_file(snakemake@input[["dsmz_rds"]]),
    sha256_file(snakemake@input[["dsmz_meta_csv"]]),
    sprintf("%dx%d", nrow(tumour_counts), ncol(tumour_counts)),
    sprintf("%dx%d", nrow(dsmz_counts), ncol(dsmz_counts)),
    length(shared_genes),
    sprintf("%dx%d", nrow(joint_vst_post), ncol(joint_vst_post)),
    "Purity-retained tumours followed by DSMZ profiles; identical across raw, corrected, VST, and coldata",
    "sva::ComBat_seq(group=NULL,full_mod=TRUE,shrink=FALSE,shrink.disp=FALSE)",
    "Tumour versus cell-line source is perfectly confounded with configured batch",
    "DESeq2::vst(blind=TRUE) fitted jointly before and jointly after ComBat-seq",
    seed,
    qc_top_genes,
    "prcomp(t(vst_matrix[feature_ids, , drop = FALSE]), center = TRUE, scale. = FALSE)",
    "top3000 selected once from the pre-ComBat-seq joint VST matrix and reused unchanged before and after correction for standalone PCA/UMAP QC",
    seed,
    "cosine",
    min(20L, ncol(joint_vst_pre) - 1L),
    0.3,
    R.version.string,
    unname(package_versions),
    vapply(key_output_paths, sha256_file, character(1))
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(provenance, snakemake@output[["provenance_tsv"]])

cat(sprintf("[INFO] Harmonised shared genes: %d\n", length(shared_genes)))
cat(sprintf("[INFO] Purity-retained primary tumours entering ComBat-seq: %d\n", ncol(tumour_counts)))
cat(sprintf("[INFO] DSMZ BRCA profiles: %d\n", ncol(dsmz_counts)))
cat(sprintf("[INFO] Final joint VST: %d genes x %d samples\n", nrow(joint_vst_post), ncol(joint_vst_post)))
cat("[SUCCESS] BRCA batch-correction and normalisation completed\n")
