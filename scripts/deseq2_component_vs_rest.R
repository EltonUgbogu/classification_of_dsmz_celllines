#!/usr/bin/env Rscript

# Component-vs-rest DESeq2 marker prioritisation.
# Unit of analysis: one resolved graph component versus all other prepared
# same-cancer cell-line samples. Wald p-values/BH-adjusted p-values provide
# inference; apeglm posterior-mode shrunken LFCs provide effect-size filtering,
# direction assignment, and deterministic marker ranking.

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(apeglm)
  library(data.table)
})

option_list <- list(
  make_option("--counts", type = "character"),
  make_option("--meta", type = "character"),
  make_option("--component", type = "character"),
  make_option("--component_col", type = "character", default = "component"),
  make_option("--sample_id_col", type = "character", default = "sample_id"),
  make_option("--outdir", type = "character"),
  make_option("--adjusted_p_value_threshold", type = "double"),
  make_option("--minimum_absolute_shrunken_log2fc", type = "double"),
  make_option("--maximum_markers_per_direction", type = "integer"),
  make_option("--minimum_base_mean", type = "double"),
  make_option("--minimum_total_gene_count", type = "integer"),
  make_option("--dispersion_fit_type", type = "character"),
  make_option("--lfc_shrinkage_method", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("counts", "meta", "component", "outdir")
for (name in required) {
  value <- opt[[name]]
  if (is.null(value) || is.na(value) || !nzchar(trimws(value))) {
    stop(sprintf("[Component DESeq2] --%s is required", name), call. = FALSE)
  }
}

is_scalar_finite_number <- function(value) {
  length(value) == 1L && !is.na(value) && is.finite(value)
}

if (!is_scalar_finite_number(opt$adjusted_p_value_threshold) ||
    opt$adjusted_p_value_threshold <= 0 || opt$adjusted_p_value_threshold > 1) {
  stop("--adjusted_p_value_threshold must be in (0, 1]", call. = FALSE)
}
if (!is_scalar_finite_number(opt$minimum_absolute_shrunken_log2fc) ||
    opt$minimum_absolute_shrunken_log2fc < 0) {
  stop("--minimum_absolute_shrunken_log2fc must be finite and non-negative", call. = FALSE)
}
if (!is_scalar_finite_number(opt$minimum_base_mean) || opt$minimum_base_mean < 0) {
  stop("--minimum_base_mean must be finite and non-negative", call. = FALSE)
}
if (!is_scalar_finite_number(opt$maximum_markers_per_direction) ||
    opt$maximum_markers_per_direction < 1 || opt$maximum_markers_per_direction != floor(opt$maximum_markers_per_direction)) {
  stop("--maximum_markers_per_direction must be a positive integer", call. = FALSE)
}
if (!is_scalar_finite_number(opt$minimum_total_gene_count) || opt$minimum_total_gene_count < 1) {
  stop("--minimum_total_gene_count must be a positive integer", call. = FALSE)
}
if (length(opt$dispersion_fit_type) != 1L || is.na(opt$dispersion_fit_type) ||
    !opt$dispersion_fit_type %in% c("parametric", "local", "mean")) {
  stop("--dispersion_fit_type must be one of: parametric, local, mean", call. = FALSE)
}
if (length(opt$lfc_shrinkage_method) != 1L || is.na(opt$lfc_shrinkage_method) ||
    opt$lfc_shrinkage_method != "apeglm") {
  stop("--lfc_shrinkage_method must be apeglm for the active component method", call. = FALSE)
}

dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "markers"), recursive = TRUE, showWarnings = FALSE)

validate_count_matrix <- function(matrix_values) {
  dim_names <- dimnames(matrix_values)
  suppressWarnings(storage.mode(matrix_values) <- "numeric")
  dimnames(matrix_values) <- dim_names
  if (any(is.na(matrix_values))) stop("Counts contain NA values", call. = FALSE)
  if (any(!is.finite(matrix_values))) stop("Counts contain non-finite values", call. = FALSE)
  if (any(matrix_values < 0)) stop("Counts contain negative values", call. = FALSE)
  if (any(abs(matrix_values - round(matrix_values)) > 1e-8)) {
    stop("DESeq2 input must contain raw integer counts", call. = FALSE)
  }
  matrix_values <- round(matrix_values)
  storage.mode(matrix_values) <- "integer"
  dimnames(matrix_values) <- dim_names
  matrix_values
}

counts_table <- fread(opt$counts, sep = "\t", data.table = FALSE, check.names = FALSE)
if (ncol(counts_table) < 2L || names(counts_table)[1L] != "gene_id") {
  stop("[Component DESeq2] Count table must have gene_id followed by sample columns", call. = FALSE)
}
gene_ids <- as.character(counts_table$gene_id)
if (any(is.na(gene_ids)) || any(!nzchar(trimws(gene_ids))) || anyDuplicated(gene_ids)) {
  stop("[Component DESeq2] Gene identifiers must be present and unique", call. = FALSE)
}
count_matrix <- as.matrix(counts_table[, -1L, drop = FALSE])
rownames(count_matrix) <- gene_ids
count_matrix <- validate_count_matrix(count_matrix)
if (anyDuplicated(colnames(count_matrix))) {
  stop("[Component DESeq2] Count sample identifiers must be unique", call. = FALSE)
}

metadata <- fread(opt$meta, sep = "\t", data.table = FALSE, check.names = FALSE)
required_metadata <- c(opt$sample_id_col, opt$component_col)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0L) {
  stop("[Component DESeq2] Metadata missing required column(s): ",
       paste(missing_metadata, collapse = ", "), call. = FALSE)
}
metadata[[opt$sample_id_col]] <- trimws(as.character(metadata[[opt$sample_id_col]]))
if (anyDuplicated(metadata[[opt$sample_id_col]])) {
  stop("[Component DESeq2] Metadata sample identifiers are duplicated", call. = FALSE)
}
metadata <- metadata[match(colnames(count_matrix), metadata[[opt$sample_id_col]]), , drop = FALSE]
if (any(is.na(metadata[[opt$sample_id_col]]))) {
  stop("[Component DESeq2] Metadata does not cover every count sample", call. = FALSE)
}
rownames(metadata) <- metadata[[opt$sample_id_col]]
if (!identical(colnames(count_matrix), rownames(metadata))) {
  stop("[Component DESeq2] Count and metadata sample orders do not agree", call. = FALSE)
}

component_values <- trimws(as.character(metadata[[opt$component_col]]))
target_component <- trimws(as.character(opt$component))
group_values <- ifelse(component_values == target_component, "FOCAL_COMPONENT", "REST")
metadata$contrast_group <- factor(group_values, levels = c("REST", "FOCAL_COMPONENT"))
n_focal <- sum(metadata$contrast_group == "FOCAL_COMPONENT")
n_rest <- sum(metadata$contrast_group == "REST")
if (n_focal < 1L || n_rest < 1L) {
  stop("[Component DESeq2] Component contrast has an empty focal or rest group", call. = FALSE)
}
if ((n_focal + n_rest) <= 2L) {
  stop("[Component DESeq2] Component contrast lacks residual degrees of freedom", call. = FALSE)
}

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ contrast_group)
keep_gene <- rowSums(counts(dds)) >= opt$minimum_total_gene_count
dds <- dds[keep_gene, ]
if (nrow(dds) == 0L) {
  stop("[Component DESeq2] No genes remain after the configured total-count prefilter", call. = FALSE)
}
dds <- DESeq(
  dds,
  quiet = TRUE,
  fitType = opt$dispersion_fit_type,
  minReplicatesForReplace = Inf
)
wald_result <- results(
  dds,
  contrast = c("contrast_group", "FOCAL_COMPONENT", "REST"),
  alpha = opt$adjusted_p_value_threshold
)

# apeglm provides a posterior mode for the shrunken LFC and a posterior
# standard deviation for its uncertainty.
shrunken_result <- lfcShrink(
  dds,
  coef = "contrast_group_FOCAL_COMPONENT_vs_REST",
  type = opt$lfc_shrinkage_method,
  quiet = TRUE
)

result_table <- as.data.frame(wald_result)
result_table$gene_id <- rownames(result_table)
shrunken_table <- as.data.frame(shrunken_result)
result_table$baseMean <- result_table$baseMean
result_table$wald_statistic <- result_table$stat
result_table$p_value <- result_table$pvalue
result_table$adjusted_p_value <- result_table$padj
result_table$log2_fold_change_unshrunk <- result_table$log2FoldChange
result_table$log2_fold_change_standard_error_unshrunk <- result_table$lfcSE
result_table$log2_fold_change_shrunken <- shrunken_table[rownames(result_table), "log2FoldChange"]
result_table$log2_fold_change_posterior_sd <- shrunken_table[rownames(result_table), "lfcSE"]
result_table$absolute_shrunken_log2_fold_change <- abs(result_table$log2_fold_change_shrunken)
result_table$effect_direction <- ifelse(
  result_table$log2_fold_change_shrunken > 0,
  "upregulated_in_focal_component_vs_rest",
  ifelse(result_table$log2_fold_change_shrunken < 0,
         "downregulated_in_focal_component_vs_rest",
         "zero_shrunken_effect")
)
result_table$contrast_id <- paste0("component_", target_component, "_vs_rest")
result_table$contrast_type <- "component_focal_vs_rest"
result_table$focal_component_id <- target_component
result_table$reference_definition <- "all_prepared_same_cancer_profiles_outside_focal_component"
result_table$lfc_shrinkage_method <- opt$lfc_shrinkage_method

result_columns <- c(
  "gene_id", "contrast_id", "contrast_type", "focal_component_id",
  "reference_definition", "baseMean", "wald_statistic", "p_value",
  "adjusted_p_value", "log2_fold_change_unshrunk",
  "log2_fold_change_standard_error_unshrunk", "log2_fold_change_shrunken",
  "log2_fold_change_posterior_sd", "absolute_shrunken_log2_fold_change",
  "effect_direction", "lfc_shrinkage_method", "stat", "pvalue", "padj",
  "log2FoldChange", "lfcSE"
)
result_table <- result_table[, intersect(result_columns, names(result_table)), drop = FALSE]
result_table <- result_table[order(result_table$adjusted_p_value, result_table$p_value, result_table$gene_id), , drop = FALSE]

result_path <- file.path(opt$outdir, "tables", paste0("component_", target_component, "_vs_rest.tsv"))
fwrite(result_table, result_path, sep = "\t")

passes_significance <- !is.na(result_table$adjusted_p_value) &
  result_table$adjusted_p_value <= opt$adjusted_p_value_threshold
passes_marker <- passes_significance &
  !is.na(result_table$baseMean) &
  result_table$baseMean >= opt$minimum_base_mean &
  !is.na(result_table$log2_fold_change_shrunken) &
  result_table$absolute_shrunken_log2_fold_change >= opt$minimum_absolute_shrunken_log2fc

marker_table <- result_table[passes_marker, , drop = FALSE]
marker_table$contrast_marker_rank <- seq_len(nrow(marker_table))
marker_columns <- c(
  "gene_id", "baseMean", "wald_statistic", "p_value", "adjusted_p_value",
  "log2_fold_change_unshrunk", "log2_fold_change_shrunken",
  "log2_fold_change_posterior_sd", "absolute_shrunken_log2_fold_change",
  "effect_direction", "contrast_marker_rank"
)

write_direction <- function(direction_name, rows, sort_columns) {
  if (nrow(rows) > 0L) {
    rows <- rows[do.call(order, sort_columns(rows)), , drop = FALSE]
    rows <- head(rows, opt$maximum_markers_per_direction)
    rows$contrast_marker_rank <- seq_len(nrow(rows))
  }
  out_tsv <- file.path(
    opt$outdir,
    "markers",
    paste0("component_", target_component, "_vs_rest_", direction_name, "_top", opt$maximum_markers_per_direction, ".tsv")
  )
  out_txt <- file.path(
    opt$outdir,
    "markers",
    paste0("component_", target_component, "_vs_rest_", direction_name, "_top", opt$maximum_markers_per_direction, ".txt")
  )
  fwrite(rows[, intersect(marker_columns, names(rows)), drop = FALSE], out_tsv, sep = "\t")
  writeLines(as.character(rows$gene_id), out_txt)
}

up_rows <- marker_table[marker_table$log2_fold_change_shrunken > 0, , drop = FALSE]
down_rows <- marker_table[marker_table$log2_fold_change_shrunken < 0, , drop = FALSE]
zero_rows <- marker_table[marker_table$log2_fold_change_shrunken == 0, , drop = FALSE]
if (nrow(zero_rows) > 0L) {
  stop("[Component DESeq2] Retained component marker has zero shrunken LFC despite effect-size filtering", call. = FALSE)
}

write_direction("UP", up_rows, function(rows) list(-rows$absolute_shrunken_log2_fold_change, rows$adjusted_p_value, rows$gene_id))
write_direction("DOWN", down_rows, function(rows) list(-rows$absolute_shrunken_log2_fold_change, rows$adjusted_p_value, rows$gene_id))

summary_table <- data.frame(
  component = target_component,
  n_focal_samples = n_focal,
  n_rest_samples = n_rest,
  n_result_genes = nrow(result_table),
  n_pvalue_nonmissing = sum(!is.na(result_table$p_value)),
  n_padj_nonmissing = sum(!is.na(result_table$adjusted_p_value)),
  n_significant_before_effect_filter = sum(passes_significance, na.rm = TRUE),
  n_markers_before_cap = nrow(marker_table),
  n_up_markers_after_cap = min(nrow(up_rows), opt$maximum_markers_per_direction),
  n_down_markers_after_cap = min(nrow(down_rows), opt$maximum_markers_per_direction),
  adjusted_p_value_threshold = opt$adjusted_p_value_threshold,
  minimum_base_mean = opt$minimum_base_mean,
  minimum_absolute_shrunken_log2fc = opt$minimum_absolute_shrunken_log2fc,
  maximum_markers_per_direction = opt$maximum_markers_per_direction,
  stringsAsFactors = FALSE
)
fwrite(summary_table, file.path(opt$outdir, "tables", paste0("component_", target_component, "_vs_rest_summary.tsv")), sep = "\t")
message(sprintf(
  "[Component DESeq2] PASS component=%s n_result_genes=%d up=%d down=%d",
  target_component,
  nrow(result_table),
  min(nrow(up_rows), opt$maximum_markers_per_direction),
  min(nrow(down_rows), opt$maximum_markers_per_direction)
))
