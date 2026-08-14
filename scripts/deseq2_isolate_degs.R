#!/usr/bin/env Rscript

# =============================================================================
# deseq2_isolate_degs.R
# =============================================================================
#
# Performs graph-derived, within-cancer-type differential expression and
# contrast-level marker prioritisation for:
#   1. focal isolate versus all other cell-line profiles from the same cancer
#      type; and
#   2. focal anchor versus same-cancer-type cell-line profiles outside the
#      focal anchor's resolved component.
#
# Size factors are estimated once from the filtered staged count matrix and
# transferred to contrast-specific DESeqDataSet objects. Each contrast is fitted
# with a two-level negative-binomial GLM and Wald test. apeglm-shrunken LFCs are
# used for marker effect-size filtering and ranking, while p-values and adjusted
# p-values remain those of the DESeq2 Wald test.
#
# Marker retention thresholds are contrast-specific. Retained genes are ranked
# by decreasing absolute shrunken LFC, increasing adjusted p-value, and gene
# identifier as the deterministic final tie-breaker. This module ends at retained
# per-contrast marker tables and an explicit contrast-level marker manifest.
#
# The script expects non-negative integer count data. Numerical validation alone
# cannot establish raw-count provenance; the upstream staged-input workflow must
# guarantee that the matrix is a valid count-scale DESeq2 input.
# =============================================================================

# -----------------------------------------------------------------------------
# Package Loading
# -----------------------------------------------------------------------------
# Required packages are checked before attachment so missing dependencies fail
# with one explicit preflight error.
required_pkgs <- c("optparse", "data.table", "DESeq2", "apeglm")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
})

# -----------------------------------------------------------------------------
# Command-Line Argument Definitions
# -----------------------------------------------------------------------------
# Command-line arguments are grouped by role:
#
# INPUT/OUTPUT:
#   --counts, --meta, --outdir, --counts_sep, --meta_sep
#
# COLUMN SPECIFICATIONS:
#   --sample_id_col, --cell_line_col, --component_col
#
# CONTRAST DEFINITIONS:
#   --isolate_list, --anchor_list, --anchor_components
#
# FILTERING THRESHOLDS:
#   isolate/anchor adjusted-p-value, shrunken-LFC, and marker-cap thresholds
#   --minimum_base_mean (shared contrast baseMean threshold)
#
# Isolate and anchor thresholds are specified separately because they define
# distinct contrast families in the pipeline. Their configured values should be
# reported explicitly with the corresponding marker-selection results.

opt_list <- list(
  # Input/output specifications
  make_option("--counts", type = "character",
              help = "Raw count matrix TSV/CSV. Rows=genes, columns=samples; first column must contain gene IDs."),
  make_option("--meta", type = "character",
              help = "Sample metadata TSV/CSV with one row per sample."),
  make_option("--outdir", type = "character", default = "deseq2_out",
              help = "Output directory"),
  make_option("--cancer_type", type = "character", default = "UNKNOWN",
              help = "Cancer type label written to retained marker tables and manifest"),
  make_option("--counts_sep", type = "character", default = "\t",
              help = "Separator for counts (default tab)"),
  make_option("--meta_sep", type = "character", default = "\t",
              help = "Separator for meta (default tab)"),

  # Column name specifications
  make_option("--sample_id_col", type = "character", default = "sample_name",
              help = "Column in meta that matches count column names"),
  make_option("--cell_line_col", type = "character", default = "DSMZ_Cell_line_norm",
              help = "Column in meta containing cell line names"),
  make_option("--component_col", type = "character", default = "component",
              help = "Column in meta for component membership (optional)"),

  # Contrast specifications

make_option("--isolate_list", type = "character",
              help = "Comma-separated isolate cell line names (e.g. CAL_51,COLO_824,DU_4475)"),
  make_option("--anchor_list", type = "character", default = NULL,
              help = "Comma-separated anchor cell line names (optional)"),
  make_option("--anchor_components", type = "character", default = NULL,
              help = "TSV mapping anchor->component (optional). Columns: anchor,component"),

  # Isolate marker thresholds
  make_option("--isolate_adjusted_p_value_threshold", type = "double",
              help = "BH-adjusted p-value cutoff for isolate markers"),
  make_option("--isolate_minimum_absolute_shrunken_log2fc", type = "double",
              help = "Minimum absolute apeglm-shrunken log2 fold-change for isolate markers"),
  make_option("--isolate_maximum_markers_per_contrast", type = "integer",
              help = "Maximum isolate markers retained after deterministic ranking"),

  # Anchor marker thresholds
  make_option("--anchor_adjusted_p_value_threshold", type = "double",
              help = "BH-adjusted p-value cutoff for anchor markers"),
  make_option("--anchor_minimum_absolute_shrunken_log2fc", type = "double",
              help = "Minimum absolute apeglm-shrunken log2 fold-change for anchor markers"),
  make_option("--anchor_maximum_markers_per_contrast", type = "integer",
              help = "Maximum anchor markers retained after deterministic ranking"),

  # Additional filtering parameters
  make_option("--minimum_base_mean", type = "double",
              help = "Minimum baseMean filter for markers"),
  make_option("--minimum_total_gene_count", type = "integer",
              help = "Minimum total count across prepared samples for the DESeq2 gene prefilter"),
  make_option("--dispersion_fit_type", type = "character",
              help = "DESeq2 dispersion fitType"),
  make_option("--lfc_shrinkage_method", type = "character",
              help = "DESeq2 LFC shrinkage method")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# Validate numerical selection parameters before reading large inputs.
validate_probability <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0 || x > 1) {
    stop(sprintf("%s must be in (0, 1]", name))
  }
}
validate_nonnegative <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop(sprintf("%s must be a finite non-negative number", name))
  }
}
validate_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 || x != floor(x)) {
    stop(sprintf("%s must be a positive integer", name))
  }
}

validate_probability(opt$isolate_adjusted_p_value_threshold, "--isolate_adjusted_p_value_threshold")
validate_probability(opt$anchor_adjusted_p_value_threshold, "--anchor_adjusted_p_value_threshold")
validate_nonnegative(opt$isolate_minimum_absolute_shrunken_log2fc, "--isolate_minimum_absolute_shrunken_log2fc")
validate_nonnegative(opt$anchor_minimum_absolute_shrunken_log2fc, "--anchor_minimum_absolute_shrunken_log2fc")
validate_nonnegative(opt$minimum_base_mean, "--minimum_base_mean")
validate_positive_integer(opt$isolate_maximum_markers_per_contrast, "--isolate_maximum_markers_per_contrast")
validate_positive_integer(opt$anchor_maximum_markers_per_contrast, "--anchor_maximum_markers_per_contrast")
validate_positive_integer(opt$minimum_total_gene_count, "--minimum_total_gene_count")
if (length(opt$dispersion_fit_type) != 1L || is.na(opt$dispersion_fit_type) ||
    !opt$dispersion_fit_type %in% c("parametric", "local", "mean")) {
  stop("--dispersion_fit_type must be one of: parametric, local, mean")
}
if (length(opt$lfc_shrinkage_method) != 1L || is.na(opt$lfc_shrinkage_method) ||
    opt$lfc_shrinkage_method != "apeglm") {
  stop("--lfc_shrinkage_method must be apeglm for the active marker-prioritisation method")
}

# -----------------------------------------------------------------------------
# Output Directory Initialisation
# -----------------------------------------------------------------------------
# The script creates a structured output directory hierarchy. Using
# showWarnings = FALSE prevents error messages if directories already exist,
# enabling re-runs without manual cleanup.
#
# Directory structure rationale:
#   - tables/  : All-gene contrast tables with raw Wald inference and raw/shrunken LFCs
#   - markers/ : Filtered gene lists for direct use in pathway analysis
#   - qc/      : Saved size factors for normalization audit and reproducibility

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "markers"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "qc"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------
# The script defines several utility functions to promote code reuse and
# maintain consistent behaviour across different processing steps.

# read_table_auto ---------------------------------------------------------------
# This function provides a unified interface for reading tabular data files.
# The data.table::fread function automatically detects file format details
# and accepts the configured TSV or CSV field separator.
#
# Parameters:
#   path: File path to the input table
#   sep:  Field separator (advisory; fread may auto-detect)
#
# Returns:
#   A standard R data.frame (converted from data.table for compatibility)

read_table_auto <- function(path, sep) {
  dt <- fread(path, sep = sep, data.table = FALSE)
  return(dt)
}

# validate_count_matrix --------------------------------------------------------
# This function validates that an input matrix contains raw integer counts
# before converting it to DESeq2's required integer storage mode.
#
# Non-integer values are not accepted. Decimal expression values such as VST,
# TPM, or other normalised quantities must fail validation rather than being
# truncated into pseudo-counts.

validate_count_matrix <- function(mat) {
  mat <- as.matrix(mat)
  dim_names <- dimnames(mat)
  suppressWarnings(storage.mode(mat) <- "numeric")
  dimnames(mat) <- dim_names

  if (any(is.na(mat))) {
    stop("Counts contain NA values")
  }

  if (any(!is.finite(mat))) {
    stop("Counts contain non-finite values")
  }

  if (any(mat < 0)) {
    stop("Counts contain negative values")
  }

  if (any(abs(mat - round(mat)) > 1e-8)) {
    stop("DESeq2 input must be raw integer counts. Non-integer values were detected.")
  }

  mat <- round(mat)
  storage.mode(mat) <- "integer"
  dimnames(mat) <- dim_names
  mat
}

# as_counts_matrix --------------------------------------------------------------
# This function converts the input data.frame to a validated integer matrix for DESeq2.
#
# DESeq2 requires count data as an integer matrix with:
#   - Rows representing genes (with gene identifiers as rownames)
#   - Columns representing samples
#   - Values as non-negative integers
#
# The input contract is explicit: the first column contains gene IDs and all
# remaining columns are sample counts. Gene IDs may themselves be numeric strings.
#
# Parameters:
#   df: A data.frame containing count data
#
# Returns:
#   An integer matrix with gene IDs as rownames

as_counts_matrix <- function(df) {
  if (ncol(df) < 2L) {
    stop("Count table must contain a gene-ID column and at least one sample column")
  }

  gene_id <- as.character(df[[1]])
  if (any(is.na(gene_id)) || any(!nzchar(trimws(gene_id)))) {
    stop("Gene identifiers contain missing or empty values")
  }
  if (anyDuplicated(gene_id)) {
    stop("Duplicate gene identifiers found in count matrix")
  }

  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- gene_id
  mat <- validate_count_matrix(mat)

  if (is.null(colnames(mat)) || anyDuplicated(colnames(mat))) {
    stop("Count matrix sample columns must be present and unique")
  }

  mat
}

# stop_if_missing ---------------------------------------------------------------
# This function provides standardised validation for required parameters.
# Centralising validation logic ensures consistent error messages and reduces
# code duplication in the main workflow.
#
# Parameters:
#   x:    Value to check
#   what: Error message to display if validation fails

stop_if_missing <- function(x, what) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || x == "") stop(what)
}

# split_csv ---------------------------------------------------------------------
# This function parses comma-separated strings into character vectors.
# Leading/trailing whitespace is trimmed from each element, and empty
# strings are removed to handle malformed input gracefully.
#
# Parameters:
#   x: A comma-separated string (e.g., "CAL_51, COLO_824, DU_4475")
#
# Returns:
#   A character vector of trimmed, non-empty elements

split_csv <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(trimws(x))) {
    return(character(0))
  }
  x <- trimws(unlist(strsplit(x, ",")))
  x[x != ""]
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

find_lfc_coef <- function(dds, factor_name, numerator, denominator) {
  available <- resultsNames(dds)
  candidates <- unique(c(
    paste0(factor_name, "_", numerator, "_vs_", denominator),
    paste0(factor_name, "_", make.names(numerator), "_vs_", make.names(denominator)),
    make.names(paste0(factor_name, "_", numerator, "_vs_", denominator))
  ))
  direct <- candidates[candidates %in% available]
  if (length(direct) == 1L) {
    return(direct[[1L]])
  }

  pattern <- paste0(
    "^", escape_regex(factor_name), "_",
    escape_regex(numerator), "_vs_",
    escape_regex(denominator), "$"
  )
  matches <- available[grepl(pattern, available)]
  if (length(matches) == 1L) {
    return(matches[[1L]])
  }

  NA_character_
}

apply_lfc_shrinkage <- function(dds, res_raw, factor_name, numerator, denominator, label) {
  res_df <- as.data.frame(res_raw)
  res_df$gene_id <- rownames(res_df)
  res_df$log2_fold_change_unshrunk <- res_df$log2FoldChange
  res_df$log2_fold_change_standard_error_unshrunk <- res_df$lfcSE

  coef_name <- find_lfc_coef(dds, factor_name, numerator, denominator)
  if (is.na(coef_name)) {
    stop(sprintf(
      "apeglm coefficient not found for %s (%s vs %s). Available coefficients: %s",
      label,
      numerator,
      denominator,
      paste(resultsNames(dds), collapse = ", ")
    ))
  }

  # apeglm provides a posterior mode for the shrunken LFC and a posterior
  # standard deviation for its uncertainty.
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = opt$lfc_shrinkage_method, quiet = TRUE)
  shrunk_df <- as.data.frame(res_shrunk)

  res_df$log2_fold_change_shrunken <- shrunk_df[rownames(res_df), "log2FoldChange"]
  res_df$log2_fold_change_posterior_sd <- shrunk_df[rownames(res_df), "lfcSE"]
  res_df$lfc_shrinkage_method <- opt$lfc_shrinkage_method
  res_df$lfc_shrinkage_coef <- coef_name
  res_df$lfc_shrinkage_applied <- TRUE

  res_df$wald_statistic <- res_df$stat
  res_df$p_value <- res_df$pvalue
  res_df$adjusted_p_value <- res_df$padj

  message(sprintf("[INFO] Applied apeglm LFC shrinkage for %s using coef=%s", label, coef_name))
  res_df
}

# write_deg_outputs -------------------------------------------------------------
# This function separates Wald-test inference from contrast-level marker
# prioritisation. Adjusted p-values define statistical significance; shrunken
# log2 fold changes define effect-size filtering, direction assignment, and
# deterministic marker ranking.

write_deg_outputs <- function(res_df, prefix, outdir_tables, outdir_markers,
                              cancer_type, contrast_id, contrast_type,
                              marker_evidence_stratum, focal_profile_id,
                              focal_component_id, reference_definition,
                              fdr_threshold,
                              minimum_absolute_shrunken_log2fc,
                              maximum_markers_per_contrast,
                              minimum_base_mean) {

  if (!("log2_fold_change_shrunken" %in% colnames(res_df))) {
    stop("Shrunken LFC column is missing; marker selection requires apeglm output")
  }

  res_df$cancer_type <- cancer_type
  res_df$contrast_id <- contrast_id
  res_df$contrast_type <- contrast_type
  res_df$marker_evidence_stratum <- marker_evidence_stratum
  res_df$focal_profile_id <- focal_profile_id
  res_df$focal_component_id <- focal_component_id
  res_df$reference_definition <- reference_definition
  res_df$absolute_shrunken_log2_fold_change <- abs(res_df$log2_fold_change_shrunken)
  res_df$effect_direction <- ifelse(
    res_df$log2_fold_change_shrunken > 0,
    "upregulated_in_focal_vs_reference",
    ifelse(
      res_df$log2_fold_change_shrunken < 0,
      "downregulated_in_focal_vs_reference",
      "zero_shrunken_effect"
    )
  )

  passes_statistical_significance <- !is.na(res_df$adjusted_p_value) &
    (res_df$adjusted_p_value <= fdr_threshold)
  passes_contrast_marker_selection <- passes_statistical_significance &
    !is.na(res_df$baseMean) &
    (res_df$baseMean >= minimum_base_mean) &
    !is.na(res_df$log2_fold_change_shrunken) &
    (res_df$absolute_shrunken_log2_fold_change >= minimum_absolute_shrunken_log2fc)

  res_df$passes_statistical_significance <- passes_statistical_significance
  res_df$passes_contrast_marker_selection <- passes_contrast_marker_selection

  full_path <- file.path(outdir_tables, paste0(prefix, ".tsv"))
  fwrite(res_df, full_path, sep = "\t")

  markers_before_cap <- res_df[passes_contrast_marker_selection, , drop = FALSE]
  markers <- markers_before_cap
  markers <- markers[
    order(
      -markers$absolute_shrunken_log2_fold_change,
      markers$adjusted_p_value,
      markers$gene_id
    ),
    ,
    drop = FALSE
  ]

  if (nrow(markers) > maximum_markers_per_contrast) {
    markers <- markers[seq_len(maximum_markers_per_contrast), , drop = FALSE]
  }
  markers$contrast_marker_rank <- seq_len(nrow(markers))

  marker_genes <- as.character(markers$gene_id)
  marker_path <- file.path(
    outdir_markers,
    paste0(prefix, "_markers_top", maximum_markers_per_contrast, ".txt")
  )
  writeLines(marker_genes, marker_path)

  retained_marker_table_path <- file.path(outdir_markers, paste0(prefix, "_retained_markers.tsv"))
  retained_column_order <- c(
    "gene_id", "cancer_type", "contrast_id", "contrast_type",
    "marker_evidence_stratum", "focal_profile_id", "focal_component_id",
    "reference_definition", "baseMean", "wald_statistic", "p_value",
    "adjusted_p_value", "log2_fold_change_unshrunk",
    "log2_fold_change_shrunken",
    "log2_fold_change_standard_error_unshrunk",
    "log2_fold_change_posterior_sd",
    "absolute_shrunken_log2_fold_change", "effect_direction",
    "contrast_marker_rank", "passes_contrast_marker_selection",
    "passes_statistical_significance"
  )
  retained_marker_table <- markers[
    ,
    intersect(retained_column_order, colnames(markers)),
    drop = FALSE
  ]
  fwrite(retained_marker_table, retained_marker_table_path, sep = "\t")

  list(
    full_path = full_path,
    marker_path = marker_path,
    retained_marker_table_path = retained_marker_table_path,
    n_markers = length(marker_genes),
    marker_genes = marker_genes,
    n_result_genes = nrow(res_df),
    n_pvalue_nonmissing = sum(!is.na(res_df$p_value)),
    n_padj_nonmissing = sum(!is.na(res_df$adjusted_p_value)),
    n_significant_before_effect_filter = sum(passes_statistical_significance, na.rm = TRUE),
    n_markers_before_cap = nrow(markers_before_cap),
    n_markers_after_cap = length(marker_genes),
    adjusted_p_value_threshold = fdr_threshold,
    minimum_base_mean = minimum_base_mean,
    minimum_absolute_shrunken_log2fc = minimum_absolute_shrunken_log2fc,
    maximum_markers_per_contrast = maximum_markers_per_contrast,
    contrast_type = contrast_type,
    marker_evidence_stratum = marker_evidence_stratum,
    focal_profile_id = focal_profile_id,
    focal_component_id = focal_component_id,
    reference_definition = reference_definition
  )
}

# -----------------------------------------------------------------------------
# Data Loading and Validation
# -----------------------------------------------------------------------------
# This section loads inputs and enforces the matrix/metadata contracts used below.
# Early validation prevents cryptic errors later in the analysis pipeline.

stop_if_missing(opt$counts, "--counts is required")
stop_if_missing(opt$meta, "--meta is required")

counts_df <- read_table_auto(opt$counts, opt$counts_sep)
meta_df   <- read_table_auto(opt$meta, opt$meta_sep)

stop_if_missing(opt$sample_id_col, "--sample_id_col is required")
stop_if_missing(opt$cell_line_col, "--cell_line_col is required")

# Validate that specified columns exist in metadata
if (!(opt$sample_id_col %in% colnames(meta_df))) {
  stop(paste0("meta missing column: ", opt$sample_id_col))
}
if (!(opt$cell_line_col %in% colnames(meta_df))) {
  stop(paste0("meta missing column: ", opt$cell_line_col))
}

counts_mat <- as_counts_matrix(counts_df)

# -----------------------------------------------------------------------------
# Sample Alignment
# -----------------------------------------------------------------------------
# DESeq2 requires that metadata rows exactly match count matrix columns.
# This section ensures proper alignment and validates sample correspondence.
#
# The alignment process:
#   1. Identifies samples present in both counts and metadata
#   2. Subsets metadata to matching samples
#   3. Reorders metadata to match count matrix column order
#   4. Sets rownames on metadata (required by DESeq2)

sample_ids <- colnames(counts_mat)
meta_df[[opt$sample_id_col]] <- as.character(meta_df[[opt$sample_id_col]])
if (any(is.na(meta_df[[opt$sample_id_col]])) ||
    any(!nzchar(trimws(meta_df[[opt$sample_id_col]])))) {
  stop("Metadata sample IDs contain missing or empty values")
}
if (anyDuplicated(meta_df[[opt$sample_id_col]])) {
  stop("Duplicate sample IDs found in metadata")
}
meta_sub <- meta_df[meta_df[[opt$sample_id_col]] %in% sample_ids, , drop = FALSE]

if (nrow(meta_sub) == 0) {
  stop("No overlapping samples between counts columns and meta sample_id_col")
}

# Reorder metadata to match count matrix column order
meta_sub <- meta_sub[match(sample_ids, meta_sub[[opt$sample_id_col]]), , drop = FALSE]
if (any(is.na(meta_sub[[opt$sample_id_col]]))) {
  missing <- sample_ids[is.na(meta_sub[[opt$sample_id_col]])]
  stop(paste0("Meta missing rows for samples: ", paste(missing, collapse = ", ")))
}

# Set rownames as required by DESeq2 and verify exact alignment
rownames(meta_sub) <- meta_sub[[opt$sample_id_col]]
if (!identical(colnames(counts_mat), rownames(meta_sub))) {
  stop("Count matrix columns and metadata rownames are not aligned")
}

# -----------------------------------------------------------------------------
# DESeqDataSet Creation and Preprocessing
# -----------------------------------------------------------------------------
# This section creates the base DESeq2 object and performs preprocessing steps
# that are common to all contrasts.
#
# Size factors are estimated once from the staged cohort and reused in all
# isolate and anchor contrast datasets. Contrast-specific designs use an
# explicit two-level focal-versus-reference factor named contrast_group.

cell_line_values <- as.character(meta_sub[[opt$cell_line_col]])
if (any(is.na(cell_line_values)) || any(!nzchar(trimws(cell_line_values)))) {
  stop("Cell-line identifiers contain missing or empty values")
}
meta_sub[[opt$cell_line_col]] <- factor(cell_line_values)

# Create a base DESeqDataSet for common gene filtering and size-factor estimation.
dds_base <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = meta_sub,
  design    = ~ 1
)

# -----------------------------------------------------------------------------
# Low-Count Gene Filtering
# -----------------------------------------------------------------------------
# A pipeline-specific light pre-filter removes genes below the configured total
# count floor across the prepared matrix. This computational filter is distinct from
# DESeq2's results-stage independent filtering.

keep_gene <- rowSums(counts(dds_base)) >= opt$minimum_total_gene_count
dds_base <- dds_base[keep_gene, ]
if (nrow(dds_base) == 0L) {
  stop("No genes remain after the total-count pre-filter")
}

# -----------------------------------------------------------------------------
# Size Factor Estimation
# -----------------------------------------------------------------------------
# Size factors provide sample-wise scaling for library-size/composition
# differences. The standard DESeq2 median-of-ratios estimator is used here.
#
# The estimated size factors are saved for inspection and reproducibility.
# Their interpretation is as sample-specific normalization factors, not as a
# stand-alone sample-quality criterion.

dds_base <- estimateSizeFactors(dds_base)

sf <- data.frame(
  sample_id = colnames(dds_base),
  size_factor = sizeFactors(dds_base)
)
fwrite(sf, file.path(opt$outdir, "qc", "size_factors.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# Contrast Specification Parsing
# -----------------------------------------------------------------------------
# The script parses the comma-separated lists of isolates and anchors into
# character vectors for iteration.

isolates <- unique(split_csv(opt$isolate_list))
anchors <- unique(split_csv(opt$anchor_list))

if (length(isolates) == 0L && length(anchors) == 0L) {
  stop("No DESeq2 marker contrasts were requested: provide at least one isolate or anchor contrast.")
}

message(sprintf(
  "[INFO] Planned DESeq2 marker contrasts: %d isolate, %d anchor",
  length(isolates),
  length(anchors)
))

# Parse anchor-to-component mapping if provided
anchor2component <- NULL
if (!is.null(opt$anchor_components)) {
  a2c <- fread(opt$anchor_components, sep = "\t", data.table = FALSE)
  if (!all(c("anchor", "component") %in% colnames(a2c))) {
    stop("--anchor_components must have columns: anchor, component")
  }
  a2c$anchor <- trimws(as.character(a2c$anchor))
  a2c$component <- trimws(as.character(a2c$component))
  if (any(is.na(a2c$anchor)) || any(!nzchar(a2c$anchor)) ||
      any(is.na(a2c$component)) || any(!nzchar(a2c$component))) {
    stop("Anchor-to-component mapping contains missing or empty values")
  }
  component_sets <- split(a2c$component, a2c$anchor)
  conflicting <- names(component_sets)[
    vapply(component_sets, function(x) length(unique(x)) > 1L, logical(1))
  ]
  if (length(conflicting) > 0L) {
    stop("Conflicting component mappings for anchors: ", paste(conflicting, collapse = ", "))
  }
  a2c <- unique(a2c[, c("anchor", "component"), drop = FALSE])
  anchor2component <- setNames(a2c$component, a2c$anchor)
}

available_cell_lines <- levels(meta_sub[[opt$cell_line_col]])
missing_isolates <- setdiff(isolates, available_cell_lines)
if (length(missing_isolates) > 0L) {
  stop("Requested isolate cell lines are absent from the staged data: ",
       paste(missing_isolates, collapse = ", "))
}

missing_anchors <- setdiff(anchors, available_cell_lines)
if (length(missing_anchors) > 0L) {
  stop("Requested anchor cell lines are absent from the staged data: ",
       paste(missing_anchors, collapse = ", "))
}

if (length(anchors) > 0L) {
  if (is.null(anchor2component)) {
    stop("Anchor contrasts were requested but --anchor_components was not provided")
  }
  missing_anchor_mappings <- setdiff(anchors, names(anchor2component))
  if (length(missing_anchor_mappings) > 0L) {
    stop("Requested anchors are missing from the anchor-to-component mapping: ",
         paste(missing_anchor_mappings, collapse = ", "))
  }
}

# -----------------------------------------------------------------------------
# Isolate Marker Analysis
# -----------------------------------------------------------------------------
# For each isolate cell line, the script performs a differential expression
# analysis comparing the focal isolate against all other same-cancer cell-line
# profiles.
#
# Statistical Approach:
#   - Binary grouping: isolate vs REST
#   - Contrast extracts the log2 fold change of isolate relative to REST
#   - Positive LFC indicates higher expression in the isolate
#   - Negative LFC indicates lower expression in the isolate
#
# Focal groups may contain one sample while the reference group contains
# multiple samples. Such contrasts have no within-cell-line replication on the
# focal side; they are used for marker prioritisation and interpreted with that
# limitation. apeglm shrinkage is applied to stabilize effect-size ranking.

all_marker_outputs <- list()

for (focal_isolate_profile_id in isolates) {
  same_cancer_reference_profile_ids <- as.character(
    meta_sub[[opt$cell_line_col]][meta_sub[[opt$cell_line_col]] != focal_isolate_profile_id]
  )
  contrast_group_values <- ifelse(
    meta_sub[[opt$cell_line_col]] == focal_isolate_profile_id,
    "FOCAL_ISOLATE",
    "OTHER_SAME_CANCER"
  )
  contrast_metadata <- meta_sub
  contrast_metadata$contrast_group <- factor(
    contrast_group_values,
    levels = c("OTHER_SAME_CANCER", "FOCAL_ISOLATE")
  )

  focal_isolate_sample_count <- sum(contrast_group_values == "FOCAL_ISOLATE")
  same_cancer_reference_sample_count <- sum(contrast_group_values == "OTHER_SAME_CANCER")

  if (focal_isolate_sample_count < 1L || same_cancer_reference_sample_count < 1L) {
    stop("Invalid isolate contrast with an empty group: ", focal_isolate_profile_id,
         " vs other same-cancer profiles")
  }
  if ((focal_isolate_sample_count + same_cancer_reference_sample_count) <= 2L) {
    stop("Isolate contrast lacks residual degrees of freedom: ",
         focal_isolate_profile_id, " vs other same-cancer profiles")
  }

  isolate_contrast_dds <- DESeqDataSetFromMatrix(
    countData = counts(dds_base),
    colData   = contrast_metadata,
    design    = ~ contrast_group
  )

  sizeFactors(isolate_contrast_dds) <- sizeFactors(dds_base)

  isolate_contrast_dds <- DESeq(
    isolate_contrast_dds,
    quiet = TRUE,
    fitType = opt$dispersion_fit_type,
    minReplicatesForReplace = Inf  # Disable automatic count replacement; results() Cook filtering remains active
  )

  res_raw <- results(
    isolate_contrast_dds,
    contrast = c("contrast_group", "FOCAL_ISOLATE", "OTHER_SAME_CANCER"),
    alpha = opt$isolate_adjusted_p_value_threshold
  )

  res_df <- apply_lfc_shrinkage(
    isolate_contrast_dds,
    res_raw,
    factor_name = "contrast_group",
    numerator = "FOCAL_ISOLATE",
    denominator = "OTHER_SAME_CANCER",
    label = focal_isolate_profile_id
  )

  col_order <- c("gene_id", "baseMean", "wald_statistic", "p_value",
                 "adjusted_p_value", "log2_fold_change_unshrunk",
                 "log2_fold_change_shrunken",
                 "log2_fold_change_standard_error_unshrunk",
                 "log2_fold_change_posterior_sd",
                 "lfc_shrinkage_method", "lfc_shrinkage_coef",
                 "lfc_shrinkage_applied", "lfcSE",
                 "stat", "pvalue", "padj")
  res_df <- res_df[, intersect(col_order, colnames(res_df)), drop = FALSE]

  prefix <- paste0("isolate_", focal_isolate_profile_id, "_vs_other_same_cancer")
  out <- write_deg_outputs(
    res_df, prefix,
    outdir_tables = file.path(opt$outdir, "tables"),
    outdir_markers = file.path(opt$outdir, "markers"),
    cancer_type = opt$cancer_type,
    contrast_id = prefix,
    contrast_type = "isolate_focal_vs_other_same_cancer",
    marker_evidence_stratum = "isolate_associated",
    focal_profile_id = focal_isolate_profile_id,
    focal_component_id = NA_character_,
    reference_definition = paste0(
      "all_other_cell_line_profiles_from_same_cancer_type;",
      "reference_profile_count=", length(unique(same_cancer_reference_profile_ids))
    ),
    fdr_threshold = opt$isolate_adjusted_p_value_threshold,
    minimum_absolute_shrunken_log2fc = opt$isolate_minimum_absolute_shrunken_log2fc,
    maximum_markers_per_contrast = opt$isolate_maximum_markers_per_contrast,
    minimum_base_mean = opt$minimum_base_mean
  )

  all_marker_outputs[[prefix]] <- out
  message(sprintf("[OK] %s: %d markers", prefix, out$n_markers))
}

# -----------------------------------------------------------------------------
# Anchor Marker Analysis (Optional)
# -----------------------------------------------------------------------------
# For each selected graph anchor, the focal cell line is compared with all
# staged cell-line samples outside its resolved component. Other members of the
# anchor's component are excluded from that contrast. Samples without a resolved
# component assignment (graph isolates) are outside the focal component and are
# therefore retained in the reference group.

if (length(anchors) > 0) {
  # Validate metadata required for anchor analysis.
  if (!(opt$component_col %in% colnames(meta_sub))) {
    stop(paste0("meta missing component column required for anchor contrasts: ", opt$component_col))
  }

  meta_sub[[opt$component_col]] <- as.character(meta_sub[[opt$component_col]])

  for (focal_anchor_profile_id in anchors) {
    if (is.na(anchor2component[[focal_anchor_profile_id]]) ||
        !nzchar(trimws(anchor2component[[focal_anchor_profile_id]]))) {
      stop("Anchor has an empty component mapping: ", focal_anchor_profile_id)
    }
    focal_component_id <- anchor2component[[focal_anchor_profile_id]]

    component_value <- meta_sub[[opt$component_col]]
    component_missing <- is.na(component_value) |
      !nzchar(trimws(component_value)) |
      toupper(trimws(component_value)) == "NA"
    anchor_rows <- meta_sub[[opt$cell_line_col]] == focal_anchor_profile_id
    if (any(component_missing[anchor_rows]) ||
        any(component_value[anchor_rows] != focal_component_id, na.rm = TRUE)) {
      stop(sprintf(
        "Anchor %s metadata component does not match anchor mapping (%s)",
        focal_anchor_profile_id,
        focal_component_id
      ))
    }
    outside_focal_component <- component_missing | component_value != focal_component_id
    keep_samples <- (meta_sub[[opt$cell_line_col]] == focal_anchor_profile_id) |
      outside_focal_component
    anchor_contrast_metadata <- meta_sub[keep_samples, , drop = FALSE]

    outside_component_reference_profile_ids <- as.character(
      anchor_contrast_metadata[[opt$cell_line_col]][
        anchor_contrast_metadata[[opt$cell_line_col]] != focal_anchor_profile_id
      ]
    )
    contrast_group_values <- ifelse(
      anchor_contrast_metadata[[opt$cell_line_col]] == focal_anchor_profile_id,
      "FOCAL_ANCHOR",
      "OUTSIDE_FOCAL_COMPONENT"
    )
    anchor_contrast_metadata$contrast_group <- factor(
      contrast_group_values,
      levels = c("OUTSIDE_FOCAL_COMPONENT", "FOCAL_ANCHOR")
    )

    focal_anchor_sample_count <- sum(contrast_group_values == "FOCAL_ANCHOR")
    outside_component_reference_sample_count <- sum(
      contrast_group_values == "OUTSIDE_FOCAL_COMPONENT"
    )

    if (focal_anchor_sample_count < 1L ||
        outside_component_reference_sample_count < 1L) {
      stop("Invalid anchor contrast with an empty group: ", focal_anchor_profile_id,
           " vs outside focal component ", focal_component_id)
    }
    if ((focal_anchor_sample_count + outside_component_reference_sample_count) <= 2L) {
      stop("Anchor contrast lacks residual degrees of freedom: ",
           focal_anchor_profile_id, " vs outside focal component ", focal_component_id)
    }

    anchor_contrast_counts <- counts(dds_base)[
      ,
      rownames(anchor_contrast_metadata),
      drop = FALSE
    ]
    anchor_contrast_dds <- DESeqDataSetFromMatrix(
      countData = anchor_contrast_counts,
      colData   = anchor_contrast_metadata,
      design    = ~ contrast_group
    )

    sizeFactors(anchor_contrast_dds) <- sizeFactors(dds_base)[rownames(anchor_contrast_metadata)]

    anchor_contrast_dds <- DESeq(
      anchor_contrast_dds,
      quiet = TRUE,
      fitType = opt$dispersion_fit_type,
      minReplicatesForReplace = Inf  # Disable automatic count replacement; results() Cook filtering remains active
    )

    res_raw <- results(
      anchor_contrast_dds,
      contrast = c("contrast_group", "FOCAL_ANCHOR", "OUTSIDE_FOCAL_COMPONENT"),
      alpha = opt$anchor_adjusted_p_value_threshold
    )

    res_df <- apply_lfc_shrinkage(
      anchor_contrast_dds,
      res_raw,
      factor_name = "contrast_group",
      numerator = "FOCAL_ANCHOR",
      denominator = "OUTSIDE_FOCAL_COMPONENT",
      label = paste0("anchor ", focal_anchor_profile_id)
    )

    col_order <- c("gene_id", "baseMean", "wald_statistic", "p_value",
                   "adjusted_p_value", "log2_fold_change_unshrunk",
                   "log2_fold_change_shrunken",
                   "log2_fold_change_standard_error_unshrunk",
                   "log2_fold_change_posterior_sd",
                   "lfc_shrinkage_method", "lfc_shrinkage_coef",
                   "lfc_shrinkage_applied", "lfcSE",
                   "stat", "pvalue", "padj")
    res_df <- res_df[, intersect(col_order, colnames(res_df)), drop = FALSE]

    prefix <- paste0("anchor_", focal_anchor_profile_id, "_vs_outside_component_", focal_component_id)
    out <- write_deg_outputs(
      res_df, prefix,
      outdir_tables = file.path(opt$outdir, "tables"),
      outdir_markers = file.path(opt$outdir, "markers"),
      cancer_type = opt$cancer_type,
      contrast_id = prefix,
      contrast_type = "anchor_focal_vs_outside_focal_component",
      marker_evidence_stratum = "anchor_associated",
      focal_profile_id = focal_anchor_profile_id,
      focal_component_id = focal_component_id,
      reference_definition = paste0(
        "same_cancer_profiles_outside_focal_resolved_component;",
        "includes_isolates_and_other_components;",
        "reference_profile_count=", length(unique(outside_component_reference_profile_ids))
      ),
      fdr_threshold = opt$anchor_adjusted_p_value_threshold,
      minimum_absolute_shrunken_log2fc = opt$anchor_minimum_absolute_shrunken_log2fc,
      maximum_markers_per_contrast = opt$anchor_maximum_markers_per_contrast,
      minimum_base_mean = opt$minimum_base_mean
    )

    all_marker_outputs[[prefix]] <- out
    message(sprintf("[OK] %s: %d markers", prefix, out$n_markers))
  }
}

# -----------------------------------------------------------------------------
# Contrast-Level Marker Manifest
# -----------------------------------------------------------------------------
# The manifest is the computational boundary between DESeq2 inference /
# contrast-level marker prioritisation and downstream graph-derived marker
# aggregation. No cross-contrast recurrence or pan-cancer feature-selection
# quantities are calculated in this module.

if (length(all_marker_outputs) == 0) {
  stop("No marker sets produced (check isolate names, inputs, or thresholds).")
}

contrast_ids <- names(all_marker_outputs)
manifest <- data.frame(
  cancer_type = opt$cancer_type,
  contrast_id = contrast_ids,
  contrast_type = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$contrast_type,
    character(1)
  ),
  marker_evidence_stratum = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$marker_evidence_stratum,
    character(1)
  ),
  focal_profile_id = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$focal_profile_id,
    character(1)
  ),
  focal_component_id = vapply(
    all_marker_outputs[contrast_ids],
    function(x) ifelse(is.na(x$focal_component_id), "", as.character(x$focal_component_id)),
    character(1)
  ),
  reference_definition = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$reference_definition,
    character(1)
  ),
  marker_table_path = vapply(
    all_marker_outputs[contrast_ids],
    function(x) file.path("markers", basename(x$retained_marker_table_path)),
    character(1)
  ),
  marker_gene_list_path = vapply(
    all_marker_outputs[contrast_ids],
    function(x) file.path("markers", basename(x$marker_path)),
    character(1)
  ),
  result_table_path = vapply(
    all_marker_outputs[contrast_ids],
    function(x) file.path("tables", basename(x$full_path)),
    character(1)
  ),
  n_result_genes = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$n_result_genes,
    integer(1)
  ),
  n_pvalue_nonmissing = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$n_pvalue_nonmissing,
    integer(1)
  ),
  n_padj_nonmissing = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$n_padj_nonmissing,
    integer(1)
  ),
  n_significant_before_effect_filter = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$n_significant_before_effect_filter,
    integer(1)
  ),
  n_markers_before_cap = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$n_markers_before_cap,
    integer(1)
  ),
  n_markers_after_cap = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$n_markers_after_cap,
    integer(1)
  ),
  adjusted_p_value_threshold = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$adjusted_p_value_threshold,
    numeric(1)
  ),
  minimum_base_mean = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$minimum_base_mean,
    numeric(1)
  ),
  minimum_absolute_shrunken_log2fc = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$minimum_absolute_shrunken_log2fc,
    numeric(1)
  ),
  maximum_markers_per_contrast = vapply(
    all_marker_outputs[contrast_ids],
    function(x) x$maximum_markers_per_contrast,
    integer(1)
  )
)
fwrite(manifest,
       file.path(opt$outdir, "markers", "contrast_level_marker_manifest.tsv"),
       sep = "\t")

# -----------------------------------------------------------------------------
# Execution Summary
# -----------------------------------------------------------------------------
# The script outputs summary statistics to facilitate validation and logging.

message(sprintf("[DONE] Wrote %d retained contrast-level marker manifests",
                nrow(manifest)))
message(sprintf("[OUT] %s",
                file.path(opt$outdir, "markers", "contrast_level_marker_manifest.tsv")))

# Write session info for reproducibility
writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))
message("[INFO] Session info written to sessionInfo.txt")
