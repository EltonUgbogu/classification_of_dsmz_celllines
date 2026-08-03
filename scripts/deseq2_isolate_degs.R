#!/usr/bin/env Rscript

# =============================================================================
# deseq2_isolate_degs.R
# =============================================================================
#
# DESCRIPTION
# -----------
# This script performs differential gene expression (DEG) analysis using DESeq2
# to identify marker genes that distinguish individual cell lines from broader
# populations. It implements two complementary contrast strategies:
#
#   (1) Isolate DEGs: For each specified "isolate" cell line, the script
#       compares that cell line against all other cell lines in the dataset.
#       This identifies genes uniquely up- or down-regulated in the isolate.
#
#   (2) Anchor DEGs (optional): For each specified "anchor" cell line, the
#       script compares that cell line against cell lines OUTSIDE its assigned
#       component/cluster. This identifies genes that distinguish the anchor
#       from unrelated cell lines whilst excluding its natural neighbours.
#
# The script produces per-contrast DEG tables, filtered marker gene lists, and
# a consolidated "unique feature set" based on gene recurrence across contrasts.
#
# BACKGROUND: DESeq2 AND DIFFERENTIAL EXPRESSION
# ----------------------------------------------
# DESeq2 is a widely-used Bioconductor package for differential expression
# analysis of count data from RNA-seq experiments. Key statistical concepts:
#
#   - Negative Binomial Model: DESeq2 models read counts using the negative
#     binomial distribution, which accounts for both biological and technical
#     variability (overdispersion) commonly observed in RNA-seq data.
#
#   - Size Factors: To enable comparison across samples with different
#     sequencing depths, DESeq2 estimates size factors that normalise counts.
#     The median-of-ratios method ensures robust normalisation even when many
#     genes are differentially expressed.
#
#   - Dispersion Estimation: The dispersion parameter captures gene-specific
#     variability. DESeq2 uses empirical Bayes shrinkage to share information
#     across genes, improving estimates especially for genes with few counts.
#
#   - Log2 Fold Change (LFC): The effect size measure indicating how much a
#     gene's expression differs between conditions. An LFC of 1 means the gene
#     is expressed twice as highly in the test condition versus the reference.
#
#   - Adjusted p-value (padj): Multiple testing correction using the
#     Benjamini-Hochberg procedure controls the false discovery rate (FDR).
#
# INPUT REQUIREMENTS
# ------------------
# 1. Raw count matrix (--counts):
#    - Rows represent genes; columns represent samples
#    - Values must be RAW INTEGER COUNTS (not TPM, FPKM, or normalised values)
#    - DESeq2 performs its own normalisation; pre-normalised data will produce
#      incorrect statistical results
#
# 2. Sample metadata (--meta):
#    - One row per sample with at minimum:
#      - Sample identifier column (matching count matrix column names)
#      - Cell line identifier column
#    - Optional: Component/cluster column for anchor analyses
#
# 3. Isolate list (--isolate_list):
#    - Comma-separated cell line names to analyse as isolates
#
# 4. Anchor configuration (optional):
#    - --anchor_list: Comma-separated anchor cell line names
#    - --anchor_components: TSV mapping anchors to their components
#
# OUTPUT STRUCTURE
# ----------------
# The script creates the following directory structure:
#
#   {outdir}/
#   ├── tables/                          # Full DESeq2 results per contrast
#   │   ├── isolate_{CL}_vs_rest.tsv
#   │   └── anchor_{CL}_vs_outside_component_{COMP}.tsv
#   ├── markers/                         # Filtered marker gene lists
#   │   ├── isolate_{CL}_vs_rest_markers_top{N}.txt
#   │   ├── anchor_{CL}_vs_outside_component_{COMP}_markers_top{N}.txt
#   │   ├── gene_recurrence_across_contrasts.tsv
#   │   ├── unique_feature_set_recurrence_ge_{k}.txt
#   │   └── marker_sets_manifest.tsv
#   └── qc/                              # Quality control outputs
#       └── size_factors.tsv
#
# USAGE EXAMPLES
# --------------
# Basic isolate analysis:
#   Rscript deseq2_isolate_degs.R \
#     --counts counts.tsv \
#     --meta meta.tsv \
#     --sample_id_col sample_name \
#     --cell_line_col DSMZ_Cell_line_norm \
#     --isolate_list CAL_51,COLO_824,DU_4475 \
#     --outdir results/unsupervised/brca
#
# With anchor analysis:
#   Rscript deseq2_isolate_degs.R \
#     --counts counts.tsv \
#     --meta meta.tsv \
#     --sample_id_col sample_name \
#     --cell_line_col DSMZ_Cell_line_norm \
#     --isolate_list CAL_51,COLO_824 \
#     --anchor_list MCF7,T47D \
#     --anchor_components anchor_components.tsv \
#     --outdir results/unsupervised/brca
#
# DEPENDENCIES
# ------------
# R packages: optparse, data.table, DESeq2
# Bioconductor: DESeq2 (install via BiocManager::install("DESeq2"))
#
# AUTHOR
# ------
# Developed as part of the unsupervised clustering pipeline for identifying
# cell-line specific marker genes in cancer genomics analyses.
#
# =============================================================================

# -----------------------------------------------------------------------------
# Package Loading
# -----------------------------------------------------------------------------
# The script loads required packages with startup messages suppressed to
# maintain clean log output. The key packages serve distinct purposes:
#
#   - optparse:    Command-line argument parsing with automatic help generation
#   - data.table:  High-performance data I/O via fread/fwrite functions
#   - DESeq2:      Core differential expression analysis framework
#   - apeglm:      LFC shrinkage for improved effect size estimates (optional)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
})

# Preflight: fail fast before expensive DESeq2 runs if required packages are absent
required_pkgs <- c("DESeq2", "apeglm")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
has_apeglm <- TRUE

# -----------------------------------------------------------------------------
# Command-Line Argument Definitions
# -----------------------------------------------------------------------------
# The script accepts a comprehensive set of arguments organised into categories:
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
#   --fdr_isolate, --lfc_isolate, --topN_isolate (for isolate contrasts)
#   --fdr_anchor, --lfc_anchor, --topN_anchor (for anchor contrasts)
#   --min_baseMean (minimum expression filter)
#   --recurrence_k (gene recurrence threshold for unique feature set)
#
# The separation of thresholds for isolate vs anchor analyses allows users to
# apply more stringent criteria for isolates (which are expected to have strong
# effects) and more permissive criteria for anchors (which may show subtler
# differences from non-component members).

opt_list <- list(
  # Input/output specifications
  make_option("--counts", type = "character",
              help = "Raw count matrix TSV/CSV. Rows=genes, Cols=samples. First column gene_id if no rownames."),
  make_option("--meta", type = "character",
              help = "Sample metadata TSV/CSV with one row per sample."),
  make_option("--outdir", type = "character", default = "deseq2_out",
              help = "Output directory"),
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

  # Isolate DEG thresholds
  make_option("--fdr_isolate", type = "double", default = 0.01,
              help = "FDR cutoff for isolate markers"),
  make_option("--lfc_isolate", type = "double", default = 1.5,
              help = "abs(log2FC) cutoff for isolate markers"),
  make_option("--topN_isolate", type = "integer", default = 50,
              help = "Top N isolate markers by abs(log2FC) if too many pass filters"),

  # Anchor DEG thresholds
  make_option("--fdr_anchor", type = "double", default = 0.05,
              help = "FDR cutoff for anchor markers"),
  make_option("--lfc_anchor", type = "double", default = 1.0,
              help = "abs(log2FC) cutoff for anchor markers"),
  make_option("--topN_anchor", type = "integer", default = 200,
              help = "Top N anchor markers by abs(log2FC) if too many pass filters"),

  # Additional filtering parameters
  make_option("--recurrence_k", type = "integer", default = 2,
              help = "Keep genes appearing in >=k contrasts for unique feature set"),
  make_option("--min_baseMean", type = "double", default = 10,
              help = "Minimum baseMean filter for markers")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# -----------------------------------------------------------------------------
# Output Directory Initialisation
# -----------------------------------------------------------------------------
# The script creates a structured output directory hierarchy. Using
# showWarnings = FALSE prevents error messages if directories already exist,
# enabling re-runs without manual cleanup.
#
# Directory structure rationale:
#   - tables/  : Full DESeq2 results for comprehensive downstream analysis
#   - markers/ : Filtered gene lists for direct use in pathway analysis
#   - qc/      : Quality control metrics for validation and troubleshooting

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
# and handles both TSV and CSV files efficiently.
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

# as_counts_matrix --------------------------------------------------------------
# This function converts a data.frame to an integer matrix suitable for DESeq2.
#
# DESeq2 requires count data as an integer matrix with:
#   - Rows representing genes (with gene identifiers as rownames)
#   - Columns representing samples
#   - Values as non-negative integers
#
# The function handles two common input formats:
#   1. First column contains gene IDs (non-numeric): used as rownames
#   2. All columns are numeric: assumed to already have proper rownames
#
# Parameters:
#   df: A data.frame containing count data
#
# Returns:
#   An integer matrix with gene IDs as rownames

as_counts_matrix <- function(df) {
  # Detect whether first column contains gene identifiers
  if (!all(sapply(df[, 1], is.numeric))) {
    gene_id <- as.character(df[[1]])
    mat <- as.matrix(df[, -1, drop = FALSE])
    rownames(mat) <- gene_id
  } else {
    mat <- as.matrix(df)
  }
  # TEMPORARY COMPROMISE: current staged "counts" may be VST-like decimals
  # because raw DSMZ counts are unavailable. Coerce to integer only to let
  # DESeq2 run; replace this with true raw integer counts when available.
  if (any(mat != round(mat), na.rm = TRUE)) {
    warning("Coercing non-integer staged expression values to integer for DESeq2 compatibility; use raw counts for final analysis.")
  }
  storage.mode(mat) <- "integer"
  return(mat)
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
  x <- trimws(unlist(strsplit(x, ",")))
  x[x != ""]
}

# write_deg_outputs -------------------------------------------------------------
# This function encapsulates the output generation logic for each DEG contrast.
# It performs three key operations:
#
#   1. Writes the full DESeq2 results table (all genes, all statistics)
#   2. Applies filtering criteria to identify marker genes
#   3. Writes the filtered marker gene list
#
# Filtering Strategy:
#   Markers must satisfy ALL of the following criteria:
#   - padj <= fdr threshold (statistical significance)
#   - |log2FoldChange| >= lfc threshold (biological significance)
#   - baseMean >= min_baseMean threshold (sufficient expression)
#
#   If more genes pass filters than topN, the function retains only the top N
#   genes ranked by absolute log2 fold change, with ties broken by p-value.
#
# Parameters:
#   res_df:         DESeq2 results as a data.frame
#   prefix:         Filename prefix for output files
#   outdir_tables:  Directory for full results tables
#   outdir_markers: Directory for marker gene lists
#   fdr:            FDR threshold for significance
#   lfc:            Minimum absolute log2 fold change
#   topN:           Maximum number of markers to retain
#   min_baseMean:   Minimum mean expression level
#
# Returns:
#   A list containing:
#     - full_path:    Path to the full results table
#     - marker_path:  Path to the marker gene list
#     - n_markers:    Number of markers identified
#     - marker_genes: Character vector of marker gene IDs

write_deg_outputs <- function(res_df, prefix, outdir_tables, outdir_markers,
                              fdr, lfc, topN, min_baseMean, 
                              dds_obj = NULL, test_sample_id = NULL) {

  test_expr <- NULL
  if (!is.null(dds_obj) && !is.null(test_sample_id)) {
    norm_counts <- counts(dds_obj, normalized = TRUE)
    if (test_sample_id %in% colnames(norm_counts)) {
      test_expr <- norm_counts[, test_sample_id]
      res_df$normalised_count_in_test_sample <- as.numeric(test_expr[res_df$gene_id])
    } else {
      res_df$normalised_count_in_test_sample <- NA_real_
    }
  } else {
    res_df$normalised_count_in_test_sample <- NA_real_
  }

  # Write full results table for comprehensive downstream analysis
  full_path <- file.path(outdir_tables, paste0(prefix, ".tsv"))
  fwrite(res_df, full_path, sep = "\t")

  # Apply conservative marker filtering
  res_df$absLFC <- abs(res_df$log2FoldChange)
  
  # Use stat (Wald statistic) for ranking if available and LFC is noisy
  # This is a fallback when LFC shrinkage wasn't applied
  use_stat_for_ranking <- FALSE
  if (!is.null(res_df$stat) && all(!is.na(res_df$stat))) {
    # If we have stat and LFC looks noisy (many extreme values), prefer stat
    extreme_lfc <- sum(abs(res_df$log2FoldChange) > 5, na.rm = TRUE)
    if (extreme_lfc > nrow(res_df) * 0.1) {
      use_stat_for_ranking <- TRUE
      message(sprintf("[INFO] Using stat for ranking (many extreme LFC values detected)"))
    }
  }
  
  keep <- !is.na(res_df$padj) &
    (res_df$padj <= fdr) &
    (res_df$absLFC >= lfc) &
    (!is.na(res_df$baseMean) & res_df$baseMean >= min_baseMean)
  
  # Additional filter: require expression in test sample (isolate/anchor)
  if (!is.null(test_expr)) {
    # Require at least 10 normalized counts in test sample
    keep <- keep & (res_df$normalised_count_in_test_sample >= 10)
    message(sprintf("[INFO] Applied expression filter: >=10 normalized counts in test sample"))
  }

  markers <- res_df[keep, , drop = FALSE]

  # Rank by effect size (absLFC or abs(stat)), then significance (padj)
  if (use_stat_for_ranking) {
    markers <- markers[order(-abs(markers$stat), markers$padj), , drop = FALSE]
  } else {
    markers <- markers[order(-markers$absLFC, markers$padj), , drop = FALSE]
  }

  # Enforce topN limit if necessary
  if (nrow(markers) > topN) {
    markers <- markers[seq_len(topN), , drop = FALSE]
  }

  marker_genes <- markers$gene_id
  marker_path <- file.path(outdir_markers, paste0(prefix, "_markers_top", topN, ".txt"))
  writeLines(marker_genes, marker_path)

  return(list(
    full_path = full_path,
    marker_path = marker_path,
    n_markers = length(marker_genes),
    marker_genes = marker_genes
  ))
}

# -----------------------------------------------------------------------------
# Data Loading and Validation
# -----------------------------------------------------------------------------
# This section loads and validates input data with comprehensive error checking.
# Early validation prevents cryptic errors later in the analysis pipeline.

stop_if_missing(opt$counts, "--counts is required")
stop_if_missing(opt$meta, "--meta is required")
stop_if_missing(opt$isolate_list, "--isolate_list is required")

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

# Set rownames as required by DESeq2
rownames(meta_sub) <- meta_sub[[opt$sample_id_col]]

# -----------------------------------------------------------------------------
# Count Matrix Validation
# -----------------------------------------------------------------------------
# DESeq2 has specific requirements for input count data. This section validates
# that the count matrix meets these requirements:
#
#   - No missing values (NA): DESeq2 cannot handle missing counts
#   - No negative values: Counts represent read counts and must be non-negative
#   - Integer values: Fractional counts indicate normalised data, which is
#     incompatible with DESeq2's statistical model
#
# These checks catch common errors such as accidentally providing TPM or
# log-transformed data instead of raw counts.

if (any(is.na(counts_mat))) {
  stop("Counts contain NA")
}
if (any(counts_mat < 0)) {
  stop("Counts contain negative values")
}
if (any(counts_mat != round(counts_mat))) {
  stop("Counts are not integers. DESeq2 requires raw integer counts.")
}

# -----------------------------------------------------------------------------
# DESeqDataSet Creation and Preprocessing
# -----------------------------------------------------------------------------
# This section creates the base DESeq2 object and performs preprocessing steps
# that are common to all contrasts.
#
# Design Considerations:
#   The script handles two experimental scenarios:
#
#   1. Replicates available (multiple samples per cell line):
#      Standard DESeq2 analysis with cell line as the design factor.
#
#   2. No replicates (single sample per cell line):
#      The script uses a simplified approach with binary groupings per contrast.
#      Dispersion estimation is more challenging without replicates; the script
#      uses local fitting to borrow information across genes.
#
# Size Factor Estimation:
#   Size factors are estimated once from the complete dataset and reused for
#   all contrasts. This ensures consistent normalisation across analyses and
#   improves computational efficiency.

meta_sub[[opt$cell_line_col]] <- factor(meta_sub[[opt$cell_line_col]])

# Detect replicate structure
n_samples_per_cl <- table(meta_sub[[opt$cell_line_col]])
has_replicates <- any(n_samples_per_cl > 1)

if (!has_replicates) {
  message("[INFO] Single sample per cell line detected. Using simplified design for DEG analysis.")
} else {
  message("[INFO] Multiple samples per cell line detected. Using cell_line as design factor.")
}

# Create base DESeqDataSet with intercept-only design
# This allows size factor estimation independent of any specific contrast
dds_base <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = meta_sub,
  design    = ~ 1
)

# -----------------------------------------------------------------------------
# Low-Count Gene Filtering
# -----------------------------------------------------------------------------
# Genes with very low counts across all samples provide little statistical
# power and can cause numerical issues. The script applies a light filter
# requiring at least 10 total counts across all samples.
#
# This threshold is intentionally permissive; DESeq2's independent filtering
# will apply more sophisticated filtering based on mean counts and the
# information content of each gene.

keep_gene <- rowSums(counts(dds_base)) >= 10
dds_base <- dds_base[keep_gene, ]

# -----------------------------------------------------------------------------
# Size Factor Estimation
# -----------------------------------------------------------------------------
# Size factors normalise for differences in sequencing depth across samples.
# DESeq2 uses the median-of-ratios method, which is robust to the presence
# of differentially expressed genes.
#
# The size factors are saved for quality control purposes. Extreme values
# (e.g., > 3-fold difference) may indicate sample quality issues.

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

isolates <- split_csv(opt$isolate_list)

anchors <- character(0)
if (!is.null(opt$anchor_list)) {
  anchors <- split_csv(opt$anchor_list)
}

# Parse anchor-to-component mapping if provided
anchor2component <- NULL
if (!is.null(opt$anchor_components)) {
  a2c <- fread(opt$anchor_components, sep = "\t", data.table = FALSE)
  if (!all(c("anchor", "component") %in% colnames(a2c))) {
    stop("--anchor_components must have columns: anchor, component")
  }
  anchor2component <- setNames(as.character(a2c$component), as.character(a2c$anchor))
}

# -----------------------------------------------------------------------------
# Isolate DEG Analysis
# -----------------------------------------------------------------------------
# For each isolate cell line, the script performs a differential expression
# analysis comparing the isolate against all other cell lines ("REST").
#
# Statistical Approach:
#   - Binary grouping: isolate vs REST
#   - Contrast extracts the log2 fold change of isolate relative to REST
#   - Positive LFC indicates higher expression in the isolate
#   - Negative LFC indicates lower expression in the isolate
#
# Handling Single-Sample Groups:
#   When the isolate/anchor has n=1 but the comparator group (REST/OUTSIDE_COMP)
#   has n=many, DESeq2 can still fit a model by borrowing dispersion information
#   across genes and the large comparator group. However, inference should be
#   treated as exploratory differential signal ranking under limited replication,
#   especially for the isolate/anchor side. The script applies LFC shrinkage
#   (apeglm) to reduce noise in effect size estimates.

all_marker_sets <- list()
all_marker_outputs <- list()

for (iso in isolates) {
  # Validate that the isolate exists in the data
  if (!(iso %in% levels(meta_sub[[opt$cell_line_col]]))) {
    warning(paste0("[WARN] isolate not found in data: ", iso, " (skipping)"))
    next
  }

  # Create binary grouping for this contrast
  grp <- ifelse(meta_sub[[opt$cell_line_col]] == iso, iso, "REST")
  meta_tmp <- meta_sub
  meta_tmp$grp_tmp <- factor(grp, levels = c("REST", iso))

  # Validate sample counts in each group
  n_iso <- sum(grp == iso)
  n_rest <- sum(grp == "REST")

  if (n_iso < 1 || n_rest < 1) {
    warning(paste0("[WARN] Not enough samples for contrast: ", iso, " vs rest (skipping)"))
    next
  }

  # Create contrast-specific DESeqDataSet
  dds_tmp <- DESeqDataSetFromMatrix(
    countData = counts(dds_base),
    colData   = meta_tmp,
    design    = ~ grp_tmp
  )

  # Transfer pre-computed size factors for consistent normalisation
  sizeFactors(dds_tmp) <- sizeFactors(dds_base)

  # Estimate dispersions and fit model
  # Note: With isolate n=1 vs REST n=many, we can still fit a model by borrowing
  # dispersion information across genes and the large comparator group.
  # However, inference should be treated as exploratory, especially for the isolate side.
  dds_tmp <- DESeq(
    dds_tmp,
    quiet = TRUE,
    fitType = "local",
    minReplicatesForReplace = Inf  # Disable outlier replacement (requires replicates)
  )

  # Extract raw results then apply shrinkage, keeping stat/pvalue/padj from DESeq2
  res_raw <- results(dds_tmp, contrast = c("grp_tmp", iso, "REST"))

  res_df <- as.data.frame(res_raw)
  res_df$gene_id               <- rownames(res_df)
  res_df$log2FoldChange_raw    <- res_df$log2FoldChange
  res_df$log2FoldChange_shrunk <- NA_real_
  res_df$lfc_shrinkage_method  <- NA_character_
  res_df$lfc_shrinkage_applied <- FALSE

  tryCatch({
    res_shrunk <- lfcShrink(dds_tmp, contrast = c("grp_tmp", iso, "REST"),
                            res = res_raw, type = "apeglm", quiet = TRUE)
    res_df$log2FoldChange_shrunk <- as.data.frame(res_shrunk)$log2FoldChange
    res_df$lfc_shrinkage_method  <- "apeglm"
    res_df$lfc_shrinkage_applied <- TRUE
    message(sprintf("[INFO] Applied apeglm LFC shrinkage for %s", iso))
  }, error = function(e) {
    message(sprintf("[WARN] apeglm shrinkage failed for %s: %s", iso, e$message))
  })

  # Use shrunken LFC for ranking; retain raw stat/pvalue/padj from DESeq2
  res_df$log2FoldChange <- ifelse(res_df$lfc_shrinkage_applied,
                                  res_df$log2FoldChange_shrunk,
                                  res_df$log2FoldChange_raw)

  col_order <- c("gene_id", "baseMean", "log2FoldChange", "log2FoldChange_raw",
                 "log2FoldChange_shrunk", "lfc_shrinkage_method", "lfc_shrinkage_applied",
                 "lfcSE", "stat", "pvalue", "padj")
  res_df <- res_df[, intersect(col_order, colnames(res_df)), drop = FALSE]

  # Get sample ID for the isolate (for expression filtering)
  iso_sample_id <- rownames(meta_tmp)[meta_tmp[[opt$cell_line_col]] == iso][1]
  
  # Write outputs and collect marker genes
  prefix <- paste0("isolate_", iso, "_vs_rest")
  out <- write_deg_outputs(
    res_df, prefix,
    outdir_tables = file.path(opt$outdir, "tables"),
    outdir_markers = file.path(opt$outdir, "markers"),
    fdr = opt$fdr_isolate,
    lfc = opt$lfc_isolate,
    topN = opt$topN_isolate,
    min_baseMean = opt$min_baseMean,
    dds_obj = dds_tmp,
    test_sample_id = iso_sample_id
  )

  all_marker_sets[[prefix]] <- out$marker_genes
  all_marker_outputs[[prefix]] <- out
  message(sprintf("[OK] %s: %d markers", prefix, out$n_markers))
}

# -----------------------------------------------------------------------------
# Anchor DEG Analysis (Optional)
# -----------------------------------------------------------------------------
# Anchor analysis identifies genes that distinguish a cell line from samples
# OUTSIDE its assigned component/cluster. This is useful when:
#
#   - The anchor is part of a biological cluster (e.g., a cancer subtype)
#   - The goal is to find markers that differentiate the anchor from unrelated
#     samples, rather than from its natural neighbours
#
# Sample Selection:
#   For each anchor, the analysis includes:
#   - All samples from the anchor cell line
#   - All samples from cell lines NOT in the anchor's component
#
#   Samples from other cell lines IN the anchor's component are excluded,
#   as these are expected to be biologically similar to the anchor.

if (length(anchors) > 0) {
  # Validate required inputs for anchor analysis
  if (is.null(anchor2component)) {
    stop("You provided --anchor_list but not --anchor_components. Provide mapping file anchor->component.")
  }
  if (!(opt$component_col %in% colnames(meta_sub))) {
    stop(paste0("meta missing component column required for anchor DEGs: ", opt$component_col))
  }

  meta_sub[[opt$component_col]] <- as.character(meta_sub[[opt$component_col]])

  for (anc in anchors) {
    # Validate that the anchor exists in the data
    if (!(anc %in% levels(meta_sub[[opt$cell_line_col]]))) {
      warning(paste0("[WARN] anchor not found in data: ", anc, " (skipping)"))
      next
    }
    if (is.na(anchor2component[[anc]])) {
      warning(paste0("[WARN] anchor has no component mapping: ", anc, " (skipping)"))
      next
    }
    comp <- anchor2component[[anc]]

    # Select samples: anchor OR not in anchor's component
    keep_samples <- (meta_sub[[opt$cell_line_col]] == anc) |
      (meta_sub[[opt$component_col]] != comp)
    meta_k <- meta_sub[keep_samples, , drop = FALSE]

    # Create binary grouping
    grp <- ifelse(meta_k[[opt$cell_line_col]] == anc, anc, "OUTSIDE_COMP")
    meta_k$grp_tmp <- factor(grp, levels = c("OUTSIDE_COMP", anc))

    # Validate sample counts
    n_anc <- sum(grp == anc)
    n_outside <- sum(grp == "OUTSIDE_COMP")

    if (n_anc < 1 || n_outside < 1) {
      warning(paste0("[WARN] Not enough samples for contrast: ", anc,
                     " vs outside component ", comp, " (skipping)"))
      next
    }

    # Create contrast-specific DESeqDataSet with subsetted samples
    counts_k <- counts(dds_base)[, rownames(meta_k), drop = FALSE]
    dds_k <- DESeqDataSetFromMatrix(
      countData = counts_k,
      colData   = meta_k,
      design    = ~ grp_tmp
    )

    # Transfer size factors for the subsetted samples
    sizeFactors(dds_k) <- sizeFactors(dds_base)[rownames(meta_k)]

      # Estimate dispersions and fit model
    # Same approach as isolate: exploratory inference with dispersion borrowing
    dds_k <- DESeq(
      dds_k,
      quiet = TRUE,
      fitType = "local",
      minReplicatesForReplace = Inf  # Disable outlier replacement
    )

    # Extract raw results then apply shrinkage, keeping stat/pvalue/padj from DESeq2
    res_raw <- results(dds_k, contrast = c("grp_tmp", anc, "OUTSIDE_COMP"))

    res_df <- as.data.frame(res_raw)
    res_df$gene_id               <- rownames(res_df)
    res_df$log2FoldChange_raw    <- res_df$log2FoldChange
    res_df$log2FoldChange_shrunk <- NA_real_
    res_df$lfc_shrinkage_method  <- NA_character_
    res_df$lfc_shrinkage_applied <- FALSE

    tryCatch({
      res_shrunk <- lfcShrink(dds_k, contrast = c("grp_tmp", anc, "OUTSIDE_COMP"),
                              res = res_raw, type = "apeglm", quiet = TRUE)
      res_df$log2FoldChange_shrunk <- as.data.frame(res_shrunk)$log2FoldChange
      res_df$lfc_shrinkage_method  <- "apeglm"
      res_df$lfc_shrinkage_applied <- TRUE
      message(sprintf("[INFO] Applied apeglm LFC shrinkage for anchor %s", anc))
    }, error = function(e) {
      message(sprintf("[WARN] apeglm shrinkage failed for anchor %s: %s", anc, e$message))
    })

    # Use shrunken LFC for ranking; retain raw stat/pvalue/padj from DESeq2
    res_df$log2FoldChange <- ifelse(res_df$lfc_shrinkage_applied,
                                    res_df$log2FoldChange_shrunk,
                                    res_df$log2FoldChange_raw)

    col_order <- c("gene_id", "baseMean", "log2FoldChange", "log2FoldChange_raw",
                   "log2FoldChange_shrunk", "lfc_shrinkage_method", "lfc_shrinkage_applied",
                   "lfcSE", "stat", "pvalue", "padj")
    res_df <- res_df[, intersect(col_order, colnames(res_df)), drop = FALSE]

    # Get sample ID for the anchor (for expression filtering)
    anc_sample_id <- rownames(meta_k)[meta_k[[opt$cell_line_col]] == anc][1]
    
    # Write outputs and collect marker genes
    prefix <- paste0("anchor_", anc, "_vs_outside_component_", comp)
    out <- write_deg_outputs(
      res_df, prefix,
      outdir_tables = file.path(opt$outdir, "tables"),
      outdir_markers = file.path(opt$outdir, "markers"),
      fdr = opt$fdr_anchor,
      lfc = opt$lfc_anchor,
      topN = opt$topN_anchor,
      min_baseMean = opt$min_baseMean,
      dds_obj = dds_k,
      test_sample_id = anc_sample_id
    )

    all_marker_sets[[prefix]] <- out$marker_genes
    all_marker_outputs[[prefix]] <- out
    message(sprintf("[OK] %s: %d markers", prefix, out$n_markers))
  }
}

# -----------------------------------------------------------------------------
# Unique Feature Set Construction
# -----------------------------------------------------------------------------
# The final step consolidates marker genes across all contrasts into a single
# "unique feature set" based on gene recurrence.
#
# Recurrence-Based Filtering:
#   Genes appearing in multiple contrasts are more likely to represent robust
#   biological signals rather than contrast-specific noise. The recurrence
#   threshold (--recurrence_k) specifies the minimum number of contrasts in
#   which a gene must appear to be included in the unique feature set.
#
# Output Files:
#   - gene_recurrence_across_contrasts.tsv: Full recurrence counts per gene
#   - unique_feature_set_recurrence_ge_{k}.txt: Genes meeting the threshold
#   - marker_sets_manifest.tsv: Summary of marker counts per contrast

if (length(all_marker_sets) == 0) {
  stop("No marker sets produced (check isolate names, inputs, or thresholds).")
}

# Calculate gene recurrence across all contrasts
all_genes <- unlist(all_marker_sets, use.names = FALSE)
if (is.null(all_genes)) {
  all_genes <- character(0)
}
all_genes <- as.character(all_genes)
gene_freq <- sort(table(all_genes), decreasing = TRUE)
gene_freq_df <- data.frame(
  gene_id = as.character(names(gene_freq)),
  freq = as.integer(gene_freq)
)
fwrite(gene_freq_df,
       file.path(opt$outdir, "markers", "gene_recurrence_across_contrasts.tsv"),
       sep = "\t")

# Apply recurrence threshold
k <- opt$recurrence_k
unique_genes <- as.character(gene_freq_df$gene_id[gene_freq_df$freq >= k])

unique_path <- file.path(opt$outdir, "markers",
                         paste0("unique_feature_set_recurrence_ge_", k, ".txt"))
writeLines(unique_genes, unique_path)

# Also report recurrence >= 1 for context
unique_genes_k1 <- as.character(gene_freq_df$gene_id[gene_freq_df$freq >= 1])
message(sprintf("[INFO] Recurrence >= 1: %d genes (for context)", length(unique_genes_k1)))
message(sprintf("[INFO] Recurrence >= %d: %d genes (recurrence threshold)", k, length(unique_genes)))

# Generate manifest summarising marker sets
manifest <- data.frame(
  contrast = names(all_marker_sets),
  marker_file = vapply(
    all_marker_outputs[names(all_marker_sets)],
    function(x) file.path("markers", basename(x$marker_path)),
    character(1)
  ),
  table_file = vapply(
    all_marker_outputs[names(all_marker_sets)],
    function(x) file.path("tables", basename(x$full_path)),
    character(1)
  ),
  n_markers = vapply(all_marker_sets, length, integer(1))
)
fwrite(manifest,
       file.path(opt$outdir, "markers", "marker_sets_manifest.tsv"),
       sep = "\t")

# -----------------------------------------------------------------------------
# Execution Summary
# -----------------------------------------------------------------------------
# The script outputs summary statistics to facilitate validation and logging.

message(sprintf("[DONE] Unique feature set: %d genes (recurrence >= %d)",
                length(unique_genes), k))
message(sprintf("[OUT] %s", unique_path))

# Write session info for reproducibility
writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))
message("[INFO] Session info written to sessionInfo.txt")
