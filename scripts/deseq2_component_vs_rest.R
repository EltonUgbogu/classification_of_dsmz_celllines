#!/usr/bin/env Rscript

# =============================================================================
# deseq2_component_vs_rest.R
# =============================================================================
#
# DESCRIPTION
# -----------
# This script performs differential gene expression analysis using DESeq2 to
# identify transcriptional signatures that distinguish samples within a
# specified graph component from all other samples in the dataset. The
# component-versus-rest design is particularly suited for characterising
# the molecular programmes underlying cluster membership identified through
# unsupervised analysis.
#
#
# BIOLOGICAL CONTEXT
# ------------------
# In cancer genomics, unsupervised clustering methods (e.g., consensus
# clustering, community detection on similarity graphs) partition samples
# into groups based on transcriptomic similarity. However, cluster membership
# alone does not reveal which genes drive the observed groupings. This script
# addresses that gap by identifying differentially expressed genes (DEGs)
# between a component of interest and the remainder of the cohort.
#
# The resulting marker genes serve multiple purposes:
#
#   (1) Biological Interpretation: Pathway enrichment analysis on DEGs can
#       reveal the biological processes, molecular functions, or cellular
#       states that characterise the component.
#
#   (2) Subtype Annotation: Comparison of DEGs against known gene signatures
#       (e.g., PAM50 for breast cancer, immunophenotype signatures) can
#       link data-driven clusters to established molecular subtypes.
#
#   (3) Biomarker Discovery: Genes with large effect sizes and robust
#       statistical support may serve as candidate biomarkers for
#       classification or prognostication.
#
#   (4) Hypothesis Generation: Unexpected DEGs may suggest novel biology
#       or previously unrecognised heterogeneity within the cohort.
#
#
# STATISTICAL FRAMEWORK
# ---------------------
# The script employs DESeq2's negative binomial generalised linear model
# framework, which is well-suited for RNA-seq count data exhibiting
# overdispersion (variance exceeding the Poisson expectation).
#
# component membership (1 = component of interest, 0 = rest)
#
# Accurate dispersion estimation is critical for reliable inference. DESeq2
# employs a three-step approach:
#
#   (1) Gene-wise maximum likelihood estimation of $\alpha_i$
#   (2) Fitting a mean-dispersion trend across all genes
#   (3) Empirical Bayes shrinkage of gene-wise estimates towards the trend
#
# This shrinkage improves power for genes with few counts whilst preserving
# sensitivity for genes with strong evidence of unusual variability.
#
# The null hypothesis $H_0: \beta_{i,1} = 0$ (no differential expression)
# is tested using the Wald statistic:
#
# Resulting $p$-values are adjusted for multiple testing using the
# Benjamini-Hochberg procedure to control the false discovery rate (FDR).
#
#
# CONTRAST DESIGN RATIONALE
# -------------------------
# The component-versus-rest design treats all non-component samples as a
# single reference group. This approach:
#
#   - Maximises statistical power by using all available data
#   - Identifies genes specific to the component relative to the cohort
#     background, rather than relative to a single comparator
#   - Is appropriate when the scientific question concerns what makes the
#     component distinctive, not how it differs from any particular group
#
# Alternative designs (e.g., pairwise component comparisons) may be
# appropriate for different scientific questions but are not implemented
# in this script.
#
#
# INPUT REQUIREMENTS
# ------------------
# 1. Raw count matrix :
#    - Tab-separated file with genes as rows and samples as columns
#    - First column contains gene identifiers
#    - Values must be raw integer counts (not normalised)
#
# 2. Sample metadata :
#    - Tab-separated file with one row per sample
#    - Must contain a sample identifier column matching count matrix columns
#    - Must contain a component membership column (integer or factor)
#
# 3. Component identifier :
#    - Integer specifying which component to compare against the rest
#
#
# OUTPUT STRUCTURE
# ----------------
# The script produces:
#
#   {outdir}/
#   +-- tables/
#   |   +-- component_{N}_vs_rest.tsv    # Full DESeq2 results
#   +-- markers/
#       +-- component_{N}_vs_rest_UP_top{M}.tsv    # Up-regulated markers
#       +-- component_{N}_vs_rest_UP_top{M}.txt    # Gene list only
#       +-- component_{N}_vs_rest_DOWN_top{M}.tsv  # Down-regulated markers
#       +-- component_{N}_vs_rest_DOWN_top{M}.txt  # Gene list only
#
#
# USAGE EXAMPLE
# -------------
#   Rscript deseq2_component_vs_rest.R \
#     --counts results/unsupervised/brca/deseq2_inputs/counts.tsv \
#     --meta results/unsupervised/brca/deseq2_inputs/metadata_with_components.tsv \
#     --component 1 \
#     --component_col component \
#     --sample_id_col sample_name \
#     --outdir results/unsupervised/brca/component_vs_rest \
#     --fdr 0.05 --lfc 1.0
#
#
# DEPENDENCIES
# ------------
# R packages: optparse, DESeq2, dplyr, readr, tibble
# Bioconductor: DESeq2 (install via BiocManager::install("DESeq2"))
#
# -----------------------------------------------------------------------------
# Package Loading
# -----------------------------------------------------------------------------
# The script loads required packages with startup messages suppressed to
# maintain clean log output in batch processing environments.
#
# Package purposes:
#   - optparse:   Command-line argument parsing
#   - DESeq2:     Differential expression analysis framework
#   - dplyr:      Data manipulation and filtering
#   - readr:      Efficient file I/O with type inference
#   - tibble:     Modern data frame operations (rownames_to_column)

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
})

# -----------------------------------------------------------------------------
# Command-Line Argument Definitions
# -----------------------------------------------------------------------------
# The script accepts arguments organised into three categories:
#
# Input Specifications:
# --counts, --meta, --component
#
# Column Mappings:
# --component_col, --sample_id_col
#
# Filtering Thresholds:
# --fdr, --lfc, --topN
#
# The FDR and LFC thresholds control the stringency of marker selection.
# Typical values for exploratory analysis are FDR $< 0.05$ and
# $|\mathrm{LFC}| > 1.0$ (2-fold change). More stringent thresholds
# (e.g., FDR $< 0.01$, $|\mathrm{LFC}| > 1.5$) may be appropriate for
# generating high-confidence marker lists.

opt_list <- list(
  make_option("--counts", type = "character",
              help = "Path to raw counts TSV"),
  make_option("--meta", type = "character",
              help = "Path to metadata TSV with component column"),
  make_option("--component", type = "integer",
              help = "Component ID to analyse"),
  make_option("--component_col", type = "character", default = "component",
              help = "Component column name"),
  make_option("--sample_id_col", type = "character", default = "sample_name",
              help = "Sample ID column name"),
  make_option("--outdir", type = "character",
              help = "Output directory"),
  make_option("--fdr", type = "numeric", default = 0.05,
              help = "FDR threshold"),
  make_option("--lfc", type = "numeric", default = 1.0,
              help = "Log2FC threshold"),
  make_option("--topN", type = "integer", default = 500,
              help = "Top N markers to extract")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# -----------------------------------------------------------------------------
# Input Validation
# -----------------------------------------------------------------------------
# Early validation of required arguments prevents cryptic errors later in
# the pipeline and provides informative feedback to the user.

if (is.null(opt$counts) || is.null(opt$meta) ||
    is.null(opt$component) || is.null(opt$outdir)) {
  stop("--counts, --meta, --component, and --outdir are required")
}

# -----------------------------------------------------------------------------
# Output Directory Initialisation
# -----------------------------------------------------------------------------
# The script creates a structured output directory hierarchy:
#
#   - tables/   : Full DESeq2 results for comprehensive downstream analysis
#   - markers/  : Filtered gene lists for pathway enrichment and annotation
#
# Using showWarnings = FALSE enables idempotent execution.

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "markers"), recursive = TRUE, showWarnings = FALSE)

cat("[INFO] Loading data...\n")

# -----------------------------------------------------------------------------
# Count Matrix Loading
# -----------------------------------------------------------------------------
# The count matrix is loaded and converted to the matrix format required by
# DESeq2. Gene identifiers from the first column are assigned as row names.
#
# IMPORTANT: DESeq2 requires raw integer counts. The script assumes that
# the input file contains unnormalised counts; providing TPM, FPKM, or
# other normalised values will produce statistically invalid results.

counts_df <- read_tsv(opt$counts, show_col_types = FALSE)

# Extract gene identifiers from the first column
gene_ids <- counts_df[[1]]

# Convert remaining columns to numeric matrix
counts_mat <- as.matrix(counts_df[, -1])
rownames(counts_mat) <- gene_ids
# TEMPORARY COMPROMISE: current staged "counts" are VST-like decimals because
# raw DSMZ counts are unavailable. Coerce to integer only to let DESeq2 run;
# replace this with true raw integer counts when available.
if (any(counts_mat != round(counts_mat), na.rm = TRUE)) {
  warning("Coercing non-integer staged expression values to integer for DESeq2 compatibility; use raw counts for final analysis.")
}
storage.mode(counts_mat) <- "integer"

# -----------------------------------------------------------------------------
# Metadata Loading and Validation
# -----------------------------------------------------------------------------
# The metadata file must contain:
#   (1) A sample identifier column matching the count matrix column names
#   (2) A component membership column indicating cluster assignment
#
# The script validates that the specified component column exists before
# proceeding to avoid downstream errors.

meta_df <- read_tsv(opt$meta, show_col_types = FALSE)

if (!opt$component_col %in% colnames(meta_df)) {
  stop(paste0("Component column '", opt$component_col, "' not found in metadata"))
}

# -----------------------------------------------------------------------------
# Binary Group Variable Construction
# -----------------------------------------------------------------------------
# The contrast design requires a binary group variable distinguishing the
# component of interest from all other samples. This transformation converts
# the multi-level component variable into a two-level factor suitable for
# the DESeq2 design formula.
#
# The naming convention "Component_N" vs "Rest" produces interpretable
# contrast names in the output and avoids issues with numeric factor levels.

meta_df$group <- ifelse(
  meta_df[[opt$component_col]] == opt$component,
  paste0("Component_", opt$component),
  "Rest"
)

# -----------------------------------------------------------------------------
# Sample Alignment
# -----------------------------------------------------------------------------
# DESeq2 requires exact correspondence between metadata rows and count matrix
# columns. This section ensures proper alignment by:
#
#   (1) Filtering metadata to samples present in the count matrix
#   (2) Subsetting the count matrix to samples with metadata
#
# This bidirectional filtering handles cases where the count matrix and
# metadata have different sample sets (e.g., due to quality control
# exclusions applied to one but not the other).

meta_df <- meta_df[meta_df[[opt$sample_id_col]] %in% colnames(counts_mat), ]
counts_mat <- counts_mat[, meta_df[[opt$sample_id_col]], drop = FALSE]

# -----------------------------------------------------------------------------
# Sample Count Reporting
# -----------------------------------------------------------------------------
# The script reports sample counts for each group to facilitate validation
# and interpretation. Imbalanced group sizes may affect statistical power
# and should be considered when interpreting results.

cat(sprintf("[INFO] Component %d samples: %d\n",
            opt$component,
            sum(meta_df$group == paste0("Component_", opt$component))))
cat(sprintf("[INFO] Rest samples: %d\n",
            sum(meta_df$group == "Rest")))
cat(sprintf("[INFO] Total genes: %d\n", nrow(counts_mat)))

# -----------------------------------------------------------------------------
# Sample Sufficiency Validation
# -----------------------------------------------------------------------------
# DESeq2 requires at least one sample in each group. Whilst the method can
# technically run with minimal sample sizes, results from severely imbalanced
# or underpowered comparisons should be interpreted with caution.
#
# The script enforces a minimum of one sample per group; users requiring
# stricter thresholds should implement additional validation.

n_comp <- sum(meta_df$group == paste0("Component_", opt$component))
n_rest <- sum(meta_df$group == "Rest")

if (n_comp < 1) {
  stop(sprintf("No samples found for component %d", opt$component))
}
if (n_rest < 1) {
  stop("No samples found for rest group")
}

# -----------------------------------------------------------------------------
# DESeqDataSet Construction
# -----------------------------------------------------------------------------
# The DESeqDataSet object encapsulates the count matrix, sample metadata,
# and experimental design in a single container. The design formula
# ~ group specifies that differential expression
# will be tested with respect to the binary group variable.
#
# DESeq2 automatically sets factor levels alphabetically, placing
# "Component_N" before "Rest". The explicit contrast specification in
# the results extraction step ensures correct directionality regardless
# of factor ordering.

cat("[INFO] Creating DESeqDataSet...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = meta_df,
  design = ~ group
)

# -----------------------------------------------------------------------------
# Low-Count Gene Filtering
# -----------------------------------------------------------------------------
# Genes with very low counts across all samples provide little statistical
# information and can cause numerical instability in dispersion estimation.
# A threshold of 10 total counts is applied:
#
#     keep gene i iff sum_{j=1}^{n} K_{ij} >= 10
#
# This threshold is intentionally permissive; DESeq2's independent filtering
# applies additional, more sophisticated filtering based on mean expression
# and information content.

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat(sprintf("[INFO] After filtering: %d genes\n", nrow(dds)))

# -----------------------------------------------------------------------------
# DESeq2 Analysis Execution
# -----------------------------------------------------------------------------
# The DESeq() function executes the complete analysis pipeline:
#
#   (1) Size factor estimation (median-of-ratios method)
#   (2) Dispersion estimation (gene-wise, trend, and shrinkage)
#   (3) Negative binomial GLM fitting
#   (4) Wald test for significance
#
# Parameter Choices:
# fitType = "local" Uses local regression for the
# mean-dispersion trend, which is more flexible than the default
# parametric fit and accommodates non-standard dispersion patterns.
# minReplicatesForReplace = Inf Disables automatic
# outlier replacement, preserving the original count values. This is appropriate
# when outliers may represent genuine biological variation rather than technical
# artefacts.

cat("[INFO] Running DESeq2...\n")
dds <- DESeq(dds, quiet = TRUE, fitType = "local", minReplicatesForReplace = Inf)

# -----------------------------------------------------------------------------
# Results Extraction
# -----------------------------------------------------------------------------
# The results() function extracts the Wald test results for the
# specified contrast. The contrast is explicitly defined as:
#
#   Component_N vs Rest (positive LFC = higher in component)
#
#Output Columns:  
# baseMean Mean normalised count across all samples
# log2FoldChange Estimated LFC (component vs rest)
# lfcSE Standard error of the LFC estimate
# stat Wald statistic (hat{beta} / SE)
# pvalue Wald test p-value
# padj Benjamini-Hochberg adjusted p-value (FDR)

contrast_name <- paste0("Component_", opt$component, "_vs_Rest")
res <- results(
  dds,
  contrast = c("group", paste0("Component_", opt$component), "Rest"),
  alpha = opt$fdr
)

# -----------------------------------------------------------------------------
# Results Formatting
# -----------------------------------------------------------------------------
# The results are converted to a data frame with gene identifiers as a
# column (rather than row names) for compatibility with downstream tools
# and file formats. Sorting by adjusted p-value facilitates rapid
# identification of the most significant genes.

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj, pvalue)

norm_counts <- counts(dds, normalized = TRUE)
component_sample_ids <- meta_df[[opt$sample_id_col]][
  meta_df$group == paste0("Component_", opt$component)
]
if (length(component_sample_ids) > 0) {
  component_expr <- rowMeans(
    norm_counts[, component_sample_ids, drop = FALSE],
    na.rm = TRUE
  )
  res_df$normalised_count_in_test_sample <- as.numeric(component_expr[res_df$gene_id])
} else {
  res_df$normalised_count_in_test_sample <- NA_real_
}

# -----------------------------------------------------------------------------
# Full Results Export
# -----------------------------------------------------------------------------
# The complete results table is written for archival and comprehensive
# downstream analysis. This file contains all tested genes regardless of
# significance, enabling:
#
#   - Custom re-filtering with different thresholds
#   - Gene set enrichment analysis (GSEA) using ranked gene lists
#   - Volcano plot generation
#   - Integration with other datasets

out_file <- file.path(opt$outdir, "tables",
                      paste0("component_", opt$component, "_vs_rest.tsv"))
write_tsv(res_df, out_file)
cat(sprintf("[OK] Full results written to: %s\n", out_file))

# -----------------------------------------------------------------------------
# Marker Gene Filtering
# -----------------------------------------------------------------------------
# Marker genes are identified by applying dual thresholds for statistical
# significance (FDR) and biological effect size (LFC). This two-criterion
# approach balances:
#
#   - Statistical rigour: FDR control limits false discoveries
#   - Biological relevance: LFC threshold excludes statistically significant
#     but biologically trivial changes
#
# Markers are separated into up-regulated (positive LFC) and down-regulated
# (negative LFC) sets to facilitate interpretation. Up-regulated genes are
# more highly expressed in the component of interest; down-regulated genes
# are more highly expressed in the rest of the cohort.
#
# Filtering Criteria:
# UP marker iff p_adj < tau_FDR and LFC > tau_LFC
# DOWN marker iff p_adj < tau_FDR and LFC < -tau_LFC

markers_up <- res_df %>%
  filter(!is.na(padj), padj < opt$fdr, log2FoldChange > opt$lfc,
         !is.na(normalised_count_in_test_sample),
         normalised_count_in_test_sample >= 10) %>%
  arrange(desc(log2FoldChange))

markers_down <- res_df %>%
  filter(!is.na(padj), padj < opt$fdr, log2FoldChange < -opt$lfc,
         !is.na(normalised_count_in_test_sample),
         normalised_count_in_test_sample >= 10) %>%
  arrange(log2FoldChange)

cat(sprintf("[INFO] UP markers: %d\n", nrow(markers_up)))
cat(sprintf("[INFO] DOWN markers: %d\n", nrow(markers_down)))

# -----------------------------------------------------------------------------
# Marker Gene Export
# -----------------------------------------------------------------------------
# Filtered marker genes are exported in two formats:
#
#   (1) TSV files containing full statistics (LFC, p-value, FDR) for
#       detailed analysis and reporting
#
#   (2) TXT files containing gene identifiers only, suitable for direct
#       input to pathway enrichment tools (e.g., Enrichr, g:Profiler, DAVID)
#
# The topN parameter limits the number of exported markers to
# prevent unwieldy gene lists. Genes are ranked by effect size (|LFC|)
# to prioritise the strongest differential expression signals.

if (nrow(markers_up) > 0) {
  markers_up_top <- head(markers_up, opt$topN)

  # Write full statistics
  out_up <- file.path(opt$outdir, "markers",
                      paste0("component_", opt$component,
                             "_vs_rest_UP_top", opt$topN, ".tsv"))
  write_tsv(markers_up_top, out_up)
  cat(sprintf("[OK] UP markers written to: %s\n", out_up))

  # Write gene list only
  out_up_genes <- file.path(opt$outdir, "markers",
                            paste0("component_", opt$component,
                                   "_vs_rest_UP_top", opt$topN, ".txt"))
  write_lines(markers_up_top$gene_id, out_up_genes)
}

if (nrow(markers_down) > 0) {
  markers_down_top <- head(markers_down, opt$topN)

  # Write full statistics
  out_down <- file.path(opt$outdir, "markers",
                        paste0("component_", opt$component,
                               "_vs_rest_DOWN_top", opt$topN, ".tsv"))
  write_tsv(markers_down_top, out_down)
  cat(sprintf("[OK] DOWN markers written to: %s\n", out_down))

  # Write gene list only
  out_down_genes <- file.path(opt$outdir, "markers",
                              paste0("component_", opt$component,
                                     "_vs_rest_DOWN_top", opt$topN, ".txt"))
  write_lines(markers_down_top$gene_id, out_down_genes)
}

# -----------------------------------------------------------------------------
# Execution Summary
# -----------------------------------------------------------------------------
# The script outputs a summary of key statistics to facilitate validation
# and provide a quick overview of results. This summary is particularly
# useful when the script is executed as part of a larger pipeline where
# detailed log review may be impractical.

cat("\n[SUMMARY]\n")
cat(sprintf("  Component %d samples: %d\n", opt$component, n_comp))
cat(sprintf("  Rest samples: %d\n", n_rest))
cat(sprintf("  UP markers (FDR<%.2f, LFC>%.2f): %d\n",
            opt$fdr, opt$lfc, nrow(markers_up)))
cat(sprintf("  DOWN markers (FDR<%.2f, LFC<%.2f): %d\n",
            opt$fdr, -opt$lfc, nrow(markers_down)))
cat("\n[OK] Analysis complete!\n")
