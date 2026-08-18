#!/usr/bin/env Rscript

# =============================================================================
# inspect_cellline_mapping.R
# =============================================================================
#
# DESCRIPTION
# -----------
# This script performs comprehensive inspection of cell-line-to-tumour
# mapping results, serving as the reciprocal analysis to inspect_tumour_mapping.R.
# While the tumour-centric script asks "which cell line best represents this
# tumour?", this script asks "which tumours does this cell line most closely
# resemble?". Both perspectives are essential for validating cell line models.
#
#
# BIOLOGICAL CONTEXT: CELL LINE MODEL VALIDATION
# ----------------------------------------------
# Cell lines are immortalised cultures derived from primary tumours that serve
# as tractable experimental models for cancer research. However, cell lines
# undergo substantial molecular evolution during establishment and long-term
# passage, raising questions about their fidelity to the originating tumour
# type. This script addresses the fundamental question: does each cell line
# retain transcriptomic features characteristic of its annotated cancer type?
#
# The Reciprocal Mapping Perspective:
# Tumour-to-cell-line mapping identifies which cell lines best model specific
# tumours, useful for selecting models for patient-derived xenograft validation
# or drug sensitivity prediction. Cell-line-to-tumour mapping, conversely,
# validates whether cell lines maintain lineage identity by asking whether
# their expression profiles most closely resemble tumours of the expected type.
#
# A cell line that maps to tumours of a different lineage may indicate:
#   - Misannotation in cell line databases (a recognised problem in the field)
#   - Substantial transcriptomic drift during long-term culture
#   - Cross-contamination with cells from another lineage
#   - Genuine biological ambiguity (e.g., tumours with mixed lineage features)
#
# Key Differences from Tumour-Centric Analysis:
#   - Unit of analysis: cell line (not tumour)
#   - Correctness criterion: best tumour lineage matches cell line lineage
#   - Confidence delta: score separation between rank-1 and rank-2 tumours
#   - Top-k consistency: composition of top-10 tumours for each cell line
#
#
# INSPECTION ANALYSES PERFORMED
# -----------------------------
# Incorrect Mapping Analysis:
# Identifies cell lines whose best-matching tumours belong to a different
# cancer type. These cases may reveal cell line misannotation, contamination,
# or significant transcriptomic drift. The analysis examines the lineage of
# incorrectly matched tumours and the presence of correct-lineage tumours
# in the top-k rankings.
#
# Confidence Distribution Analysis:
# Quantifies mapping certainty using confidence delta (score separation) and
# same-lineage fraction metrics. Per-lineage statistics reveal whether certain
# cell line collections have systematically lower confidence, potentially
# indicating heterogeneous derivation conditions or annotation quality issues.
#
# Low Confidence Case Identification:
# Flags cell lines with confidence delta below a configurable threshold
# (default: 0.03) for manual review. Low confidence may indicate multiple
# similarly-matching tumour types or atypical expression profiles.
#
# Component-Level Mapping (optional):
# When tumour graph components are available, assesses whether cell lines
# map coherently to biologically meaningful tumour clusters.
#
# Top-k Lineage Consistency:
# Evaluates how consistently the top-k tumour matches belong to the expected
# cancer type. Cell lines with mixed-lineage top-10 rankings warrant closer
# examination for potential annotation errors or biological ambiguity.
#
# Score Distribution Analysis:
# Analyses best-match score distributions across cell line lineages,
# revealing whether certain cancer types have systematically better or
# worse tumour representation in the reference cohort.
#
#
# INPUT REQUIREMENTS
# ------------------
# 1. Cell line mapping summary (--summary):
#    Per-cell-line best match, scores, and confidence metrics with columns:
#    cell_line, cell_lineage, best_tumour, best_tumour_lineage, best_score,
#    confidence_delta, confidence_frac_same_lineage, is_correct
#
# 2. Full rankings (--rankings):
#    Top-k tumours per cell line with columns:
#    cell_line, tumour, rank, tumour_lineage
#
# 3. Tumour components (--components, optional):
#    Tumour-to-component mappings for component-level analysis
#
#
# OUTPUT FILES
# ------------
#   incorrect_mappings_analysis.tsv
#       Detailed analysis of each incorrect mapping including top-k breakdown.
#
#   confidence_statistics.tsv
#       Per-lineage summary statistics for confidence metrics.
#
#   low_confidence_cases.tsv
#       Cell lines with confidence delta below threshold, flagged for review.
#
#   component_usage.tsv (if components provided)
#       Summary of cell line mappings per tumour graph component.
#
#   component_by_lineage.tsv (if components provided)
#       Cross-tabulation of cell line lineage and component assignment.
#
#   topk_lineage_consistency.tsv
#       Per-cell-line assessment of lineage consistency in top-k rankings.
#
#   score_statistics.tsv
#       Per-lineage summary of best-match score distributions.
#
#   INSPECTION_SUMMARY.md
#       Human-readable report summarising all inspection findings.
#
#
# USAGE EXAMPLE
# -------------
#   Rscript inspect_cellline_mapping.R \
#     --summary cellline_mapping_summary.tsv \
#     --rankings cellline_to_tumour_rankings.tsv \
#     --components tumour_components.tsv \
#     --output-dir inspection_cellline \
#     --top-k 10 \
#     --low-conf-thr 0.03
#
#
# DEPENDENCIES
# ------------
# R packages: data.table, ggplot2, optparse
#
# =============================================================================


# =============================================================================
# PACKAGE LOADING
# =============================================================================
# The script requires three R packages for its core functionality:
#
#   data.table: Provides efficient data manipulation through reference
#               semantics and optimised grouping operations. Essential for
#               handling large-scale mapping results where cell lines are
#               compared against thousands of tumour samples. Key functions
#               include fread() for fast file reading, := for in-place
#               column assignment, and by= for grouped aggregations.
#
#   ggplot2:    The grammar of graphics implementation in R. Loaded for
#               potential inspection visualisations, though the current
#               implementation focuses on tabular outputs.
#
#   optparse:   Provides GNU-style command-line argument parsing with
#               support for default values, type checking, and automatic
#               help generation. Enables flexible parameterisation of
#               analysis thresholds and file paths.
#
# The suppressPackageStartupMessages() wrapper prevents verbose loading
# messages that would clutter pipeline logs.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})


# =============================================================================
# COMMAND-LINE ARGUMENT DEFINITIONS
# =============================================================================
# Command-line arguments enable flexible script execution with configurable
# thresholds and input/output paths. This design facilitates both interactive
# use and integration into automated analysis pipelines.
#
# Required Arguments:
#   --summary:    Path to cell line mapping summary file containing per-cell-
#                 line best matches and confidence metrics
#   --rankings:   Path to complete rankings file with top-k tumours per cell
#                 line
#   --output-dir: Directory for all output files (created if non-existent)
#
# Optional Arguments:
#   --components:   Path to tumour component assignments. Unlike the tumour
#                   inspection script which uses cell line components, this
#                   script optionally uses tumour components to assess whether
#                   cell lines map to coherent tumour clusters.
#   --top-k:        Number of top-ranked tumours to consider for consistency
#                   analysis (default: 10). Higher values provide more robust
#                   consistency estimates but may include less relevant matches.
#   --low-conf-thr: Threshold for flagging low-confidence mappings based on
#                   confidence delta (default: 0.03). This threshold should
#                   be calibrated based on the score distribution of the
#                   specific dataset.
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(c("--summary"), type = "character", default = NULL,
              help = "Path to cell line mapping summary TSV"),
  make_option(c("--rankings"), type = "character", default = NULL,
              help = "Path to cell line-to-tumour rankings TSV"),
  make_option(c("--components"), type = "character", default = NULL,
              help = "Optional: tumour components TSV (tumour -> component)"),
  make_option(c("--output-dir"), type = "character", default = NULL,
              help = "Output directory for inspection results"),
  make_option(c("--top-k"), type = "integer", default = 10,
              help = "Top-k to inspect (default: 10)"),
  make_option(c("--low-conf-thr"), type = "double", default = 0.03,
              help = "Low-confidence threshold on confidence_delta (default: 0.03)")
)

# -----------------------------------------------------------------------------
# Parse arguments from command line
# The OptionParser object processes command-line inputs according to the
# defined option_list specification, returning a named list of argument values.
# -----------------------------------------------------------------------------

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# =============================================================================
# INPUT VALIDATION
# =============================================================================
# Validates that all required arguments are provided before proceeding with
# analysis. Early termination with informative error messages prevents
# confusing downstream failures from missing inputs.
#
# The output directory is created with recursive = TRUE to handle nested
# directory structures, and showWarnings = FALSE to suppress messages when
# the directory already exists.
# -----------------------------------------------------------------------------

if (is.null(opt$summary) || is.null(opt$rankings) || is.null(opt$`output-dir`)) {
  stop("--summary, --rankings, and --output-dir are required")
}
dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("CELLLINE -> TUMOUR MAPPING INSPECTION\n")
cat("========================================\n\n")


# =============================================================================
# STEP 1: LOAD INPUT DATA
# =============================================================================
# This step loads the pre-computed mapping results required for inspection
# analysis. The data.table fread() function provides fast reading of
# tab-separated files with automatic type inference and memory-efficient
# storage.
#
# Summary Data Structure:
# The summary file contains one row per cell line with columns:
#   - cell_line:                   Unique identifier for each cell line
#   - cell_lineage:                Annotated cancer type of the cell line
#   - best_tumour:                 Identifier of the top-ranked tumour match
#   - best_tumour_lineage:         Cancer type of the best-matched tumour
#   - best_score:                  Similarity score for the top-1 match
#   - is_correct:                  Boolean indicating lineage-concordant mapping
#   - confidence_delta:            Score difference between rank-1 and rank-2
#   - confidence_frac_same_lineage: Fraction of top-k tumours from same lineage
#
# Rankings Data Structure:
# The rankings file contains multiple rows per cell line (one per ranked
# tumour) with columns:
#   - cell_line:      Cell line identifier matching summary$cell_line
#   - tumour:         Tumour sample identifier
#   - rank:           Position in the ranking (1 = best match)
#   - tumour_lineage: Cancer type annotation of the tumour
#
# Components Data Structure (optional):
# When provided, contains tumour to graph component mappings:
#   - tumour (or sample): Tumour identifier
#   - component:          Integer component assignment from graph clustering
#
# Column Validation:
# The script performs explicit validation of required columns to provide
# informative error messages rather than cryptic downstream failures.
# This is particularly important for pipeline integration where input
# files may be generated by different upstream processes.
# -----------------------------------------------------------------------------

cat("[1] Loading data...\n")

# Load the per-cell-line mapping summary using fast file reading
summary <- fread(opt$summary)

# Load the complete rankings with all top-k matches per cell line
rankings <- fread(opt$rankings)

# Conditionally load tumour component assignments if file path is provided
# Note: Unlike tumour inspection which uses cell line components, this script
# uses tumour components to assess which tumour clusters cell lines map to
if (!is.null(opt$components) && file.exists(opt$components)) {
  components <- fread(opt$components)
  cat("  Loaded components\n")
} else {
  components <- NULL
  cat("  No components file provided\n")
}

# Report data dimensions for verification
cat("  Summary: ", nrow(summary), " cell lines\n", sep = "")
cat("  Rankings: ", nrow(rankings), " entries\n", sep = "")

# -----------------------------------------------------------------------------
# Column Validation for Summary File
# Verifies that all required columns are present in the summary file.
# Missing columns would cause cryptic errors in downstream analysis steps,
# so explicit validation with informative error messages is preferred.
# -----------------------------------------------------------------------------

req_sum <- c("cell_line", "cell_lineage", "best_tumour", "best_tumour_lineage",
             "best_score", "confidence_delta", "confidence_frac_same_lineage", "is_correct")
missing_sum <- setdiff(req_sum, names(summary))
if (length(missing_sum) > 0) {
  stop("Summary is missing columns: ", paste(missing_sum, collapse = ", "))
}

# -----------------------------------------------------------------------------
# Column Validation for Rankings File
# The rankings file requires fewer columns as it primarily provides the
# ordered list of tumour matches for each cell line.
# -----------------------------------------------------------------------------

req_rank <- c("cell_line", "tumour", "rank", "tumour_lineage")
missing_rank <- setdiff(req_rank, names(rankings))
if (length(missing_rank) > 0) {
  stop("Rankings is missing columns: ", paste(missing_rank, collapse = ", "))
}

# -----------------------------------------------------------------------------
# Data Quality Checks
# Missing lineage annotations can compromise the accuracy of consistency
# statistics. Warnings alert the user to potential data quality issues
# without terminating execution, allowing the analysis to proceed with
# available data.
# -----------------------------------------------------------------------------

# Check for missing cell line lineage annotations
if (any(is.na(summary$cell_lineage))) {
  warning("Summary has NA cell_lineage for some cell lines.")
}

# Check for missing tumour lineage annotations in rankings
if (any(is.na(rankings$tumour_lineage))) {
  warning("Rankings has NA tumour_lineage for some entries. Top-k stats may be affected.")
}

# -----------------------------------------------------------------------------
# Column Name Harmonisation
# The rankings file may use either 'cell_lineage' or 'cl_lineage' to denote
# the cell line's cancer type. This section standardises to 'cl_lineage'
# and adds the column by merging with summary data if not present.
#
# This harmonisation is necessary because:
#   - Different upstream scripts may use different naming conventions
#   - The top-k consistency analysis requires comparing tumour_lineage to
#     the cell line's lineage for each ranking entry
# -----------------------------------------------------------------------------

# Rename cell_lineage to cl_lineage if present (standardise naming)
if ("cell_lineage" %in% names(rankings) && !"cl_lineage" %in% names(rankings)) {
  setnames(rankings, "cell_lineage", "cl_lineage")
}

# If cl_lineage is still missing, merge from summary table
# This ensures every ranking entry has the cell line's lineage for comparison
if (!"cl_lineage" %in% names(rankings)) {
  rankings <- merge(rankings,
                    unique(summary[, .(cell_line, cl_lineage = cell_lineage)]),
                    by = "cell_line", all.x = TRUE)
}


# =============================================================================
# STEP 2: INCORRECT MAPPING ANALYSIS
# =============================================================================
# Incorrect mappings occur when a cell line's best-matching tumour originates
# from a different cancer type than the cell line's annotated lineage. These
# cases are particularly important for cell line validation as they may
# indicate fundamental problems with the cell line model.
#
# Biological and Technical Causes of Incorrect Mappings:
#
# 1. Cell Line Misannotation:
#    Cell line databases contain documented annotation errors where cell lines
#    are attributed to incorrect cancer types. The ICLAC (International Cell
#    Line Authentication Committee) maintains a register of misidentified
#    cell lines, and transcriptomic mapping can help identify additional cases.
#
# 2. Cross-Contamination:
#    Cell line cross-contamination is a well-documented problem in research
#    laboratories. A cell line contaminated with faster-growing cells from
#    another lineage will exhibit the transcriptomic profile of the
#    contaminant rather than the original cells.
#
# 3. Transcriptomic Drift:
#    Long-term culture can induce substantial changes in gene expression
#    profiles. Cell lines may lose lineage-specific features and acquire
#    more generic proliferative or stress-response signatures.
#
# 4. Biological Ambiguity:
#    Some tumours have genuinely ambiguous lineage features, particularly
#    poorly differentiated cancers or those arising at anatomical boundaries
#    between tissue types.
#
# Interpretation:
# The analysis examines whether correct-lineage tumours appear in the top-k
# rankings. If correct-lineage tumours are present but not ranked first,
# the mapping error may be marginal. If correct-lineage tumours are entirely
# absent from the top-k, this suggests a more fundamental mismatch.
# -----------------------------------------------------------------------------

cat("\n[2] Analysing incorrect mappings...\n")

# Subset to only incorrectly mapped cell lines using data.table filtering
# The is_correct column indicates whether the cell line and best tumour
# share the same cancer type annotation
incorrect <- summary[is_correct == FALSE]
cat("  Incorrect mappings: ", nrow(incorrect), "\n", sep = "")

# Extract key inspection fields for each incorrect mapping
# This creates a focused data.table with columns most relevant for
# understanding why each mapping failed
incorrect_analysis <- incorrect[, .(
  cell_line,
  cell_lineage,
  best_tumour,
  best_tumour_lineage,
  best_score,
  confidence_delta,
  confidence_frac_same_lineage
)]

# Examine top-k ranking composition for each incorrect mapping
# This reveals whether correct-lineage tumours were close alternatives
for (i in seq_len(nrow(incorrect_analysis))) {
  # Get the cell line identifier and expected lineage
  cl_id <- incorrect_analysis$cell_line[i]
  cl_lineage <- incorrect_analysis$cell_lineage[i]

  # Extract the top-k ranked tumours for this cell line
  # Uses the configurable top-k parameter (default: 10)
  cl_ranks <- rankings[cell_line == cl_id & rank <= opt$`top-k`]

  # Handle cases where rankings may be incomplete
  if (nrow(cl_ranks) == 0) {
    incorrect_analysis[i, n_correct_in_top10 := NA_integer_]
    incorrect_analysis[i, top10_lineages := NA_character_]
  } else {
    # Count tumours matching the cell line's expected lineage
    n_correct <- sum(cl_ranks$tumour_lineage == cl_lineage, na.rm = TRUE)
    incorrect_analysis[i, n_correct_in_top10 := n_correct]
    
    # Record all lineages present in top-k for pattern analysis
    incorrect_analysis[i, top10_lineages := paste(cl_ranks$tumour_lineage, collapse = ";")]
  }
}

# Write the incorrect mapping analysis to output file
incorrect_file <- file.path(opt$`output-dir`, "incorrect_mappings_analysis.tsv")
fwrite(incorrect_analysis, file = incorrect_file, sep = "\t")
cat("  Saved to:", incorrect_file, "\n")

# Display cross-tabulation of cell line lineage versus best tumour lineage
# Rows represent cell line lineages, columns represent tumour lineages
# This reveals systematic patterns in incorrect mappings
cat("\n  Breakdown by cell-line lineage:\n")
print(table(incorrect$cell_lineage, incorrect$best_tumour_lineage))


# =============================================================================
# STEP 3: CONFIDENCE DISTRIBUTION ANALYSIS
# =============================================================================
# Confidence metrics quantify the certainty of each cell-line-to-tumour
# mapping. Systematic analysis across cell line lineages reveals whether
# certain cancer types have inherently lower mapping confidence.
#
# Confidence Metrics Explained:
#
#   confidence_delta: The difference in similarity score between the rank-1
#   and rank-2 tumour matches. This measures how distinctly the best match
#   stands out from alternatives.
#     - High delta (>0.10): Strong preference for top match; the cell line
#       clearly resembles one tumour more than others
#     - Medium delta (0.03-0.10): Moderate confidence; top match is preferred
#       but alternatives are relatively close
#     - Low delta (<0.03): Multiple similarly-scored tumours; the cell line
#       may represent a common transcriptional state shared by many tumours
#
#   confidence_frac_same_lineage: The fraction of top-k tumours that share
#   the same cancer type as the cell line. This measures lineage-specific
#   affinity.
#     - High fraction (>0.8): Cell line profile is typical for its lineage
#     - Medium fraction (0.4-0.8): Some cross-lineage similarity exists
#     - Low fraction (<0.4): Cell line profile is atypical; may indicate
#       drift, contamination, or misannotation
#
# Per-Lineage Interpretation:
# Cell line collections with systematically lower confidence may indicate:
#   - Heterogeneous derivation conditions across laboratories
#   - Variable passage numbers affecting expression profiles
#   - Incomplete tumour representation for that cancer type
#   - Greater molecular diversity within the lineage
# -----------------------------------------------------------------------------

cat("\n[3] Analysing confidence distributions...\n")

# Compute per-lineage summary statistics for confidence metrics
# The by = .(lineage = cell_lineage) groups by cell line cancer type
conf_summary <- summary[, .(
  # Sample size for this lineage
  n = .N,
  
  # Central tendency measures for confidence delta
  mean_delta = mean(confidence_delta, na.rm = TRUE),
  median_delta = median(confidence_delta, na.rm = TRUE),
  
  # Interquartile range bounds capturing central 50% of distribution
  q25_delta = quantile(confidence_delta, 0.25, na.rm = TRUE),
  q75_delta = quantile(confidence_delta, 0.75, na.rm = TRUE),
  
  # Mean same-lineage fraction across cell lines of this type
  mean_frac_same = mean(confidence_frac_same_lineage, na.rm = TRUE),
  
  # Count of low-confidence cases using configurable threshold
  n_low_conf = sum(confidence_delta < opt$`low-conf-thr`, na.rm = TRUE)
), by = .(lineage = cell_lineage)]

# Write confidence statistics to output file
conf_file <- file.path(opt$`output-dir`, "confidence_statistics.tsv")
fwrite(conf_summary, file = conf_file, sep = "\t")
cat("  Saved to:", conf_file, "\n")


# =============================================================================
# STEP 4: LOW CONFIDENCE CASE IDENTIFICATION
# =============================================================================
# This step identifies cell lines with confidence delta below a configurable
# threshold, flagging them for manual review. Low-confidence mappings warrant
# investigation as they may indicate problematic cell line models.
#
# Reasons for Low Confidence in Cell Line Mappings:
#
# 1. Generic Proliferative Signatures:
#    Cell lines adapted to culture conditions may lose lineage-specific
#    gene expression and acquire generic proliferative signatures. These
#    "pan-cancer" profiles match many tumour types equally well.
#
# 2. Dedifferentiation:
#    Long-term culture can cause loss of differentiation markers that
#    distinguish cancer types. Dedifferentiated cell lines may resemble
#    stem-like or mesenchymal states common across multiple lineages.
#
# 3. Incomplete Tumour Reference:
#    If the tumour reference cohort lacks adequate representation of the
#    cell line's specific molecular subtype, no tumour will match well,
#    resulting in low confidence scores across all candidates.
#
# 4. Technical Batch Effects:
#    Systematic differences between cell line and tumour expression
#    platforms can reduce matching scores globally, compressing the
#    score distribution and reducing confidence deltas.
#
# Manual Review Considerations:
# Low-confidence cell lines should be evaluated for:
#   - Authentication status (STR profiling, mycoplasma testing)
#   - Passage number and culture conditions
#   - Consistency with published characterisation data
#   - Potential utility despite transcriptomic ambiguity
# -----------------------------------------------------------------------------

cat("\n[4] Identifying low confidence cases...\n")

# Filter to cell lines with confidence delta below the threshold
# Uses the configurable --low-conf-thr parameter (default: 0.03)
low_conf <- summary[confidence_delta < opt$`low-conf-thr`]
cat("  Low confidence cases: ", nrow(low_conf), "\n", sep = "")

# Write low-confidence cases with selected inspection columns
low_conf_file <- file.path(opt$`output-dir`, "low_confidence_cases.tsv")
fwrite(low_conf[, .(cell_line, cell_lineage, best_tumour_lineage,
                    is_correct, confidence_delta, confidence_frac_same_lineage)],
       file = low_conf_file, sep = "\t")
cat("  Saved to:", low_conf_file, "\n")


# =============================================================================
# STEP 5: COMPONENT-LEVEL MAPPING ANALYSIS (OPTIONAL)
# =============================================================================
# When tumour graph components are available, this step analyses how cell
# lines distribute across tumour clusters. This is the reciprocal of the
# tumour inspection script's component analysis, which examines tumour
# mappings to cell line components.
#
# Biological Rationale:
# Tumour similarity graphs constructed from expression data partition into
# communities (components) representing distinct molecular states. These
# may correspond to:
#   - Cancer type (lineage-specific components)
#   - Molecular subtypes within a lineage
#   - Shared biological programmes (e.g., immune-infiltrated, proliferative)
#
# Questions Addressed:
#
#   Component Usage: Which tumour components receive the most cell line
#   mappings? Components with many mappings contain tumours that closely
#   resemble available cell line models.
#
#   Per-Lineage Preferences: Do cell lines of a given lineage preferentially
#   map to specific tumour components? Ideally, breast cancer cell lines
#   should map to breast tumour-dominated components.
#
#   Accuracy by Component: Are certain tumour components associated with
#   higher cell line mapping accuracy?
#
# Note on Data Requirements:
# This analysis requires tumour (not cell line) component assignments.
# The components file must contain columns mapping tumour identifiers to
# component numbers. If the file uses 'sample' instead of 'tumour' for the
# identifier column, the script handles the renaming automatically.
# -----------------------------------------------------------------------------

if (!is.null(components)) {
  cat("\n[5] Analysing component-level mapping...\n")

  # Handle column naming variations in components file
  # Some pipelines use 'sample' instead of 'tumour' as the identifier column
  if ("sample" %in% names(components) && !"tumour" %in% names(components)) {
    setnames(components, "sample", "tumour")
  }
  
  # Validate required columns are present after any renaming
  if (!all(c("tumour", "component") %in% names(components))) {
    stop("Components file must have columns: tumour, component (or sample, component)")
  }

  # Create a named vector mapping tumour identifiers to component numbers
  # This enables fast lookup of component assignment for any tumour
  comp_map <- setNames(components$component, components$tumour)
  
  # Assign component to each cell line based on its top-1 tumour match
  summary[, top1_component := comp_map[best_tumour]]

  # Summarise component usage across all cell lines
  comp_usage <- summary[, .(
    # Number of cell lines mapping to this tumour component
    n_cell_lines = .N,
    
    # Number of correct (lineage-concordant) mappings within this component
    n_correct = sum(is_correct),
    
    # Mean confidence delta for mappings to this component
    mean_confidence = mean(confidence_delta, na.rm = TRUE)
  ), by = top1_component]

  # Write component usage summary
  comp_file <- file.path(opt$`output-dir`, "component_usage.tsv")
  fwrite(comp_usage, file = comp_file, sep = "\t")
  cat("  Saved to:", comp_file, "\n")

  # Cross-tabulate cell line lineage and tumour component assignment
  # Reveals whether cell lines map to lineage-appropriate tumour clusters
  comp_by_lineage <- summary[, .(n = .N), by = .(cell_lineage, top1_component)]

  comp_by_lineage_file <- file.path(opt$`output-dir`, "component_by_lineage.tsv")
  fwrite(comp_by_lineage, file = comp_by_lineage_file, sep = "\t")
  cat("  Saved to:", comp_by_lineage_file, "\n")
}


# =============================================================================
# STEP 6: TOP-K LINEAGE CONSISTENCY ANALYSIS
# =============================================================================
# This step evaluates how consistently the top-k tumour matches belong to
# the expected cancer type for each cell line. Lineage consistency provides
# a robust measure of whether cell lines maintain transcriptomic identity
# with their annotated cancer type.
#
# Interpretation of Consistency Levels:
#
#   High Consistency (all top-k same lineage):
#   The cell line's transcriptomic profile is clearly and unambiguously
#   associated with its annotated cancer type. Tumours of that type form
#   a distinct cluster that the cell line matches preferentially.
#
#   Mixed Consistency (some top-k from different lineages):
#   Multiple cancer types are represented in the top-k tumour matches.
#   This pattern may indicate:
#     - Biological similarity between cancer types (shared developmental
#       origin or convergent transcriptional programmes)
#     - Cell line with atypical features for its annotated lineage
#     - Partial transcriptomic drift towards a generic cancer state
#     - Potential annotation uncertainty
#
#   Low Consistency (most top-k from different lineages):
#   The cell line profile is highly atypical for its annotated lineage.
#   Strong candidates for re-evaluation of annotation accuracy.
#
# Implications for Model Selection:
# Cell lines with high lineage consistency are better validated as
# representative models of their cancer type. Those with mixed consistency
# should be used with awareness of their ambiguous identity, or reserved
# for studies where lineage specificity is less critical.
# -----------------------------------------------------------------------------

cat("\n[6] Analysing top-k lineage consistency...\n")

# Compute lineage consistency metrics for each cell line
# The analysis uses the configurable top-k parameter (default: 10)
topk_consistency <- rankings[rank <= opt$`top-k`, .(
  # Total number of tumours in the top-k (should equal top-k if complete)
  n_total = .N,
  
  # Count of tumours matching the cell line's lineage
  # Compares tumour_lineage to cl_lineage (cell line's annotated type)
  n_same_lineage = sum(tumour_lineage == cl_lineage, na.rm = TRUE),
  
  # Fraction of top-k from the same lineage
  frac_same = mean(tumour_lineage == cl_lineage, na.rm = TRUE)
), by = .(cell_line, cl_lineage)]

# Flag cell lines with perfect consistency (all top-k from same lineage)
# These represent the most confidently validated cell line models
topk_consistency[, all_same := n_same_lineage == n_total]

# Write lineage consistency results
consistency_file <- file.path(opt$`output-dir`, "topk_lineage_consistency.tsv")
fwrite(topk_consistency, file = consistency_file, sep = "\t")
cat("  Saved to:", consistency_file, "\n")

# Report count of perfectly consistent cell lines
cat("\n  Cell lines with all top-10 tumours from same lineage: ",
    sum(topk_consistency$all_same), " / ", nrow(topk_consistency), "\n", sep = "")


# =============================================================================
# STEP 7: SCORE DISTRIBUTION ANALYSIS
# =============================================================================
# This step analyses the distribution of best-match similarity scores
# across cell line lineages. Score distributions reveal how well cell lines
# match their best tumour candidates and whether certain cancer types have
# systematically better tumour representation.
#
# Interpreting Score Distributions:
#
#   High Mean/Median Scores:
#   Cell lines of this lineage tend to have tumours that closely match
#   their transcriptomic profiles. This suggests good tumour representation
#   in the reference cohort and/or well-preserved lineage identity in the
#   cell lines.
#
#   Low Mean/Median Scores:
#   Cell lines poorly match even their best tumour candidates. This may
#   indicate:
#     - Insufficient tumour representation for this cancer type in the
#       reference cohort
#     - Systematic transcriptomic drift in cell lines of this lineage
#     - Technical batch effects between platforms
#     - Cell lines derived from unusual tumour subtypes
#
#   Wide Score Range:
#   High variability in matching quality within the lineage. Some cell
#   lines have excellent tumour matches while others do not, suggesting
#   heterogeneity in cell line quality or derivation conditions.
#
# Relationship to Accuracy:
# The mean_correct column reports the fraction of cell lines with lineage-
# concordant top-1 mappings. Comparing this to score statistics reveals
# whether low scores predict incorrect mappings within each lineage.
# -----------------------------------------------------------------------------

cat("\n[7] Analysing score distributions...\n")

# Compute per-lineage summary statistics for best-match scores
score_summary <- summary[, .(
  # Sample size
  n = .N,
  
  # Central tendency of best-match scores
  mean_score = mean(best_score, na.rm = TRUE),
  median_score = median(best_score, na.rm = TRUE),
  
  # Score range indicating variability in matching quality
  min_score = min(best_score, na.rm = TRUE),
  max_score = max(best_score, na.rm = TRUE),
  
  # Top-1 accuracy (fraction of correct mappings)
  mean_correct = mean(is_correct, na.rm = TRUE)
), by = .(lineage = cell_lineage)]

# Write score statistics to output file
score_file <- file.path(opt$`output-dir`, "score_statistics.tsv")
fwrite(score_summary, file = score_file, sep = "\t")
cat("  Saved to:", score_file, "\n")


# =============================================================================
# STEP 8: GENERATE SUMMARY REPORT
# =============================================================================
# This step produces a human-readable Markdown report that consolidates
# all inspection findings for rapid assessment of cell line mapping quality.
# The report serves as the primary entry point for quality control review.
#
# Report Sections:
#
#   Overall Performance: Global metrics including total cell line count,
#   top-1 accuracy, incorrect mapping count, and mean confidence measures.
#   Provides immediate assessment of cell line panel quality.
#
#   Per-Lineage Confidence: The confidence statistics table showing
#   sample size, confidence delta distribution, and low-confidence counts
#   for each cancer type. Identifies lineages requiring attention.
#
#   Incorrect Mappings: Count and cross-tabulation of incorrect mappings
#   by cell line lineage and best tumour lineage. Reveals systematic
#   patterns that may indicate annotation issues or biological relationships.
#
#   Low Confidence Cases: Count and threshold used for flagging cases
#   requiring manual review.
# -----------------------------------------------------------------------------

cat("\n[8] Creating inspection summary...\n")

# Construct the Markdown-formatted summary report
summary_text <- paste0(
  "# Cell line -> tumour Mapping Inspection Summary\n\n",
  
  # Overall performance metrics section
  "## Overall Performance\n",
  "- Total cell lines: ", nrow(summary), "\n",
  "- Top-1 accuracy: ", round(mean(summary$is_correct) * 100, 2), "%\n",
  "- Incorrect mappings: ", sum(!summary$is_correct), "\n",
  "- Mean confidence delta: ", round(mean(summary$confidence_delta, na.rm = TRUE), 3), "\n",
  "- Mean top-", opt$`top-k`, " same-lineage fraction: ",
  round(mean(summary$confidence_frac_same_lineage, na.rm = TRUE), 3), "\n\n",
  
  # Per-lineage confidence statistics
  "## Per-Lineage Confidence\n",
  paste(capture.output(print(conf_summary)), collapse = "\n"), "\n\n",
  
  # Incorrect mapping breakdown
  "## Incorrect Mappings\n",
  ifelse(nrow(incorrect) > 0,
         paste0("- Count: ", nrow(incorrect), "\n",
                "- Breakdown:\n",
                paste(capture.output(print(table(incorrect$cell_lineage,
                                                 incorrect$best_tumour_lineage))),
                      collapse = "\n")),
         "- None\n"),
  
  # Low confidence case summary
  "\n\n## Low Confidence Cases\n",
  "- Count: ", nrow(low_conf), "\n",
  "- Threshold: confidence_delta < ", opt$`low-conf-thr`, "\n"
)

# Write the summary report to a Markdown file
report_file <- file.path(opt$`output-dir`, "INSPECTION_SUMMARY.md")
writeLines(summary_text, report_file)
cat("  Saved to:", report_file, "\n")


# =============================================================================
# EXECUTION SUMMARY
# =============================================================================
# Final output confirming successful completion and indicating where output
# files have been saved. This structured termination message facilitates
# integration with pipeline logging systems and provides clear confirmation
# that the analysis completed without error.
# -----------------------------------------------------------------------------

cat("\n========================================\n")
cat("INSPECTION COMPLETE\n")
cat("========================================\n")
cat("\nOutput files saved to:", opt$`output-dir`, "\n")