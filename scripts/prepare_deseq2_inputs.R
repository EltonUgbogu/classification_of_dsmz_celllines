#!/usr/bin/env Rscript

# ============================================================
# Prepare DESeq2 inputs: extract counts and metadata per disease
# ============================================================
#
# This script:
# 1) Reads DSMZ_count_gene.rds (all cell lines)
# 2) Filters to disease-specific cell lines (from node stats)
# 3) Creates metadata with sample_id, cell_line columns
# 4) Writes counts TSV and metadata TSV for DESeq2
#
# Usage:
# Rscript scripts/prepare_deseq2_inputs.R \
#   --profile brca \
#   --dsmz_counts data/dsmz/DSMZ_count_gene.rds \
#   --dsmz_meta data/dsmz/DSMZ_metadata.csv \
#   --node_stats results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all/dsmz_cellline_graph_node_stats.tsv \
#   --outdir results/unsupervised/brca/deseq2_inputs
#
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(readr)
})

# -----------------------------------------------------------------------------
# Cell Line Canonicalization Function
# -----------------------------------------------------------------------------
# This function normalizes cell line names by:
#   - Converting to uppercase
#   - Removing all non-alphanumeric characters (underscores, hyphens, spaces, etc.)
#
# This handles naming inconsistencies like:
#   - LAN_5, LAN-5, LAN5 → all become LAN5
#   - CHP_126, CHP-126, CHP126 → all become CHP126
#
# Parameters:
#   x: Character vector of cell line names
#
# Returns:
#   Character vector of canonicalized names

canon_cl <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[^A-Z0-9]", "", x)   # remove underscores, hyphens, spaces, etc.
  x
}

opt_list <- list(
  make_option("--profile", type="character", help="Profile name: brca, nbl, or rbl"),
  make_option("--dsmz_counts", type="character", default="data/dsmz/DSMZ_count_gene.rds", help="Path to DSMZ count RDS"),
  make_option("--dsmz_meta", type="character", default="data/dsmz/DSMZ_metadata.csv", help="Path to DSMZ metadata CSV"),
  make_option("--node_stats", type="character", help="Path to node stats TSV with cell_line column"),
  make_option("--outdir", type="character", help="Output directory for counts.tsv and metadata.tsv")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# Read node stats to get cell lines for this disease
node_stats <- read_tsv(opt$node_stats, show_col_types = FALSE)
target_cell_lines <- unique(node_stats$cell_line)
target_canon <- canon_cl(target_cell_lines)
message(sprintf("[INFO] Found %d cell lines for %s", length(target_cell_lines), opt$profile))

# Read DSMZ counts
counts_mat <- readRDS(opt$dsmz_counts)
counts_df <- as.data.frame(counts_mat, check.names = FALSE)
counts_df <- data.frame(gene_id = rownames(counts_mat), counts_df, check.names = FALSE)
message(sprintf("[INFO] DSMZ counts: %d genes, %d samples", nrow(counts_df), ncol(counts_df) - 1))

# Read DSMZ metadata
meta_df <- read_csv(opt$dsmz_meta, show_col_types = FALSE)
message(sprintf("[INFO] DSMZ metadata: %d rows", nrow(meta_df)))

# Normalize cell line names in metadata (handle variations like CHP-126 vs CHP_126)
# The metadata has both Cell_Line and DSMZ_Cell_line_norm columns
# We create both a "clean" version (preferring normalized) and canonical versions
meta_df <- meta_df %>%
  mutate(
    cell_line_clean = ifelse(!is.na(DSMZ_Cell_line_norm) & DSMZ_Cell_line_norm != "", 
                             DSMZ_Cell_line_norm, 
                             Cell_Line),
    cell_line_canon = canon_cl(cell_line_clean),
    Cell_Line_canon = canon_cl(Cell_Line)
  )

# Find matching samples using canonicalized names
# This handles naming inconsistencies like LAN_5 vs LAN5 vs LAN-5
matched_samples <- character()
matched_cell_lines <- character()

for (i in seq_along(target_cell_lines)) {
  cl <- target_cell_lines[i]
  cl_can <- target_canon[i]
  
  # Special case: if cl directly matches a sample column, use it (RBL case)
  if (cl %in% colnames(counts_df)) {
    matched_samples <- c(matched_samples, cl)
    matched_cell_lines <- c(matched_cell_lines, cl)
    next
  }
  
  # Try matching via metadata sample_name (direct match)
  if (cl %in% meta_df$sample_name) {
    if (cl %in% colnames(counts_df)) {
      matched_samples <- c(matched_samples, cl)
      matched_cell_lines <- c(matched_cell_lines, cl)
      next
    }
  }
  
  # Match by canonicalized cell line name
  # This handles LAN_5 → LAN5, CHP_126 → CHP126, etc.
  meta_matches <- meta_df %>%
    filter(cell_line_canon == cl_can | Cell_Line_canon == cl_can)
  
  if (nrow(meta_matches) > 0) {
    # Retain all count-bearing libraries for the biological cell-line group.
    # This keeps RBL_15 and RBL_20 replicate profiles under one grouped focal identity.
    sample_names <- unique(meta_matches$sample_name[!is.na(meta_matches$sample_name) & meta_matches$sample_name != ""])
    for (sn in sample_names) {
      if (sn %in% colnames(counts_df)) {
        matched_samples <- c(matched_samples, sn)
        matched_cell_lines <- c(matched_cell_lines, cl)  # Keep node-stats label as cell_line
      }
    }
  }
}

# Deduplicate matched samples (keep first mapping)
keep_idx <- !duplicated(matched_samples)
matched_samples <- matched_samples[keep_idx]
matched_cell_lines <- matched_cell_lines[keep_idx]

message(sprintf("[INFO] Matched %d samples for %d target cell lines",
                length(matched_samples), length(target_cell_lines)))

# Report unmatched cell lines for debugging
unmatched <- target_cell_lines[!(target_cell_lines %in% matched_cell_lines)]
if (length(unmatched) > 0) {
  message("[WARN] Unmatched cell lines from node stats:")
  message(paste("  ", unmatched, collapse = "\n"))
  message("[WARN] These may not exist in DSMZ metadata/counts, or have naming mismatches")
}

if (length(matched_samples) == 0) {
  stop("No samples matched! Check cell line name mappings.")
}

# Filter counts to matched samples
counts_subset <- counts_df[, c("gene_id", matched_samples), drop = FALSE]

# Create metadata for matched samples
meta_subset <- data.frame(
  sample_id = matched_samples,
  cell_line = matched_cell_lines,
  stringsAsFactors = FALSE
) %>%
  left_join(
    meta_df %>% select(sample_name, Cell_Line, Disease, Project_id),
    by = c("sample_id" = "sample_name")
  ) %>%
  mutate(cell_line_display = dplyr::coalesce(as.character(Cell_Line), cell_line)) %>%
  select(sample_id, cell_line, cell_line_display, everything())

# Write outputs
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Write counts TSV (gene_id as first column, then samples)
counts_out <- counts_subset[, c("gene_id", matched_samples), drop = FALSE]

write_tsv(counts_out, file.path(opt$outdir, "counts.tsv"))
message(sprintf("[OK] Wrote counts: %s", file.path(opt$outdir, "counts.tsv")))

# Write metadata TSV
write_tsv(meta_subset, file.path(opt$outdir, "metadata.tsv"))
message(sprintf("[OK] Wrote metadata: %s", file.path(opt$outdir, "metadata.tsv")))

message(sprintf("[DONE] Prepared inputs for %s: %d samples, %d genes", 
                opt$profile, length(matched_samples), nrow(counts_out)))
