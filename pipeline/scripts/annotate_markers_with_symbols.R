#!/usr/bin/env Rscript

# ============================================================
# Annotate marker gene lists with gene symbols
# ============================================================
#
# This script adds gene symbols to Ensembl ID-based marker lists
# using org.Hs.eg.db (Bioconductor annotation package).
#
# Usage:
# Rscript scripts/annotate_markers_with_symbols.R \
#   --input results/unsupervised/component_consensus_markers_v2/brca/component_0/component_consensus.UP.core.tsv \
#   --output results/unsupervised/component_consensus_markers_v2/brca/component_0/component_consensus.UP.core.annotated.tsv
#
# Or for batch processing:
# find results/unsupervised/component_consensus_markers_v2 -name "component_consensus.*.tsv" | \
#   while read f; do
#     out="${f%.tsv}.annotated.tsv"
#     Rscript scripts/annotate_markers_with_symbols.R --input "$f" --output "$out"
#   done
#
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
})

opt_list <- list(
  make_option("--input", type="character", help="Input TSV file with gene_id column (Ensembl IDs)"),
  make_option("--output", type="character", help="Output TSV file with added symbol column"),
  make_option("--gene_id_col", type="character", default="gene_id", help="Column name containing Ensembl IDs")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("--input and --output are required")
}

# Read marker file
df <- read.delim(opt$input, sep="\t", stringsAsFactors = FALSE)

if (!(opt$gene_id_col %in% colnames(df))) {
  stop(paste0("Column '", opt$gene_id_col, "' not found in input file"))
}

# Strip version from Ensembl IDs
df$ensembl_id <- sub("\\..*$", "", df[[opt$gene_id_col]])

# Get unique Ensembl IDs
unique_ids <- unique(df$ensembl_id[!is.na(df$ensembl_id)])

# Map to symbols
map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique_ids,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
)

# Merge
df_annotated <- merge(
  df,
  map,
  by.x = "ensembl_id",
  by.y = "ENSEMBL",
  all.x = TRUE
)

# Rename SYMBOL to symbol for consistency
if ("SYMBOL" %in% colnames(df_annotated)) {
  df_annotated <- df_annotated %>%
    rename(symbol = SYMBOL)
}

# Reorder columns: put symbol early, keep original order
if ("symbol" %in% colnames(df_annotated)) {
  other_cols <- setdiff(colnames(df_annotated), c("ensembl_id", "symbol", "GENENAME"))
  df_annotated <- df_annotated[, c(opt$gene_id_col, "ensembl_id", "symbol", "GENENAME", other_cols)]
  df_annotated <- df_annotated[, !is.na(colnames(df_annotated))]
}

# Write output
write.table(df_annotated, opt$output, sep="\t", quote=FALSE, row.names=FALSE)

n_mapped <- sum(!is.na(df_annotated$symbol))
n_total <- nrow(df_annotated)
message(sprintf("[OK] Annotated %d/%d genes (%.1f%%)", n_mapped, n_total, 100*n_mapped/n_total))
message(sprintf("[OUT] %s", opt$output))
