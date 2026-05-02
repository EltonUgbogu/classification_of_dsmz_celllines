#!/usr/bin/env Rscript

# ============================================================
# Add component assignments to metadata from node stats
# ============================================================
#
# This script reads the dsmz_cellline_graph_node_stats.tsv file
# and adds component assignments to your metadata file.
#
# Usage:
# Rscript scripts/add_component_to_metadata.R \
#   --node_stats results/unsupervised/brca/tumour_neighbourhoods/final_consensus_all/dsmz_cellline_graph_node_stats.tsv \
#   --meta data/brca/metadata.tsv \
#   --cell_line_col DSMZ_Cell_line_norm \
#   --output data/brca/metadata_with_components.tsv
#
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

opt_list <- list(
  make_option("--node_stats", type="character", help="Path to dsmz_cellline_graph_node_stats.tsv"),
  make_option("--meta", type="character", help="Input metadata TSV/CSV"),
  make_option("--cell_line_col", type="character", default="DSMZ_Cell_line_norm", help="Column in meta containing cell line names"),
  make_option("--output", type="character", help="Output metadata with component column"),
  make_option("--meta_sep", type="character", default="\t", help="Separator for input meta (default tab)"),
  make_option("--output_sep", type="character", default="\t", help="Separator for output (default tab)")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# Read node stats
node_stats <- fread(opt$node_stats, sep = "\t", data.table = FALSE)

# Read metadata
meta <- fread(opt$meta, sep = opt$meta_sep, data.table = FALSE)

# Check cell_line_col exists
if (!(opt$cell_line_col %in% colnames(meta))) {
  stop(paste0("Metadata missing column: ", opt$cell_line_col))
}

# Merge component assignments
meta_with_comp <- meta %>%
  left_join(
    node_stats %>% select(cell_line, component, is_isolate),
    by = setNames("cell_line", opt$cell_line_col)
  )

# Check for missing assignments
missing <- meta_with_comp %>%
  filter(is.na(component)) %>%
  distinct(!!sym(opt$cell_line_col))

if (nrow(missing) > 0) {
  warning(paste0("Some cell lines have no component assignment:\n",
                 paste(missing[[opt$cell_line_col]], collapse = ", ")))
}

# Write output
fwrite(meta_with_comp, opt$output, sep = opt$output_sep)

message(sprintf("[OK] Added component assignments to metadata"))
message(sprintf("[OUT] %s", opt$output))
message(sprintf("Added columns: component, is_isolate"))
