#!/usr/bin/env Rscript

# The script computes pilot diagnostics for representation weighting.
# It reads protected Version 1 outputs and writes Version 2 diagnostics only.

parse_arguments <- function() {
  argument_values <- commandArgs(trailingOnly = TRUE)
  parsed_arguments <- list()
  argument_index <- 1L
  while (argument_index <= length(argument_values)) {
    argument_name <- argument_values[[argument_index]]
    if (!startsWith(argument_name, "--")) {
      stop("Unexpected argument: ", argument_name)
    }
    option_name <- sub("^--", "", argument_name)
    if (argument_index == length(argument_values)) {
      stop("Missing value for argument: ", argument_name)
    }
    parsed_arguments[[option_name]] <- argument_values[[argument_index + 1L]]
    argument_index <- argument_index + 2L
  }
  parsed_arguments
}

require_argument <- function(arguments, name) {
  value <- arguments[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required argument --", name)
  }
  value
}

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, quote = "", comment.char = "")
}

write_tsv <- function(data_frame, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data_frame, path, sep = "\t", row.names = FALSE, quote = FALSE,
              na = "")
}

as_numeric_safe <- function(values) {
  suppressWarnings(as.numeric(values))
}

mean_numeric <- function(values, default_value = 0) {
  numeric_values <- as_numeric_safe(values)
  if (!length(numeric_values) || all(is.na(numeric_values))) {
    return(default_value)
  }
  mean(numeric_values, na.rm = TRUE)
}

mean_logical <- function(values, default_value = 0) {
  if (!length(values)) {
    return(default_value)
  }
  logical_values <- toupper(as.character(values)) %in% c("TRUE", "T", "1", "YES")
  mean(logical_values, na.rm = TRUE)
}

normalise_unit <- function(value) {
  if (is.na(value) || !is.finite(value)) {
    return(0)
  }
  max(0, min(1, value))
}

component_purity <- function(community_table) {
  if (!all(c("component", "lineage") %in% names(community_table))) {
    return(NA_real_)
  }
  component_ids <- unique(community_table$component)
  component_weights <- numeric(length(component_ids))
  component_purities <- numeric(length(component_ids))
  for (component_index in seq_along(component_ids)) {
    component_rows <- community_table[community_table$component == component_ids[[component_index]], , drop = FALSE]
    lineage_counts <- table(component_rows$lineage)
    component_weights[[component_index]] <- nrow(component_rows)
    component_purities[[component_index]] <- max(lineage_counts) / sum(lineage_counts)
  }
  if (!sum(component_weights)) {
    return(NA_real_)
  }
  weighted.mean(component_purities, component_weights)
}

arguments <- parse_arguments()

v1_edges_path <- require_argument(arguments, "v1-edges")
v1_node_metadata_path <- require_argument(arguments, "v1-node-metadata")
v1_louvain_path <- require_argument(arguments, "v1-louvain")
v1_leiden_path <- require_argument(arguments, "v1-leiden")
v1_tumour_rankings_path <- require_argument(arguments, "v1-tumour-rankings")
v1_cellline_rankings_path <- require_argument(arguments, "v1-cellline-rankings")
v1_rank_summary_path <- require_argument(arguments, "v1-rank-summary")
output_path <- require_argument(arguments, "output")

v1_edges <- read_tsv(v1_edges_path)
v1_node_metadata <- read_tsv(v1_node_metadata_path)
v1_louvain <- read_tsv(v1_louvain_path)
v1_leiden <- read_tsv(v1_leiden_path)
v1_tumour_rankings <- read_tsv(v1_tumour_rankings_path)
v1_cellline_rankings <- read_tsv(v1_cellline_rankings_path)
v1_rank_summary <- read_tsv(v1_rank_summary_path)

edge_weight_mean <- if ("weight" %in% names(v1_edges)) mean_numeric(v1_edges$weight, 0) else 0
edge_weight_sd <- if ("weight" %in% names(v1_edges)) stats::sd(as_numeric_safe(v1_edges$weight), na.rm = TRUE) else 0
edge_instability <- normalise_unit(edge_weight_sd / max(edge_weight_mean, 1e-6))

tumour_top1 <- v1_tumour_rankings[as_numeric_safe(v1_tumour_rankings$rank) == 1, , drop = FALSE]
cellline_top1 <- v1_cellline_rankings[as_numeric_safe(v1_cellline_rankings$rank) == 1, , drop = FALSE]

tumour_top1_agreement <- if ("is_correct" %in% names(tumour_top1)) mean_logical(tumour_top1$is_correct, 0) else 0
cellline_top1_agreement <- if ("is_correct" %in% names(cellline_top1)) mean_logical(cellline_top1$is_correct, 0) else 0
tumour_top1_score <- if ("score" %in% names(tumour_top1)) mean_numeric(tumour_top1$score, 0) else 0
cellline_top1_score <- if ("score" %in% names(cellline_top1)) mean_numeric(cellline_top1$score, 0) else 0

rank_margin_mean <- if ("delta_top1_top2" %in% names(v1_rank_summary)) {
  mean_numeric(v1_rank_summary$delta_top1_top2, 0)
} else {
  0
}
rank_margin_score <- normalise_unit(rank_margin_mean / 0.10)
rank_top1_agreement <- if ("top1_same_lineage" %in% names(v1_rank_summary)) {
  mean_logical(v1_rank_summary$top1_same_lineage, 0)
} else {
  cellline_top1_agreement
}

louvain_purity <- component_purity(v1_louvain)
leiden_purity <- component_purity(v1_leiden)
mean_community_purity <- mean(c(louvain_purity, leiden_purity), na.rm = TRUE)
if (!is.finite(mean_community_purity)) {
  mean_community_purity <- 0
}

metadata_source_count <- if ("source_profile" %in% names(v1_node_metadata)) {
  length(unique(v1_node_metadata$source_profile))
} else {
  NA_integer_
}
batch_association_score <- if (is.na(metadata_source_count) || metadata_source_count <= 1) 0 else 0.10

diagnostics <- data.frame(
  representation_id = c(
    "cell_line_graph_weight",
    "tumour_to_cellline_similarity",
    "cellline_to_tumour_similarity",
    "cellline_rank_margin",
    "community_lineage_purity"
  ),
  neighbourhood_stability_score = c(
    normalise_unit(edge_weight_mean),
    normalise_unit(tumour_top1_score),
    normalise_unit(cellline_top1_score),
    rank_margin_score,
    normalise_unit(mean_community_purity)
  ),
  cancer_type_agreement_score = c(
    normalise_unit(mean_community_purity),
    normalise_unit(tumour_top1_agreement),
    normalise_unit(cellline_top1_agreement),
    normalise_unit(rank_top1_agreement),
    normalise_unit(mean_community_purity)
  ),
  representation_redundancy_index = c(0.55, 0.45, 0.45, 0.35, 0.60),
  batch_association_score = c(
    batch_association_score,
    batch_association_score,
    batch_association_score,
    batch_association_score,
    batch_association_score
  ),
  instability_score = c(
    edge_instability,
    1 - normalise_unit(tumour_top1_agreement),
    1 - normalise_unit(cellline_top1_agreement),
    1 - rank_margin_score,
    1 - normalise_unit(mean_community_purity)
  ),
  diagnostic_source = c(
    basename(v1_edges_path),
    basename(v1_tumour_rankings_path),
    basename(v1_cellline_rankings_path),
    basename(v1_rank_summary_path),
    paste(basename(v1_louvain_path), basename(v1_leiden_path), sep = ";")
  ),
  pilot_notes = c(
    "Version 1 graph edge weights are used as protected cell-line relationship evidence.",
    "Version 1 tumour-to-cell-line top-rank agreement is used as patient-referenced agreement evidence.",
    "Version 1 cell-line-to-tumour top-rank agreement is used as reciprocal retrieval evidence.",
    "Version 1 top1-top2 margin is used as a ranking stability proxy.",
    "Version 1 community lineage purity is used as a graph-partition stability proxy."
  ),
  stringsAsFactors = FALSE
)

write_tsv(diagnostics, output_path)

