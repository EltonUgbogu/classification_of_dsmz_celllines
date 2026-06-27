#!/usr/bin/env Rscript

# The script builds weighted Version 2 consensus evidence from protected
# Version 1 graph and ranking outputs.

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

numeric_argument <- function(arguments, name, default_value) {
  value <- arguments[[name]]
  if (is.null(value) || !nzchar(value)) {
    return(default_value)
  }
  numeric_value <- suppressWarnings(as.numeric(value))
  if (!is.finite(numeric_value)) {
    stop("Argument --", name, " must be numeric")
  }
  numeric_value
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

normalise_vector <- function(values) {
  numeric_values <- as_numeric_safe(values)
  numeric_values[!is.finite(numeric_values)] <- NA_real_
  if (all(is.na(numeric_values))) {
    return(rep(0, length(values)))
  }
  value_min <- min(numeric_values, na.rm = TRUE)
  value_max <- max(numeric_values, na.rm = TRUE)
  if (!is.finite(value_min) || !is.finite(value_max) || value_max == value_min) {
    return(rep(1, length(values)))
  }
  normalised_values <- (numeric_values - value_min) / (value_max - value_min)
  normalised_values[is.na(normalised_values)] <- 0
  pmax(0, pmin(1, normalised_values))
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

map_by_sample <- function(table_data, value_column) {
  if (!all(c("sample", value_column) %in% names(table_data))) {
    return(setNames(character(0), character(0)))
  }
  mapped_values <- table_data[[value_column]]
  names(mapped_values) <- table_data$sample
  mapped_values
}

arguments <- parse_arguments()

v1_edges <- read_tsv(require_argument(arguments, "v1-edges"))
v1_louvain <- read_tsv(require_argument(arguments, "v1-louvain"))
v1_leiden <- read_tsv(require_argument(arguments, "v1-leiden"))
v1_tumour_rankings <- read_tsv(require_argument(arguments, "v1-tumour-rankings"))
v1_cellline_rankings <- read_tsv(require_argument(arguments, "v1-cellline-rankings"))
v1_rank_summary <- read_tsv(require_argument(arguments, "v1-rank-summary"))
weights <- read_tsv(require_argument(arguments, "weights"))
top_rank_count <- numeric_argument(arguments, "top-rank-count", 10)
output_path <- require_argument(arguments, "output")

required_edge_columns <- c("from", "to", "weight")
missing_edge_columns <- setdiff(required_edge_columns, names(v1_edges))
if (length(missing_edge_columns)) {
  stop("Version 1 edge table is missing required columns: ", paste(missing_edge_columns, collapse = ", "))
}
if (!all(c("representation_id", "representation_weight") %in% names(weights))) {
  stop("Weight table must contain representation_id and representation_weight")
}

weight_vector <- as.numeric(weights$representation_weight)
names(weight_vector) <- weights$representation_id
weight_for <- function(representation_id) {
  value <- weight_vector[[representation_id]]
  if (is.null(value) || is.na(value) || !is.finite(value)) {
    return(0)
  }
  value
}

source_node <- pmin(v1_edges$from, v1_edges$to)
target_node <- pmax(v1_edges$from, v1_edges$to)
edge_key <- paste(source_node, target_node, sep = "||")
unique_edge_index <- !duplicated(edge_key)
edge_table <- data.frame(
  source_node = source_node[unique_edge_index],
  target_node = target_node[unique_edge_index],
  v1_edge_weight = as_numeric_safe(v1_edges$weight[unique_edge_index]),
  stringsAsFactors = FALSE
)

cell_line_graph_weight_score <- normalise_vector(edge_table$v1_edge_weight)

louvain_component <- map_by_sample(v1_louvain, "component")
louvain_lineage <- map_by_sample(v1_louvain, "lineage")
leiden_component <- map_by_sample(v1_leiden, "component")

source_louvain_component <- louvain_component[edge_table$source_node]
target_louvain_component <- louvain_component[edge_table$target_node]
source_leiden_component <- leiden_component[edge_table$source_node]
target_leiden_component <- leiden_component[edge_table$target_node]
source_lineage <- louvain_lineage[edge_table$source_node]
target_lineage <- louvain_lineage[edge_table$target_node]

same_louvain_component <- !is.na(source_louvain_component) & source_louvain_component == target_louvain_component
same_leiden_component <- !is.na(source_leiden_component) & source_leiden_component == target_leiden_component
same_lineage <- !is.na(source_lineage) & source_lineage == target_lineage

community_lineage_purity_score <- ifelse(same_louvain_component & same_leiden_component, 1,
                                         ifelse(same_louvain_component | same_leiden_component, 0.75,
                                                ifelse(same_lineage, 0.50, 0)))

rank_margin_global_score <- if ("delta_top1_top2" %in% names(v1_rank_summary)) {
  pmax(0, pmin(1, mean_numeric(v1_rank_summary$delta_top1_top2, 0) / 0.10))
} else {
  0
}
cellline_rank_margin_score <- rep(rank_margin_global_score, nrow(edge_table))

tumour_top_rankings <- v1_tumour_rankings[as_numeric_safe(v1_tumour_rankings$rank) <= top_rank_count, , drop = FALSE]
tumour_to_cellline_similarity_score <- if ("is_correct" %in% names(tumour_top_rankings)) {
  rep(mean_logical(tumour_top_rankings$is_correct, 0), nrow(edge_table))
} else {
  rep(0, nrow(edge_table))
}

cellline_top_rankings <- v1_cellline_rankings[as_numeric_safe(v1_cellline_rankings$rank) <= top_rank_count, , drop = FALSE]
cellline_to_tumour_similarity_score <- if ("is_correct" %in% names(cellline_top_rankings)) {
  rep(mean_logical(cellline_top_rankings$is_correct, 0), nrow(edge_table))
} else {
  rep(0, nrow(edge_table))
}

weighted_consensus_score <- (
  weight_for("cell_line_graph_weight") * cell_line_graph_weight_score +
    weight_for("tumour_to_cellline_similarity") * tumour_to_cellline_similarity_score +
    weight_for("cellline_to_tumour_similarity") * cellline_to_tumour_similarity_score +
    weight_for("cellline_rank_margin") * cellline_rank_margin_score +
    weight_for("community_lineage_purity") * community_lineage_purity_score
)

edge_output <- data.frame(
  source_node = edge_table$source_node,
  target_node = edge_table$target_node,
  v1_edge_weight = edge_table$v1_edge_weight,
  cell_line_graph_weight_score = cell_line_graph_weight_score,
  tumour_to_cellline_similarity_score = tumour_to_cellline_similarity_score,
  cellline_to_tumour_similarity_score = cellline_to_tumour_similarity_score,
  cellline_rank_margin_score = cellline_rank_margin_score,
  community_lineage_purity_score = community_lineage_purity_score,
  weighted_consensus_score = weighted_consensus_score,
  same_louvain_component = same_louvain_component,
  same_leiden_component = same_leiden_component,
  same_lineage = same_lineage,
  edge_resolution_status = "unclassified_prebootstrap",
  boundary_status = "not_classified",
  boundary_reason = "",
  weighting_method = "deterministic_simplex_normalised_pilot",
  pilot_mapping_limitation = "Version 1 ranking tables are group-level whereas graph nodes are profile-level; ranking evidence is used as global calibration until profile-level mapping is implemented.",
  stringsAsFactors = FALSE
)

write_tsv(edge_output, output_path)

