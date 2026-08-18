#!/usr/bin/env Rscript

# The script compares protected Version 1 outputs with Version 2 development
# outputs. It writes comparison tables only in the Version 2 namespace.

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

write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
}

as_numeric_safe <- function(values) {
  suppressWarnings(as.numeric(values))
}

arguments <- parse_arguments()
v1_edges <- read_tsv(require_argument(arguments, "v1-edges"))
v1_rank_summary <- read_tsv(require_argument(arguments, "v1-rank-summary"))
weighted_edges <- read_tsv(require_argument(arguments, "weighted-edges"))
boundary_rankings <- read_tsv(require_argument(arguments, "boundary-rankings"))
pdg_nodes <- read_tsv(require_argument(arguments, "pdg-nodes"))
edge_comparison_out <- require_argument(arguments, "edge-comparison-out")
ranking_comparison_out <- require_argument(arguments, "ranking-comparison-out")
node_status_comparison_out <- require_argument(arguments, "node-status-comparison-out")
report_out <- require_argument(arguments, "report-out")

v1_edge_table <- data.frame(
  source_node = pmin(v1_edges$from, v1_edges$to),
  target_node = pmax(v1_edges$from, v1_edges$to),
  v1_present = TRUE,
  v1_edge_weight = as_numeric_safe(v1_edges$weight),
  stringsAsFactors = FALSE
)
v1_edge_table <- v1_edge_table[!duplicated(paste(v1_edge_table$source_node, v1_edge_table$target_node, sep = "||")), , drop = FALSE]

v2_edge_table <- data.frame(
  source_node = weighted_edges$source_node,
  target_node = weighted_edges$target_node,
  v2_present = TRUE,
  v2_weighted_consensus_score = as_numeric_safe(weighted_edges$weighted_consensus_score),
  v2_edge_resolution_status = weighted_edges$edge_resolution_status,
  v2_boundary_status = weighted_edges$boundary_status,
  v2_boundary_reason = weighted_edges$boundary_reason,
  stringsAsFactors = FALSE
)
edge_comparison <- merge(v1_edge_table, v2_edge_table, by = c("source_node", "target_node"), all = TRUE)
edge_comparison$v1_present[is.na(edge_comparison$v1_present)] <- FALSE
edge_comparison$v2_present[is.na(edge_comparison$v2_present)] <- FALSE
edge_comparison$score_delta_v2_minus_v1 <- edge_comparison$v2_weighted_consensus_score - edge_comparison$v1_edge_weight
edge_comparison$comparison_status <- ifelse(
  edge_comparison$v1_present & edge_comparison$v2_present,
  "shared_edge",
  ifelse(edge_comparison$v1_present, "v1_only_edge", "v2_only_edge")
)
write_tsv(edge_comparison, edge_comparison_out)

ranking_identifier <- if ("cell_line_group" %in% names(v1_rank_summary)) {
  v1_rank_summary$cell_line_group
} else if ("cell_line" %in% names(v1_rank_summary)) {
  v1_rank_summary$cell_line
} else {
  paste0("ranking_", seq_len(nrow(v1_rank_summary)))
}
v1_ranking_table <- data.frame(
  ranking_id = ranking_identifier,
  v1_top1_score = if ("top1_score" %in% names(v1_rank_summary)) as_numeric_safe(v1_rank_summary$top1_score) else NA_real_,
  v1_delta_top1_top2 = if ("delta_top1_top2" %in% names(v1_rank_summary)) as_numeric_safe(v1_rank_summary$delta_top1_top2) else NA_real_,
  v1_top1_same_lineage = if ("top1_same_lineage" %in% names(v1_rank_summary)) v1_rank_summary$top1_same_lineage else NA,
  stringsAsFactors = FALSE
)
ranking_comparison <- merge(v1_ranking_table, boundary_rankings, by = "ranking_id", all = TRUE)
ranking_comparison$comparison_status <- ifelse(
  !is.na(ranking_comparison$v1_top1_score) & !is.na(ranking_comparison$node_resolution_status),
  "shared_ranking",
  ifelse(!is.na(ranking_comparison$v1_top1_score), "v1_only_ranking", "v2_only_ranking")
)
write_tsv(ranking_comparison, ranking_comparison_out)

node_status_comparison <- data.frame(
  node_id = pdg_nodes$node_id,
  v2_node_resolution_status = pdg_nodes$node_resolution_status,
  v2_boundary_status = pdg_nodes$boundary_status,
  bridge_like_anchor_score = if ("bridge_like_anchor_score" %in% names(pdg_nodes)) pdg_nodes$bridge_like_anchor_score else NA,
  v1_reference_status = "present_in_protected_v1_graph_input",
  comparison_status = "v1_node_with_v2_status",
  stringsAsFactors = FALSE
)
write_tsv(node_status_comparison, node_status_comparison_out)

report_lines <- c(
  "# Version 1 versus Version 2 comparison report",
  "",
  "This report compares protected Version 1 inputs with Version 2 development outputs.",
  "",
  "## Comparison tables",
  "",
  paste0("- Edge comparison: `", edge_comparison_out, "`"),
  paste0("- Ranking comparison: `", ranking_comparison_out, "`"),
  paste0("- Node status comparison: `", node_status_comparison_out, "`"),
  "",
  "## Interpretation",
  "",
  "The comparison is a development check. It does not validate Version 2 scientifically and does not modify Version 1 outputs.",
  "",
  "## Counts",
  "",
  paste0("- Edge comparison rows: ", nrow(edge_comparison)),
  paste0("- Ranking comparison rows: ", nrow(ranking_comparison)),
  paste0("- Node status comparison rows: ", nrow(node_status_comparison))
)
write_lines(report_lines, report_out)

