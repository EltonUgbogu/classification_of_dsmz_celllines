#!/usr/bin/env Rscript

# The script classifies Version 2 boundary cases without forcing binary calls.
# It writes edge, ranking, node, and summary tables in the neutral Version 2 namespace.

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

write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
}

as_numeric_safe <- function(values) {
  numeric_values <- suppressWarnings(as.numeric(values))
  numeric_values[!is.finite(numeric_values)] <- NA_real_
  numeric_values
}

as_logical_safe <- function(values) {
  toupper(as.character(values)) %in% c("TRUE", "T", "1", "YES")
}

arguments <- parse_arguments()
weighted_edges <- read_tsv(require_argument(arguments, "weighted-edges"))
v1_louvain <- read_tsv(require_argument(arguments, "v1-louvain"))
v1_leiden <- read_tsv(require_argument(arguments, "v1-leiden"))
v1_rank_summary <- read_tsv(require_argument(arguments, "v1-rank-summary"))
summary_out <- require_argument(arguments, "summary-out")
boundary_edges_out <- require_argument(arguments, "boundary-edges-out")
boundary_rankings_out <- require_argument(arguments, "boundary-rankings-out")
bridge_nodes_out <- require_argument(arguments, "bridge-nodes-out")
report_out <- require_argument(arguments, "report-out")
ranking_margin <- numeric_argument(arguments, "ranking-margin", 0.02)
bridge_fraction <- numeric_argument(arguments, "bridge-fraction", 0.35)

required_edge_columns <- c(
  "source_node", "target_node", "weighted_consensus_score",
  "edge_resolution_status", "boundary_status", "boundary_reason"
)
missing_edge_columns <- setdiff(required_edge_columns, names(weighted_edges))
if (length(missing_edge_columns)) {
  stop("Weighted edge table is missing required columns: ", paste(missing_edge_columns, collapse = ", "))
}

boundary_edges <- weighted_edges[weighted_edges$edge_resolution_status == "boundary", , drop = FALSE]
if (!nrow(boundary_edges)) {
  boundary_edges <- weighted_edges[FALSE, , drop = FALSE]
}
write_tsv(boundary_edges, boundary_edges_out)

rank_identifier <- if ("cell_line_group" %in% names(v1_rank_summary)) {
  v1_rank_summary$cell_line_group
} else if ("cell_line" %in% names(v1_rank_summary)) {
  v1_rank_summary$cell_line
} else {
  paste0("ranking_", seq_len(nrow(v1_rank_summary)))
}

top1_score <- if ("top1_score" %in% names(v1_rank_summary)) as_numeric_safe(v1_rank_summary$top1_score) else rep(NA_real_, nrow(v1_rank_summary))
score_rank2 <- if ("score_rank2" %in% names(v1_rank_summary)) as_numeric_safe(v1_rank_summary$score_rank2) else top1_score
decision_margin <- if ("delta_top1_top2" %in% names(v1_rank_summary)) {
  as_numeric_safe(v1_rank_summary$delta_top1_top2)
} else {
  top1_score - score_rank2
}
top1_same_lineage <- if ("top1_same_lineage" %in% names(v1_rank_summary)) {
  as_logical_safe(v1_rank_summary$top1_same_lineage)
} else {
  rep(NA, nrow(v1_rank_summary))
}

node_resolution_status <- ifelse(
  is.na(decision_margin),
  "unresolved",
  ifelse(decision_margin < ranking_margin, "boundary",
         ifelse(!is.na(top1_same_lineage) & !top1_same_lineage, "unresolved", "resolved"))
)
ranking_boundary_reason <- ifelse(
  node_resolution_status == "boundary",
  "low_rank_margin",
  ifelse(node_resolution_status == "unresolved", "weak_weighted_consensus", "")
)

boundary_rankings <- data.frame(
  ranking_id = rank_identifier,
  top_1_similarity = top1_score,
  top_2_similarity = score_rank2,
  decision_margin = decision_margin,
  node_resolution_status = node_resolution_status,
  boundary_status = ifelse(node_resolution_status == "boundary", "boundary", "not_boundary"),
  boundary_reason = ranking_boundary_reason,
  top1_same_lineage = top1_same_lineage,
  classification_rule = "decision_margin = top_1_similarity - top_2_similarity",
  stringsAsFactors = FALSE
)
write_tsv(boundary_rankings, boundary_rankings_out)

all_nodes <- sort(unique(c(weighted_edges$source_node, weighted_edges$target_node)))
edge_status_is_boundary <- weighted_edges$edge_resolution_status == "boundary"
bridge_rows <- lapply(all_nodes, function(node_id) {
  node_edge_index <- weighted_edges$source_node == node_id | weighted_edges$target_node == node_id
  node_edges <- weighted_edges[node_edge_index, , drop = FALSE]
  total_edge_count <- nrow(node_edges)
  boundary_edge_count <- sum(node_edges$edge_resolution_status == "boundary", na.rm = TRUE)
  retained_edge_count <- sum(node_edges$edge_resolution_status == "retained", na.rm = TRUE)
  boundary_fraction <- if (total_edge_count > 0) boundary_edge_count / total_edge_count else 0
  bridge_like_anchor_score <- boundary_fraction
  node_status <- ifelse(
    boundary_fraction >= bridge_fraction & retained_edge_count > 0,
    "bridge_like",
    ifelse(boundary_edge_count > 0, "boundary", ifelse(retained_edge_count > 0, "resolved", "unresolved"))
  )
  data.frame(
    node_id = node_id,
    total_edge_count = total_edge_count,
    retained_edge_count = retained_edge_count,
    boundary_edge_count = boundary_edge_count,
    rejected_edge_count = sum(node_edges$edge_resolution_status == "rejected", na.rm = TRUE),
    bridge_like_anchor_score = bridge_like_anchor_score,
    node_resolution_status = node_status,
    boundary_status = ifelse(node_status %in% c("bridge_like", "boundary"), "boundary_related", "not_boundary"),
    boundary_reason = ifelse(node_status == "bridge_like", "mixed_component_membership",
                             ifelse(node_status == "boundary", "unstable_representation_support", "")),
    stringsAsFactors = FALSE
  )
})
bridge_like_nodes <- do.call(rbind, bridge_rows)

louvain_lookup <- v1_louvain[, intersect(c("sample", "component", "lineage"), names(v1_louvain)), drop = FALSE]
if (all(c("sample", "component") %in% names(louvain_lookup))) {
  names(louvain_lookup)[names(louvain_lookup) == "sample"] <- "node_id"
  names(louvain_lookup)[names(louvain_lookup) == "component"] <- "v1_louvain_component"
  bridge_like_nodes <- merge(bridge_like_nodes, louvain_lookup, by = "node_id", all.x = TRUE)
}
leiden_lookup <- v1_leiden[, intersect(c("sample", "component"), names(v1_leiden)), drop = FALSE]
if (all(c("sample", "component") %in% names(leiden_lookup))) {
  names(leiden_lookup) <- c("node_id", "v1_leiden_component")
  bridge_like_nodes <- merge(bridge_like_nodes, leiden_lookup, by = "node_id", all.x = TRUE)
}
write_tsv(bridge_like_nodes, bridge_nodes_out)

edge_summary <- as.data.frame(table(
  case_type = "edge",
  status = weighted_edges$edge_resolution_status
), stringsAsFactors = FALSE)
ranking_summary <- as.data.frame(table(
  case_type = "ranking",
  status = boundary_rankings$node_resolution_status
), stringsAsFactors = FALSE)
node_summary <- as.data.frame(table(
  case_type = "node",
  status = bridge_like_nodes$node_resolution_status
), stringsAsFactors = FALSE)
summary_table <- rbind(edge_summary, ranking_summary, node_summary)
names(summary_table) <- c("case_type", "status", "count")
summary_table$classification_status <- "pilot_boundary_classification"
write_tsv(summary_table, summary_out)

report_lines <- c(
  "# Boundary-case report",
  "",
  "This Version 2 report classifies uncertain cases without forcing every edge or node into a binary retained/rejected decision.",
  "",
  "## Inputs",
  "",
  "- Protected Version 1 Louvain and Leiden community assignments were read as inputs.",
  "- Protected Version 1 cell-line ranking summaries were read as inputs.",
  "- Version 2 weighted consensus edges were read from the neutral `results/version2/` namespace.",
  "",
  "## Boundary rules",
  "",
  paste0("- Ranking decision margin threshold: ", ranking_margin),
  paste0("- Bridge-like boundary fraction threshold: ", bridge_fraction),
  "- Graph edge boundary calls use `edge_resolution_status` and `boundary_reason` from the weighted consensus table.",
  "",
  "## Scientific status",
  "",
  "This is a pilot boundary-classification implementation. It is suitable for development inspection but is not yet a scientifically validated decision layer."
)
write_lines(report_lines, report_out)

