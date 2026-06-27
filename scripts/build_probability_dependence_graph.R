#!/usr/bin/env Rscript

# The script builds a Version 2 probability-dependence graph layer.
# The graph is a development representation of weighted patient-referenced
# support and must not be described as a causal or validated biological network.

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

xml_escape <- function(values) {
  escaped_values <- gsub("&", "&amp;", as.character(values), fixed = TRUE)
  escaped_values <- gsub("<", "&lt;", escaped_values, fixed = TRUE)
  escaped_values <- gsub(">", "&gt;", escaped_values, fixed = TRUE)
  escaped_values <- gsub("\"", "&quot;", escaped_values, fixed = TRUE)
  escaped_values
}

arguments <- parse_arguments()
v1_edges <- read_tsv(require_argument(arguments, "v1-edges"))
weighted_edges <- read_tsv(require_argument(arguments, "weighted-edges"))
boundary_edges <- read_tsv(require_argument(arguments, "boundary-edges"))
bridge_nodes <- read_tsv(require_argument(arguments, "bridge-nodes"))
edges_out <- require_argument(arguments, "edges-out")
nodes_out <- require_argument(arguments, "nodes-out")
graphml_out <- require_argument(arguments, "graphml-out")
summary_out <- require_argument(arguments, "summary-out")
report_out <- require_argument(arguments, "report-out")

required_columns <- c(
  "source_node", "target_node", "weighted_consensus_score",
  "bootstrap_support_lower", "bootstrap_support_upper",
  "edge_resolution_status", "boundary_status", "boundary_reason"
)
missing_columns <- setdiff(required_columns, names(weighted_edges))
if (length(missing_columns)) {
  stop("Weighted edge table is missing required columns: ", paste(missing_columns, collapse = ", "))
}

probability_dependence_edges <- data.frame(
  source_node = weighted_edges$source_node,
  target_node = weighted_edges$target_node,
  probability_dependence_weight = weighted_edges$weighted_consensus_score,
  weighted_consensus_score = weighted_edges$weighted_consensus_score,
  bootstrap_support_lower = weighted_edges$bootstrap_support_lower,
  bootstrap_support_upper = weighted_edges$bootstrap_support_upper,
  edge_resolution_status = weighted_edges$edge_resolution_status,
  boundary_status = weighted_edges$boundary_status,
  boundary_reason = weighted_edges$boundary_reason,
  v1_edge_weight = if ("v1_edge_weight" %in% names(weighted_edges)) weighted_edges$v1_edge_weight else NA,
  interpretation_label = ifelse(
    weighted_edges$edge_resolution_status == "retained",
    "resolved edge",
    ifelse(weighted_edges$edge_resolution_status == "boundary", "uncertain dependence", "rejected edge")
  ),
  terminology_guardrail = "probability-dependence graph; not a causal graph or validated biological interaction network",
  stringsAsFactors = FALSE
)
write_tsv(probability_dependence_edges, edges_out)

all_nodes <- sort(unique(c(probability_dependence_edges$source_node, probability_dependence_edges$target_node)))
node_degree <- vapply(all_nodes, function(node_id) {
  sum(probability_dependence_edges$source_node == node_id | probability_dependence_edges$target_node == node_id)
}, numeric(1))
node_boundary_degree <- vapply(all_nodes, function(node_id) {
  node_edge_index <- probability_dependence_edges$source_node == node_id | probability_dependence_edges$target_node == node_id
  sum(probability_dependence_edges$edge_resolution_status[node_edge_index] == "boundary")
}, numeric(1))
probability_dependence_nodes <- data.frame(
  node_id = all_nodes,
  probability_dependence_degree = node_degree,
  boundary_edge_count = node_boundary_degree,
  node_resolution_status = ifelse(node_boundary_degree > 0, "boundary", "resolved"),
  boundary_status = ifelse(node_boundary_degree > 0, "boundary_related", "not_boundary"),
  bridge_like_anchor_score = NA_real_,
  stringsAsFactors = FALSE
)
if ("node_id" %in% names(bridge_nodes)) {
  bridge_columns <- intersect(
    c("node_id", "node_resolution_status", "boundary_status", "bridge_like_anchor_score", "boundary_reason"),
    names(bridge_nodes)
  )
  bridge_subset <- bridge_nodes[, bridge_columns, drop = FALSE]
  probability_dependence_nodes <- merge(
    probability_dependence_nodes,
    bridge_subset,
    by = "node_id",
    all.x = TRUE,
    suffixes = c("", "_bridge")
  )
  if ("node_resolution_status_bridge" %in% names(probability_dependence_nodes)) {
    replacement_index <- !is.na(probability_dependence_nodes$node_resolution_status_bridge)
    probability_dependence_nodes$node_resolution_status[replacement_index] <- probability_dependence_nodes$node_resolution_status_bridge[replacement_index]
  }
  if ("boundary_status_bridge" %in% names(probability_dependence_nodes)) {
    replacement_index <- !is.na(probability_dependence_nodes$boundary_status_bridge)
    probability_dependence_nodes$boundary_status[replacement_index] <- probability_dependence_nodes$boundary_status_bridge[replacement_index]
  }
  if ("bridge_like_anchor_score_bridge" %in% names(probability_dependence_nodes)) {
    probability_dependence_nodes$bridge_like_anchor_score <- probability_dependence_nodes$bridge_like_anchor_score_bridge
  }
  probability_dependence_nodes <- probability_dependence_nodes[, !grepl("_bridge$", names(probability_dependence_nodes)), drop = FALSE]
}
write_tsv(probability_dependence_nodes, nodes_out)

graphml_lines <- c(
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
  "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">",
  "  <key id=\"d0\" for=\"edge\" attr.name=\"probability_dependence_weight\" attr.type=\"double\"/>",
  "  <key id=\"d1\" for=\"edge\" attr.name=\"edge_resolution_status\" attr.type=\"string\"/>",
  "  <key id=\"d2\" for=\"edge\" attr.name=\"boundary_status\" attr.type=\"string\"/>",
  "  <key id=\"d3\" for=\"node\" attr.name=\"node_resolution_status\" attr.type=\"string\"/>",
  "  <key id=\"d4\" for=\"node\" attr.name=\"bridge_like_anchor_score\" attr.type=\"double\"/>",
  "  <graph id=\"version2_graph\" edgedefault=\"undirected\">"
)
for (node_index in seq_len(nrow(probability_dependence_nodes))) {
  node_id <- xml_escape(probability_dependence_nodes$node_id[[node_index]])
  graphml_lines <- c(
    graphml_lines,
    paste0("    <node id=\"", node_id, "\">"),
    paste0("      <data key=\"d3\">", xml_escape(probability_dependence_nodes$node_resolution_status[[node_index]]), "</data>"),
    paste0("      <data key=\"d4\">", probability_dependence_nodes$bridge_like_anchor_score[[node_index]], "</data>"),
    "    </node>"
  )
}
for (edge_index in seq_len(nrow(probability_dependence_edges))) {
  edge_id <- paste0("e", edge_index)
  graphml_lines <- c(
    graphml_lines,
    paste0("    <edge id=\"", edge_id, "\" source=\"", xml_escape(probability_dependence_edges$source_node[[edge_index]]), "\" target=\"", xml_escape(probability_dependence_edges$target_node[[edge_index]]), "\">"),
    paste0("      <data key=\"d0\">", probability_dependence_edges$probability_dependence_weight[[edge_index]], "</data>"),
    paste0("      <data key=\"d1\">", xml_escape(probability_dependence_edges$edge_resolution_status[[edge_index]]), "</data>"),
    paste0("      <data key=\"d2\">", xml_escape(probability_dependence_edges$boundary_status[[edge_index]]), "</data>"),
    "    </edge>"
  )
}
graphml_lines <- c(graphml_lines, "  </graph>", "</graphml>")
write_lines(graphml_lines, graphml_out)

summary_table <- data.frame(
  metric = c(
    "v1_input_edge_rows",
    "version2_probability_dependence_edges",
    "version2_probability_dependence_nodes",
    "retained_edges",
    "boundary_edges",
    "rejected_edges",
    "boundary_edges_reported"
  ),
  value = c(
    nrow(v1_edges),
    nrow(probability_dependence_edges),
    nrow(probability_dependence_nodes),
    sum(probability_dependence_edges$edge_resolution_status == "retained"),
    sum(probability_dependence_edges$edge_resolution_status == "boundary"),
    sum(probability_dependence_edges$edge_resolution_status == "rejected"),
    nrow(boundary_edges)
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary_table, summary_out)

report_lines <- c(
  "# Probability-dependence graph report",
  "",
  "This Version 2 graph layer represents probability-dependence development evidence derived from weighted patient-referenced support.",
  "",
  "## Terminology boundary",
  "",
  "Use: probability-dependence graph, weighted patient-referenced support, boundary edge, resolved edge, uncertain dependence, bridge-like anchor.",
  "",
  "Avoid: causal graph, true biological network, validated interaction network, diagnostic network.",
  "",
  "## Scientific status",
  "",
  "The current graph is a development artefact. It is not a final biological interaction network and has not been scientifically validated.",
  "",
  "## Output paths",
  "",
  paste0("- Edges: `", edges_out, "`"),
  paste0("- Nodes: `", nodes_out, "`"),
  paste0("- GraphML: `", graphml_out, "`"),
  paste0("- Summary: `", summary_out, "`")
)
write_lines(report_lines, report_out)
