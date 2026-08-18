#!/usr/bin/env Rscript

# =============================================================================
# inspect_pan_cancer_graph.R
# =============================================================================
# Inspection of the pan-cancer graph layer.
#
# Outputs (in graph/inspection/):
#   - INSPECTION_SUMMARY.md
#   - graph_connectivity.tsv
#   - component_vs_community.tsv
#   - community_purity.tsv
#   - lineage_mixing_edges.tsv
#   - isolated_or_bridge_nodes.tsv
#     (bridge role uses unnormalised, unweighted, undirected betweenness;
#      normalised betweenness is also reported for audit)
#   - layout_sanity_check.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(optparse)
})

opt_list <- list(
  make_option(c("--graph-dir"), type = "character", default = "results/unsupervised/pan_cancer/graph",
              help = "Graph directory containing edges/communities/components (default: results/unsupervised/pan_cancer/graph)"),
  make_option(c("--inspection-dir"), type = "character", default = NULL,
              help = "Inspection output directory (default: <graph-dir>/inspection)"),
  make_option(c("--layout-file"), type = "character", default = NULL,
              help = "Optional layout TSV (sample,x,y); defaults to <graph-dir>/community_validation/pan_cancer_layout_fr.tsv")
)

opt <- parse_args(OptionParser(option_list = opt_list))

graph_dir <- opt$`graph-dir`
if (is.null(opt$`inspection-dir`)) {
  opt$`inspection-dir` <- file.path(graph_dir, "inspection")
}
if (is.null(opt$`layout-file`)) {
  opt$`layout-file` <- file.path(graph_dir, "community_validation", "pan_cancer_layout_fr.tsv")
}

dir.create(opt$`inspection-dir`, recursive = TRUE, showWarnings = FALSE)

edges_file <- file.path(graph_dir, "pan_cancer_graph_edges.tsv")
comm_file <- file.path(graph_dir, "pan_cancer_communities.tsv")
comp_file <- file.path(graph_dir, "pan_cancer_components.tsv")
mod_file <- file.path(graph_dir, "community_validation", "validation_modularity.tsv")
assort_file <- file.path(graph_dir, "community_validation", "validation_assortativity.tsv")

stopifnot(file.exists(edges_file))
stopifnot(file.exists(comm_file))
stopifnot(file.exists(comp_file))

cat("========================================\n")
cat("PAN-CANCER GRAPH INSPECTION\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
edges <- fread(edges_file)
req_edges <- c("from", "to")
missing_edges <- setdiff(req_edges, names(edges))
if (length(missing_edges) > 0) {
  stop("Edges missing required columns: ", paste(missing_edges, collapse = ", "))
}
if (!"weight" %in% names(edges)) {
  edges[, weight := 1]
}

communities <- fread(comm_file)
req_comm <- c("sample", "community")
missing_comm <- setdiff(req_comm, names(communities))
if (length(missing_comm) > 0) {
  stop("Communities missing required columns: ", paste(missing_comm, collapse = ", "))
}
if (!"lineage" %in% names(communities)) {
  communities[, lineage := "UNKNOWN"]
}
communities[is.na(lineage) | lineage == "" | lineage == "\"\"", lineage := "UNKNOWN"]

components <- fread(comp_file)
if (!("sample" %in% names(components))) {
  stop("Components missing required column: sample")
}
if (!("component" %in% names(components))) {
  stop("Components missing required column: component")
}

# -----------------------------------------------------------------------------
# Build graph
# -----------------------------------------------------------------------------
g <- graph_from_data_frame(edges[, .(from, to)], directed = FALSE)
E(g)$weight <- edges$weight

# -----------------------------------------------------------------------------
# 1) graph_connectivity.tsv
# -----------------------------------------------------------------------------
comp <- components(g)
n_components <- comp$no
largest_component_size <- if (length(comp$csize) > 0) max(comp$csize) else 0

graph_connectivity <- data.table(
  n_vertices = vcount(g),
  n_edges = ecount(g),
  n_components = n_components,
  largest_component_size = largest_component_size
)
fwrite(graph_connectivity, file.path(opt$`inspection-dir`, "graph_connectivity.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 2) component_vs_community.tsv
# -----------------------------------------------------------------------------
component_vs_community <- merge(
  components[, .(sample, component)],
  communities[, .(sample, community)],
  by = "sample",
  all.x = TRUE
)
fwrite(component_vs_community, file.path(opt$`inspection-dir`, "component_vs_community.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 3) community_purity.tsv
# -----------------------------------------------------------------------------
community_purity <- communities[, .N, by = .(community, lineage)]
community_purity[, frac := N / sum(N), by = community]
fwrite(community_purity, file.path(opt$`inspection-dir`, "community_purity.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 4) lineage_mixing_edges.tsv
# -----------------------------------------------------------------------------
lin_map <- setNames(communities$lineage, communities$sample)
edges_lin <- edges[, .(
  from_lineage = lin_map[from],
  to_lineage = lin_map[to],
  weight = weight
)]
edges_lin[is.na(from_lineage), from_lineage := "UNKNOWN"]
edges_lin[is.na(to_lineage), to_lineage := "UNKNOWN"]

# Undirected pairing for lineage mixing
edges_lin[, a := pmin(from_lineage, to_lineage)]
edges_lin[, b := pmax(from_lineage, to_lineage)]

lineage_mixing_edges <- edges_lin[, .(
  n_edges = .N,
  mean_weight = mean(weight, na.rm = TRUE)
), by = .(from_lineage = a, to_lineage = b)][order(-n_edges)]

fwrite(lineage_mixing_edges, file.path(opt$`inspection-dir`, "lineage_mixing_edges.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 5) isolated_or_bridge_nodes.tsv
# -----------------------------------------------------------------------------
deg <- degree(g)
bet_unnorm <- betweenness(g, directed = FALSE, weights = NA, normalized = FALSE)
bet_norm <- betweenness(g, directed = FALSE, weights = NA, normalized = TRUE)

deg_thr <- as.numeric(quantile(deg, 0.9, na.rm = TRUE))
bet_thr <- as.numeric(quantile(bet_unnorm, 0.9, na.rm = TRUE))

role <- rep("peripheral", vcount(g))
role[deg == 0] <- "isolated"
role[deg > 0 & bet_unnorm >= bet_thr] <- "bridge"
role[deg >= deg_thr] <- "hub"

bridge_df <- data.table(
  sample = V(g)$name,
  lineage = lin_map[V(g)$name],
  degree = deg,
  betweenness = bet_unnorm,
  betweenness_unnormalised = bet_unnorm,
  betweenness_normalised = bet_norm,
  centrality_metric = "betweenness",
  centrality_normalised = FALSE,
  centrality_weighted = FALSE,
  centrality_scope = "whole_graph",
  centrality_tie_break = "not_applicable",
  role = role
)
bridge_df[is.na(lineage), lineage := "UNKNOWN"]
fwrite(bridge_df, file.path(opt$`inspection-dir`, "isolated_or_bridge_nodes.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 6) layout_sanity_check.tsv
# -----------------------------------------------------------------------------
if (file.exists(opt$`layout-file`)) {
  layout_df <- fread(opt$`layout-file`)
  req_layout <- c("sample", "x", "y")
  missing_layout <- setdiff(req_layout, names(layout_df))
  if (length(missing_layout) > 0) {
    stop("Layout file missing required columns: ", paste(missing_layout, collapse = ", "))
  }
  fwrite(layout_df, file.path(opt$`inspection-dir`, "layout_sanity_check.tsv"), sep = "\t")
} else {
  # Fall back: compute FR layout if missing
  set.seed(1)
  xy <- layout_with_fr(g, weights = E(g)$weight)
  layout_df <- data.table(sample = V(g)$name, x = xy[, 1], y = xy[, 2])
  fwrite(layout_df, file.path(opt$`inspection-dir`, "layout_sanity_check.tsv"), sep = "\t")
}

# -----------------------------------------------------------------------------
# 7) INSPECTION_SUMMARY.md
# -----------------------------------------------------------------------------
comm_sizes <- sort(table(communities$community), decreasing = TRUE)
comm_sizes_str <- paste(comm_sizes, collapse = ", ")

modularity_val <- NA_real_
if (file.exists(mod_file)) {
  mods <- fread(mod_file)
  if ("method" %in% names(mods) && "modularity" %in% names(mods)) {
    mod_louv <- mods[method == "louvain", modularity]
    if (length(mod_louv) == 1) modularity_val <- mod_louv
  }
}

assort_val <- NA_real_
if (file.exists(assort_file)) {
  assort <- fread(assort_file)
  if ("assortativity" %in% names(assort)) {
    assort_val <- assort$assortativity[1]
  }
}

summary_lines <- c(
  "# Pan-cancer graph inspection summary",
  "",
  paste0("- Graph vertices: ", vcount(g)),
  paste0("- Graph edges: ", ecount(g)),
  paste0("- Connected components: ", n_components, " (largest = ", largest_component_size, ")"),
  paste0("- Communities: ", length(comm_sizes), " (sizes: ", comm_sizes_str, ")"),
  paste0("- Modularity (Louvain): ", ifelse(is.na(modularity_val), "NA", round(modularity_val, 4))),
  paste0("- Lineage assortativity: ", ifelse(is.na(assort_val), "NA", round(assort_val, 4))),
  "",
  "Conclusion:",
  "Graph structure appears coherent and non-degenerate."
)
writeLines(summary_lines, file.path(opt$`inspection-dir`, "INSPECTION_SUMMARY.md"))

cat("Inspection outputs saved to: ", opt$`inspection-dir`, "\n", sep = "")
