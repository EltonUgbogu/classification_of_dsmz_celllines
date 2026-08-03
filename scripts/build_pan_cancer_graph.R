#!/usr/bin/env Rscript
# =============================================================================
# build_pan_cancer_graph.R
# =============================================================================
#
# Build validation kNN graph for pan-cancer data.
#
# NOTE: This is a simplified kNN graph for validation purposes (negative
# control testing). For production analysis, this should be replaced with
# the intersection consensus policy (best_overall_dir ∩ winner_dir) used in
# disease-specific analysis.
#
# This script computes kNN graphs from correlation matrix to test hematologic
# rejection. The graph structure is used to verify that HEME samples form
# separate components and do not mix with solid cancers.
#
# USAGE:
# Rscript build_pan_cancer_graph.R \
#   --cor-matrix results/unsupervised/pan_cancer/pan_cancer_cor.rds \
#   --k 20 \
#   --output-dir results/unsupervised/pan_cancer/graph
#
# =============================================================================

# Suppress package startup messages for cleaner console output
suppressPackageStartupMessages({
  library(data.table)  # Fast data manipulation and file I/O
  library(optparse)    # Command-line argument parsing
})

# -----------------------------------------------------------------------------
# Connected components via depth-first search (DFS)
# -----------------------------------------------------------------------------
# igraph is intentionally avoided here to keep dependencies minimal for this
# validation script. The function implements a standard DFS-based labelling:
# each node gets an integer component ID; nodes sharing an ID are reachable
# from one another through undirected edges.
find_components <- function(edges) {
  
  # Collect all unique node names from both ends of every edge
  nodes <- unique(c(edges$from, edges$to))
  
  # Build a forward adjacency list (from → list of to-nodes)
  adj_list <- split(edges$to, edges$from)
  
  # Extend adjacency list with reverse edges to make the graph undirected
  # (i.e. if A→B exists, also record B→A)
  for (to_node in unique(edges$to)) {
    from_nodes <- edges[to == to_node, from]
    for (from_node in from_nodes) {
      if (!from_node %in% names(adj_list)) {
        adj_list[[from_node]] <- character(0)
      }
      if (!to_node %in% adj_list[[from_node]]) {
        adj_list[[from_node]] <- c(adj_list[[from_node]], to_node)
      }
    }
  }
  
  # Initialise visited flags and component assignment vectors, keyed by node name
  visited    <- setNames(rep(FALSE, length(nodes)), nodes)
  components <- setNames(rep(0L,    length(nodes)), nodes)
  comp_id    <- 1L
  
  # Recursive DFS: marks all nodes reachable from `node` with `comp` ID.
  # <<- is used to modify visited/components in the enclosing environment.
  dfs <- function(node, comp) {
    if (visited[node]) return          # Skip already-visited nodes
    visited[node]    <<- TRUE
    components[node] <<- comp
    neighbors <- adj_list[[node]]
    if (length(neighbors) > 0) {
      for (neighbor in neighbors) {
        if (!visited[neighbor]) {
          dfs(neighbor, comp)          # Recurse into unvisited neighbours
        }
      }
    }
  }
  
  # Iterate over all nodes; start a new DFS (and increment component ID)
  # whenever an unvisited node is encountered
  for (node in nodes) {
    if (!visited[node]) {
      dfs(node, comp_id)
      comp_id <- comp_id + 1L
    }
  }
  
  return(components)  # Named integer vector: node → component ID
}

# -----------------------------------------------------------------------------
# Command-line options
# -----------------------------------------------------------------------------
option_list <- list(
  make_option(c("--cor-matrix"), type="character", default=NULL,
              help="Path to correlation matrix RDS file"),
  make_option(c("--k"), type="integer", default=20,
              help="Number of nearest neighbours for kNN (default: 20)"),
  make_option(c("--output-dir"), type="character", default=NULL,
              help="Output directory for graph files"),
  make_option(c("--min-cor"), type="double", default=0.0,
              help="Minimum correlation threshold for edges (default: 0.0)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Both --cor-matrix and --output-dir are mandatory
if (is.null(opt$`cor-matrix`) || is.null(opt$`output-dir`)) {
  stop("--cor-matrix and --output-dir are required")
}

# Create output directory (including any missing parents) if it doesn't exist
dir.create(opt$`output-dir`, recursive=TRUE, showWarnings=FALSE)

# -----------------------------------------------------------------------------
# Load correlation matrix
# -----------------------------------------------------------------------------
cat("Loading correlation matrix...\n")
cor_mat <- readRDS(opt$`cor-matrix`)
cat("  Correlation matrix: ", nrow(cor_mat), " x ", ncol(cor_mat), "\n", sep="")

# -----------------------------------------------------------------------------
# Build kNN graph
# -----------------------------------------------------------------------------
# For each sample, find its k most correlated neighbours (excluding self).
# Each qualifying pair becomes an undirected weighted edge.
cat("Building kNN graph (k=", opt$k, ")...\n", sep="")

edges_list <- list()
edge_id    <- 1

for (sample in colnames(cor_mat)) {
  
  # Extract row of correlations and remove the self-correlation entry
  cor_vec <- cor_mat[sample, ]
  cor_vec <- cor_vec[names(cor_vec) != sample]
  
  # Drop any neighbours below the minimum correlation threshold
  cor_vec <- cor_vec[cor_vec >= opt$`min-cor`]
  
  # Sort descending and take the top k (or fewer if not enough neighbours pass the filter)
  nn_sorted <- sort(cor_vec, decreasing=TRUE)
  nn_ids    <- names(nn_sorted)[1:min(opt$k, length(nn_sorted))]
  
  # Record one edge per selected neighbour
  for (neighbor in nn_ids) {
    edges_list[[edge_id]] <- data.table(
      from      = sample,
      to        = neighbor,
      weight    = cor_mat[sample, neighbor],
      direction = "both"   # Flagged as undirected for downstream consumers
    )
    edge_id <- edge_id + 1
  }
}

edges <- rbindlist(edges_list)

# Deduplicate: since the loop adds both A→B and B→A, collapse to a canonical
# undirected edge by sorting node names and keeping unique pairs
edges[, key := paste(pmin(from, to), pmax(from, to), sep="|")]
edges <- unique(edges, by="key")
edges[, key := NULL]

cat("  Created ", nrow(edges), " edges\n", sep="")

# -----------------------------------------------------------------------------
# Find connected components
# -----------------------------------------------------------------------------
cat("Finding connected components...\n")
comp_membership <- find_components(edges)
comp_sizes      <- table(comp_membership)

# Assemble a per-node data.table with component ID and the size of that component
comp_df <- data.table(
  sample    = names(comp_membership),
  component = as.integer(comp_membership),
  comp_size = as.integer(comp_sizes[as.character(comp_membership)])
)

cat("  Graph: ", length(unique(c(edges$from, edges$to))), " vertices, ", nrow(edges), " edges\n", sep="")
cat("  Number of components: ", length(unique(comp_df$component)), "\n", sep="")
cat("  Component sizes:\n")
print(table(comp_df$comp_size))  # Distribution of how many nodes each component contains

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------
cat("Saving graph outputs...\n")

# Edge list — consumed by plot_pan_cancer_graph.R via --edges
edges_file <- file.path(opt$`output-dir`, "pan_cancer_graph_edges.tsv")
fwrite(edges, file=edges_file, sep="\t")
cat("  Edges saved to:", edges_file, "\n")

# Component assignments — consumed by plot_pan_cancer_graph.R via --components
comp_file <- file.path(opt$`output-dir`, "pan_cancer_components.tsv")
fwrite(comp_df, file=comp_file, sep="\t")
cat("  Components saved to:", comp_file, "\n")

# Note: no igraph object is saved; igraph is not required for this validation script

cat("[OK] Pan-cancer graph built successfully\n")
