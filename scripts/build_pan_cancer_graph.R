#!/usr/bin/env Rscript
# =============================================================================
# build_pan_cancer_graph.R
# =============================================================================
#
# Build a simplified k-nearest-neighbour graph for pan-cancer data.
#
# This script builds a simplified k-nearest-neighbour graph from a sample-sample
# correlation matrix for pan-cancer feature-space inspection. Depending on the
# input matrix, it may support auxiliary negative-control inspection, but it is
# not the primary patient-referenced support-threshold consensus graph.
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
# Connected components via depth-first search
# -----------------------------------------------------------------------------
# igraph is intentionally avoided here to keep dependencies minimal for this
# graph-inspection script. The function labels every supplied node, including
# isolates, after constructing an undirected adjacency list.
get_components <- function(edges, nodes) {
  nodes <- unique(as.character(nodes))
  adj_list <- setNames(vector("list", length(nodes)), nodes)

  if (!is.null(edges) && nrow(edges) > 0) {
    edges[, from := as.character(from)]
    edges[, to := as.character(to)]

    for (i in seq_len(nrow(edges))) {
      a <- edges$from[i]
      b <- edges$to[i]

      if (!a %in% names(adj_list)) adj_list[[a]] <- character(0)
      if (!b %in% names(adj_list)) adj_list[[b]] <- character(0)

      adj_list[[a]] <- unique(c(adj_list[[a]], b))
      adj_list[[b]] <- unique(c(adj_list[[b]], a))
    }
  }

  visited <- setNames(rep(FALSE, length(adj_list)), names(adj_list))
  comp_id <- integer(length(adj_list))
  names(comp_id) <- names(adj_list)

  current_comp <- 0L

  for (node in names(adj_list)) {
    if (!visited[[node]]) {
      current_comp <- current_comp + 1L
      stack <- node

      while (length(stack) > 0) {
        current <- stack[[length(stack)]]
        stack <- stack[-length(stack)]

        if (!visited[[current]]) {
          visited[[current]] <- TRUE
          comp_id[[current]] <- current_comp

          neighbours <- adj_list[[current]]
          neighbours <- neighbours[neighbours %in% names(visited)]
          unvisited <- neighbours[!visited[neighbours]]

          stack <- c(stack, unvisited)
        }
      }
    }
  }

  data.table(
    sample = names(comp_id),
    component = as.integer(comp_id)
  )
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
if (is.na(opt$k) || opt$k < 1) {
  stop("--k must be an integer >= 1")
}

# Create output directory (including any missing parents) if it doesn't exist
dir.create(opt$`output-dir`, recursive=TRUE, showWarnings=FALSE)

# -----------------------------------------------------------------------------
# Load correlation matrix
# -----------------------------------------------------------------------------
cat("Loading correlation matrix...\n")
cor_mat <- readRDS(opt$`cor-matrix`)
if (!is.matrix(cor_mat) && !is.data.frame(cor_mat)) {
  stop("Correlation input must be a matrix-like object")
}

cor_mat <- as.matrix(cor_mat)

if (nrow(cor_mat) != ncol(cor_mat)) {
  stop("Correlation matrix must be square")
}

if (is.null(rownames(cor_mat)) || is.null(colnames(cor_mat))) {
  stop("Correlation matrix must have rownames and colnames")
}

if (!setequal(rownames(cor_mat), colnames(cor_mat))) {
  stop("Correlation matrix rownames and colnames do not match")
}

cor_mat <- cor_mat[rownames(cor_mat), rownames(cor_mat), drop=FALSE]

if (nrow(cor_mat) < 2) {
  stop("At least two samples are required to build a graph")
}
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

  # Drop missing and below-threshold neighbours, then use sample ID as the
  # deterministic tie-breaker for equal correlations.
  nn_dt <- data.table(
    id = names(cor_vec),
    similarity = as.numeric(cor_vec)
  )
  nn_dt <- nn_dt[!is.na(similarity) & similarity >= opt$`min-cor`]
  setorder(nn_dt, -similarity, id)

  n_keep <- min(opt$k, nrow(nn_dt))
  if (n_keep == 0) next
  nn_ids <- nn_dt$id[seq_len(n_keep)]

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

if (length(edges_list) > 0) {
  edges <- rbindlist(edges_list)
} else {
  edges <- data.table(from=character(), to=character(), weight=numeric(), direction=character())
}

# Deduplicate: since the loop adds both A→B and B→A, collapse to a canonical
# undirected edge by sorting node names and keeping unique pairs
if (nrow(edges) > 0) {
  edges[, key := paste(pmin(from, to), pmax(from, to), sep="|")]
  edges <- unique(edges, by="key")
  edges[, key := NULL]
}

cat("  Created ", nrow(edges), " edges\n", sep="")

# -----------------------------------------------------------------------------
# Find connected components
# -----------------------------------------------------------------------------
cat("Finding connected components...\n")
comp_df <- get_components(edges, nodes=colnames(cor_mat))
comp_sizes <- table(comp_df$component)

# Assemble a per-node data.table with component ID and the size of that component
comp_df[, comp_size := as.integer(comp_sizes[as.character(component)])]

cat("  Graph: ", nrow(cor_mat), " vertices, ", nrow(edges), " edges\n", sep="")
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
