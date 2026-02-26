#!/usr/bin/env Rscript
# =============================================================================
# plot_dsmz_communities.R
# =============================================================================
# DSMZ-only (NG- ↔ NG-) community graph from the mixed pan_cancer edges file.
# Filters edges to cell-line ↔ cell-line (NG- only), builds undirected graph,
# runs Louvain, plots communities, writes filtered edges + node table.
#
# Input:  pan_cancer_graph_edges.tsv (mixed; may contain non-NG- nodes)
# Output: Fig_pan_cancer_DSMZ_cellline_communities.pdf
#         pan_cancer_DSMZ_cellline_edges.tsv
#         pan_cancer_DSMZ_cellline_communities.tsv
#
# Usage (from .../pan_cancer/graph/ or with paths):
#   Rscript plot_dsmz_communities.R --edges pan_cancer_graph_edges.tsv [--outdir .]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(optparse)
})

option_list <- list(
  make_option(c("--edges"), type = "character", default = NULL,
              help = "Path to pan_cancer_graph_edges.tsv (mixed edges)"),
  make_option(c("--outdir"), type = "character", default = ".",
              help = "Output directory for PDF and TSVs"),
  make_option(c("--seed"), type = "integer", default = 1,
              help = "Random seed for layout")
)

opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$edges), file.exists(opt$edges))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

is_ng <- function(x) !is.na(x) & grepl("^NG-", x)

# -----------------------------------------------------------------------------
# Load edges and normalize column names
# -----------------------------------------------------------------------------
edges <- fread(opt$edges)
edge_cols <- names(edges)
from_col <- intersect(edge_cols, c("from", "source", "cell_line_1", "cell1", "u", "A"))[1]
to_col   <- intersect(edge_cols, c("to", "target", "cell_line_2", "cell2", "v", "B"))[1]
if (is.na(from_col) || is.na(to_col)) {
  stop("Could not infer edge columns. Expected from/to or cell_line_1/cell_line_2.")
}
setnames(edges, c(from_col, to_col), c("from", "to"))

# Optional weight column
weight_col <- intersect(edge_cols, c("weight", "w", "similarity", "cor"))[1]
if (!is.na(weight_col)) {
  setnames(edges, weight_col, "weight")
} else {
  edges[, weight := 1.0]
}

# Keep ONLY DSMZ cell-line ↔ cell-line edges (NG- ↔ NG-)
edges <- edges[is_ng(from) & is_ng(to)]
if (nrow(edges) == 0) stop("No NG-↔NG- edges found after filtering.")

# Build undirected graph (keep weights)
g <- graph_from_data_frame(edges[, .(from, to, weight)], directed = FALSE)

# Community detection (Louvain)
comm <- cluster_louvain(g, weights = E(g)$weight)
V(g)$community <- membership(comm)
V(g)$degree <- degree(g)

# Plot
set.seed(opt$seed)
lay <- layout_with_fr(g)
pal <- grDevices::hcl.colors(length(unique(V(g)$community)), "Spectral")
vcol <- pal[as.integer(factor(V(g)$community))]

pdf(file.path(opt$outdir, "Fig_pan_cancer_DSMZ_cellline_communities.pdf"), width = 11, height = 8.5)
par(mar = c(0, 0, 2, 0))
plot(
  g, layout = lay,
  vertex.color = vcol,
  vertex.size = pmax(3, pmin(10, 3 + sqrt(V(g)$degree))),
  vertex.label = NA,
  vertex.frame.color = "white",
  edge.color = adjustcolor("black", alpha.f = 0.12),
  edge.width = 0.6,
  main = sprintf("DSMZ cell line communities (NG- only) | n=%d, m=%d, k=%d",
                 vcount(g), ecount(g), length(unique(V(g)$community)))
)
dev.off()

# Save filtered edges + node communities
edges_out <- as.data.table(get.data.frame(g, what = "edges"))
fwrite(edges_out, file.path(opt$outdir, "pan_cancer_DSMZ_cellline_edges.tsv"), sep = "\t")

nodes_out <- data.table(
  sample_id = V(g)$name,
  community = V(g)$community,
  degree    = V(g)$degree
)
fwrite(nodes_out, file.path(opt$outdir, "pan_cancer_DSMZ_cellline_communities.tsv"), sep = "\t")

cat("DSMZ-only graph summary:\n")
cat("  Nodes:", vcount(g), "\n")
cat("  Edges:", ecount(g), "\n")
cat("  Communities:", length(unique(V(g)$community)), "\n")
cat("Wrote:\n")
cat("  - ", file.path(opt$outdir, "Fig_pan_cancer_DSMZ_cellline_communities.pdf"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "pan_cancer_DSMZ_cellline_edges.tsv"), "\n", sep = "")
cat("  - ", file.path(opt$outdir, "pan_cancer_DSMZ_cellline_communities.tsv"), "\n", sep = "")
