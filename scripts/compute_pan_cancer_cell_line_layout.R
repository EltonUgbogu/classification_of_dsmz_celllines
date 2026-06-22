#!/usr/bin/env Rscript
# =============================================================================
# compute_pan_cancer_cell_line_layout.R
# =============================================================================
#
# Computes Fruchterman-Reingold layout for the cell-line-only similarity graph.
# Uses a fixed seed for reproducibility.
#
# Output:
#   pan_cancer_cell_line_layout.tsv  (columns: sample, x, y)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(optparse)
})

option_list <- list(
  make_option("--edges",  type="character", default=NULL,
    help="[REQUIRED] Edge list TSV (from, to, weight)"),
  make_option("--out",    type="character", default=NULL,
    help="[REQUIRED] Output layout TSV path"),
  make_option("--seed",   type="integer",   default=42,
    help="Random seed (default: 42)"),
  make_option("--layout", type="character", default="fr",
    help="Layout algorithm: fr (Fruchterman-Reingold) or kk (Kamada-Kawai) (default: fr)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$edges) || is.null(opt$out)) {
  stop("--edges and --out are required")
}
dir.create(dirname(opt$out), recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

# ---------------------------------------------------------------------------
# Load data and build graph
# ---------------------------------------------------------------------------
cat("[1] Loading edges:", opt$edges, "\n")
edges <- fread(opt$edges)
E_und <- edges[, .(weight=max(weight, na.rm=TRUE)), by=.(a=pmin(from,to), b=pmax(from,to))]
setnames(E_und, c("a","b"), c("from","to"))

g <- graph_from_data_frame(E_und, directed=FALSE)
E(g)$weight <- E_und$weight
cat("  Nodes:", vcount(g), "  Edges:", ecount(g), "\n")

# ---------------------------------------------------------------------------
# Compute layout
# ---------------------------------------------------------------------------
cat("[2] Computing layout (method=", opt$layout, ", seed=", opt$seed, ")...\n")
set.seed(opt$seed)
if (opt$layout == "kk") {
  layout_m <- layout_with_kk(g, weights=E(g)$weight)
} else {
  layout_m <- layout_with_fr(g, weights=E(g)$weight)
}

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
layout_dt <- data.table(sample=V(g)$name, x=layout_m[,1], y=layout_m[,2])
fwrite(layout_dt, opt$out, sep="\t")
cat("[OK] Layout written to:", opt$out, "\n")
cat("  X range: [", round(min(layout_dt$x),3), ",", round(max(layout_dt$x),3), "]\n")
cat("  Y range: [", round(min(layout_dt$y),3), ",", round(max(layout_dt$y),3), "]\n")
