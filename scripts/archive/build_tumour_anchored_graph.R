#!/usr/bin/env Rscript
# =============================================================================
# build_tumour_anchored_graph.R
# =============================================================================
# Tumour-anchored cell line communities:
#   - Keeps only NG- cell lines and tumours from pan_cancer_expr.rds
#   - Cell↔tumour Spearman → cell↔cell similarity (correlation of tumour profiles)
#   - kNN weighted graph → Louvain communities → PDF + TSVs in graph/
#
# Usage (from .../pan_cancer/):
#   Rscript build_tumour_anchored_graph.R [--rds=FILE] [--outdir=DIR] [--k=15] [--min_edge=0.2] [--method=spearman] [--seed=1]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

# Optional args (override via command line if desired)
args <- commandArgs(trailingOnly = TRUE)
in_rds    <- "pan_cancer_expr.rds"
outdir    <- "graph"
k         <- 15L
min_edge  <- 0.20
cor_method <- "spearman"
seed      <- 1L
for (a in args) {
  if (grepl("^--rds=", a))      in_rds    <- sub("^--rds=", "", a)
  if (grepl("^--outdir=", a))   outdir    <- sub("^--outdir=", "", a)
  if (grepl("^--k=", a))        k         <- as.integer(sub("^--k=", "", a))
  if (grepl("^--min_edge=", a)) min_edge  <- as.numeric(sub("^--min_edge=", "", a))
  if (grepl("^--method=", a))   cor_method <- sub("^--method=", "", a)
  if (grepl("^--seed=", a))     seed      <- as.integer(sub("^--seed=", "", a))
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

x <- readRDS(in_rds)
meta <- as.data.table(x$meta)
expr <- as.matrix(x$expr)

stopifnot(all(c("sample_id", "type") %in% names(meta)))
stopifnot(!is.null(colnames(expr)))

# --- orient expr to genes x samples (should be 258 x N) ---
if (!is.null(rownames(expr)) && any(grepl("^NG-", rownames(expr))) && !any(grepl("^NG-", colnames(expr)))) {
  expr <- t(expr)
}
if (ncol(expr) == 258 && nrow(expr) != 258) expr <- t(expr)
if (nrow(expr) != 258) stop("Expected 258 genes in rows; got dim(expr) = ", paste(dim(expr), collapse = " x "))

# keep only samples present in expr
meta <- meta[sample_id %in% colnames(expr)]
expr <- expr[, meta$sample_id, drop = FALSE]

# --- split samples ---
cell_ids <- meta[type == "cell_line", sample_id]
cell_ids <- cell_ids[grepl("^NG-", cell_ids)]   # enforce NG- only
tum_ids  <- meta[type == "tumour", sample_id]

cat("Cell lines (NG-):", length(cell_ids), "\n")
cat("Tumours:", length(tum_ids), "\n")
if (length(cell_ids) < 3) stop("Too few NG- cell lines after filtering.")
if (length(tum_ids) < 10) stop("Too few tumours; cannot do tumour-anchored similarity.")

# --- build tumour-anchored similarities ---
# expr is genes x samples
Xc <- expr[, cell_ids, drop = FALSE]  # 258 x n_cells
Xt <- expr[, tum_ids,  drop = FALSE]  # 258 x n_tumours

stopifnot(nrow(Xc) == nrow(Xt))

# (1) cell ↔ tumour (columns are samples; rows are genes)
sim_ct <- cor(Xc, Xt, method = cor_method, use = "pairwise.complete.obs")
rownames(sim_ct) <- colnames(Xc)
colnames(sim_ct) <- colnames(Xt)
# sim_ct is cell_lines x tumours (e.g. 58 x 1128)

# (2) cell ↔ cell similarity via tumour profiles (correlate rows of sim_ct)
sim_cc <- cor(t(sim_ct), method = cor_method, use = "pairwise.complete.obs")
diag(sim_cc) <- NA_real_

# 3) build kNN graph (weighted, undirected)
k <- min(k, nrow(sim_cc) - 1)

edges <- list()
for (i in seq_len(nrow(sim_cc))) {
  v <- sim_cc[i, ]
  o <- order(v, decreasing = TRUE, na.last = NA)
  if (length(o) == 0) next
  o <- o[seq_len(min(k, length(o)))]
  keep <- o[v[o] >= min_edge]
  if (length(keep) == 0) next
  for (j in keep) {
    edges[[length(edges) + 1]] <- list(from = rownames(sim_cc)[i], to = rownames(sim_cc)[j], weight = v[j])
  }
}
edges_dt <- rbindlist(edges)
if (nrow(edges_dt) == 0) stop("No edges kept. Lower min_edge or raise k.")

# undirected unique
edges_dt[, a := pmin(from, to)]
edges_dt[, b := pmax(from, to)]
edges_dt <- unique(edges_dt[, .(from = a, to = b, weight = weight)])

g <- graph_from_data_frame(edges_dt, directed = FALSE)

# 4) communities (Louvain)
comm <- cluster_louvain(g, weights = E(g)$weight)
V(g)$community <- membership(comm)
V(g)$degree <- degree(g)

# attach lineage for later coloring/inspection
lin <- meta[, .(sample_id, lineage)]
lin <- lin[sample_id %in% V(g)$name]
V(g)$lineage <- lin$lineage[match(V(g)$name, lin$sample_id)]

cat("Graph nodes:", vcount(g), " edges:", ecount(g),
    " communities:", length(unique(V(g)$community)), "\n")

# 5) plot communities
set.seed(seed)
lay <- layout_with_fr(g)

pal <- grDevices::hcl.colors(length(unique(V(g)$community)), "Spectral")
vcol <- pal[as.integer(factor(V(g)$community))]

pdf(file.path(outdir, "Fig_pan_cancer_tumour_anchored_communities.pdf"), width = 11, height = 8.5)
par(mar = c(0, 0, 2, 0))
plot(
  g, layout = lay,
  vertex.color = vcol,
  vertex.size = pmax(3, pmin(10, 3 + sqrt(V(g)$degree))),
  vertex.label = NA,
  vertex.frame.color = "white",
  edge.color = adjustcolor("black", alpha.f = 0.12),
  edge.width = 0.6,
  main = sprintf("DSMZ cell line communities (tumour-anchored) | NG- only | n=%d, m=%d, k=%d",
                 vcount(g), ecount(g), length(unique(V(g)$community)))
)
dev.off()

# 6) save outputs
fwrite(edges_dt, file.path(outdir, "pan_cancer_tumour_anchored_cellline_edges.tsv"), sep = "\t")
fwrite(data.table(sample_id = V(g)$name,
                  lineage = V(g)$lineage,
                  community = V(g)$community,
                  degree = V(g)$degree),
       file.path(outdir, "pan_cancer_tumour_anchored_cellline_communities.tsv"), sep = "\t")
saveRDS(g, file.path(outdir, "pan_cancer_tumour_anchored_graph.rds"))

cat("Wrote:\n",
    " - ", file.path(outdir, "Fig_pan_cancer_tumour_anchored_communities.pdf"), "\n",
    " - ", file.path(outdir, "pan_cancer_tumour_anchored_cellline_edges.tsv"), "\n",
    " - ", file.path(outdir, "pan_cancer_tumour_anchored_cellline_communities.tsv"), "\n",
    " - ", file.path(outdir, "pan_cancer_tumour_anchored_graph.rds"), "\n", sep = "")
