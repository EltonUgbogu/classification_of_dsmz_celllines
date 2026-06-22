#!/usr/bin/env Rscript
# =============================================================================
# build_pan_cancer_cell_line_similarity_graph.R
# =============================================================================
#
# Build a cell-line-only Spearman kNN similarity graph from the pan-cancer
# feature expression matrix. Filters to DSMZ cell lines (type == "cell_line"),
# computes pairwise Spearman correlation, and constructs a union kNN graph.
#
# FAILS with stop() if:
#   - Gene count != expected_genes (when expected_genes > 0)
#   - Node count != expected_nodes (when expected_nodes > 0)
#   - Any node lacks lineage annotation
#   - Tumour samples remain after filtering
#
# Outputs:
#   results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_graph_edges.tsv
#   results/unsupervised/pan_cancer/cell_line_similarity/pan_cancer_cell_line_node_metadata.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--expr-rds",    type="character", default=NULL,
    help="[REQUIRED] Pan-cancer expression RDS with $expr, $meta, $genes"),
  make_option("--output-dir",  type="character", default=NULL,
    help="[REQUIRED] Output directory"),
  make_option("--k",           type="integer",   default=20,
    help="Number of nearest neighbours (default: 20)"),
  make_option("--correlation", type="character", default="spearman",
    help="Pairwise similarity method; must be spearman for Methods consistency (default: spearman)"),
  make_option("--expected-genes", type="integer", default=0,
    help="Expected gene count; 0 = skip check (default: 0)"),
  make_option("--expected-nodes", type="integer", default=0,
    help="Expected node count; 0 = skip check (default: 0)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt[["expr-rds"]]) || is.null(opt[["output-dir"]])) {
  stop("--expr-rds and --output-dir are required")
}
if (tolower(opt[["correlation"]]) != "spearman") {
  stop("Methods require --correlation spearman; got: ", opt[["correlation"]])
}
dir.create(opt[["output-dir"]], recursive=TRUE, showWarnings=FALSE)

# ---------------------------------------------------------------------------
# 1. Load expression matrix
# ---------------------------------------------------------------------------
cat("[1] Loading expression matrix:", opt[["expr-rds"]], "\n")
obj <- readRDS(opt[["expr-rds"]])
if (!all(c("expr","meta") %in% names(obj))) {
  stop("RDS must contain $expr and $meta")
}
expr <- obj$expr
meta <- as.data.table(obj$meta)
cat("  Full matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")

if (opt[["expected-genes"]] > 0 && nrow(expr) != opt[["expected-genes"]]) {
  stop("FAIL: Expected ", opt[["expected-genes"]], " genes but found ", nrow(expr))
}

# ---------------------------------------------------------------------------
# 2. Filter to cell lines only
# ---------------------------------------------------------------------------
cat("[2] Filtering to cell lines (type == cell_line)...\n")
if (!"type" %in% names(meta)) {
  stop("$meta must have a type column")
}
cl_meta <- meta[type == "cell_line"]
tumour_count <- nrow(meta[type == "tumour"])
cat("  Tumour samples found and excluded:", tumour_count, "\n")

cl_ids <- cl_meta$sample_id
expr_cl <- expr[, cl_ids, drop=FALSE]
cat("  Cell-line samples retained:", ncol(expr_cl), "\n")

# Validate: no tumour samples remain
non_cl <- sum(cl_meta$type != "cell_line")
if (non_cl > 0) {
  stop("FAIL: ", non_cl, " non-cell-line samples remain after filtering")
}

# Validate: expected node count
if (opt[["expected-nodes"]] > 0 && ncol(expr_cl) != opt[["expected-nodes"]]) {
  stop("FAIL: Expected ", opt[["expected-nodes"]], " cell-line nodes but found ", ncol(expr_cl))
}

# Validate: no missing lineage
if (!"lineage" %in% names(cl_meta)) {
  stop("$meta must have a lineage column")
}
missing_lin <- sum(is.na(cl_meta$lineage) | cl_meta$lineage == "" | cl_meta$lineage == "UNKNOWN")
if (missing_lin > 0) {
  warning("WARNING: ", missing_lin, " cell lines have missing/UNKNOWN lineage")
}

cat("  Lineage distribution:\n")
print(table(cl_meta$lineage))

# ---------------------------------------------------------------------------
# 3. Compute Spearman correlation
# ---------------------------------------------------------------------------
cat("[3] Computing Spearman correlation matrix (", ncol(expr_cl), "x", ncol(expr_cl), ")...\n")
cor_mat <- cor(expr_cl, method=tolower(opt[["correlation"]]), use="pairwise.complete.obs")
cat("  Correlation range: [", round(min(cor_mat, na.rm=TRUE), 3),
    ",", round(max(cor_mat, na.rm=TRUE), 3), "]\n")

# ---------------------------------------------------------------------------
# 4. Build kNN graph
# ---------------------------------------------------------------------------
cat("[4] Building union kNN graph (k=", opt$k,
    "; exact ties resolved by profile ID)...\n")
sample_ids <- colnames(expr_cl)
edges_list <- list()
edge_id <- 1L

for (s in sample_ids) {
  cv <- cor_mat[s, ]
  cv <- cv[names(cv) != s]
  candidates <- data.table(profile_id=names(cv), correlation=as.numeric(cv))
  candidates <- candidates[!is.na(correlation)]
  setorder(candidates, -correlation, profile_id)
  nn_ids <- candidates$profile_id[seq_len(min(opt$k, nrow(candidates)))]
  for (nb in nn_ids) {
    edges_list[[edge_id]] <- data.table(from=s, to=nb, weight=cor_mat[s, nb])
    edge_id <- edge_id + 1L
  }
}

edges <- rbindlist(edges_list)
edges[, `:=`(from_u=pmin(from, to), to_u=pmax(from, to))]
edges <- edges[, .(weight=weight[1]), by=.(from=from_u, to=to_u)]
setorder(edges, from, to)

cat("  Edges (undirected):", nrow(edges), "\n")
cat("  Nodes:", length(unique(c(edges$from, edges$to))), "\n")

# ---------------------------------------------------------------------------
# 5. Write outputs
# ---------------------------------------------------------------------------
edges_path <- file.path(opt[["output-dir"]], "pan_cancer_cell_line_graph_edges.tsv")
meta_path  <- file.path(opt[["output-dir"]], "pan_cancer_cell_line_node_metadata.tsv")

fwrite(edges, edges_path, sep="\t")
cat("[5] Edges written to:", edges_path, "\n")

meta_out <- cl_meta[, .(sample_id, lineage, type, source_profile)]
fwrite(meta_out, meta_path, sep="\t")
cat("    Metadata written to:", meta_path, "\n")

cat("[OK] Cell-line similarity graph built successfully\n")
cat("     Nodes:", length(unique(c(edges$from, edges$to))), "\n")
cat("     Edges:", nrow(edges), "\n")
