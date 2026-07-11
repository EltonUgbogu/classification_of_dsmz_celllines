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
#   - The input RDS lacks $expr or $meta
#   - Metadata lacks sample_id, type, or cancer_type annotations
#   - Gene count != expected_genes (when expected_genes > 0)
#   - Node count != expected_nodes (when expected_nodes > 0)
#   - Any node lacks cancer_type annotation
#   - Tumour samples remain after filtering
#   - No graph edges are selected under the configured k and correlation values
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
    help="Expected node count; 0 = skip check (default: 0)"),
  make_option("--allow-empty-graph", action="store_true", default=FALSE,
    help="Allow an empty graph edge table when no valid neighbours are selected")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt[["expr-rds"]]) || is.null(opt[["output-dir"]])) {
  stop("--expr-rds and --output-dir are required")
}
if (tolower(opt[["correlation"]]) != "spearman") {
  stop("Methods require --correlation spearman; got: ", opt[["correlation"]])
}
if (is.na(opt$k) || opt$k < 1) {
  stop("--k must be an integer >= 1")
}
dir.create(opt[["output-dir"]], recursive=TRUE, showWarnings=FALSE)

# ---------------------------------------------------------------------------
# 1. Load expression matrix
# ---------------------------------------------------------------------------
cat("[1] Loading expression matrix:", opt[["expr-rds"]], "\n")
obj <- readRDS(opt[["expr-rds"]])
if (!"expr" %in% names(obj)) {
  stop("Input RDS object must contain $expr")
}

if (!"meta" %in% names(obj)) {
  stop("Input RDS object must contain $meta")
}
expr <- obj$expr
meta <- as.data.table(obj$meta)
cat("  Full matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")

if (is.null(colnames(expr))) {
  stop("$expr must have sample IDs in colnames")
}

if (!"sample_id" %in% names(meta)) {
  stop("$meta must contain a sample_id column")
}

if (anyDuplicated(meta$sample_id)) {
  stop("Duplicate sample_id values found in metadata")
}

if (!"type" %in% names(meta)) {
  stop("$meta must contain a type column")
}

if (!"cancer_type" %in% names(meta)) {
  if ("lineage" %in% names(meta)) {
    warning("Using legacy metadata column 'lineage' to populate 'cancer_type'")
    meta[, cancer_type := lineage]
  } else {
    stop("$meta must contain cancer_type. Legacy lineage column was also not found.")
  }
}

if (opt[["expected-genes"]] > 0 && nrow(expr) != opt[["expected-genes"]]) {
  stop("FAIL: Expected ", opt[["expected-genes"]], " genes but found ", nrow(expr))
}

# ---------------------------------------------------------------------------
# 2. Filter to cell lines only
# ---------------------------------------------------------------------------
cat("[2] Filtering to cell lines (type == cell_line)...\n")
cl_meta <- meta[type == "cell_line"]
tumour_count <- nrow(meta[type == "tumour"])
cat("  Tumour samples found and excluded:", tumour_count, "\n")

if (nrow(cl_meta) < 2) {
  stop("At least two cell-line samples are required to build a similarity graph")
}

missing_expr <- setdiff(cl_meta$sample_id, colnames(expr))
if (length(missing_expr) > 0) {
  stop("Cell-line sample IDs missing from expression matrix: ",
       paste(head(missing_expr, 20), collapse=", "))
}

expr_cl <- expr[, cl_meta$sample_id, drop=FALSE]
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

# The script validates that cell-line metadata contain cancer-type annotations.
missing_ct <- sum(is.na(cl_meta$cancer_type) | cl_meta$cancer_type == "" | cl_meta$cancer_type == "UNKNOWN")
if (missing_ct > 0) {
  stop("FAIL: ", missing_ct, " cell lines have missing/UNKNOWN cancer_type annotation")
}

cat("  Cancer-type distribution:\n")
print(table(cl_meta$cancer_type))

# ---------------------------------------------------------------------------
# 3. Compute Spearman correlation
# ---------------------------------------------------------------------------
cat("[3] Computing Spearman correlation matrix (", ncol(expr_cl), "x", ncol(expr_cl), ")...\n")
cor_mat <- cor(expr_cl, method=tolower(opt[["correlation"]]), use="pairwise.complete.obs")
if (all(is.na(cor_mat[upper.tri(cor_mat)]))) {
  stop("All pairwise Spearman correlations are NA; check for constant profiles or invalid expression matrix")
}
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
  nn_dt <- data.table(
    id = names(cv),
    similarity = as.numeric(cv)
  )
  nn_dt <- nn_dt[!is.na(similarity)]
  setorder(nn_dt, -similarity, id)

  n_keep <- min(opt$k, nrow(nn_dt))
  if (n_keep == 0) next
  nn_ids <- nn_dt$id[seq_len(n_keep)]

  for (nb in nn_ids) {
    edges_list[[edge_id]] <- data.table(from=s, to=nb, weight=cor_mat[s, nb])
    edge_id <- edge_id + 1L
  }
}

if (length(edges_list) > 0) {
  edges <- rbindlist(edges_list)
  edges[, `:=`(from_u=pmin(from, to), to_u=pmax(from, to))]
  edges <- edges[, .(weight=weight[1]), by=.(from=from_u, to=to_u)]
  setorder(edges, from, to)
} else {
  edges <- data.table(from=character(), to=character(), weight=numeric())
}

if (nrow(edges) == 0 && !isTRUE(opt[["allow-empty-graph"]])) {
  stop("No graph edges were selected; check expression matrix, correlations, and --k")
}

cat("  Edges (undirected):", nrow(edges), "\n")
cat("  Nodes:", length(sample_ids), "\n")

# ---------------------------------------------------------------------------
# 5. Write outputs
# ---------------------------------------------------------------------------
edges_path <- file.path(opt[["output-dir"]], "pan_cancer_cell_line_graph_edges.tsv")
meta_path  <- file.path(opt[["output-dir"]], "pan_cancer_cell_line_node_metadata.tsv")

fwrite(edges, edges_path, sep="\t")
cat("[5] Edges written to:", edges_path, "\n")

required_cols <- c("sample_id", "cancer_type", "type", "source_profile")
missing_cols <- setdiff(required_cols, names(cl_meta))

if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse=", "))
}

meta_out <- cl_meta[, .(sample_id, cancer_type, type, source_profile)]
fwrite(meta_out, meta_path, sep="\t")
cat("    Metadata written to:", meta_path, "\n")

cat("[OK] Cell-line similarity graph built successfully\n")
cat("     Nodes:", length(sample_ids), "\n")
cat("     Edges:", nrow(edges), "\n")
