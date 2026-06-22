#!/usr/bin/env Rscript
# =============================================================================
# compute_pan_cancer_cell_line_validation.R
# =============================================================================
#
# Validates the cell-line-only similarity graph and community assignments:
#   - Node count
#   - Edge count
#   - Community count
#   - Weighted Newman-Girvan modularity for Louvain and Leiden
#   - Unweighted nominal lineage assortativity
#   - Tumour count
#   - Missing lineage count
#
# Outputs:
#   validation_modularity.tsv
#   validation_assortativity.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(optparse)
})

option_list <- list(
  make_option("--edges", type="character", default=NULL,
    help="[REQUIRED] Edge list TSV (from, to, weight)"),
  make_option("--communities", type="character", default=NULL,
    help="[REQUIRED] Louvain community TSV (sample, component, lineage, ...)"),
  make_option("--leiden-communities", type="character", default=NULL,
    help="[OPTIONAL] Leiden community TSV (sample, component, lineage, ...)"),
  make_option("--out-dir", type="character", default=NULL,
    help="[REQUIRED] Output directory"),
  make_option("--metadata", type="character", default=NULL,
    help="[OPTIONAL] Node metadata TSV (sample_id, type, ...) for accurate tumour_count")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (any(sapply(c("edges", "communities", "out-dir"), function(k) is.null(opt[[k]])))) {
  stop("--edges, --communities, and --out-dir are required")
}
dir.create(opt[["out-dir"]], recursive=TRUE, showWarnings=FALSE)

load_communities <- function(path, algorithm) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  dt <- fread(path)
  if (!all(c("sample", "component", "lineage") %in% names(dt))) {
    stop(algorithm, " communities must contain sample, component, and lineage columns")
  }
  dt[, algorithm := algorithm]
  dt
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat("Loading edges:", opt$edges, "\n")
edges <- fread(opt$edges)
if (!all(c("from", "to", "weight") %in% names(edges))) {
  stop("edges must contain from, to, and weight columns")
}
E_und <- edges[, .(weight=max(as.numeric(weight), na.rm=TRUE)),
               by=.(from=pmin(from,to), to=pmax(from,to))]
setorder(E_und, from, to)

cat("Loading Louvain communities:", opt$communities, "\n")
comm_louvain <- load_communities(opt$communities, "Louvain")
comm_leiden <- load_communities(opt[["leiden-communities"]], "Leiden")
comm_all <- rbindlist(Filter(Negate(is.null), list(comm_louvain, comm_leiden)),
                      use.names=TRUE, fill=TRUE)

# Optional metadata for accurate type-based tumour count
node_meta <- NULL
meta_types <- NULL
if (!is.null(opt$metadata) && nzchar(opt$metadata) && file.exists(opt$metadata)) {
  cat("Loading metadata:", opt$metadata, "\n")
  node_meta <- fread(opt$metadata)
  if ("type" %in% names(node_meta) && "sample_id" %in% names(node_meta)) {
    meta_types <- setNames(node_meta$type, node_meta$sample_id)
  }
}

# ---------------------------------------------------------------------------
# Build graph
# ---------------------------------------------------------------------------
vertex_ids <- sort(unique(c(E_und$from, E_und$to, comm_all$sample)))
if (!is.null(node_meta) && "sample_id" %in% names(node_meta)) {
  vertex_ids <- sort(unique(c(vertex_ids, node_meta$sample_id)))
}
g <- graph_from_data_frame(
  E_und,
  directed=FALSE,
  vertices=data.table(name=vertex_ids)
)
E(g)$weight <- E_und$weight

lineage_map <- setNames(comm_louvain$lineage, comm_louvain$sample)
V(g)$lineage <- lineage_map[V(g)$name]
V(g)$lineage[is.na(V(g)$lineage) | V(g)$lineage == ""] <- "UNKNOWN"

# ---------------------------------------------------------------------------
# Compute metrics
# ---------------------------------------------------------------------------
n_nodes <- vcount(g)
n_edges <- ecount(g)

if (!is.null(meta_types)) {
  node_types <- meta_types[V(g)$name]
  tumour_count <- sum(!is.na(node_types) & node_types == "tumour", na.rm=TRUE)
  unknown_type <- sum(is.na(node_types))
  if (unknown_type > 0) {
    warning(unknown_type, " nodes have no type information in metadata; treated as cell lines")
  }
} else {
  tumour_count <- sum(!grepl("^NG-", V(g)$name) &
                        (is.na(V(g)$lineage) | V(g)$lineage != "HEME"))
}

miss_lin <- sum(is.na(V(g)$lineage) | V(g)$lineage == "UNKNOWN")

modularity_rows <- comm_all[, {
  membership_map <- setNames(component, sample)
  membership_vec <- as.integer(membership_map[V(g)$name])
  missing_membership <- sum(is.na(membership_vec))
  keep <- which(!is.na(membership_vec))
  g_sub <- induced_subgraph(g, vids=keep)
  membership_sub <- membership_vec[keep]
  .(
    metric = "weighted_newman_girvan_modularity",
    value = modularity(g_sub, membership_sub, weights=E(g_sub)$weight),
    n_communities = length(unique(membership_sub)),
    n_nodes = vcount(g_sub),
    n_edges = ecount(g_sub),
    missing_membership = missing_membership
  )
}, by=algorithm]

lin_levels <- sort(unique(na.omit(V(g)$lineage[V(g)$lineage != "UNKNOWN"])))
lin_id <- match(V(g)$lineage, lin_levels)
keep_v <- which(!is.na(lin_id))
g_lin <- induced_subgraph(g, vids=keep_v)
lin_id2 <- match(V(g_lin)$lineage, lin_levels)
assort_r <- assortativity_nominal(g_lin, types=lin_id2, directed=FALSE)

# ---------------------------------------------------------------------------
# Print summary
# ---------------------------------------------------------------------------
cat("\n=== VALIDATION SUMMARY ===\n")
cat("Node count:           ", n_nodes, "\n")
cat("Edge count:           ", n_edges, "\n")
for (i in seq_len(nrow(modularity_rows))) {
  cat(modularity_rows$algorithm[i], " communities: ",
      modularity_rows$n_communities[i], "\n", sep="")
  cat(modularity_rows$algorithm[i], " weighted modularity Q: ",
      round(modularity_rows$value[i], 4), "\n", sep="")
}
cat("Assortativity r:      ", round(assort_r, 4), "\n")
cat("Tumour samples:       ", tumour_count, "\n")
cat("Missing lineage:      ", miss_lin, "\n")
cat("==========================\n\n")

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
mod_file <- file.path(opt[["out-dir"]], "validation_modularity.tsv")
fwrite(modularity_rows[order(algorithm)], mod_file, sep="\t")
cat("Saved:", mod_file, "\n")

assort_file <- file.path(opt[["out-dir"]], "validation_assortativity.tsv")
fwrite(data.table(
  metric = "unweighted_nominal_lineage_assortativity",
  value = assort_r,
  n_nodes_used = length(keep_v),
  n_edges_used = ecount(g_lin),
  edge_weights_used = FALSE,
  tumour_count = tumour_count,
  missing_lineage = miss_lin
), assort_file, sep="\t")
cat("Saved:", assort_file, "\n")

cat("[OK] Validation complete\n")
