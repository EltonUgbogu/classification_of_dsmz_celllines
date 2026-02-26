#!/usr/bin/env Rscript

# =============================================================================
# subcluster_community.R
# =============================================================================
# Sub-cluster a specific community to resolve mixed lineages
#
# INPUT:
#   --edges      pan_cancer_graph_edges.tsv
#   --communities pan_cancer_communities.tsv (with component column)
#   --meta       pan_cancer_cellline_metadata.tsv
#   --target-comm integer (which community to subcluster)
#   --out        output TSV for sub-communities
#   --method     louvain|leiden (default: louvain)
#   --resolution numeric (for leiden, default: 1.5)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(optparse)
})

opt_list <- list(
  make_option("--edges", type="character", default=NULL,
              help="Path to pan_cancer_graph_edges.tsv"),
  make_option("--communities", type="character", default=NULL,
              help="Path to pan_cancer_communities.tsv"),
  make_option("--meta", type="character", default=NULL,
              help="Path to metadata TSV with sample_id + lineage"),
  make_option("--target-comm", type="integer", default=NULL,
              help="Community ID to subcluster"),
  make_option("--out", type="character", default="subcommunities.tsv",
              help="Output TSV for sub-communities"),
  make_option("--method", type="character", default="louvain",
              help="Method: louvain|leiden (default: louvain)"),
  make_option("--resolution", type="numeric", default=1.5,
              help="Resolution parameter for Leiden (default: 1.5)"),
  make_option("--seed", type="integer", default=1,
              help="Random seed (default: 1)")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$edges) || !file.exists(opt$edges)) {
  stop("--edges is required and must exist")
}
if (is.null(opt$communities) || !file.exists(opt$communities)) {
  stop("--communities is required and must exist")
}
if (is.null(opt$`target-comm`)) {
  stop("--target-comm is required (which community to subcluster)")
}

set.seed(opt$seed)

cat("========================================\n")
cat("SUB-COMMUNITY ANALYSIS\n")
cat("========================================\n\n")

# Load data
cat("[1] Loading data...\n")
E <- fread(opt$edges)
comms <- fread(opt$communities)

# Handle flexible column names
if ("component" %in% names(comms)) {
  comm_col <- "component"
} else if ("community" %in% names(comms)) {
  comm_col <- "community"
} else {
  stop("communities file must have 'component' or 'community' column")
}

if ("sample" %in% names(comms)) {
  sample_col <- "sample"
} else if ("sample_id" %in% names(comms)) {
  sample_col <- "sample_id"
} else {
  stop("communities file must have 'sample' or 'sample_id' column")
}

# Build graph
E[, a := pmin(from, to)]
E[, b := pmax(from, to)]
E_und <- E[, .(weight = max(weight, na.rm=TRUE)), by = .(a,b)]
setnames(E_und, c("a","b"), c("from","to"))

g <- graph_from_data_frame(E_und, directed = FALSE)
E(g)$weight <- E_und$weight

# Get target community nodes
target_nodes <- comms[get(comm_col) == opt$`target-comm`, get(sample_col)]
cat("  Target community: ", opt$`target-comm`, "\n", sep="")
cat("  Nodes in community: ", length(target_nodes), "\n", sep="")

# Extract subgraph
g_sub <- induced_subgraph(g, vids = target_nodes)
cat("  Subgraph nodes: ", vcount(g_sub), "\n", sep="")
cat("  Subgraph edges: ", ecount(g_sub), "\n\n", sep="")

# Load metadata if provided
if (!is.null(opt$meta) && file.exists(opt$meta)) {
  meta <- fread(opt$meta)
  meta_map <- setNames(meta$lineage, meta$sample_id)
  V(g_sub)$lineage <- meta_map[V(g_sub)$name]
  V(g_sub)$lineage[is.na(V(g_sub)$lineage)] <- "UNKNOWN"
  cat("[2] Metadata loaded. Lineage distribution:\n")
  print(table(V(g_sub)$lineage))
  cat("\n")
} else {
  cat("[2] No metadata provided...\n\n")
}

# Sub-cluster
cat("[3] Computing sub-communities...\n")
if (opt$method == "louvain") {
  cl_sub <- cluster_louvain(g_sub, weights = E(g_sub)$weight)
  cat("  Sub-communities: ", length(unique(membership(cl_sub))), "\n", sep="")
  cat("  Modularity: ", round(modularity(cl_sub), 4), "\n\n", sep="")
} else if (opt$method == "leiden") {
  if (!"cluster_leiden" %in% ls("package:igraph")) {
    stop("cluster_leiden() not available. Use --method louvain or upgrade igraph.")
  }
  cl_sub <- cluster_leiden(g_sub, weights = E(g_sub)$weight, 
                           resolution_parameter = opt$resolution)
  cat("  Resolution parameter: ", opt$resolution, "\n", sep="")
  cat("  Sub-communities: ", length(unique(membership(cl_sub))), "\n", sep="")
  cat("  Modularity: ", round(modularity(cl_sub), 4), "\n\n", sep="")
} else {
  stop("--method must be louvain or leiden")
}

# Create output
subcomms <- data.table(
  sample = V(g_sub)$name,
  parent_community = opt$`target-comm`,
  sub_community = as.integer(membership(cl_sub))
)
subcomms[, sub_community_size := as.integer(table(sub_community)[as.character(sub_community)])]

if (!is.null(opt$meta) && file.exists(opt$meta)) {
  subcomms[, lineage := V(g_sub)$lineage]
  
  # Lineage enrichment per sub-community
  cat("[4] Sub-community lineage enrichment:\n")
  enrich <- subcomms[lineage != "UNKNOWN", .N, by=.(sub_community, lineage)]
  enrich_tot <- subcomms[lineage != "UNKNOWN", .N, by=sub_community]
  setnames(enrich_tot, "N", "n_total")
  enrich <- merge(enrich, enrich_tot, by="sub_community", all.x=TRUE)
  enrich[, fraction := N / n_total]
  enrich <- enrich[order(sub_community, -fraction)]
  print(enrich)
  cat("\n")
}

# Save
cat("[5] Writing sub-communities...\n")
fwrite(subcomms[order(sub_community, sample)], opt$out, sep="\t")
cat("  Saved: ", opt$out, "\n\n", sep="")

cat("========================================\n")
cat("DONE\n")
cat("========================================\n")
