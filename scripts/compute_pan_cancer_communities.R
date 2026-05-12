#!/usr/bin/env Rscript

# =============================================================================
# compute_pan_cancer_communities.R
# =============================================================================
#
# DESCRIPTION
# -----------
# This script performs community detection on a pan-cancer similarity graph
# to identify groups of samples (cell lines or tumours) that share similar
# transcriptomic profiles. Community detection partitions the graph into
# densely connected subgroups, revealing biologically meaningful clusters
# that often correspond to cancer types, molecular subtypes, or shared
# transcriptional programmes.
#
#
# BIOLOGICAL CONTEXT: COMMUNITY STRUCTURE IN CANCER TRANSCRIPTOMICS
# -----------------------------------------------------------------
# Gene expression profiles from cancer samples can be represented as a
# similarity network where nodes represent samples and edges represent
# transcriptomic similarity (typically Spearman correlation). This network
# representation enables the application of graph-theoretic methods to
# discover structure in the data.
#
# Community detection algorithms identify groups of nodes that are more
# densely connected to each other than to the rest of the network. In the
# context of cancer transcriptomics, these communities typically represent:
#
#   Cancer Lineages:
#   Samples from the same cancer type tend to cluster together due to
#   tissue-of-origin gene expression signatures. Breast cancers share
#   epithelial programmes distinct from neuroblastomas or haematologic
#   malignancies.
#
#   Molecular Subtypes:
#   Within a cancer type, community structure may resolve known molecular
#   subtypes. For example, breast cancer communities might separate into
#   luminal, HER2-enriched, and basal-like groups reflecting distinct
#   transcriptional programmes and clinical behaviour.
#
#   Shared Biological Programmes:
#   Some communities may span multiple cancer types, grouping samples
#   that share transcriptional features such as epithelial-mesenchymal
#   transition (EMT), proliferation signatures, or immune infiltration
#   patterns.
#
#
# COMMUNITY DETECTION ALGORITHMS
# ------------------------------
# The script implements two widely-used community detection methods:
#
# Louvain Algorithm:
# A greedy modularity optimisation algorithm that iteratively merges nodes
# into communities to maximise modularity (the fraction of edges within
# communities compared to a random network). Louvain is computationally
# efficient and scalable, making it suitable for large graphs. However,
# it may produce slightly different results across runs due to the order
# of node processing, and can suffer from the resolution limit problem
# where small communities may be inappropriately merged.
#
# Leiden Algorithm:
# An improved version of Louvain that addresses several limitations of
# its predecessor. Leiden guarantees that all communities are well-connected
# (no disconnected subcommunities) and provides more stable partitions.
# The algorithm includes a refinement phase that can split poorly-connected
# communities, often producing more biologically interpretable results.
#
# Both algorithms utilise edge weights to inform community structure, giving
# stronger influence to highly correlated sample pairs. The resolution
# parameter (default: 1.0) controls the granularity of detected communities,
# with higher values producing more numerous, smaller communities.
#
#
# VALIDATION METRICS
# ------------------
# The script computes several metrics to assess community quality:
#
# Modularity:
# Measures the strength of community structure by comparing the density
# of within-community edges to the expected density under a null model
# (configuration model preserving degree distribution). Values range from
# -0.5 to 1.0, with higher values indicating stronger community structure.
# Typical biological networks exhibit modularity in the range of 0.3 to 0.7.
# Values above 0.3 are generally considered indicative of significant
# community structure.
#
# Assortativity:
# Measures the tendency of nodes to connect preferentially with others
# of the same type (in this case, cancer lineage). Assortativity ranges
# from -1 to 1, where positive values indicate assortative mixing (like
# connects to like) and negative values indicate disassortative mixing.
# High positive assortativity validates that the similarity metric captures
# biologically meaningful relationships where samples of the same cancer
# type are more similar to each other than to samples of different types.
#
# Community Purity:
# For each community, purity is defined as the fraction of samples belonging
# to the most common (dominant) lineage within that community. High purity
# (approaching 1.0) indicates that communities correspond to single cancer
# types, validating lineage-specific clustering. Lower purity values suggest
# cross-lineage groupings that may reflect shared biological programmes,
# technical batch effects, or suboptimal community detection parameters.
#
#
# INPUT REQUIREMENTS
# ------------------
# Required:
#   --edges: Path to edge list TSV containing the similarity graph with
#            columns:
#     - from:   Source node identifier (sample ID)
#     - to:     Target node identifier (sample ID)
#     - weight: Edge weight representing similarity (e.g., Spearman rho)
#
# Optional:
#   --meta:   Path to metadata TSV with sample annotations containing:
#     - sample_id: Sample identifier matching graph node names
#     - lineage:   Cancer type annotation for validation analyses
#
#   --out:    Output path for community assignments TSV
#             (default: pan_cancer_communities.tsv)
#
#   --outdir: Output directory for validation tables
#             (default: same directory as --out)
#
#   --method: Community detection algorithm to use
#             Options: louvain, leiden, or both (default: both)
#
#   --seed:   Random seed for reproducibility (default: 1)
#
#
# OUTPUT FILES
# ------------
# Primary Output:
#   pan_cancer_communities.tsv
#       Per-sample community assignments with columns:
#       - sample:    Sample identifier
#       - component: Community assignment (integer)
#       - comp_size: Number of samples in the assigned community
#       - lineage:   Cancer type annotation (if metadata provided)
#
# Validation Tables:
#   validation_community_summary.tsv
#       Per-community statistics including size, internal edge counts,
#       mean within-community edge weight, and degree distribution metrics.
#
#   validation_modularity.tsv
#       Modularity scores and community counts for each algorithm run.
#
#   validation_lineage_enrichment.tsv (if metadata provided)
#       Community-by-lineage contingency table with purity scores,
#       enabling assessment of lineage coherence within communities.
#
#   validation_assortativity.tsv (if metadata provided)
#       Graph-level assortativity coefficient measuring the tendency
#       of samples to connect with others of the same lineage.
#
# Graph Files:
#   pan_cancer_graph.rds
#       R-native serialised igraph object for downstream analysis.
#
#   pan_cancer_graph.graphml
#       GraphML format for interoperability with other tools (e.g., Cytoscape).
#
#   pan_cancer_graph.gml
#       GML format for compatibility with legacy graph analysis software.
#
#
# USAGE EXAMPLE
# -------------
#   Rscript compute_pan_cancer_communities.R \
#     --edges pan_cancer_graph_edges.tsv \
#     --meta sample_metadata.tsv \
#     --out pan_cancer_communities.tsv \
#     --outdir validation \
#     --method both \
#     --seed 42
#
#
# DEPENDENCIES
# ------------
# R packages: data.table, igraph, optparse
#
# Note: The Leiden algorithm requires igraph version 1.2.5 or later.
# If cluster_leiden() is unavailable, the script will fall back to
# Louvain or terminate with an informative error message.
#
# =============================================================================


# =============================================================================
# PACKAGE LOADING
# =============================================================================
# The script requires three packages for its core functionality:
#
#   data.table: Provides efficient data manipulation with memory-optimised
#               storage and fast aggregation operations. Essential for
#               handling large edge lists and computing community statistics.
#
#   igraph:     Comprehensive graph analysis library providing network
#               construction, community detection algorithms (Louvain, Leiden),
#               and graph-theoretic metrics (modularity, assortativity, degree).
#               The igraph package is the de facto standard for network
#               analysis in R.
#
#   optparse:   Enables GNU-style command-line argument parsing with support
#               for default values, type checking, and automatic help generation.
#
# The suppressPackageStartupMessages() wrapper prevents verbose loading
# messages that would clutter pipeline output logs.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({

  library(data.table)
  library(igraph)
  library(optparse)
})


# =============================================================================
# COMMAND-LINE ARGUMENT DEFINITIONS
# =============================================================================
# Command-line arguments enable flexible script execution with configurable
# algorithm selection, output paths, and reproducibility settings.
#
# The script supports running both Louvain and Leiden algorithms for
# comparison, with Leiden preferred as the primary output when available
# due to its improved theoretical guarantees.
# -----------------------------------------------------------------------------

opt_list <- list(
  make_option("--edges", type = "character", default = NULL,
              help = "Path to pan_cancer_graph_edges.tsv"),
  make_option("--meta", type = "character", default = NULL,
              help = "Optional metadata TSV with sample_id + lineage"),
  make_option("--meta-fallback", type = "character", default = NULL,
              help = "Optional fallback metadata CSV/TSV used to fill blank lineage labels"),
  make_option("--out", type = "character", default = "pan_cancer_communities.tsv",
              help = "Output TSV for community assignments"),
  make_option("--outdir", type = "character", default = NULL,
              help = "Output directory for validation tables (default: dirname(out))"),
  make_option("--method", type = "character", default = "both",
              help = "Community method: louvain|leiden|both (default: both)"),
  make_option("--seed", type = "integer", default = 1,
              help = "Random seed (default: 1)")
)

opt <- parse_args(OptionParser(option_list = opt_list))


# =============================================================================
# INPUT VALIDATION AND SETUP
# =============================================================================
# The script validates that required inputs exist before proceeding.
# Early validation prevents confusing downstream errors and provides
# clear feedback about missing dependencies.
#
# The random seed is set early to ensure reproducibility of community
# detection results, as both Louvain and Leiden involve stochastic
# elements in their optimisation procedures.
# -----------------------------------------------------------------------------

if (is.null(opt$edges) || !file.exists(opt$edges)) {
  stop("--edges is required and must exist")
}

# Default output directory to the same location as the main output file
if (is.null(opt$outdir)) opt$outdir <- dirname(opt$out)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Set random seed for reproducibility of stochastic community detection
set.seed(opt$seed)

lineage_from_values <- function(values) {
  out <- rep(NA_character_, length(values))
  if (is.null(values)) {
    return(out)
  }
  out[grepl("BRCA|Breast", values, ignore.case = TRUE)] <- "BRCA"
  out[grepl("NBL|Neuroblastoma", values, ignore.case = TRUE)] <- "NBL"
  out[grepl("RBL|Retinoblastoma", values, ignore.case = TRUE)] <- "RBL"
  out[grepl("HEME|Hema|LL-100|Leukemia|Lymphoma", values, ignore.case = TRUE)] <- "HEME"
  out
}

derive_lineage <- function(dt, default_lineage = "UNKNOWN") {
  lineage <- if ("lineage" %in% names(dt)) as.character(dt$lineage) else rep(NA_character_, nrow(dt))
  lineage[is.na(lineage) | lineage == ""] <- NA_character_
  if ("lineage_fbk" %in% names(dt)) {
    idx <- is.na(lineage) | lineage == ""
    lineage[idx] <- as.character(dt$lineage_fbk[idx])
  }
  for (col in c("cancer_type", "cancer_type_fbk", "Disease", "Disease_fbk",
                "group2", "group2_fbk", "organ", "organ_fbk")) {
    if (col %in% names(dt)) {
      guess <- lineage_from_values(dt[[col]])
      idx <- is.na(lineage) | lineage == ""
      lineage[idx] <- guess[idx]
    }
  }
  lineage[is.na(lineage) | lineage == ""] <- default_lineage
  toupper(lineage)
}

cat("========================================\n")
cat("PAN-CANCER COMMUNITY DETECTION\n")
cat("========================================\n\n")


# =============================================================================
# STEP 1: LOAD EDGES AND BUILD GRAPH
# =============================================================================
# This step loads the edge list and constructs an undirected, weighted graph
# suitable for community detection. The graph representation enables
# application of network analysis algorithms to identify sample clusters.
#
# Edge List Processing:
# The input edge list may contain directed edges (A->B and B->A recorded
# separately) or duplicate edges from different analysis steps. The script
# collapses these to undirected edges by:
#   1. Canonicalising edge orientation (smaller ID first)
#   2. Taking the maximum weight for duplicate edges
#
# This approach preserves the strongest observed similarity between any
# pair of samples while ensuring the graph is symmetric.
#
# Graph Construction:
# The igraph library constructs a graph object from the edge data frame,
# with edge weights stored as an edge attribute. Weighted community
# detection algorithms use these weights to preferentially group samples
# with stronger transcriptomic similarity.
# -----------------------------------------------------------------------------

cat("[1] Loading edges...\n")

# Load edge list using fast file reading
E <- fread(opt$edges)

# Validate required columns are present
req <- c("from", "to", "weight")
miss <- setdiff(req, names(E))
if (length(miss) > 0) {
  stop("Edge file missing required columns: ", paste(miss, collapse = ", "))
}

# -----------------------------------------------------------------------------
# Convert to undirected edges by canonicalising node order
# For each edge, place the lexicographically smaller node ID first.
# This ensures that edges A-B and B-A are treated as the same edge.
# -----------------------------------------------------------------------------

E[, a := pmin(from, to)]
E[, b := pmax(from, to)]

# Collapse duplicate edges by taking the maximum weight
# This preserves the strongest observed similarity between sample pairs
E_und <- E[, .(weight = max(weight, na.rm = TRUE)), by = .(a, b)]
setnames(E_und, c("a", "b"), c("from", "to"))

cat("  Edges (raw):       ", nrow(E), "\n", sep = "")
cat("  Edges (undirected):", nrow(E_und), "\n", sep = "")

# -----------------------------------------------------------------------------
# Construct igraph object from the edge data frame
# The directed = FALSE argument creates an undirected graph where edges
# can be traversed in either direction.
# -----------------------------------------------------------------------------

g <- graph_from_data_frame(E_und, directed = FALSE)

# Assign edge weights as a graph attribute
# Community detection algorithms will use these weights to inform clustering
E(g)$weight <- E_und$weight

cat("  Nodes: ", vcount(g), "\n", sep = "")
cat("  Edges: ", ecount(g), "\n\n", sep = "")

# -----------------------------------------------------------------------------
# Save Graph in Multiple Formats
# The graph is saved in three formats to support different downstream uses:
#
#   RDS: Native R serialisation preserving all igraph object attributes.
#        Ideal for continued analysis in R pipelines.
#
#   GraphML: XML-based format widely supported by graph visualisation tools
#            including Cytoscape, Gephi, and yEd. Preserves node and edge
#            attributes with explicit typing.
#
#   GML: Graph Modelling Language format for compatibility with legacy
#        software and certain network analysis packages.
# -----------------------------------------------------------------------------

graph_rds <- file.path(opt$outdir, "pan_cancer_graph.rds")
graph_graphml <- file.path(opt$outdir, "pan_cancer_graph.graphml")
graph_gml <- file.path(opt$outdir, "pan_cancer_graph.gml")

saveRDS(g, graph_rds)
write_graph(g, graph_graphml, format = "graphml")
write_graph(g, graph_gml, format = "gml")

cat("  Saved graph object: ", graph_rds, "\n", sep = "")
cat("  Saved GraphML:      ", graph_graphml, "\n", sep = "")
cat("  Saved GML:          ", graph_gml, "\n\n", sep = "")


# =============================================================================
# STEP 2: LOAD METADATA AND ANNOTATE VERTICES
# =============================================================================
# When metadata is provided, this step annotates graph vertices with cancer
# lineage information. This annotation enables downstream validation analyses
# that assess whether detected communities correspond to known cancer types.
#
# Metadata Alignment:
# The script creates a mapping from sample identifiers to lineage annotations,
# then applies this mapping to graph vertices. Samples without metadata
# annotations receive the "UNKNOWN" label, allowing the analysis to proceed
# while flagging incomplete annotation coverage.
#
# Biological Rationale:
# Lineage annotations provide ground truth for evaluating community detection
# quality. If communities correspond to cancer types, the analysis validates
# that transcriptomic similarity captures biologically meaningful relationships.
# Communities spanning multiple lineages may indicate shared biology (e.g.,
# EMT programmes) or potential data quality issues requiring investigation.
# -----------------------------------------------------------------------------

meta <- NULL
if (!is.null(opt$meta) && file.exists(opt$meta)) {
  cat("[2] Loading metadata...\n")
  meta <- fread(opt$meta)
  fallback_meta <- NULL
  if (!is.null(opt$`meta-fallback`) && file.exists(opt$`meta-fallback`)) {
    fallback_meta <- fread(opt$`meta-fallback`)
    if (!("sample_id" %in% names(fallback_meta))) {
      stop("Fallback metadata must contain sample_id")
    }
    keep_cols <- intersect(c("sample_id", "lineage", "cancer_type", "Disease", "group2", "organ"), names(fallback_meta))
    fallback_meta <- unique(fallback_meta[, ..keep_cols])
    fallback_cols <- setdiff(names(fallback_meta), "sample_id")
    setnames(fallback_meta, fallback_cols, paste0(fallback_cols, "_fbk"))
  }
  
  # Normalize lineage column if provided as cancer_type
  if (!("lineage" %in% names(meta)) && ("cancer_type" %in% names(meta))) {
    meta[, lineage := fcase(
      cancer_type %like% "BRCA|Breast", "BRCA",
      cancer_type %like% "NBL|Neuroblastoma", "NBL",
      cancer_type %like% "RBL|Retinoblastoma", "RBL",
      cancer_type %like% "HEME|Hematologic", "HEME",
      default = NA_character_
    )]
  }

  # Validate required metadata columns
  if (!all(c("sample_id", "lineage") %in% names(meta))) {
    stop("meta must contain columns: sample_id, lineage")
  }
  
  meta <- unique(meta, by = "sample_id")
  if (!is.null(fallback_meta)) {
    meta <- merge(meta, fallback_meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  }
  meta[, lineage := derive_lineage(meta)]
  meta[lineage == "", lineage := "UNKNOWN"]
  
  # Create unique sample-to-lineage mapping
  # Deduplication handles cases where samples appear multiple times in metadata
  meta <- unique(meta[, .(sample_id, lineage)])
  meta_map <- setNames(meta$lineage, meta$sample_id)

  # Apply lineage annotations to graph vertices
  V(g)$lineage <- meta_map[V(g)$name]
  
  # Assign "UNKNOWN" to samples lacking metadata
  # This explicit labelling enables downstream filtering and prevents
  # NA-related issues in validation analyses
  V(g)$lineage[is.na(V(g)$lineage)] <- "UNKNOWN"
  
  # Report annotation coverage
  na_lin <- sum(V(g)$lineage == "UNKNOWN")
  cat("  Lineage annotated for: ", vcount(g) - na_lin, " / ", vcount(g), "\n", sep = "")
  cat("  UNKNOWN lineage nodes: ", na_lin, "\n\n", sep = "")
} else {
  cat("[2] No metadata provided (skipping lineage validation)...\n\n")
}


# =============================================================================
# STEP 3: COMMUNITY DETECTION
# =============================================================================
# This step applies community detection algorithms to partition the graph
# into densely connected subgroups. The script supports two algorithms
# with different characteristics and guarantees.
#
# Algorithm Selection Rationale:
# Running both algorithms enables comparison of results and selection of
# the most appropriate partition for downstream analysis. When both are
# computed, Leiden is preferred as the primary output due to its stronger
# theoretical guarantees (connected communities, no resolution limit).
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Louvain Algorithm Implementation
# The Louvain algorithm performs greedy modularity optimisation through
# two iterative phases:
#   1. Local moving: Each node is moved to the community that maximises
#      modularity gain
#   2. Aggregation: Communities are collapsed into super-nodes and the
#      process repeats
#
# The algorithm terminates when no further modularity improvement is possible.
# Edge weights influence the modularity calculation, giving stronger weight
# to edges between highly similar samples.
# -----------------------------------------------------------------------------

do_louvain <- function(g) {
  # cluster_louvain uses the weight edge attribute automatically when present
  cluster_louvain(g, weights = E(g)$weight)
}

# -----------------------------------------------------------------------------
# Leiden Algorithm Implementation
# The Leiden algorithm improves upon Louvain with three key enhancements:
#   1. Guaranteed connectivity: All detected communities are connected subgraphs
#   2. Refinement phase: Communities can be split if poorly connected internally
#   3. Faster convergence: More efficient local moving heuristics
#
# The resolution parameter controls community granularity:
#   - Lower resolution (< 1.0): Fewer, larger communities
#   - Higher resolution (> 1.0): More, smaller communities
#   - Default (1.0): Standard modularity optimisation
#
# The function includes version compatibility handling for the resolution
# parameter, which was renamed between igraph versions.
# -----------------------------------------------------------------------------

do_leiden <- function(g, resolution = 1.0) {
  # Check if Leiden is available in the installed igraph version
  if ("cluster_leiden" %in% ls("package:igraph")) {
    # Attempt with current parameter name, fall back to deprecated name
    tryCatch({
      cluster_leiden(g, weights = E(g)$weight, resolution = resolution)
    }, error = function(e) {
      # Older igraph versions used 'resolution_parameter' instead of 'resolution'
      cluster_leiden(g, weights = E(g)$weight, resolution_parameter = resolution)
    })
  } else {
    stop("cluster_leiden() not available in your igraph build. ",
         "Upgrade igraph or use --method louvain.")
  }
}

# Validate method argument
methods <- tolower(opt$method)
if (!methods %in% c("louvain", "leiden", "both")) {
  stop("--method must be louvain|leiden|both")
}

cat("[3] Computing communities...\n")

# Initialise result objects
cl_louv <- NULL
cl_leid <- NULL

# -----------------------------------------------------------------------------
# Execute Louvain Algorithm
# Report the number of detected communities and the achieved modularity score.
# Higher modularity indicates stronger community structure.
# -----------------------------------------------------------------------------

if (methods %in% c("louvain", "both")) {
  cat("  - Louvain...\n")
  cl_louv <- do_louvain(g)
  cat("    Louvain communities: ", length(unique(membership(cl_louv))), "\n", sep = "")
  cat("    Louvain modularity:  ", round(modularity(cl_louv), 4), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# Execute Leiden Algorithm
# Includes error handling for igraph version compatibility and sanity checking
# for degenerate results (one community per node indicates algorithm failure).
# -----------------------------------------------------------------------------

if (methods %in% c("leiden", "both")) {
  cat("  - Leiden...\n")
  cl_leid <- tryCatch(do_leiden(g, resolution = 1.0), error = function(e) {
    message("    Leiden failed: ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(cl_leid)) {
    n_comms <- length(unique(membership(cl_leid)))
    cat("    Leiden communities: ", n_comms, "\n", sep = "")
    
    # Modularity calculation may fail in some igraph versions
    mod_leid <- tryCatch(modularity(cl_leid), error = function(e) NA_real_)
    if (!is.na(mod_leid)) {
      cat("    Leiden modularity:  ", round(mod_leid, 4), "\n", sep = "")
    } else {
      cat("    Leiden modularity:  (not available)\n", sep = "")
    }
    
    # Sanity check: One community per node indicates algorithm failure
    # This can occur with excessively high resolution or numerical issues
    if (n_comms == vcount(g)) {
      warning("    WARNING: Leiden produced one community per node. ",
              "This may indicate resolution is too high or algorithm issue.")
    }
  }
}

# -----------------------------------------------------------------------------
# Select Primary Community Assignment
# When both algorithms are run, Leiden is preferred due to its theoretical
# advantages. The primary assignment is used for the main output file and
# downstream validation analyses.
# -----------------------------------------------------------------------------

primary <- "louvain"
if (!is.null(cl_leid)) {
  n_comms_leid <- length(unique(membership(cl_leid)))
  if (n_comms_leid < vcount(g)) {
    primary <- "leiden"
  }
}
cat("\n  Primary assignment: ", primary, "\n\n", sep = "")

# Extract membership vector from the selected algorithm
mem_primary <- if (primary == "leiden") membership(cl_leid) else membership(cl_louv)

# -----------------------------------------------------------------------------
# FAIL-FAST SANITY CHECKS (prevents silent regressions)
# -----------------------------------------------------------------------------
stopifnot(!is.null(mem_primary))
stopifnot(length(mem_primary) > 1)

comm_sizes <- sort(table(mem_primary), decreasing = TRUE)

# 1) Degenerate check: all singleton communities => invalid
if (max(comm_sizes) <= 1) {
  stop(sprintf(
    "Degenerate community assignment (%s): all communities are singletons (max size = %d).",
    primary, max(comm_sizes)
  ))
}

# 2) Require at least one meaningful community
if (max(comm_sizes) < 5) {
  stop(sprintf(
    "Suspicious community assignment (%s): max community size = %d (< 5). Refusing to write outputs.",
    primary, max(comm_sizes)
  ))
}

# 3) Basic graph consistency: graph must have edges
if (ecount(g) == 0) {
  stop("Graph has 0 edges. Community detection is undefined; refusing to proceed.")
}

# 4) Membership must cover exactly the graph vertices
vnames <- V(g)$name
if (!setequal(names(mem_primary), vnames)) {
  missing <- setdiff(vnames, names(mem_primary))
  extra <- setdiff(names(mem_primary), vnames)
  stop(sprintf(
    "Membership/graph mismatch: missing=%d extra=%d (check joins / naming).",
    length(missing), length(extra)
  ))
}

# 5) Within-community edges must be non-zero
ends_mat <- ends(g, E(g), names = TRUE)
same_comm <- mem_primary[ends_mat[, 1]] == mem_primary[ends_mat[, 2]]
within_edges_n <- sum(same_comm, na.rm = TRUE)
if (within_edges_n == 0) {
  stop(sprintf(
    "Invalid community assignment (%s): 0 within-community edges. Refusing to write outputs.",
    primary
  ))
}

message(sprintf(
  "Fail-fast checks passed: primary=%s; n_communities=%d; max_size=%d; within_edges=%d",
  primary, length(comm_sizes), max(comm_sizes), within_edges_n
))


# =============================================================================
# STEP 4: CREATE AND SAVE COMMUNITY ASSIGNMENTS
# =============================================================================
# This step constructs the primary output table containing per-sample
# community assignments. The output format uses "community" terminology
# to reflect true Louvain/Leiden memberships.
#
# Output Columns:
#   sample:    Sample identifier (matches graph vertex names)
#   community: Integer community assignment from the primary algorithm
#   community_size: Number of samples in the assigned community
#   lineage:   Cancer type annotation (if metadata was provided)
# -----------------------------------------------------------------------------

# Create community assignment table
communities <- data.table(
  sample = V(g)$name,
  community = as.integer(mem_primary)
)

# Calculate community sizes and attach to each sample
# This enables downstream filtering by community size without additional joins
communities[, community_size := as.integer(table(community)[as.character(community)])]

# Attach lineage annotations if metadata was provided
if (!is.null(meta)) {
  communities[, lineage := V(g)$lineage]
}

# Write primary output file sorted by community and sample for readability
cat("[4] Writing pan_cancer_communities.tsv...\n")
fwrite(communities[order(community, sample)], opt$out, sep = "\t")
cat("  Saved: ", opt$out, "\n\n", sep = "")


# =============================================================================
# STEP 5: VALIDATION TABLES
# =============================================================================
# This section generates diagnostic tables for assessing community detection
# quality and biological coherence. These validation outputs enable quality
# control and inform interpretation of downstream analyses.
# -----------------------------------------------------------------------------

cat("[5] Writing validation tables...\n")

# Store community assignments as vertex attribute for edge analysis
V(g)$community <- mem_primary


# =============================================================================
# VALIDATION 5.1: COMMUNITY SUMMARY STATISTICS
# =============================================================================
# This table provides per-community structural statistics for assessing
# the quality and characteristics of detected communities.
#
# Metrics Computed:
#   n_nodes:            Number of samples in the community
#   n_within_edges:     Number of edges connecting samples within the community
#   mean_within_weight: Average edge weight for within-community edges
#   mean_degree:        Average node degree within the community
#   median_degree:      Median node degree (robust to outliers)
#   max_degree:         Maximum node degree (identifies hub samples)
#
# Interpretation:
# Communities with high internal edge counts and weights relative to their
# size indicate tightly-connected clusters with strong transcriptomic
# similarity. Degree statistics reveal whether communities contain hub
# samples (high max_degree) or have uniform connectivity (similar mean
# and median degree).
# -----------------------------------------------------------------------------

# Convert graph edges to data.table for efficient manipulation
ed <- as.data.table(as_data_frame(g, what = "edges"))

# Assign community membership to each edge endpoint
ed[, c_from := V(g)$community[match(from, V(g)$name)]]
ed[, c_to := V(g)$community[match(to, V(g)$name)]]

# Flag within-community edges (both endpoints in same community)
ed[, is_within := (c_from == c_to)]

# Calculate community sizes
comm_sizes <- communities[, .(n_nodes = .N), by = community]

# Calculate within-community edge statistics
# These metrics indicate how densely connected each community is internally
comm_within <- ed[is_within == TRUE, .(
  n_within_edges = .N,
  mean_within_weight = mean(weight)
), by = c_from]
setnames(comm_within, "c_from", "community")

# Merge size and edge statistics
comm_summary <- merge(comm_sizes, comm_within, by = "community", all.x = TRUE)

# Handle communities with no internal edges (singletons or disconnected)
comm_summary[is.na(n_within_edges), `:=`(n_within_edges = 0L, mean_within_weight = NA_real_)]

# -----------------------------------------------------------------------------
# Calculate Degree Distribution Statistics per Community
# Node degree (number of connections) characterises sample connectivity.
# Per-community degree statistics reveal structural heterogeneity.
# -----------------------------------------------------------------------------

deg <- degree(g, mode = "all")
deg_dt <- data.table(
  sample = names(deg),
  degree = as.integer(deg),
  community = as.integer(mem_primary[names(deg)])
)

deg_comm <- deg_dt[, .(
  mean_degree = mean(degree),
  median_degree = median(degree),
  max_degree = max(degree)
), by = community]

# Merge degree statistics into community summary
comm_summary <- merge(comm_summary, deg_comm, by = "community", all.x = TRUE)

# Write community summary sorted by size (largest first)
comm_summary_file <- file.path(opt$outdir, "validation_community_summary.tsv")
fwrite(comm_summary[order(-n_nodes)], comm_summary_file, sep = "\t")
cat("  Saved: ", comm_summary_file, "\n", sep = "")


# =============================================================================
# VALIDATION 5.2: MODULARITY COMPARISON
# =============================================================================
# This table records modularity scores for each algorithm run, enabling
# comparison of community detection quality across methods.
#
# Modularity Interpretation:
# Modularity measures the fraction of edges falling within communities
# minus the expected fraction if edges were distributed randomly (while
# preserving node degrees). The value ranges from -0.5 to 1.0:
#
#   Q > 0.3:  Significant community structure
#   Q > 0.5:  Strong community structure
#   Q > 0.7:  Very strong community structure (rare in biological networks)
#
# Comparing modularity between Louvain and Leiden indicates whether the
# algorithms find similar quality partitions. Large differences may
# suggest the community structure is sensitive to algorithmic choices.
# -----------------------------------------------------------------------------

mods <- data.table(method = character(), modularity = numeric(), n_communities = integer())

if (!is.null(cl_louv)) {
  mods <- rbind(mods, data.table(
    method = "louvain",
    modularity = modularity(cl_louv),
    n_communities = length(unique(membership(cl_louv)))
  ))
}

if (!is.null(cl_leid)) {
  mod_leid <- tryCatch(modularity(cl_leid), error = function(e) NA_real_)
  mods <- rbind(mods, data.table(
    method = "leiden",
    modularity = mod_leid,
    n_communities = length(unique(membership(cl_leid)))
  ))
}

mods_file <- file.path(opt$outdir, "validation_modularity.tsv")
fwrite(mods, mods_file, sep = "\t")
cat("  Saved: ", mods_file, "\n", sep = "")


# =============================================================================
# VALIDATION 5.3: LINEAGE COHERENCE ANALYSIS
# =============================================================================
# When metadata is available, this section assesses whether detected
# communities correspond to known cancer lineages. Two complementary
# metrics are computed: graph-level assortativity and per-community purity.
#
# These validation analyses address the fundamental question: does the
# transcriptomic similarity network capture biologically meaningful
# relationships where samples of the same cancer type cluster together?
# -----------------------------------------------------------------------------

if (!is.null(meta)) {
  
  # ---------------------------------------------------------------------------
  # Assortativity Analysis
  # Assortativity measures the tendency of nodes to connect preferentially
  # with others of the same type. For nominal (categorical) attributes like
  # cancer lineage, assortativity ranges from -1 to 1:
  #
  #   r > 0:  Assortative mixing (like connects to like)
  #   r = 0:  Random mixing (no preference)
  #   r < 0:  Disassortative mixing (like avoids like)
  #
  # High positive assortativity indicates that the similarity metric
  # successfully captures lineage-specific expression patterns, with
  # samples of the same cancer type being more similar to each other
  # than to samples of different types.
  # ---------------------------------------------------------------------------
  
  # Convert lineage labels to integer codes required by assortativity_nominal()
  lin <- V(g)$lineage
  lin_levels <- sort(unique(na.omit(lin)))
  lin_id <- match(lin, lin_levels)

  # Restrict to vertices with valid lineage annotations
  # assortativity_nominal() handles NAs poorly, so exclusion is necessary
  keep_v <- which(!is.na(lin_id))
  g_lin <- induced_subgraph(g, vids = keep_v)
  lin_id2 <- match(V(g_lin)$lineage, lin_levels)

  # Compute nominal assortativity coefficient
  assort <- assortativity_nominal(g_lin, types = lin_id2, directed = FALSE)
  
  assort_file <- file.path(opt$outdir, "validation_assortativity.tsv")
  fwrite(data.table(metric = "assortativity_nominal_lineage", value = assort), assort_file, sep = "\t")
  cat("  Saved: ", assort_file, "\n", sep = "")

  # ---------------------------------------------------------------------------
  # Community Purity Analysis
  # Purity measures the lineage homogeneity of each community, defined as
  # the fraction of samples belonging to the most common lineage within
  # that community.
  #
  #   Purity = (count of dominant lineage) / (total samples in community)
  #
  # Interpretation:
  #   Purity = 1.0:  All samples are from a single lineage (perfect coherence)
  #   Purity > 0.8:  Community is dominated by one lineage with minor mixing
  #   Purity < 0.5:  Community contains substantial multi-lineage membership
  #
  # Low-purity communities may represent:
  #   - Shared transcriptional programmes across cancer types
  #   - Batch effects causing spurious cross-lineage clustering
  #   - Suboptimal community detection parameters
  #   - Genuine biological ambiguity (e.g., poorly differentiated tumours)
  # ---------------------------------------------------------------------------
  
  # Compute community-by-lineage contingency counts (excluding UNKNOWN)
  comm_lin <- communities[lineage != "UNKNOWN", .N, by = .(community, lineage)]
  
  # Compute total samples per community (for purity denominator)
  comm_tot <- communities[lineage != "UNKNOWN", .N, by = community]
  setnames(comm_tot, "N", "n_in_comm")

  # Identify dominant lineage and compute purity for each community
  purity <- comm_lin[, .(
    max_n = max(N),
    top_lineage = lineage[which.max(N)]
  ), by = community]
  
  purity <- merge(purity, comm_tot, by = "community", all.x = TRUE)
  purity[, purity := max_n / n_in_comm]

  # Create enrichment table with full contingency information
  # This enables downstream visualisation and detailed analysis
  enrichment_file <- file.path(opt$outdir, "validation_lineage_enrichment.tsv")
  enrichment_out <- merge(
    comm_lin,
    purity[, .(community, n_in_comm, top_lineage, purity)],
    by = "community",
    all.x = TRUE
  )
  fwrite(enrichment_out[order(-purity, -n_in_comm)], enrichment_file, sep = "\t")
  cat("  Saved: ", enrichment_file, "\n", sep = "")

  # Display summary of top communities for quick quality assessment
  cat("\nLineage coherence quick view (top 10 communities by purity):\n")
  print(purity[order(-purity, -n_in_comm)][1:min(10, .N)])
}


# =============================================================================
# EXECUTION SUMMARY
# =============================================================================
# Final output summarising the completed analysis and locations of all
# generated files. This structured summary facilitates verification that
# all expected outputs were produced successfully.
# -----------------------------------------------------------------------------

cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")
cat("Primary communities written to: ", opt$out, "\n", sep = "")
cat("Validation tables written to:   ", opt$outdir, "\n", sep = "")
