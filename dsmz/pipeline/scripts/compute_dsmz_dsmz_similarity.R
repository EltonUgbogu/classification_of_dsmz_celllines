# ==============================================================================
# compute_dsmz_dsmz_similarity.R
# Cell Line Similarity Network Construction from Consensus Neighbourhoods
# ==============================================================================
#
# PURPOSE:
# This script constructs a similarity network among DSMZ cell lines based on
# the overlap of their tumour neighbourhood profiles. The underlying principle
# is that cell lines which consistently associate with similar tumours
# across multiple clustering methods likely share biological characteristics.
#
# BIOLOGICAL RATIONALE:
# If two cell lines have highly correlated p_consensus vectors, they tend to
# be neighbours of the same tumours. This suggests they occupy similar regions
# of the molecular landscape and may represent similar cancer subtypes
# or share functional properties.

# OUTPUTS:
# - Similarity matrix (Pearson correlation of p_consensus profiles)
# - Graph with edges connecting similar cell lines
# - Community detection results (Louvain and Leiden algorithms)
# - visualisations
#
# ==============================================================================

# ------------------------------------------------------------------------------
# PACKAGE LOADING
# ------------------------------------------------------------------------------
# suppressPackageStartupMessages() prevents verbose startup messages from
# cluttering the console, which is useful in automated pipelines.

suppressPackageStartupMessages({
  library(dplyr)      # Data manipulation (filter, mutate, summarise, etc.)
  library(readr)      # Fast and consistent file I/O (read_csv, write_tsv)
  library(tidyr)      # Data reshaping (pivot_wider, pivot_longer)
  library(ggplot2)    # Grammar of graphics plotting system
  library(ggrepel)    # Automatic label placement to avoid overlaps
  library(purrr)      # Functional programming tools (map, reduce)
  library(igraph)     # Core graph/network analysis functions
  library(tidygraph)  # Tidy interface to igraph (dplyr-like syntax for graphs)
  library(ggraph)     # ggplot2 extension for graph visualisation
  library(viridis)    # Colourblind-friendly colour palettes
  library(pheatmap)   # heatmaps with clustering
  library(grid)       # Low-level graphics (used for fine-tuning plots)
  library(optparse)   # Command-line argument parsing
  library(yaml)       # YAML configuration file parsing
})

# ------------------------------------------------------------------------------
# COMMAND-LINE ARGUMENT PARSING
# ------------------------------------------------------------------------------
option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml",
              help = "Path to config.yaml [default %default]"),
  make_option("--direction", type = "character", default = NULL,
              help = "Direction identifier (pam50_euc, pam50_corr, hvg_euc, hvg_corr)")
)

opt <- parse_args(OptionParser(option_list = option_list))
cfg <- yaml::read_yaml(opt$config)

# Validate required arguments
if (is.null(opt$direction)) {
  stop("Please supply --direction.")
}

direction <- opt$direction

# Validate direction against allowed values from configuration
directions <- cfg$tumour_neighbourhoods$directions
if (is.null(directions) || length(directions) == 0) {
  directions <- c("hvg_euc", "hvg_corr")
}

if (!direction %in% directions) {
  stop("Invalid --direction: ", direction,
       "\nAllowed (from config tumour_neighbourhoods$directions): ",
       paste(directions, collapse = ", "))
}

# Construct file paths from configuration
unsup_root <- cfg$paths$unsup_root
base_dir <- file.path(unsup_root, "tumour_neighbourhoods", direction)
cons_dir <- file.path(base_dir, "final_consensus")

in_rds <- file.path(cons_dir, sprintf("Final_consensus_tumour_neighbourhoods_%s.rds", direction))

cat("=== DSMZ-DSMZ similarity from p_consensus ===\n")
cat("Using config file:\n  ", opt$config, "\n")
cat("Direction: ", direction, "\n", sep = "")
cat("Input consensus file:\n  ", in_rds, "\n")
cat("Edge threshold will be determined from data (90th percentile)\n\n")

# Validate input file exists
if (!file.exists(in_rds)) {
  stop("Consensus RDS not found: ", in_rds,
       "\nRun tumour_neighbourhood_p_consensus for direction ", direction, " first.")
}

# Load the consensus data (cell_line, tumor_id, p_consensus triplets)
consensus_pairs <- readRDS(in_rds)

cat("Summary of consensus_pairs:\n")
print(dplyr::glimpse(consensus_pairs))

# ------------------------------------------------------------------------------
# SECTION 1: CONSTRUCT P_CONSENSUS MATRIX
# ------------------------------------------------------------------------------
# The consensus_pairs data is in long format (one row per cell line-tumour pair).
# This section reshapes it into a wide matrix where:
#   - Rows = cell lines
#   - Columns = Tumours
#   - Values = p_consensus(c, t)
#
# This matrix representation allows each cell line to be treated as a vector
# in "tumour space", where the coordinates are the consensus probabilities
# with each tumour.

mat_wide <- consensus_pairs %>%
  select(cell_line, tumor_id, p_consensus) %>%
  # pivot_wider() transforms long data to wide format
  # names_from: column whose values become new column names
  # values_from: column whose values fill the new columns
  # values_fill: value to use for missing combinations (tumours not in neighbourhood)
  tidyr::pivot_wider(
    names_from  = tumor_id,
    values_from = p_consensus,
    values_fill = 0  # Cell lines with no association to a tumour get 0
  )

# Convert tibble to matrix for correlation computation
mat <- mat_wide %>% as.data.frame()
rownames(mat) <- mat$cell_line
mat$cell_line <- NULL
mat <- as.matrix(mat)

cat("\nMatrix dimensions (DSMZ x Tumour):\n")
print(dim(mat))

# ------------------------------------------------------------------------------
# SECTION 2: COMPUTE CELL LINE SIMILARITY
# ------------------------------------------------------------------------------
# Similarity between two cell lines is defined as the Pearson correlation of
# their p_consensus vectors across all Tumours.
#
# INTERPRETATION:
#   r ≈ 1:  Cell lines have very similar tumour neighbourhood profiles
#   r ≈ 0:  Cell lines have uncorrelated (independent) profiles
#   r ≈ -1: Cell lines have opposite profiles (rare in this context)
#
# The cor() function computes pairwise correlations.
# t(mat) transposes the matrix so that cor() computes correlations between
# rows (cell lines) rather than columns (tumours).

sim_mat <- cor(t(mat), method = "pearson", use = "pairwise.complete.obs")

cat("\nSimilarity matrix dimensions (DSMZ x DSMZ):\n")
print(dim(sim_mat))

# ------------------------------------------------------------------------------
# SECTION 3: EXPORT SIMILARITY DATA
# ------------------------------------------------------------------------------

out_rds_mat  <- file.path(cons_dir, sprintf("DSMZ_DSMZ_similarity_matrix_%s.rds", direction))
out_tsv_long <- file.path(cons_dir, sprintf("DSMZ_DSMZ_similarity_pairs_%s.tsv", direction))

saveRDS(sim_mat, out_rds_mat)

# Convert symmetric matrix to long format for easier inspection and plotting
# as.table() converts matrix to table, as.data.frame() to data frame
# filter(cell_line1 < cell_line2) keeps only upper triangle (avoids duplicates)
sim_long <- as.data.frame(as.table(sim_mat)) %>%
  mutate(across(c(Var1, Var2), as.character)) %>%
  rename(cell_line1 = Var1, cell_line2 = Var2, similarity = Freq) %>%
  filter(cell_line1 < cell_line2)

readr::write_tsv(sim_long, out_tsv_long)

cat("\nSaved DSMZ-DSMZ similarity matrix to:\n  ", out_rds_mat, "\n")
cat("Saved DSMZ-DSMZ pairwise similarities to:\n  ", out_tsv_long, "\n")

# ------------------------------------------------------------------------------
# SECTION 3B: DATA-DRIVEN EDGE THRESHOLD
# ------------------------------------------------------------------------------
# Rather than using an arbitrary threshold for graph edges, a data-driven
# approach selects the 90th percentile of similarity values. This ensures
# that only the strongest associations are represented as edges, while
# adapting to the specific distribution of the data.

cat("\n=== DSMZ-DSMZ similarity summary ===\n")
print(summary(sim_long$similarity))

cat("\nTop 10 highest DSMZ-DSMZ similarities:\n")
sim_long %>% arrange(desc(similarity)) %>% slice(1:10) %>% print()

cat("\nCount of pairs by threshold (0.7):\n")
sim_long %>% mutate(ge_0.7 = similarity >= 0.7) %>% count(ge_0.7) %>% print()

# quantile() computes the specified percentile of the distribution
edge_threshold <- quantile(sim_long$similarity, 0.9, na.rm = TRUE)
cat("\nData-driven edge threshold (90th percentile):", edge_threshold, "\n")

# ------------------------------------------------------------------------------
# SECTION 4: SIMILARITY DISTRIBUTION HISTOGRAM
# ------------------------------------------------------------------------------
# Visualising the distribution of pairwise similarities helps assess:
#   - Whether most cell lines are weakly or strongly correlated
#   - The appropriateness of the chosen edge threshold
#   - Potential outliers or bimodal structure

hist_pdf <- file.path(cons_dir, sprintf("Fig_DSMZ_DSMZ_similarity_histogram_%s.pdf", direction))

p_hist <- ggplot(sim_long, aes(x = similarity)) +
  # geom_histogram() creates the histogram
  # binwidth controls the width of each bar
  geom_histogram(binwidth = 0.05, colour = "white", fill = "#2c7bb6") +
  # Add vertical line at the threshold for reference
  geom_vline(xintercept = edge_threshold, linetype = "dashed", colour = "red", linewidth = 0.7) +
  # coord_cartesian() sets axis limits without removing data
  coord_cartesian(xlim = c(-1, 1)) +
  theme_minimal(base_size = 14) +
  labs(
    title    = sprintf("Distribution of DSMZ-DSMZ similarity (%s)", direction),
    subtitle = "Pearson correlation of tumour neighbourhood p_consensus(c, t)",
    x = "Similarity (Pearson r)",
    y = "Number of DSMZ-DSMZ pairs"
  ) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(hist_pdf, p_hist, width = 7, height = 5)
cat("\nHistogram saved to:\n  ", hist_pdf, "\n")

# ------------------------------------------------------------------------------
# SECTION 4B: PER-CELL-LINE ANCHORING STRENGTH
# ------------------------------------------------------------------------------
# This scatter plot characterises how reliably each cell line maps to
# tumour space. Cell lines in the upper-right quadrant have both strong
# maximum consensus and many strongly-anchored neighbours.

scatter_pdf <- file.path(cons_dir, sprintf("Fig_DSMZ_p_consensus_cell_scatter_%s.pdf", direction))

# Compute summary statistics per cell line
per_cell <- consensus_pairs %>%
  group_by(cell_line) %>%
  summarise(
    max_p_consensus = max(p_consensus, na.rm = TRUE),
    n_tumours       = n(),
    frac_ge_0_7     = mean(p_consensus >= 0.7, na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  mutate(strong_anchor = max_p_consensus >= 0.7)

cat("\nCell-line anchoring summary:\n")
print(per_cell %>% arrange(desc(max_p_consensus)))

p_scatter <- ggplot(per_cell, aes(x = max_p_consensus, y = frac_ge_0_7)) +
  geom_point(aes(size = n_tumours, colour = strong_anchor), alpha = 0.9) +
  # geom_label_repel() adds labels that automatically avoid overlapping
  geom_label_repel(
    aes(label = cell_line), size = 3, label.size = 0,
    label.padding = grid::unit(0.15, "lines"), fill = "white", alpha = 0.9,
    max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25
  ) +
  # Reference lines at 0.7 threshold
  geom_vline(xintercept = 0.7, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0.7, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(
    name = "Max p_consensus",
    values = c(`FALSE` = "red3", `TRUE` = "#2c7bb6"),
    labels = c(`FALSE` = "< 0.7", `TRUE` = ">= 0.7")
  ) +
  scale_size_continuous(name = "Number of tumour neighbours", range = c(3, 7)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "right",
        plot.title = element_text(face = "bold")) +
  labs(
    title = sprintf("Anchoring strength of DSMZ cancer cell lines (%s)", direction),
    subtitle = "x: max p_consensus(c,t); y: fraction of tumour neighbours with p_consensus >= 0.7",
    x = "Max p_consensus per cell line",
    y = "Tumour neighbours with strong consensus (%)"
  )

ggsave(scatter_pdf, p_scatter, width = 8, height = 6)
cat("\nPer-cell scatter saved to:\n  ", scatter_pdf, "\n")

# ------------------------------------------------------------------------------
# SECTION 5: CONSTRUCT SIMILARITY GRAPH
# ------------------------------------------------------------------------------
# A graph is constructed where:
#   - Nodes = cell lines
#   - Edges = similarity above threshold
#   - Edge weights = similarity values
#
# This representation enables network analysis techniques such as community
# detection to identify groups of functionally related cell lines.

edges_tsv <- file.path(cons_dir, sprintf("DSMZ_DSMZ_graph_edges_%s.tsv", direction))

# Filter to pairs above threshold
graph_edges <- sim_long %>%
  filter(similarity >= edge_threshold) %>%
  arrange(desc(similarity))

# Handle edge case of no edges above threshold
if (nrow(graph_edges) == 0) {
  cat("\n[WARN] No DSMZ-DSMZ pairs with similarity >= ", edge_threshold, ".\n", sep = "")
  cat("Skipping graph construction.\n")
  readr::write_tsv(graph_edges, edges_tsv)
  quit(save = "no", status = 0)
}

readr::write_tsv(graph_edges, edges_tsv)
cat("\nGraph edges (similarity >= ", edge_threshold, ") saved to:\n  ", edges_tsv, "\n", sep = "")

# ------------------------------------------------------------------------------
# SECTION 5B: NODE-LEVEL SUMMARY
# ------------------------------------------------------------------------------
# For each cell line, compute:
#   - degree: number of edges (connected neighbours)
#   - mean_edge_sim: average similarity to connected neighbours
#   - is_outlier: whether the cell line has no connections

all_cell_lines <- rownames(sim_mat)

# Compute statistics for nodes that have edges
node_summary_edge <- bind_rows(
  graph_edges %>% select(cell_line = cell_line1, similarity),
  graph_edges %>% select(cell_line = cell_line2, similarity)
) %>%
  group_by(cell_line) %>%
  summarise(degree = n(), mean_edge_sim = mean(similarity),
            max_edge_sim = max(similarity), .groups = "drop")

# Join with all cell lines to include isolated nodes (degree 0)
node_summary <- tibble(cell_line = all_cell_lines) %>%
  left_join(node_summary_edge, by = "cell_line") %>%
  # coalesce() replaces NA with specified value (0 for isolated nodes)
  mutate(degree = dplyr::coalesce(degree, 0L), is_outlier = degree == 0L) %>%
  arrange(desc(degree), desc(mean_edge_sim))

nodes_tsv <- file.path(cons_dir, sprintf("DSMZ_DSMZ_graph_node_summary_%s.tsv", direction))
readr::write_tsv(node_summary, nodes_tsv)
cat("Node summary saved to:\n  ", nodes_tsv, "\n")

# ------------------------------------------------------------------------------
# SECTION 6: GRAPH CONSTRUCTION WITH COMMUNITY DETECTION
# ------------------------------------------------------------------------------
# tidygraph provides a tidy interface to igraph, allowing dplyr-like syntax
# for graph manipulation. Community detection algorithms identify densely
# connected subgroups within the network.
#
# LOUVAIN ALGORITHM:
# A greedy modularity optimisation algorithm. Fast and widely used.
#
# LEIDEN ALGORITHM:
# An improved version of Louvain that guarantees well-connected communities.
# Generally preferred for biological networks.

nodes_annot_tsv <- file.path(cons_dir, sprintf("DSMZ_DSMZ_graph_node_annotations_%s.tsv", direction))
comm_summary_tsv <- file.path(cons_dir, sprintf("DSMZ_DSMZ_graph_community_summary_%s.tsv", direction))

cat("\n=== Building DSMZ-DSMZ graph ===\n")

# tbl_graph() creates a tidygraph object from nodes and edges data frames
graph_tbl <- tidygraph::tbl_graph(
  nodes = node_summary %>% rename(name = cell_line),
  edges = graph_edges %>% rename(from = cell_line1, to = cell_line2),
  directed = FALSE
) %>%
  # mutate() on a tbl_graph operates on the active element (nodes by default)
  mutate(
    # centrality_degree(): number of edges connected to each node
    degree          = centrality_degree(),
    # centrality_betweenness(): how often a node lies on shortest paths
    betweenness     = centrality_betweenness(),
    # group_components(): identifies disconnected subgraphs
    component       = as.factor(group_components()),
    # group_louvain(): Louvain community detection
    community_louv  = as.factor(group_louvain(weights = similarity)),
    # group_leiden(): Leiden community detection (resolution controls granularity)
    community_leid  = as.factor(group_leiden(weights = similarity, resolution = 1.0))
  )

# Convert to igraph for additional summary statistics
ig <- as.igraph(graph_tbl)
comp <- igraph::components(ig)

cat("\nGraph summary:\n")
cat("  # nodes:               ", igraph::gorder(ig), "\n")
cat("  # edges:               ", igraph::gsize(ig), "\n")
# edge_density(): proportion of possible edges that exist
cat("  Density:               ", igraph::edge_density(ig), "\n")
cat("  # connected components:", comp$no, "\n")

cat("\nCommunity (Louvain) sizes:\n")
graph_tbl %>% as_tibble() %>% count(community_louv, name = "n") %>%
  arrange(desc(n)) %>% print()

cat("\nCommunity (Leiden) sizes:\n")
graph_tbl %>% as_tibble() %>% count(community_leid, name = "n") %>%
  arrange(desc(n)) %>% print()

# Export detailed node annotations
node_annotations <- graph_tbl %>%
  as_tibble() %>%
  transmute(
    cell_line = name, degree, betweenness,
    component = as.character(component),
    community_louv = as.character(community_louv),
    community_leid = as.character(community_leid),
    mean_edge_sim, max_edge_sim, is_outlier
  ) %>%
  arrange(is_outlier, community_leid, desc(degree))

readr::write_tsv(node_annotations, nodes_annot_tsv)
cat("Node annotations saved to:\n  ", nodes_annot_tsv, "\n")

# Export community-level summary
community_summary <- node_annotations %>%
  group_by(community_leid) %>%
  summarise(
    n_members = n(),
    # paste() with collapse concatenates member names into single string
    members = paste(sort(cell_line), collapse = ";"),
    mean_degree = mean(degree),
    mean_mean_edge_sim = mean(mean_edge_sim, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_members)) %>%
  rename(community_leiden = community_leid)

readr::write_tsv(community_summary, comm_summary_tsv)
cat("Community summary saved to:\n  ", comm_summary_tsv, "\n")

# ------------------------------------------------------------------------------
# SECTION 7: LOUVAIN VS LEIDEN COMPARISON
# ------------------------------------------------------------------------------
# Comparing community assignments from different algorithms assesses the
# robustness of the detected structure. High agreement suggests stable
# community structure; disagreement may indicate ambiguous boundaries.

cat("\n=== Louvain vs Leiden comparison ===\n")

# table() creates a contingency table of community assignments
comm_table <- table(
  Louvain = node_annotations$community_louv,
  Leiden  = node_annotations$community_leid
)

comm_table_tsv <- file.path(cons_dir,
  sprintf("DSMZ_DSMZ_Louvain_vs_Leiden_community_table_%s.tsv", direction))
readr::write_tsv(as.data.frame(comm_table), comm_table_tsv)
cat("Louvain vs Leiden contingency table saved to:\n  ", comm_table_tsv, "\n")

# Create heatmap of community overlap
comm_mat <- as.matrix(comm_table)
heatmap_pdf <- file.path(cons_dir,
  sprintf("Fig_DSMZ_DSMZ_Louvain_vs_Leiden_heatmap_%s.pdf", direction))

if (all(dim(comm_mat) > 0)) {
  pdf(heatmap_pdf, width = 7, height = 6)
  # pheatmap() creates annotated heatmaps with optional clustering
  pheatmap::pheatmap(
    comm_mat, cluster_rows = FALSE, cluster_cols = FALSE,
    display_numbers = TRUE, number_format = "%.0f", fontsize_number = 10,
    main = "Overlap of Louvain vs Leiden communities", angle_col = 45
  )
  dev.off()
  cat("Community overlap heatmap saved to:\n  ", heatmap_pdf, "\n")
}

# ------------------------------------------------------------------------------
# SECTION 8: GRAPH VISUALISATIONS
# ------------------------------------------------------------------------------
# ggraph extends ggplot2 for graph visualisation. The Fruchterman-Reingold
# layout ("fr") positions nodes to minimise edge crossings while keeping
# connected nodes close together.

graph_leiden_pdf <- file.path(cons_dir, sprintf("Fig_DSMZ_DSMZ_graph_Leiden_%s.pdf", direction))
graph_louvain_pdf <- file.path(cons_dir, sprintf("Fig_DSMZ_DSMZ_graph_Louvain_%s.pdf", direction))
graph_minimal_pdf <- file.path(cons_dir, sprintf("Fig_DSMZ_DSMZ_graph_minimal_%s.pdf", direction))

# Set seed for reproducible layout
set.seed(123)

# --------------------------------------------------------------------------
# LEIDEN-COLOURED GRAPH
# --------------------------------------------------------------------------
p_graph_leiden <- ggraph(graph_tbl, layout = "fr") +
  # geom_edge_link() draws edges between nodes
  # edge_alpha maps similarity to transparency
  geom_edge_link(aes(edge_alpha = similarity), show.legend = TRUE) +
  scale_edge_alpha(name = "Similarity (r)", range = c(0.15, 0.9)) +
  # geom_node_point() draws nodes
  # size = degree, colour = community, shape distinguishes outliers
  geom_node_point(aes(size = degree, colour = community_leid, shape = is_outlier), alpha = 0.95) +
  scale_shape_manual(name = "Outlier (no edges)", values = c(`FALSE` = 19, `TRUE` = 17)) +
  scale_size_continuous(name = "Degree", range = c(2.5, 8)) +
  # viridis colour scale is colourblind-friendly
  scale_colour_viridis_d(name = "Leiden community", option = "D") +
  # geom_node_text() with repel = TRUE avoids label overlaps
  geom_node_text(aes(label = name), size = 2.8, repel = TRUE, colour = "black") +
  theme_void(base_size = 14) +
  theme(legend.position = "right", plot.title = element_text(face = "bold")) +
  labs(
    title = sprintf("DSMZ-DSMZ similarity graph (Leiden, %s)", direction),
    subtitle = "Edges connect cell lines with similar tumour neighbourhood profiles"
  )

ggsave(graph_leiden_pdf, p_graph_leiden, width = 8, height = 6)
cat("\nLeiden graph saved to:\n  ", graph_leiden_pdf, "\n")

# --------------------------------------------------------------------------
# LOUVAIN-COLOURED GRAPH
# --------------------------------------------------------------------------
p_graph_louvain <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity), show.legend = TRUE) +
  scale_edge_alpha(name = "Similarity (r)", range = c(0.15, 0.9)) +
  geom_node_point(aes(size = degree, colour = community_louv, shape = is_outlier), alpha = 0.95) +
  scale_shape_manual(name = "Outlier (no edges)", values = c(`FALSE` = 19, `TRUE` = 17)) +
  scale_size_continuous(name = "Degree", range = c(2.5, 8)) +
  scale_colour_viridis_d(name = "Louvain community", option = "D") +
  geom_node_text(aes(label = name), size = 2.8, repel = TRUE, colour = "black") +
  theme_void(base_size = 14) +
  theme(legend.position = "right", plot.title = element_text(face = "bold")) +
  labs(
    title = sprintf("DSMZ-DSMZ similarity graph (Louvain, %s)", direction),
    subtitle = "Edges connect cell lines with similar tumour neighbourhood profiles"
  )

ggsave(graph_louvain_pdf, p_graph_louvain, width = 8, height = 6)
cat("\nLouvain graph saved to:\n  ", graph_louvain_pdf, "\n")

# --------------------------------------------------------------------------
# MINIMAL GRAPH (CLEAN VIEW WITHOUT LEGENDS)
# --------------------------------------------------------------------------
p_graph_minimal <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(aes(edge_alpha = similarity), show.legend = FALSE) +
  geom_node_point(aes(size = degree, colour = community_leid), alpha = 0.9) +
  geom_node_text(aes(label = name), size = 2.5, repel = TRUE, colour = "black") +
  scale_size_continuous(range = c(2.5, 8)) +
  scale_colour_viridis_d() +
  scale_edge_alpha(range = c(0.2, 0.9)) +
  theme_void(base_size = 12) +
  labs(
    title = sprintf("DSMZ-DSMZ similarity graph (minimal, %s)", direction),
    subtitle = "Edge alpha proportional to Pearson r"
  )

ggsave(graph_minimal_pdf, p_graph_minimal, width = 8, height = 6)
cat("\nMinimal graph saved to:\n  ", graph_minimal_pdf, "\n")

cat("\n=== Done: DSMZ-DSMZ graph from p_consensus (", direction, ") ===\n", sep = "")