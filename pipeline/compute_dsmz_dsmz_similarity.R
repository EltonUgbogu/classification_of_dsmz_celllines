#!/usr/bin/env Rscript

# compute_dsmz_dsmz_similarity.R
# Build DSMZ–DSMZ similarity graph from p_consensus(c, t) for a given direction
#
# Outputs (under {unsup_root}/tumour_neighbourhoods/{direction}/final_consensus):
#   DSMZ_DSMZ_similarity_matrix_{direction}.rds
#   DSMZ_DSMZ_similarity_pairs_{direction}.tsv
#   DSMZ_DSMZ_graph_edges_{direction}.tsv
#   DSMZ_DSMZ_graph_node_summary_{direction}.tsv
#   DSMZ_DSMZ_graph_node_annotations_{direction}.tsv
#   DSMZ_DSMZ_graph_community_summary_{direction}.tsv
#   DSMZ_DSMZ_Louvain_vs_Leiden_community_table_{direction}.tsv
#   Fig_DSMZ_DSMZ_similarity_histogram_{direction}.pdf
#   Fig_DSMZ_p_consensus_cell_scatter_{direction}.pdf
#   Fig_DSMZ_DSMZ_Louvain_vs_Leiden_heatmap_{direction}.pdf
#   Fig_DSMZ_DSMZ_graph_Leiden_{direction}.pdf
#   Fig_DSMZ_DSMZ_graph_Louvain_{direction}.pdf
#   Fig_DSMZ_DSMZ_graph_minimal_{direction}.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(viridis)
  library(pheatmap)
  library(grid)
  library(optparse)
  library(yaml)
})

# ----------------------------------------------------------------------
# 0) Config, direction, paths
# ----------------------------------------------------------------------
option_list <- list(
  make_option(
    "--config",
    type    = "character",
    default = "config/config.yaml",
    help    = "Path to config.yaml [default %default]"
  ),
  make_option(
    "--direction",
    type    = "character",
    default = NULL,
    help    = "Direction identifier (pam50_euc, pam50_corr, hvg_euc, hvg_corr)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
cfg <- yaml::read_yaml(opt$config)

if (is.null(opt$direction)) {
  stop("Please supply --direction.")
}

direction <- opt$direction

# Directions are config-driven (BRCA may have PAM50; NBL/RBL won't)
directions <- cfg$tumour_neighbourhoods$directions
if (is.null(directions) || length(directions) == 0) {
  directions <- c("hvg_euc", "hvg_corr")  # fallback
}

if (!direction %in% directions) {
  stop("Invalid --direction: ", direction,
       "\nAllowed (from config tumour_neighbourhoods$directions): ",
       paste(directions, collapse = ", "))
}

unsup_root <- cfg$paths$unsup_root

base_dir <- file.path(unsup_root, "tumour_neighbourhoods", direction)
cons_dir <- file.path(base_dir, "final_consensus")

in_rds <- file.path(
  cons_dir,
  sprintf("Final_consensus_tumour_neighbourhoods_%s.rds", direction)
)

cat("=== DSMZ–DSMZ similarity from p_consensus ===\n")
cat("Using config file:\n  ", opt$config, "\n")
cat("Direction: ", direction, "\n", sep = "")
cat("Input consensus file:\n  ", in_rds, "\n")
cat("Edge threshold will be determined from data (90th percentile)\n\n")

if (!file.exists(in_rds)) {
  stop("Consensus RDS not found: ", in_rds,
       "\nRun tumour_neighbourhood_p_consensus for direction ", direction, " first.")
}

consensus_pairs <- readRDS(in_rds)

cat("Summary of consensus_pairs:\n")
print(dplyr::glimpse(consensus_pairs))

# ----------------------------------------------------------------------
# 1) Build wide matrix: rows = cell lines, columns = TCGA tumours
#    values = p_consensus(c, t), NA -> 0
# ----------------------------------------------------------------------
mat_wide <- consensus_pairs %>%
  select(cell_line, tumor_id, p_consensus) %>%
  tidyr::pivot_wider(
    names_from  = tumor_id,
    values_from = p_consensus,
    values_fill = 0
  )

mat <- mat_wide %>% as.data.frame()
rownames(mat) <- mat$cell_line
mat$cell_line <- NULL
mat <- as.matrix(mat)

cat("\nMatrix dimensions (DSMZ x TCGA):\n")
print(dim(mat))

# ----------------------------------------------------------------------
# 2) DSMZ–DSMZ similarity = Pearson correlation over p_consensus vectors
# ----------------------------------------------------------------------
sim_mat <- cor(t(mat), method = "pearson", use = "pairwise.complete.obs")

cat("\nSimilarity matrix dimensions (DSMZ x DSMZ):\n")
print(dim(sim_mat))

# ----------------------------------------------------------------------
# 3) Export similarity matrix + long table
# ----------------------------------------------------------------------
out_rds_mat  <- file.path(
  cons_dir, sprintf("DSMZ_DSMZ_similarity_matrix_%s.rds", direction)
)
out_tsv_long <- file.path(
  cons_dir, sprintf("DSMZ_DSMZ_similarity_pairs_%s.tsv", direction)
)

saveRDS(sim_mat, out_rds_mat)

sim_long <- as.data.frame(as.table(sim_mat)) %>%
  mutate(across(c(Var1, Var2), as.character)) %>%
  rename(
    cell_line1 = Var1,
    cell_line2 = Var2,
    similarity = Freq
  ) %>%
  filter(cell_line1 < cell_line2)  # upper triangle only

readr::write_tsv(sim_long, out_tsv_long)

cat("\nSaved DSMZ–DSMZ similarity matrix to:\n  ", out_rds_mat, "\n")
cat("Saved DSMZ–DSMZ pairwise similarities to:\n  ", out_tsv_long, "\n")

# ----------------------------------------------------------------------
# 3b) Inspect similarity distribution and set data-driven threshold
# ----------------------------------------------------------------------
cat("\n=== DSMZ–DSMZ similarity summary ===\n")
print(summary(sim_long$similarity))

cat("\nTop 10 highest DSMZ–DSMZ similarities:\n")
sim_long %>%
  arrange(desc(similarity)) %>%
  slice(1:10) %>%
  print()

cat("\nCount of pairs by threshold (0.7):\n")
sim_long %>%
  mutate(ge_0.7 = similarity >= 0.7) %>%
  count(ge_0.7) %>%
  print()

edge_threshold <- quantile(sim_long$similarity, 0.9, na.rm = TRUE)
cat("\nData-driven edge threshold (90th percentile):", edge_threshold, "\n")
cat("Using this threshold for graph construction.\n")

# ----------------------------------------------------------------------
# 4) QC plot: similarity distribution
# ----------------------------------------------------------------------
hist_pdf <- file.path(
  cons_dir, sprintf("Fig_DSMZ_DSMZ_similarity_histogram_%s.pdf", direction)
)

p_hist <- ggplot(sim_long, aes(x = similarity)) +
  geom_histogram(
    binwidth = 0.05,
    colour   = "white",
    fill     = "#2c7bb6"
  ) +
  geom_vline(
    xintercept = edge_threshold,
    linetype   = "dashed",
    colour     = "red",
    linewidth  = 0.7
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_minimal(base_size = 14) +
  labs(
    title    = sprintf("Distribution of DSMZ-DSMZ similarity (%s)", direction),
    subtitle = "Pearson correlation of tumour neighbourhood p_consensus(c, t)",
    x        = "Similarity (Pearson r)",
    y        = "Number of DSMZ-DSMZ pairs"
  ) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(hist_pdf, p_hist, width = 7, height = 5)
cat("\nHistogram of DSMZ–DSMZ similarities saved to:\n  ", hist_pdf, "\n")

# ----------------------------------------------------------------------
# 4b) Per-cell-line anchoring strength scatter
# ----------------------------------------------------------------------
cat("\n=== Computing per-cell-line anchoring strength ===\n")

scatter_pdf <- file.path(
  cons_dir, sprintf("Fig_DSMZ_p_consensus_cell_scatter_%s.pdf", direction)
)

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

p_scatter <- ggplot(
  per_cell,
  aes(x = max_p_consensus, y = frac_ge_0_7)
) +
  geom_point(
    aes(size = n_tumours, colour = strong_anchor),
    alpha = 0.9
  ) +
  geom_label_repel(
    aes(label = cell_line),
    size          = 3,
    label.size    = 0,
    label.padding = grid::unit(0.15, "lines"),
    fill          = "white",
    alpha         = 0.9,
    max.overlaps  = Inf,
    box.padding   = 0.4,
    point.padding = 0.25
  ) +
  geom_vline(xintercept = 0.7, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0.7, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = 0.02)
  ) +
  scale_colour_manual(
    name   = "Max p_consensus",
    values = c(`FALSE` = "red3", `TRUE` = "#2c7bb6"),
    labels = c(`FALSE` = "< 0.7", `TRUE` = "≥ 0.7")
  ) +
  scale_size_continuous(
    name  = "Number of tumour neighbours",
    range = c(3, 7)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold"),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  labs(
    title    = sprintf("Anchoring strength of DSMZ breast cancer cell lines (%s)", direction),
    subtitle = "x: max p_consensus(c,t); y: fraction of tumour neighbours with p_consensus ≥ 0.7",
    x        = "Max p_consensus per cell line",
    y        = "Tumour neighbours with strong consensus (%)"
  )

ggsave(scatter_pdf, p_scatter, width = 8, height = 6)
cat("\nPer-cell p_consensus scatter saved to:\n  ", scatter_pdf, "\n")

# ----------------------------------------------------------------------
# 5) Build DSMZ–DSMZ graph
# ----------------------------------------------------------------------
edges_tsv <- file.path(
  cons_dir, sprintf("DSMZ_DSMZ_graph_edges_%s.tsv", direction)
)

graph_edges <- sim_long %>%
  filter(similarity >= edge_threshold) %>%
  arrange(desc(similarity))

if (nrow(graph_edges) == 0) {
  cat("\n[WARN] No DSMZ-DSMZ pairs with similarity >= ", edge_threshold, ".\n", sep = "")
  cat("Skipping graph construction. Check similarity summary / lower threshold.\n")

  readr::write_tsv(graph_edges, edges_tsv)
  cat("Empty edge list saved to:\n  ", edges_tsv, "\n")
  quit(save = "no", status = 0)
}

readr::write_tsv(graph_edges, edges_tsv)
cat("\nDSMZ–DSMZ graph edges (similarity >= ", edge_threshold, ") saved to:\n  ",
    edges_tsv, "\n", sep = "")

# ----------------------------------------------------------------------
# 5b) Node-level summary
# ----------------------------------------------------------------------
all_cell_lines <- rownames(sim_mat)

node_summary_edge <- bind_rows(
  graph_edges %>%
    select(cell_line = cell_line1, similarity),
  graph_edges %>%
    select(cell_line = cell_line2, similarity)
) %>%
  group_by(cell_line) %>%
  summarise(
    degree         = n(),
    mean_edge_sim  = mean(similarity),
    max_edge_sim   = max(similarity),
    .groups        = "drop"
  )

node_summary <- tibble(cell_line = all_cell_lines) %>%
  left_join(node_summary_edge, by = "cell_line") %>%
  mutate(
    degree        = dplyr::coalesce(degree, 0L),
    is_outlier    = degree == 0L
  ) %>%
  arrange(desc(degree), desc(mean_edge_sim))

nodes_tsv <- file.path(
  cons_dir, sprintf("DSMZ_DSMZ_graph_node_summary_%s.tsv", direction)
)
readr::write_tsv(node_summary, nodes_tsv)
cat("DSMZ-DSMZ graph node summary (including degree-0 outliers) saved to:\n  ",
    nodes_tsv, "\n")

# ----------------------------------------------------------------------
# 6) Graph object + annotations + community summary
# ----------------------------------------------------------------------
nodes_annot_tsv <- file.path(
  cons_dir, sprintf("DSMZ_DSMZ_graph_node_annotations_%s.tsv", direction)
)
comm_summary_tsv <- file.path(
  cons_dir, sprintf("DSMZ_DSMZ_graph_community_summary_%s.tsv", direction)
)

cat("\n=== Building DSMZ-DSMZ graph visualisation ===\n")

graph_tbl <- tidygraph::tbl_graph(
  nodes = node_summary %>% rename(name = cell_line),
  edges = graph_edges %>%
    rename(from = cell_line1, to = cell_line2),
  directed = FALSE
) %>%
  mutate(
    degree          = centrality_degree(),
    betweenness     = centrality_betweenness(),
    component       = as.factor(group_components()),
    community_louv  = as.factor(group_louvain(weights = similarity)),
    community_leid  = as.factor(group_leiden(weights = similarity, resolution = 1.0))
  )

ig <- as.igraph(graph_tbl)

comp <- igraph::components(ig)
cat("\nGraph summary:\n")
cat("  # nodes:               ", igraph::gorder(ig), "\n")
cat("  # edges:               ", igraph::gsize(ig),  "\n")
cat("  Density:               ", igraph::edge_density(ig), "\n")
cat("  # connected components:", comp$no, "\n")
cat("  Component sizes:       ", paste(sort(comp$csize, decreasing = TRUE), collapse = ", "), "\n\n")

cat("Community (Louvain) sizes:\n")
community_sizes_louv <- graph_tbl %>%
  as_tibble() %>%
  count(community_louv, name = "n") %>%
  arrange(desc(n))
print(community_sizes_louv)

cat("\nCommunity (Leiden) sizes:\n")
community_sizes_leid <- graph_tbl %>%
  as_tibble() %>%
  count(community_leid, name = "n") %>%
  arrange(desc(n))
print(community_sizes_leid)

node_annotations <- graph_tbl %>%
  as_tibble() %>%
  transmute(
    cell_line        = name,
    degree,
    betweenness,
    component        = as.character(component),
    community_louv   = as.character(community_louv),
    community_leid   = as.character(community_leid),
    mean_edge_sim    = mean_edge_sim,
    max_edge_sim     = max_edge_sim,
    is_outlier       = is_outlier
  ) %>%
  arrange(is_outlier, community_leid, desc(degree))

readr::write_tsv(node_annotations, nodes_annot_tsv)
cat("Full node annotations saved to:\n  ", nodes_annot_tsv, "\n")

community_summary <- node_annotations %>%
  group_by(community_leid) %>%
  summarise(
    n_members          = n(),
    members            = paste(sort(cell_line), collapse = ";"),
    mean_degree        = mean(degree),
    mean_mean_edge_sim = mean(mean_edge_sim, na.rm = TRUE),
    .groups            = "drop"
  ) %>%
  arrange(desc(n_members)) %>%
  rename(community_leiden = community_leid)

readr::write_tsv(community_summary, comm_summary_tsv)
cat("Community summary saved to:\n  ", comm_summary_tsv, "\n")

# ----------------------------------------------------------------------
# 7) Louvain vs Leiden overlap
# ----------------------------------------------------------------------
cat("\n=== Louvain vs Leiden community comparison ===\n")

comm_table <- table(
  Louvain = node_annotations$community_louv,
  Leiden  = node_annotations$community_leid
)

comm_table_tsv <- file.path(
  cons_dir, sprintf("DSMZ_DSMZ_Louvain_vs_Leiden_community_table_%s.tsv", direction)
)
readr::write_tsv(as.data.frame(comm_table), comm_table_tsv)
cat("Louvain vs Leiden community contingency table saved to:\n  ",
    comm_table_tsv, "\n")

comm_mat <- as.matrix(comm_table)

heatmap_pdf <- file.path(
  cons_dir, sprintf("Fig_DSMZ_DSMZ_Louvain_vs_Leiden_heatmap_%s.pdf", direction)
)

if (all(dim(comm_mat) > 0)) {
  pdf(heatmap_pdf, width = 7, height = 6)
  pheatmap::pheatmap(
    comm_mat,
    cluster_rows    = FALSE,
    cluster_cols    = FALSE,
    display_numbers = TRUE,
    number_format   = "%.0f",
    fontsize_number = 10,
    main            = "Overlap of Louvain vs Leiden communities",
    angle_col       = 45
  )
  dev.off()
  cat("Louvain vs Leiden community overlap heatmap saved to:\n  ",
      heatmap_pdf, "\n")
} else {
  cat("[WARN] Louvain vs Leiden table is empty – skipping heatmap.\n")
}

# ----------------------------------------------------------------------
# 8) Graph plots (Leiden, Louvain, minimal)
# ----------------------------------------------------------------------
graph_leiden_pdf <- file.path(
  cons_dir, sprintf("Fig_DSMZ_DSMZ_graph_Leiden_%s.pdf", direction)
)
graph_louvain_pdf <- file.path(
  cons_dir, sprintf("Fig_DSMZ_DSMZ_graph_Louvain_%s.pdf", direction)
)
graph_minimal_pdf <- file.path(
  cons_dir, sprintf("Fig_DSMZ_DSMZ_graph_minimal_%s.pdf", direction)
)

set.seed(123)

p_graph_leiden <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(
    aes(edge_alpha = similarity),
    show.legend = TRUE
  ) +
  scale_edge_alpha(
    name  = "Similarity (r)",
    range = c(0.15, 0.9)
  ) +
  geom_node_point(
    aes(size = degree, colour = community_leid, shape = is_outlier),
    alpha = 0.95
  ) +
  scale_shape_manual(
    name   = "Outlier (no edges)",
    values = c(`FALSE` = 19, `TRUE` = 17)
  ) +
  scale_size_continuous(
    name  = "Degree\n(# neighbours)",
    range = c(2.5, 8)
  ) +
  scale_colour_viridis_d(
    name      = "Leiden community",
    option    = "D",
    direction = 1
  ) +
  geom_node_text(
    aes(label = name),
    size          = 2.8,
    repel         = TRUE,
    colour        = "black",
    fontface      = "plain",
    point.padding = grid::unit(0.2, "lines")
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position   = "right",
    plot.title        = element_text(face = "bold", hjust = 0),
    plot.subtitle     = element_text(hjust = 0),
    plot.margin       = margin(10, 10, 10, 10),
    legend.box        = "vertical",
    legend.box.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    title    = sprintf("DSMZ-DSMZ similarity graph (Leiden, %s)", direction),
    subtitle = "Nodes = cell lines; edges = high similarity of tumour neighbourhood p_consensus(c,t)",
    caption  = "Node size: degree; node colour: Leiden community; triangles = outliers (degree 0); edge alpha: Pearson r"
  )

ggsave(graph_leiden_pdf, p_graph_leiden, width = 8, height = 6)
cat("\nLeiden graph figure saved to:\n  ", graph_leiden_pdf, "\n")

p_graph_louvain <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(
    aes(edge_alpha = similarity),
    show.legend = TRUE
  ) +
  scale_edge_alpha(
    name  = "Similarity (r)",
    range = c(0.15, 0.9)
  ) +
  geom_node_point(
    aes(size = degree, colour = community_louv, shape = is_outlier),
    alpha = 0.95
  ) +
  scale_shape_manual(
    name   = "Outlier (no edges)",
    values = c(`FALSE` = 19, `TRUE` = 17)
  ) +
  scale_size_continuous(
    name  = "Degree\n(# neighbours)",
    range = c(2.5, 8)
  ) +
  scale_colour_viridis_d(
    name      = "Louvain community",
    option    = "D",
    direction = 1
  ) +
  geom_node_text(
    aes(label = name),
    size          = 2.8,
    repel         = TRUE,
    colour        = "black",
    fontface      = "plain",
    point.padding = grid::unit(0.2, "lines")
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position   = "right",
    plot.title        = element_text(face = "bold", hjust = 0),
    plot.subtitle     = element_text(hjust = 0),
    plot.margin       = margin(10, 10, 10, 10),
    legend.box        = "vertical",
    legend.box.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    title    = sprintf("DSMZ-DSMZ similarity graph (Louvain, %s)", direction),
    subtitle = "Nodes = cell lines; edges = high similarity of tumour neighbourhood p_consensus(c,t)",
    caption  = "Node size: degree; node colour: Louvain community; triangles = outliers (degree 0); edge alpha: Pearson r"
  )

ggsave(graph_louvain_pdf, p_graph_louvain, width = 8, height = 6)
cat("\nLouvain graph figure saved to:\n  ", graph_louvain_pdf, "\n")

p_graph_minimal <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(
    aes(edge_alpha = similarity),
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(size = degree, colour = community_leid),
    alpha = 0.9
  ) +
  geom_node_text(
    aes(label = name),
    size          = 2.5,
    repel         = TRUE,
    colour        = "black",
    point.padding = grid::unit(0.2, "lines")
  ) +
  scale_size_continuous(range = c(2.5, 8)) +
  scale_colour_viridis_d() +
  scale_edge_alpha(range = c(0.2, 0.9)) +
  theme_void(base_size = 12) +
  labs(
    title    = sprintf("DSMZ-DSMZ similarity graph (minimal view, %s)", direction),
    subtitle = "Edge alpha: proportional to Pearson r between tumour neighbourhood profiles"
  )

ggsave(graph_minimal_pdf, p_graph_minimal, width = 8, height = 6)
cat("\nQC graph (minimal) saved to:\n  ", graph_minimal_pdf, "\n")

cat("\n=== Done: DSMZ-DSMZ graph from p_consensus (", direction, ") ===\n", sep = "")
