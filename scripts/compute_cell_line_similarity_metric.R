#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(igraph)
  library(optparse)
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml"),
  make_option("--profile", type = "character", default = NULL),
  make_option("--direction", type = "character", default = NULL),
  make_option("--out_base", type = "character", default = NULL),
  make_option("--similarity_metric", type = "character", default = "pearson"),
  make_option("--consensus_threshold", type = "double", default = NA_real_),
  make_option("--similarity_quantile", type = "double", default = NA_real_)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$direction) || !nzchar(opt$direction)) {
  stop("Please supply --direction.", call. = FALSE)
}
if (is.null(opt$out_base) || !nzchar(opt$out_base)) {
  stop("Please supply --out_base.", call. = FALSE)
}
if (!opt$similarity_metric %in% c("pearson", "jaccard")) {
  stop("--similarity_metric must be 'pearson' or 'jaccard'", call. = FALSE)
}

config_dir <- dirname(normalizePath(opt$config))
script_dir <- file.path(config_dir, "..", "scripts")
lib_config_path <- file.path(script_dir, "lib_config.R")
similarity_utils_path <- file.path(config_dir, "..", "R", "patient_referenced_similarity_utils.R")
if (file.exists(lib_config_path)) {
  source(lib_config_path)
  profile_override <- if (!is.null(opt$profile)) opt$profile else Sys.getenv("SNAKEMAKE_PROFILE", unset = "")
  if (identical(profile_override, "")) profile_override <- NULL
  cfg <- read_profiled_config(opt$config, profile_override = profile_override)
} else {
  cfg <- yaml::read_yaml(opt$config)
}
if (!file.exists(similarity_utils_path)) {
  stop("Missing required helper: ", similarity_utils_path, call. = FALSE)
}
source(similarity_utils_path)

strip_prefix <- function(x) sub("^(CELL:|TUMOUR:|TUMOR:)", "", as.character(x))

collapse_semicolon_ids <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  ids <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  ids <- trimws(ids)
  ids <- sort(unique(ids[nzchar(ids)]))
  if (!length(ids)) NA_character_ else paste(ids, collapse = ";")
}

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) NA_character_ else x[[1]]
}

direction <- opt$direction
graph_cfg <- cfg$patient_referenced_graph
if (is.null(graph_cfg)) {
  stop("Missing required config section: patient_referenced_graph", call. = FALSE)
}
required_graph_keys <- c("p_consensus_threshold", "similarity_quantile")
missing_graph_keys <- required_graph_keys[!required_graph_keys %in% names(graph_cfg)]
if (length(missing_graph_keys) > 0) {
  stop(
    "Missing required patient_referenced_graph key(s): ",
    paste(missing_graph_keys, collapse = ", "),
    call. = FALSE
  )
}
threshold <- if (is.finite(opt$consensus_threshold)) {
  opt$consensus_threshold
} else {
  graph_cfg$p_consensus_threshold
}
similarity_quantile <- if (is.finite(opt$similarity_quantile)) {
  opt$similarity_quantile
} else {
  graph_cfg$similarity_quantile
}
if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) || threshold <= 0 || threshold > 1) {
  stop("p_consensus_threshold must satisfy 0 < threshold <= 1.", call. = FALSE)
}
if (!is.numeric(similarity_quantile) || length(similarity_quantile) != 1L || is.na(similarity_quantile) ||
    similarity_quantile <= 0 || similarity_quantile >= 1) {
  stop("similarity_quantile must satisfy 0 < q < 1.", call. = FALSE)
}

pipe_root <- normalizePath(file.path(dirname(opt$config), ".."))
abs_from_root <- function(p) {
  if (is.null(p) || is.na(p) || !nzchar(p)) return(p)
  if (grepl("^/", p) || grepl("^[A-Za-z]:[\\\\/]", p)) return(p)
  file.path(pipe_root, p)
}

unsup_root <- abs_from_root(cfg$paths$unsup_root)
consensus_rds <- file.path(
  unsup_root,
  "tumour_neighbourhoods",
  direction,
  "final_consensus",
  sprintf("Final_consensus_tumour_neighbourhoods_%s.rds", direction)
)
if (!file.exists(consensus_rds)) {
  stop("Consensus RDS not found: ", consensus_rds, call. = FALSE)
}

out_base <- opt$out_base
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

consensus_pairs <- readRDS(consensus_rds) %>%
  {
    df <- .
    if (!"tumour_id" %in% colnames(df) && "tumor_id" %in% colnames(df)) {
      df <- df %>% rename(tumour_id = tumor_id)
    }
    if (!"cell_line" %in% colnames(df) && "cell_tech_id" %in% colnames(df)) {
      df <- df %>% mutate(cell_line = cell_tech_id)
    }
    df
  } %>%
  mutate(
    cell_line = strip_prefix(cell_line),
    tumour_id = strip_prefix(tumour_id),
    sample_id = if ("sample_id" %in% colnames(.)) strip_prefix(sample_id) else strip_prefix(cell_tech_id),
    cell_tech_id = if ("cell_tech_id" %in% colnames(.)) strip_prefix(cell_tech_id) else as.character(cell_line),
    cell_line_display = if ("cell_line_display" %in% colnames(.)) as.character(cell_line_display) else as.character(cell_line)
  )

cp <- consensus_pairs %>%
  mutate(
    profile_key = dplyr::coalesce(
      as.character(.data$sample_id),
      as.character(.data$cell_tech_id),
      as.character(.data$cell_line)
    )
  )

cell_label_map <- cp %>%
  group_by(cell_line) %>%
  summarise(
    sample_id = collapse_semicolon_ids(sample_id),
    cell_tech_id = collapse_semicolon_ids(profile_key),
    cell_line_display = first_non_missing(cell_line_display),
    .groups = "drop"
  )

profile_long <- cp %>%
  select(profile_key, cell_line, tumour_id, p_consensus, n_methods) %>%
  group_by(profile_key, cell_line, tumour_id) %>%
  summarise(
    p_consensus = mean(p_consensus, na.rm = TRUE),
    n_methods = max(n_methods, na.rm = TRUE),
    .groups = "drop"
  )

profile_wide <- profile_long %>%
  select(profile_key, cell_line, tumour_id, p_consensus) %>%
  pivot_wider(names_from = tumour_id, values_from = p_consensus, values_fill = 0)

profile_meta <- profile_wide %>% select(profile_key, cell_line)
profile_mat <- profile_wide %>% select(-profile_key, -cell_line) %>% as.data.frame() %>% as.matrix()
cell_line_counts <- table(profile_meta$cell_line)
mean_mat <- rowsum(profile_mat, group = profile_meta$cell_line, reorder = TRUE)
mean_mat <- sweep(mean_mat, 1, as.numeric(cell_line_counts[rownames(mean_mat)]), "/")

thresholded_mat <- threshold_restrict_mean_profiles(mean_mat, threshold)
binary_mat <- ifelse(thresholded_mat > 0, 1, 0)
selected_tumour_count <- rowSums(binary_mat, na.rm = TRUE)

selected_tumour_stats <- tibble(
  cohort = cfg$profile %||% opt$profile %||% "",
  direction = direction,
  similarity_metric = opt$similarity_metric,
  cell_line = rownames(mean_mat),
  n_clustering_methods = profile_long %>%
    group_by(cell_line) %>%
    summarise(n_clustering_methods = max(n_methods, na.rm = TRUE), .groups = "drop") %>%
    right_join(tibble(cell_line = rownames(mean_mat)), by = "cell_line") %>%
    pull(n_clustering_methods),
  p_consensus_threshold = threshold,
  minimum_method_count = ceiling(threshold * n_clustering_methods),
  n_candidate_tumours = ncol(mean_mat),
  n_tumours_ge_threshold = selected_tumour_count,
  fraction_tumours_ge_threshold = rowMeans(binary_mat, na.rm = TRUE)
)

nodes <- rownames(binary_mat)
sim_mat <- matrix(NA_real_, nrow = nrow(binary_mat), ncol = nrow(binary_mat), dimnames = list(nodes, nodes))
diag(sim_mat) <- 1
pair_detail_rows <- vector("list", length = 0L)

if (nrow(binary_mat) >= 2) {
  for (idx in combn(seq_len(nrow(binary_mat)), 2, simplify = FALSE)) {
    i <- idx[[1]]
    j <- idx[[2]]
    pair_result <- compute_pairwise_active_union_similarity(
      thresholded_mat[i, ],
      thresholded_mat[j, ],
      metric = opt$similarity_metric
    )

    sim_mat[i, j] <- pair_result$similarity
    sim_mat[j, i] <- pair_result$similarity
    pair_detail_rows[[length(pair_detail_rows) + 1L]] <- tibble(
      cell_line1 = nodes[[i]],
      cell_line2 = nodes[[j]],
      similarity = pair_result$similarity,
      n_pairwise_active_tumours = pair_result$n_pairwise_active_tumours,
      n_shared_selected_tumours = pair_result$n_shared_selected_tumours,
      undefined_similarity_reason = pair_result$undefined_similarity_reason
    )
  }
}

if (opt$similarity_metric == "pearson") {
  similarity_definition <- sprintf(
    "Pearson correlation of threshold-restricted p_consensus fractions over the pairwise active-tumour union defined by p_consensus >= %.2f after biological-cell-line mean pooling.",
    threshold
  )
  similarity_question <- sprintf(
    "Among patient tumours selected by at least one of the two biological cell lines under p_consensus >= %.2f, how concordant are the threshold-restricted tumour-wise p_consensus fractions?",
    threshold
  )
  joint_zero_contribution <- "Excluded by restricting Pearson correlation to the pairwise active-tumour union."
  fraction_magnitude_used <- TRUE
  vector_semantics <- "thresholded_mean_p_consensus_profile_over_pairwise_active_tumour_union"
} else {
  similarity_definition <- "Binary Jaccard overlap of threshold-selected tumour memberships after biological-cell-line mean pooling."
  similarity_question <- sprintf(
    "Do two biological cell lines select the same strong patient tumours under p_consensus >= %.2f?",
    threshold
  )
  joint_zero_contribution <- "Excluded by restricting Jaccard similarity to the pairwise active-tumour union."
  fraction_magnitude_used <- FALSE
  vector_semantics <- "thresholded_binary_mean_p_consensus_profile_over_pairwise_active_tumour_union"
}

sim_long <- if (length(pair_detail_rows) > 0) {
  bind_rows(pair_detail_rows)
} else {
  tibble(
    cell_line1 = character(),
    cell_line2 = character(),
    similarity = numeric(),
    n_pairwise_active_tumours = integer(),
    n_shared_selected_tumours = integer(),
    undefined_similarity_reason = character()
  )
} %>%
  left_join(cell_label_map %>% select(cell_line1 = cell_line, sample_id1 = sample_id, cell_line1_display = cell_line_display), by = "cell_line1") %>%
  left_join(cell_label_map %>% select(cell_line2 = cell_line, sample_id2 = sample_id, cell_line2_display = cell_line_display), by = "cell_line2") %>%
  mutate(similarity_metric = opt$similarity_metric)

selected_count_map <- tibble(cell_line = names(selected_tumour_count), selected_tumour_count = as.integer(selected_tumour_count))

sim_long <- sim_long %>%
  left_join(selected_count_map %>% rename(cell_line1 = cell_line, cell_line1_selected_tumours = selected_tumour_count), by = "cell_line1") %>%
  left_join(selected_count_map %>% rename(cell_line2 = cell_line, cell_line2_selected_tumours = selected_tumour_count), by = "cell_line2")

defined_sims <- sim_long %>% filter(!is.na(similarity))
undefined_pairs <- sim_long %>% filter(is.na(similarity))
edge_threshold <- if (nrow(defined_sims) > 0) {
  as.numeric(quantile(defined_sims$similarity, similarity_quantile, na.rm = TRUE, names = FALSE))
} else {
  NA_real_
}
graph_edges <- if (is.finite(edge_threshold)) {
  sim_long %>% filter(!is.na(similarity), similarity >= edge_threshold) %>% arrange(desc(similarity))
} else {
  sim_long %>% slice(0)
}

all_cell_lines <- rownames(sim_mat)
node_summary_edge <- bind_rows(
  graph_edges %>% select(cell_line = cell_line1, similarity),
  graph_edges %>% select(cell_line = cell_line2, similarity)
) %>%
  group_by(cell_line) %>%
  summarise(
    degree = n(),
    mean_edge_sim = mean(similarity),
    max_edge_sim = max(similarity),
    weighted_strength = sum(similarity, na.rm = TRUE),
    .groups = "drop"
  )

node_summary <- tibble(cell_line = all_cell_lines) %>%
  left_join(node_summary_edge, by = "cell_line") %>%
  mutate(
    degree = coalesce(degree, 0L),
    weighted_strength = coalesce(weighted_strength, 0),
    is_outlier = degree == 0L
  ) %>%
  left_join(cell_label_map, by = "cell_line") %>%
  arrange(desc(degree), desc(weighted_strength), desc(mean_edge_sim))

n_nodes <- length(all_cell_lines)
n_candidate_pairs <- nrow(sim_long)
n_defined_pairs <- nrow(defined_sims)
n_na_pairs <- n_candidate_pairs - n_defined_pairs
n_pairs_gt_similarity_threshold <- if (is.finite(edge_threshold)) sum(defined_sims$similarity > edge_threshold, na.rm = TRUE) else 0L
n_pairs_eq_similarity_threshold <- if (is.finite(edge_threshold)) sum(defined_sims$similarity == edge_threshold, na.rm = TRUE) else 0L
n_edges <- nrow(graph_edges)
density <- if (n_nodes > 1) (2 * n_edges) / (n_nodes * (n_nodes - 1)) else 0

provenance_tbl <- tibble(
  cohort = cfg$profile %||% opt$profile %||% "",
  direction = direction,
  similarity_metric = opt$similarity_metric,
  similarity_definition = similarity_definition,
  similarity_question = similarity_question,
  vector_semantics = vector_semantics,
  binary_tumour_membership = opt$similarity_metric == "jaccard",
  pairwise_union_definition = "pairwise_active_tumour_union: tumours selected by at least one biological cell line after thresholding",
  joint_zero_contribution = joint_zero_contribution,
  fraction_magnitude_used = fraction_magnitude_used,
  zero_fill_rule = "absent tumour entries are filled with 0 before biological-cell-line collapse",
  replicate_collapse = "arithmetic mean across profile-level vectors within biological cell line",
  biological_cell_line_pooling = "arithmetic_mean",
  p_consensus_threshold = threshold,
  similarity_quantile = similarity_quantile,
  n_nodes = n_nodes,
  n_candidate_pairs = n_candidate_pairs,
  n_defined_pairs = n_defined_pairs,
  n_na_pairs = n_na_pairs,
  min_pairwise_active_tumours = if (n_candidate_pairs > 0) min(sim_long$n_pairwise_active_tumours, na.rm = TRUE) else NA_real_,
  pairwise_active_tumour_count_q1 = if (n_candidate_pairs > 0) as.numeric(stats::quantile(sim_long$n_pairwise_active_tumours, 0.25, na.rm = TRUE, names = FALSE)) else NA_real_,
  median_pairwise_active_tumours = if (n_candidate_pairs > 0) stats::median(sim_long$n_pairwise_active_tumours, na.rm = TRUE) else NA_real_,
  pairwise_active_tumour_count_q3 = if (n_candidate_pairs > 0) as.numeric(stats::quantile(sim_long$n_pairwise_active_tumours, 0.75, na.rm = TRUE, names = FALSE)) else NA_real_,
  max_pairwise_active_tumours = if (n_candidate_pairs > 0) max(sim_long$n_pairwise_active_tumours, na.rm = TRUE) else NA_real_,
  n_pairs_lt_3_active_tumours = if (n_candidate_pairs > 0) sum(sim_long$n_pairwise_active_tumours < 3, na.rm = TRUE) else 0L,
  n_pairs_lt_5_active_tumours = if (n_candidate_pairs > 0) sum(sim_long$n_pairwise_active_tumours < 5, na.rm = TRUE) else 0L,
  n_pairs_lt_10_active_tumours = if (n_candidate_pairs > 0) sum(sim_long$n_pairwise_active_tumours < 10, na.rm = TRUE) else 0L,
  undefined_similarity_reasons = if (nrow(undefined_pairs) > 0) {
    paste(
      undefined_pairs %>%
        count(undefined_similarity_reason, name = "n") %>%
        transmute(label = paste0(undefined_similarity_reason, "=", n)) %>%
        pull(label),
      collapse = ";"
    )
  } else {
    ""
  },
  computed_similarity_threshold = edge_threshold,
  n_pairs_gt_similarity_threshold = n_pairs_gt_similarity_threshold,
  n_pairs_eq_similarity_threshold = n_pairs_eq_similarity_threshold,
  threshold_tie_count = n_pairs_eq_similarity_threshold,
  selected_edge_count = n_edges,
  edge_fraction = if (n_candidate_pairs > 0) n_edges / n_candidate_pairs else NA_real_,
  selected_edge_fraction = if (n_candidate_pairs > 0) n_edges / n_candidate_pairs else NA_real_,
  graph_density = density
)

write_tsv(as.data.frame(sim_mat), file.path(out_base, sprintf("cell_line_similarity_matrix_%s.tsv", direction)))
saveRDS(sim_mat, file.path(out_base, sprintf("cell_line_similarity_matrix_%s.rds", direction)))
write_tsv(sim_long, file.path(out_base, sprintf("cell_line_similarity_pairs_%s.tsv", direction)))
write_tsv(undefined_pairs, file.path(out_base, sprintf("cell_line_similarity_undefined_pairs_%s.tsv", direction)))
write_tsv(graph_edges, file.path(out_base, sprintf("cell_line_similarity_graph_edges_%s.tsv", direction)))
write_tsv(node_summary, file.path(out_base, sprintf("cell_line_similarity_graph_node_summary_%s.tsv", direction)))
write_tsv(
  node_summary %>%
    transmute(
      cell_line,
      degree,
      weighted_strength,
      betweenness = NA_real_,
      component = NA_character_,
      community_louv = NA_character_,
      community_leid = NA_character_,
      mean_edge_sim,
      max_edge_sim,
      is_outlier
    ),
  file.path(out_base, sprintf("cell_line_similarity_graph_node_annotations_%s.tsv", direction))
)
write_tsv(selected_tumour_stats, file.path(out_base, sprintf("cell_line_similarity_selected_tumours_%s.tsv", direction)))
write_tsv(provenance_tbl, file.path(out_base, sprintf("cell_line_similarity_graph_provenance_%s.tsv", direction)))

graph_obj <- igraph::graph_from_data_frame(
  graph_edges %>% select(cell_line1, cell_line2, similarity),
  directed = FALSE,
  vertices = tibble(name = all_cell_lines)
)
igraph::write_graph(graph_obj, file.path(out_base, sprintf("cell_line_similarity_graph_%s.graphml", direction)), format = "graphml")

cat("[OK] Wrote metric-specific graph products for ", direction, " (", opt$similarity_metric, ").\n", sep = "")
