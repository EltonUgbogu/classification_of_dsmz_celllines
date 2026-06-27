command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- command_arguments[grepl("^--file=", command_arguments)]
current_script_path <- if (length(file_argument) > 0) sub("^--file=", "", file_argument[[1]]) else file.path("research_framework", "optimisation", "scripts", "01_score_representation_redundancy.R")
source(file.path(dirname(normalizePath(current_script_path, mustWork = TRUE)), "optimisation_utils.R"))
load_required_packages(c("data.table"))

adjusted_rand_index <- function(labels_a, labels_b) {
  complete_index <- !is.na(labels_a) & !is.na(labels_b)
  labels_a <- labels_a[complete_index]
  labels_b <- labels_b[complete_index]
  if (length(labels_a) < 2) {
    return(NA_real_)
  }
  contingency <- table(labels_a, labels_b)
  choose2 <- function(x) x * (x - 1) / 2
  sum_comb <- sum(choose2(contingency))
  row_comb <- sum(choose2(rowSums(contingency)))
  column_comb <- sum(choose2(colSums(contingency)))
  total_comb <- choose2(sum(contingency))
  expected_index <- row_comb * column_comb / total_comb
  max_index <- (row_comb + column_comb) / 2
  if (isTRUE(all.equal(max_index, expected_index))) {
    return(NA_real_)
  }
  (sum_comb - expected_index) / (max_index - expected_index)
}

read_consensus_vector <- function(path) {
  representation_metadata <- parse_representation_from_path(path)
  table <- read_table_if_exists(path)
  if (is.null(table)) {
    return(NULL)
  }
  cell_column <- if ("cell_line" %in% names(table)) "cell_line" else if ("cell_tech_id" %in% names(table)) "cell_tech_id" else NA_character_
  tumour_column <- if ("tumour_id" %in% names(table)) "tumour_id" else if ("tumour" %in% names(table)) "tumour" else NA_character_
  if (is.na(cell_column) || is.na(tumour_column) || !"p_consensus" %in% names(table)) {
    warning("Skipping final consensus table with unsupported columns: ", path)
    return(NULL)
  }
  data.frame(
    cohort_id = representation_metadata$cohort_id,
    representation_id = representation_metadata$representation_id,
    feature_set_id = representation_metadata$feature_set_id,
    distance_metric = representation_metadata$distance_metric,
    transformation_id = representation_metadata$transformation_id,
    cell_line_id = as.character(table[[cell_column]]),
    tumour_sample_id = as.character(table[[tumour_column]]),
    p_consensus_observed = safe_numeric(table$p_consensus, default = 0),
    source_file = path,
    stringsAsFactors = FALSE
  )
}

read_node_clusters <- function(paths) {
  cluster_tables <- lapply(paths, function(path) {
    representation_metadata <- parse_representation_from_path(path)
    table <- read_table_if_exists(path)
    if (is.null(table) || !"cell_line" %in% names(table)) {
      return(NULL)
    }
    cluster_column <- c("community_leid", "community_louv", "component")[c("community_leid", "community_louv", "component") %in% names(table)][1]
    if (is.na(cluster_column)) {
      return(NULL)
    }
    data.frame(
      cohort_id = representation_metadata$cohort_id,
      representation_id = representation_metadata$representation_id,
      cell_line_id = as.character(table$cell_line),
      cluster_label = as.character(table[[cluster_column]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, cluster_tables[!vapply(cluster_tables, is.null, logical(1))])
}

arguments <- parse_key_value_arguments()
config <- load_optimisation_config(arguments$config)

output_root <- config_value(config, c("outputs", "optimisation_results_root"), file.path("research_framework", "optimisation", "results"))
docs_root <- config_value(config, c("outputs", "optimisation_docs_root"), file.path("research_framework", "optimisation", "docs"))
threshold <- as.numeric(config_value(config, c("representation_diagnostics", "strong_neighbourhood_threshold"), 0.7))

redundancy_output <- file.path(output_root, "representation_diagnostics", "representation_redundancy.tsv")
matrix_output <- file.path(output_root, "representation_diagnostics", "representation_distance_matrix_correlations.tsv")
nearest_output <- file.path(output_root, "representation_diagnostics", "nearest_neighbour_overlap.tsv")
report_output <- file.path(docs_root, "02_representation_diagnostics_report.md")

consensus_paths <- final_consensus_paths(config)
if (length(consensus_paths) == 0) {
  stop("No final consensus neighbourhood tables were found. The representation redundancy step requires thesis final consensus tables.", call. = FALSE)
}

consensus_tables <- lapply(consensus_paths, read_consensus_vector)
consensus_data <- do.call(rbind, consensus_tables[!vapply(consensus_tables, is.null, logical(1))])
if (is.null(consensus_data) || nrow(consensus_data) == 0) {
  stop("Final consensus neighbourhood tables were found, but no usable p_consensus rows could be read.", call. = FALSE)
}

representation_metadata <- unique(consensus_data[, c("cohort_id", "representation_id", "feature_set_id", "distance_metric", "transformation_id", "source_file")])
representation_metadata$cohort_representation_id <- paste(representation_metadata$cohort_id, representation_metadata$representation_id, sep = "::")
consensus_data$cohort_representation_id <- paste(consensus_data$cohort_id, consensus_data$representation_id, sep = "::")
consensus_data$pair_key <- paste(consensus_data$cell_line_id, consensus_data$tumour_sample_id, sep = "||")

representation_ids <- sort(unique(consensus_data$cohort_representation_id))
pairwise_rows <- list()
nearest_rows <- list()

strong_sets <- split(consensus_data[consensus_data$p_consensus_observed >= threshold, ], consensus_data$cohort_representation_id[consensus_data$p_consensus_observed >= threshold])
vector_sets <- split(consensus_data[, c("pair_key", "p_consensus_observed")], consensus_data$cohort_representation_id)

row_counter <- 1
nearest_counter <- 1
for (first_index in seq_along(representation_ids)) {
  if (first_index >= length(representation_ids)) {
    next
  }
  for (second_index in seq.int(first_index + 1, length(representation_ids))) {
    first_id <- representation_ids[[first_index]]
    second_id <- representation_ids[[second_index]]
    first_meta <- representation_metadata[representation_metadata$cohort_representation_id == first_id, ][1, ]
    second_meta <- representation_metadata[representation_metadata$cohort_representation_id == second_id, ][1, ]
    if (!identical(first_meta$cohort_id, second_meta$cohort_id)) {
      next
    }
    first_vector <- vector_sets[[first_id]]
    second_vector <- vector_sets[[second_id]]
    merged_vector <- merge(first_vector, second_vector, by = "pair_key", all = TRUE, suffixes = c("_1", "_2"))
    merged_vector$p_consensus_observed_1[is.na(merged_vector$p_consensus_observed_1)] <- 0
    merged_vector$p_consensus_observed_2[is.na(merged_vector$p_consensus_observed_2)] <- 0
    matrix_correlation <- if (stats::sd(merged_vector$p_consensus_observed_1) > 0 && stats::sd(merged_vector$p_consensus_observed_2) > 0) {
      stats::cor(merged_vector$p_consensus_observed_1, merged_vector$p_consensus_observed_2, method = "spearman")
    } else {
      NA_real_
    }
    first_strong <- if (!is.null(strong_sets[[first_id]])) unique(strong_sets[[first_id]]$pair_key) else character(0)
    second_strong <- if (!is.null(strong_sets[[second_id]])) unique(strong_sets[[second_id]]$pair_key) else character(0)
    union_size <- length(union(first_strong, second_strong))
    nearest_neighbour_jaccard <- if (union_size == 0) NA_real_ else length(intersect(first_strong, second_strong)) / union_size
    feature_overlap_jaccard <- if (identical(first_meta$feature_set_id, second_meta$feature_set_id)) 1 else NA_real_
    pairwise_rows[[row_counter]] <- data.frame(
      cohort_id = first_meta$cohort_id,
      representation_id_1 = first_meta$representation_id,
      representation_id_2 = second_meta$representation_id,
      feature_set_id_1 = first_meta$feature_set_id,
      feature_set_id_2 = second_meta$feature_set_id,
      distance_metric_1 = first_meta$distance_metric,
      distance_metric_2 = second_meta$distance_metric,
      matrix_correlation = matrix_correlation,
      feature_overlap_jaccard = feature_overlap_jaccard,
      n_pairs_compared = nrow(merged_vector),
      stringsAsFactors = FALSE
    )
    nearest_rows[[nearest_counter]] <- data.frame(
      cohort_id = first_meta$cohort_id,
      representation_id_1 = first_meta$representation_id,
      representation_id_2 = second_meta$representation_id,
      nearest_neighbour_jaccard = nearest_neighbour_jaccard,
      strong_threshold = threshold,
      n_union_pairs = union_size,
      stringsAsFactors = FALSE
    )
    row_counter <- row_counter + 1
    nearest_counter <- nearest_counter + 1
  }
}

matrix_table <- if (length(pairwise_rows) == 0) data.frame() else do.call(rbind, pairwise_rows)
nearest_table <- if (length(nearest_rows) == 0) data.frame() else do.call(rbind, nearest_rows)

node_data <- read_node_clusters(node_annotation_paths(config))
cluster_rows <- list()
if (!is.null(node_data) && nrow(node_data) > 0 && nrow(matrix_table) > 0) {
  for (row_index in seq_len(nrow(matrix_table))) {
    row <- matrix_table[row_index, ]
    first_clusters <- node_data[node_data$cohort_id == row$cohort_id & node_data$representation_id == row$representation_id_1, c("cell_line_id", "cluster_label")]
    second_clusters <- node_data[node_data$cohort_id == row$cohort_id & node_data$representation_id == row$representation_id_2, c("cell_line_id", "cluster_label")]
    merged_clusters <- merge(first_clusters, second_clusters, by = "cell_line_id", suffixes = c("_1", "_2"))
    cluster_rows[[row_index]] <- adjusted_rand_index(merged_clusters$cluster_label_1, merged_clusters$cluster_label_2)
  }
  matrix_table$cluster_ari <- unlist(cluster_rows)
} else if (nrow(matrix_table) > 0) {
  matrix_table$cluster_ari <- NA_real_
}

summary_rows <- lapply(seq_len(nrow(representation_metadata)), function(row_index) {
  metadata_row <- representation_metadata[row_index, ]
  pair_rows <- matrix_table[matrix_table$cohort_id == metadata_row$cohort_id & (matrix_table$representation_id_1 == metadata_row$representation_id | matrix_table$representation_id_2 == metadata_row$representation_id), , drop = FALSE]
  nearest_rows_for_representation <- nearest_table[nearest_table$cohort_id == metadata_row$cohort_id & (nearest_table$representation_id_1 == metadata_row$representation_id | nearest_table$representation_id_2 == metadata_row$representation_id), , drop = FALSE]
  mean_abs_correlation <- if (nrow(pair_rows) == 0) NA_real_ else mean(abs(pair_rows$matrix_correlation), na.rm = TRUE)
  mean_nearest_jaccard <- if (nrow(nearest_rows_for_representation) == 0) NA_real_ else mean(nearest_rows_for_representation$nearest_neighbour_jaccard, na.rm = TRUE)
  mean_cluster_ari <- if (nrow(pair_rows) == 0 || !"cluster_ari" %in% names(pair_rows)) NA_real_ else mean(pair_rows$cluster_ari, na.rm = TRUE)
  redundancy_score <- max(c(mean_abs_correlation, mean_nearest_jaccard, mean_cluster_ari), na.rm = TRUE)
  if (!is.finite(redundancy_score)) {
    redundancy_score <- NA_real_
  }
  redundancy_class <- if (is.na(redundancy_score)) {
    "not_estimated"
  } else if (redundancy_score >= 0.85) {
    "highly_redundant"
  } else if (redundancy_score >= 0.60) {
    "partly_redundant"
  } else {
    "complementary"
  }
  data.frame(
    cohort_id = metadata_row$cohort_id,
    representation_id = metadata_row$representation_id,
    feature_set_id = metadata_row$feature_set_id,
    distance_metric = metadata_row$distance_metric,
    transformation_id = metadata_row$transformation_id,
    matrix_correlation = mean_abs_correlation,
    nearest_neighbour_jaccard = mean_nearest_jaccard,
    cluster_ari = mean_cluster_ari,
    feature_overlap_jaccard = NA_real_,
    source_association_score = NA_real_,
    cancer_type_preservation_score = ifelse(metadata_row$cohort_id == "multicohort_cancer", NA_real_, 1),
    redundancy_class = redundancy_class,
    stringsAsFactors = FALSE
  )
})

redundancy_table <- do.call(rbind, summary_rows)
write_tsv_table(redundancy_table, redundancy_output)
write_tsv_table(matrix_table, matrix_output)
write_tsv_table(nearest_table, nearest_output)

class_counts <- as.data.frame(table(redundancy_table$redundancy_class), stringsAsFactors = FALSE)
names(class_counts) <- c("redundancy_class", "n_representations")

report_lines <- c(
  "# Representation diagnostics report",
  "",
  "## Redundancy diagnostics",
  "",
  "The redundancy step compared representation-specific p-consensus value matrices within each cohort. Pairwise matrix agreement was estimated by Spearman correlation after missing pair values were set to zero. Nearest-neighbour overlap was estimated by Jaccard overlap of cell-line--tumour pairs with p-consensus values at or above the configured threshold.",
  "",
  "## Assumptions",
  "",
  paste0("- Strong neighbourhood threshold: ", threshold, "."),
  "- Missing cell-line--tumour pairs were treated as observed p-consensus value 0 for matrix comparisons.",
  "- Feature-set overlap was only recorded where direct feature membership tables were available; otherwise it remains unestimated.",
  "",
  "## Redundancy classes",
  "",
  paste0("- ", class_counts$redundancy_class, ": ", class_counts$n_representations, " representation(s)."),
  "",
  "## Outputs",
  "",
  paste0("- `", redundancy_output, "`"),
  paste0("- `", matrix_output, "`"),
  paste0("- `", nearest_output, "`"),
  "",
  "## Manual inspection required",
  "",
  "- Inspect highly redundant representations before treating them as independent evidence sources.",
  "- Inspect complementary representations with low matrix correlation and low nearest-neighbour overlap because they may encode distinct tumour-neighbourhood structure or noisy evidence."
)
write_markdown(report_lines, report_output)
cat("Representation redundancy diagnostics completed: ", redundancy_output, "\n", sep = "")
