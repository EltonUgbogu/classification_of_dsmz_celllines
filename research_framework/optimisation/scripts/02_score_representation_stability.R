command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- command_arguments[grepl("^--file=", command_arguments)]
current_script_path <- if (length(file_argument) > 0) sub("^--file=", "", file_argument[[1]]) else file.path("research_framework", "optimisation", "scripts", "02_score_representation_stability.R")
source(file.path(dirname(normalizePath(current_script_path, mustWork = TRUE)), "optimisation_utils.R"))
load_required_packages(c("data.table"))

arguments <- parse_key_value_arguments()
config <- load_optimisation_config(arguments$config)

output_root <- config_value(config, c("outputs", "optimisation_results_root"), file.path("research_framework", "optimisation", "results"))
docs_root <- config_value(config, c("outputs", "optimisation_docs_root"), file.path("research_framework", "optimisation", "docs"))
threshold <- as.numeric(config_value(config, c("representation_diagnostics", "strong_neighbourhood_threshold"), 0.7))
sensitivity_window <- as.numeric(config_value(config, c("representation_stability", "threshold_sensitivity_window"), 0.05))
bootstrap_iterations <- as.integer(config_value(config, c("representation_stability", "bootstrap_iterations"), 100))
enable_lightweight_bootstrap <- isTRUE(config_value(config, c("representation_stability", "enable_lightweight_bootstrap"), FALSE))

stability_output <- file.path(output_root, "representation_diagnostics", "representation_stability.tsv")
report_output <- file.path(docs_root, "02_representation_diagnostics_report.md")

consensus_paths <- final_consensus_paths(config)
if (length(consensus_paths) == 0) {
  stop("No final consensus neighbourhood tables were found. The stability step requires thesis final consensus tables.", call. = FALSE)
}

stability_rows <- lapply(consensus_paths, function(path) {
  representation_metadata <- parse_representation_from_path(path)
  table <- read_table_if_exists(path)
  if (is.null(table) || !"p_consensus" %in% names(table)) {
    warning("Skipping stability scoring for unsupported final consensus table: ", path)
    return(NULL)
  }
  p_consensus_values <- safe_numeric(table$p_consensus, default = 0)
  support_indicator <- p_consensus_values >= threshold
  mean_neighbourhood_recurrence <- mean(p_consensus_values, na.rm = TRUE)
  median_neighbourhood_recurrence <- stats::median(p_consensus_values, na.rm = TRUE)
  threshold_sensitivity <- mean(abs(p_consensus_values - threshold) <= sensitivity_window, na.rm = TRUE)
  bootstrap_recurrence <- NA_real_
  if (enable_lightweight_bootstrap && length(p_consensus_values) > 1) {
    set.seed(20260627)
    bootstrap_means <- replicate(bootstrap_iterations, mean(sample(p_consensus_values, replace = TRUE), na.rm = TRUE))
    bootstrap_recurrence <- mean(bootstrap_means, na.rm = TRUE)
  }
  graph_edge_path <- sub("Final_consensus_tumour_neighbourhoods_.*\\.tsv$", paste0("cell_line_similarity_graph_edges_", representation_metadata$representation_id, ".tsv"), path)
  edge_stability_score <- NA_real_
  if (file.exists(graph_edge_path)) {
    edge_table <- read_table_if_exists(graph_edge_path)
    similarity_column <- c("similarity", "weight", "edge_weight")[c("similarity", "weight", "edge_weight") %in% names(edge_table)][1]
    if (!is.na(similarity_column)) {
      edge_stability_score <- mean(safe_numeric(edge_table[[similarity_column]], default = NA_real_), na.rm = TRUE)
    }
  }
  cluster_stability_score <- ifelse(is.na(edge_stability_score), mean_neighbourhood_recurrence, mean(c(edge_stability_score, mean_neighbourhood_recurrence), na.rm = TRUE))
  ranking_stability_score <- 1 - threshold_sensitivity
  stability_score <- mean(c(mean_neighbourhood_recurrence, cluster_stability_score, edge_stability_score, ranking_stability_score), na.rm = TRUE)
  stability_class <- if (is.na(stability_score)) {
    "not_estimated"
  } else if (stability_score >= 0.75) {
    "stable"
  } else if (stability_score >= 0.50) {
    "moderately_stable"
  } else {
    "unstable"
  }
  data.frame(
    cohort_id = representation_metadata$cohort_id,
    representation_id = representation_metadata$representation_id,
    feature_set_id = representation_metadata$feature_set_id,
    distance_metric = representation_metadata$distance_metric,
    transformation_id = representation_metadata$transformation_id,
    bootstrap_iterations = ifelse(enable_lightweight_bootstrap, bootstrap_iterations, 0L),
    mean_neighbourhood_recurrence = mean_neighbourhood_recurrence,
    median_neighbourhood_recurrence = median_neighbourhood_recurrence,
    bootstrap_neighbourhood_recurrence = bootstrap_recurrence,
    cluster_stability_score = cluster_stability_score,
    edge_stability_score = edge_stability_score,
    ranking_stability_score = ranking_stability_score,
    threshold_sensitivity_score = 1 - threshold_sensitivity,
    stability_score = stability_score,
    stability_class = stability_class,
    source_file = path,
    stringsAsFactors = FALSE
  )
})

stability_table <- do.call(rbind, stability_rows[!vapply(stability_rows, is.null, logical(1))])
if (is.null(stability_table) || nrow(stability_table) == 0) {
  stop("No usable representation stability rows could be computed.", call. = FALSE)
}
write_tsv_table(stability_table, stability_output)

class_counts <- as.data.frame(table(stability_table$stability_class), stringsAsFactors = FALSE)
names(class_counts) <- c("stability_class", "n_representations")

report_lines <- c(
  "",
  "## Stability diagnostics",
  "",
  "The stability step estimated representation-level stability from observed p-consensus values, edge similarity when graph edges were available, and sensitivity to the configured strong-neighbourhood threshold.",
  "",
  "## Assumptions",
  "",
  paste0("- Strong neighbourhood threshold: ", threshold, "."),
  paste0("- Threshold sensitivity window: ", sensitivity_window, "."),
  paste0("- Lightweight bootstrap enabled: ", enable_lightweight_bootstrap, "."),
  "- If no bootstrap artefacts were available and lightweight bootstrap was disabled, stability was estimated from existing thesis outputs only.",
  "",
  "## Stability classes",
  "",
  paste0("- ", class_counts$stability_class, ": ", class_counts$n_representations, " representation(s)."),
  "",
  "## Output",
  "",
  paste0("- `", stability_output, "`"),
  "",
  "## Manual inspection required",
  "",
  "- Inspect representations with low ranking stability scores because they have many values close to the threshold.",
  "- Inspect representations with missing edge stability scores because graph-edge evidence was unavailable or unreadable."
)
append_markdown(report_lines, report_output)
cat("Representation stability diagnostics completed: ", stability_output, "\n", sep = "")

