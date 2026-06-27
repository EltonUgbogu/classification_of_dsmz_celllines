command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- command_arguments[grepl("^--file=", command_arguments)]
current_script_path <- if (length(file_argument) > 0) sub("^--file=", "", file_argument[[1]]) else file.path("research_framework", "optimisation", "scripts", "05_infer_probabilistic_graph_edges.R")
source(file.path(dirname(normalizePath(current_script_path, mustWork = TRUE)), "optimisation_utils.R"))
load_required_packages(c("data.table"))

arguments <- parse_key_value_arguments()
config <- load_optimisation_config(arguments$config)

output_root <- config_value(config, c("outputs", "optimisation_results_root"), file.path("research_framework", "optimisation", "results"))
docs_root <- config_value(config, c("outputs", "optimisation_docs_root"), file.path("research_framework", "optimisation", "docs"))
null_iterations <- as.integer(config_value(config, c("probabilistic_graph", "null_iterations"), 1000))
null_model <- config_value(config, c("probabilistic_graph", "null_model"), "within_cancer_tumour_permutation")
random_seed <- as.integer(config_value(config, c("probabilistic_graph", "random_seed"), 20260627))

posterior_path <- file.path(output_root, "probabilistic_neighbourhoods", "posterior_neighbourhood_probabilities.tsv.gz")
edges_output <- file.path(output_root, "probabilistic_graphs", "probabilistic_cellline_edges.tsv")
null_output <- file.path(output_root, "probabilistic_graphs", "probabilistic_edge_null_summary.tsv")
report_output <- file.path(docs_root, "05_probabilistic_graph_report.md")

if (!file.exists(posterior_path)) {
  stop("Posterior neighbourhood probabilities are required before probabilistic graph edge inference.", call. = FALSE)
}
posterior_table <- data.table::fread(posterior_path, sep = "\t", data.table = FALSE, showProgress = FALSE)
if (nrow(posterior_table) == 0) {
  stop("Posterior neighbourhood table is empty.", call. = FALSE)
}

set.seed(random_seed)
cancer_types <- sort(unique(posterior_table$cancer_type))
edge_rows <- list()
null_rows <- list()
edge_counter <- 1
null_counter <- 1

for (cancer_type in cancer_types) {
  cohort_table <- posterior_table[posterior_table$cancer_type == cancer_type, ]
  wide_table <- reshape(
    cohort_table[, c("cell_line_id", "tumour_sample_id", "posterior_mean_neighbourhood_probability")],
    idvar = "cell_line_id",
    timevar = "tumour_sample_id",
    direction = "wide"
  )
  if (nrow(wide_table) < 2) {
    next
  }
  cell_lines <- wide_table$cell_line_id
  probability_matrix <- as.matrix(wide_table[, setdiff(names(wide_table), "cell_line_id"), drop = FALSE])
  probability_matrix[is.na(probability_matrix)] <- 0
  expected_shared_matrix <- probability_matrix %*% t(probability_matrix)
  pair_indices <- which(upper.tri(expected_shared_matrix), arr.ind = TRUE)
  if (nrow(pair_indices) == 0) {
    next
  }
  observed_shared <- expected_shared_matrix[pair_indices]
  null_matrix <- matrix(NA_real_, nrow = nrow(pair_indices), ncol = null_iterations)
  for (iteration in seq_len(null_iterations)) {
    permuted_matrix <- t(apply(probability_matrix, 1, sample, replace = FALSE))
    null_shared_matrix <- permuted_matrix %*% t(permuted_matrix)
    null_matrix[, iteration] <- null_shared_matrix[pair_indices]
  }
  null_mean <- rowMeans(null_matrix, na.rm = TRUE)
  null_sd <- apply(null_matrix, 1, stats::sd, na.rm = TRUE)
  null_p_value <- vapply(seq_along(observed_shared), function(index) {
    (sum(null_matrix[index, ] >= observed_shared[[index]], na.rm = TRUE) + 1) / (null_iterations + 1)
  }, numeric(1))
  posterior_edge_probability <- pmin(pmax(1 - null_p_value, 0), 1)
  edge_q_value <- bh_q_value(null_p_value)
  edge_table <- data.frame(
    cell_line_1 = cell_lines[pair_indices[, 1]],
    cell_line_2 = cell_lines[pair_indices[, 2]],
    cancer_type_1 = cancer_type,
    cancer_type_2 = cancer_type,
    expected_shared_neighbourhood_support = observed_shared,
    null_mean_shared_support = null_mean,
    null_sd_shared_support = null_sd,
    null_p_value = null_p_value,
    edge_q_value = edge_q_value,
    posterior_edge_probability = posterior_edge_probability,
    edge_probability_class = ifelse(posterior_edge_probability >= 0.95, "high_probability",
      ifelse(posterior_edge_probability >= 0.80, "moderate_probability",
        ifelse(posterior_edge_probability >= 0.50, "uncertain_probability", "low_probability"))),
    edge_source = "posterior_neighbourhood_overlap",
    notes = paste0("Null model: ", null_model),
    stringsAsFactors = FALSE
  )
  edge_rows[[edge_counter]] <- edge_table
  null_rows[[null_counter]] <- data.frame(
    cancer_type = cancer_type,
    null_model = null_model,
    null_iterations = null_iterations,
    n_cell_lines = length(cell_lines),
    n_tumours = ncol(probability_matrix),
    mean_expected_shared_neighbourhood_support = mean(observed_shared, na.rm = TRUE),
    mean_null_shared_support = mean(null_mean, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  edge_counter <- edge_counter + 1
  null_counter <- null_counter + 1
}

edge_table_all <- if (length(edge_rows) == 0) data.frame() else do.call(rbind, edge_rows)
null_summary_table <- if (length(null_rows) == 0) data.frame() else do.call(rbind, null_rows)
if (nrow(edge_table_all) == 0) {
  stop("No probabilistic graph edges could be estimated. At least two cell lines per cancer type are required.", call. = FALSE)
}

write_tsv_table(edge_table_all, edges_output)
write_tsv_table(null_summary_table, null_output)

class_counts <- as.data.frame(table(edge_table_all$edge_probability_class), stringsAsFactors = FALSE)
names(class_counts) <- c("edge_probability_class", "n_edges")

report_lines <- c(
  "# Probabilistic graph report",
  "",
  "## Scope",
  "",
  "The probabilistic graph step estimated cell-line graph edges from expected shared probabilistic tumour-neighbourhood evidence.",
  "",
  "## Model",
  "",
  "- The initial implementation focuses on latent neighbourhood events `N_ct` and graph-edge events `E_cd`.",
  "- Expected shared neighbourhood evidence was computed as the sum over tumours of products of posterior neighbourhood probabilities for two cell lines.",
  paste0("- Null model: ", null_model, "."),
  paste0("- Null iterations: ", null_iterations, "."),
  "- Posterior edge probability was estimated as one minus the empirical upper-tail null probability.",
  "",
  "## Edge probability classes",
  "",
  paste0("- ", class_counts$edge_probability_class, ": ", class_counts$n_edges, " edge(s)."),
  "",
  "## Outputs",
  "",
  paste0("- `", edges_output, "`"),
  paste0("- `", null_output, "`"),
  "",
  "## Manual inspection required",
  "",
  "- Inspect high-probability edges that are absent from deterministic thesis graphs.",
  "- Inspect deterministic edges with low posterior edge probabilities.",
  "- Inspect cancer types where the null mean is close to observed shared neighbourhood evidence because threshold selection may be unstable."
)
write_markdown(report_lines, report_output)
cat("Probabilistic graph edge inference completed: ", edges_output, "\n", sep = "")

