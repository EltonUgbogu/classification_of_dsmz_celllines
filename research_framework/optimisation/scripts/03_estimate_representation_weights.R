command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- command_arguments[grepl("^--file=", command_arguments)]
current_script_path <- if (length(file_argument) > 0) sub("^--file=", "", file_argument[[1]]) else file.path("research_framework", "optimisation", "scripts", "03_estimate_representation_weights.R")
source(file.path(dirname(normalizePath(current_script_path, mustWork = TRUE)), "optimisation_utils.R"))
load_required_packages(c("data.table"))

arguments <- parse_key_value_arguments()
config <- load_optimisation_config(arguments$config)

output_root <- config_value(config, c("outputs", "optimisation_results_root"), file.path("research_framework", "optimisation", "results"))
docs_root <- config_value(config, c("outputs", "optimisation_docs_root"), file.path("research_framework", "optimisation", "docs"))
minimum_weight <- as.numeric(config_value(config, c("representation_weights", "minimum_weight"), 0.01))
normalise_weights <- isTRUE(config_value(config, c("representation_weights", "normalise_weights"), TRUE))

redundancy_path <- file.path(output_root, "representation_diagnostics", "representation_redundancy.tsv")
stability_path <- file.path(output_root, "representation_diagnostics", "representation_stability.tsv")
weights_output <- file.path(output_root, "representation_weights", "representation_weights.tsv")
report_output <- file.path(docs_root, "02_representation_diagnostics_report.md")

redundancy_table <- read_table_if_exists(redundancy_path)
stability_table <- read_table_if_exists(stability_path)
if (is.null(redundancy_table) || is.null(stability_table)) {
  stop("Representation weights require both redundancy and stability tables.", call. = FALSE)
}

merged_table <- merge(
  redundancy_table,
  stability_table,
  by = c("cohort_id", "representation_id", "feature_set_id", "distance_metric", "transformation_id"),
  all = TRUE,
  suffixes = c("_redundancy", "_stability")
)

redundancy_signal <- pmax(
  safe_numeric(merged_table$matrix_correlation, default = NA_real_),
  safe_numeric(merged_table$nearest_neighbour_jaccard, default = NA_real_),
  safe_numeric(merged_table$cluster_ari, default = NA_real_),
  na.rm = TRUE
)
redundancy_signal[!is.finite(redundancy_signal)] <- 0.5
non_redundancy_score <- pmin(pmax(1 - redundancy_signal, 0.05), 1)

stability_score <- safe_numeric(merged_table$stability_score, default = 0.5)
batch_independence_score <- 1 - safe_numeric(merged_table$source_association_score, default = 0)
batch_independence_score <- pmin(pmax(batch_independence_score, 0.05), 1)
cancer_type_preservation_score <- safe_numeric(merged_table$cancer_type_preservation_score, default = 0.5)
cancer_type_preservation_score <- pmin(pmax(cancer_type_preservation_score, 0.05), 1)

raw_weight <- stability_score * non_redundancy_score * batch_independence_score * cancer_type_preservation_score
raw_weight[is.na(raw_weight) | !is.finite(raw_weight)] <- minimum_weight
raw_weight <- pmax(raw_weight, minimum_weight)
if (normalise_weights) {
  weight_normalised <- ave(raw_weight, merged_table$cohort_id, FUN = function(x) x / sum(x))
} else {
  weight_normalised <- raw_weight
}

weight_class <- ifelse(weight_normalised >= stats::quantile(weight_normalised, 0.75, na.rm = TRUE), "high_weight",
  ifelse(weight_normalised <= stats::quantile(weight_normalised, 0.25, na.rm = TRUE), "low_weight", "medium_weight"))

weights_table <- data.frame(
  cohort_id = merged_table$cohort_id,
  representation_id = merged_table$representation_id,
  feature_set_id = merged_table$feature_set_id,
  distance_metric = merged_table$distance_metric,
  transformation_id = merged_table$transformation_id,
  stability_score = stability_score,
  non_redundancy_score = non_redundancy_score,
  batch_independence_score = batch_independence_score,
  cancer_type_preservation_score = cancer_type_preservation_score,
  representation_reliability_weight = raw_weight,
  weight_normalised = weight_normalised,
  weight_class = weight_class,
  notes = ifelse(is.na(merged_table$source_association_score), "batch/source score not estimated; neutral independence score used", ""),
  stringsAsFactors = FALSE
)

write_tsv_table(weights_table, weights_output)

top_weights <- weights_table[order(-weights_table$weight_normalised), ]
top_weights <- utils::head(top_weights, 10)
low_weights <- weights_table[order(weights_table$weight_normalised), ]
low_weights <- utils::head(low_weights, 10)

report_lines <- c(
  "",
  "## Representation reliability weights",
  "",
  "Representation reliability weights were estimated as the product of stability, non-redundancy, batch-independence, and cancer-type preservation scores. Weights were normalised within cohort when configured.",
  "",
  "## Highest-weighted representations",
  "",
  paste0("- `", top_weights$cohort_id, "::", top_weights$representation_id, "`: weight ", signif(top_weights$weight_normalised, 4), "."),
  "",
  "## Lowest-weighted representations",
  "",
  paste0("- `", low_weights$cohort_id, "::", low_weights$representation_id, "`: weight ", signif(low_weights$weight_normalised, 4), "."),
  "",
  "## Output",
  "",
  paste0("- `", weights_output, "`"),
  "",
  "## Manual inspection required",
  "",
  "- Inspect low-weight representations before excluding them from scientific interpretation.",
  "- Inspect representations with neutral batch-independence scores because metadata associations were not estimated from available files."
)
append_markdown(report_lines, report_output)
cat("Representation weights completed: ", weights_output, "\n", sep = "")

