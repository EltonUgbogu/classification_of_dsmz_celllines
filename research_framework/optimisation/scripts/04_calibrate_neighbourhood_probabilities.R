command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- command_arguments[grepl("^--file=", command_arguments)]
current_script_path <- if (length(file_argument) > 0) sub("^--file=", "", file_argument[[1]]) else file.path("research_framework", "optimisation", "scripts", "04_calibrate_neighbourhood_probabilities.R")
source(file.path(dirname(normalizePath(current_script_path, mustWork = TRUE)), "optimisation_utils.R"))
load_required_packages(c("data.table"))

read_support_observations <- function(path) {
  representation_metadata <- parse_representation_from_path(path)
  table <- read_table_if_exists(path)
  if (is.null(table) || !"p_consensus" %in% names(table)) {
    return(NULL)
  }
  cell_column <- if ("cell_line" %in% names(table)) "cell_line" else if ("cell_tech_id" %in% names(table)) "cell_tech_id" else NA_character_
  tumour_column <- if ("tumour_id" %in% names(table)) "tumour_id" else if ("tumour" %in% names(table)) "tumour" else NA_character_
  if (is.na(cell_column) || is.na(tumour_column)) {
    return(NULL)
  }
  data.frame(
    cohort_id = representation_metadata$cohort_id,
    representation_id = representation_metadata$representation_id,
    cell_line_id = as.character(table[[cell_column]]),
    tumour_sample_id = as.character(table[[tumour_column]]),
    p_consensus_observed = safe_numeric(table$p_consensus, default = 0),
    stringsAsFactors = FALSE
  )
}

arguments <- parse_key_value_arguments()
config <- load_optimisation_config(arguments$config)

output_root <- config_value(config, c("outputs", "optimisation_results_root"), file.path("research_framework", "optimisation", "results"))
docs_root <- config_value(config, c("outputs", "optimisation_docs_root"), file.path("research_framework", "optimisation", "docs"))
prior_alpha <- as.numeric(config_value(config, c("probabilistic_neighbourhoods", "prior_alpha"), 1))
prior_beta <- as.numeric(config_value(config, c("probabilistic_neighbourhoods", "prior_beta"), 1))
credible_interval <- as.numeric(config_value(config, c("probabilistic_neighbourhoods", "credible_interval"), 0.95))
evidence_threshold <- as.numeric(config_value(config, c("probabilistic_neighbourhoods", "evidence_support_threshold"), 0.7))

weights_path <- file.path(output_root, "representation_weights", "representation_weights.tsv")
posterior_output <- file.path(output_root, "probabilistic_neighbourhoods", "posterior_neighbourhood_probabilities.tsv.gz")
summary_output <- file.path(output_root, "probabilistic_neighbourhoods", "neighbourhood_probability_summary.tsv")
report_output <- file.path(docs_root, "04_probabilistic_neighbourhood_report.md")

weights_table <- read_table_if_exists(weights_path)
if (is.null(weights_table)) {
  stop("Representation weights are required before calibrating posterior neighbourhood probabilities.", call. = FALSE)
}

consensus_paths <- final_consensus_paths(config)
support_tables <- lapply(consensus_paths, read_support_observations)
support_data <- do.call(rbind, support_tables[!vapply(support_tables, is.null, logical(1))])
if (is.null(support_data) || nrow(support_data) == 0) {
  stop("No usable final consensus observations were available for posterior neighbourhood modelling.", call. = FALSE)
}

support_data <- merge(
  support_data,
  weights_table[, c("cohort_id", "representation_id", "weight_normalised")],
  by = c("cohort_id", "representation_id"),
  all.x = TRUE
)
support_data$weight_normalised[is.na(support_data$weight_normalised)] <- 0
support_data$support_indicator <- as.integer(support_data$p_consensus_observed >= evidence_threshold)
support_data$weighted_support_contribution <- support_data$weight_normalised * support_data$support_indicator

representation_counts <- aggregate(
  weight_normalised ~ cohort_id,
  data = unique(weights_table[, c("cohort_id", "representation_id", "weight_normalised")]),
  FUN = sum
)
names(representation_counts)[names(representation_counts) == "weight_normalised"] <- "effective_representation_count"
raw_representation_counts <- aggregate(representation_id ~ cohort_id, data = unique(weights_table[, c("cohort_id", "representation_id")]), FUN = length)
names(raw_representation_counts)[names(raw_representation_counts) == "representation_id"] <- "total_representation_count"

weighted_support <- aggregate(
  weighted_support_contribution ~ cohort_id + cell_line_id + tumour_sample_id,
  data = support_data,
  FUN = sum
)
support_count <- aggregate(
  support_indicator ~ cohort_id + cell_line_id + tumour_sample_id,
  data = support_data,
  FUN = sum
)
observed_max <- aggregate(
  p_consensus_observed ~ cohort_id + cell_line_id + tumour_sample_id,
  data = support_data,
  FUN = max
)

posterior_table <- Reduce(function(x, y) merge(x, y, by = c("cohort_id", "cell_line_id", "tumour_sample_id"), all = TRUE), list(weighted_support, support_count, observed_max))
names(posterior_table)[names(posterior_table) == "weighted_support_contribution"] <- "weighted_support"
names(posterior_table)[names(posterior_table) == "support_indicator"] <- "support_count"
posterior_table <- merge(posterior_table, representation_counts, by = "cohort_id", all.x = TRUE)
posterior_table <- merge(posterior_table, raw_representation_counts, by = "cohort_id", all.x = TRUE)
posterior_table$effective_representation_count[is.na(posterior_table$effective_representation_count)] <- 0
posterior_table$total_representation_count[is.na(posterior_table$total_representation_count)] <- 0

posterior_table$prior_alpha <- prior_alpha
posterior_table$prior_beta <- prior_beta
posterior_table$posterior_alpha <- prior_alpha + posterior_table$weighted_support
posterior_table$posterior_beta <- prior_beta + posterior_table$effective_representation_count - posterior_table$weighted_support
posterior_table$posterior_beta <- pmax(posterior_table$posterior_beta, 1e-9)
posterior_table$posterior_mean_neighbourhood_probability <- posterior_table$posterior_alpha / (posterior_table$posterior_alpha + posterior_table$posterior_beta)
posterior_table$posterior_probability_ge_0_5 <- stats::pbeta(0.5, posterior_table$posterior_alpha, posterior_table$posterior_beta, lower.tail = FALSE)
posterior_table$posterior_probability_ge_0_7 <- stats::pbeta(0.7, posterior_table$posterior_alpha, posterior_table$posterior_beta, lower.tail = FALSE)
posterior_table$posterior_probability_ge_0_8 <- stats::pbeta(0.8, posterior_table$posterior_alpha, posterior_table$posterior_beta, lower.tail = FALSE)
interval_tail <- (1 - credible_interval) / 2
posterior_table$posterior_interval_lower <- stats::qbeta(interval_tail, posterior_table$posterior_alpha, posterior_table$posterior_beta)
posterior_table$posterior_interval_upper <- stats::qbeta(1 - interval_tail, posterior_table$posterior_alpha, posterior_table$posterior_beta)
posterior_table$probabilistic_neighbourhood_class <- ifelse(
  posterior_table$posterior_probability_ge_0_8 >= 0.8, "high_posterior_probability",
  ifelse(posterior_table$posterior_probability_ge_0_7 >= 0.5, "moderate_posterior_probability",
    ifelse(posterior_table$posterior_probability_ge_0_5 >= 0.5, "weak_posterior_probability", "low_posterior_probability"))
)
posterior_table$cancer_type <- posterior_table$cohort_id

posterior_table <- posterior_table[, c(
  "cell_line_id",
  "tumour_sample_id",
  "cancer_type",
  "support_count",
  "total_representation_count",
  "p_consensus_observed",
  "weighted_support",
  "effective_representation_count",
  "prior_alpha",
  "prior_beta",
  "posterior_alpha",
  "posterior_beta",
  "posterior_mean_neighbourhood_probability",
  "posterior_probability_ge_0_5",
  "posterior_probability_ge_0_7",
  "posterior_probability_ge_0_8",
  "posterior_interval_lower",
  "posterior_interval_upper",
  "probabilistic_neighbourhood_class"
)]

write_tsv_gzip(posterior_table, posterior_output)

summary_table <- aggregate(
  posterior_mean_neighbourhood_probability ~ cancer_type + probabilistic_neighbourhood_class,
  data = posterior_table,
  FUN = function(x) c(n = length(x), median = stats::median(x), mean = mean(x))
)
summary_table <- do.call(data.frame, summary_table)
names(summary_table) <- c("cancer_type", "probabilistic_neighbourhood_class", "n_pairs", "median_posterior_mean", "mean_posterior_mean")
write_tsv_table(summary_table, summary_output)

report_lines <- c(
  "# Probabilistic neighbourhood report",
  "",
  "## Scope",
  "",
  "The probabilistic neighbourhood step treated thesis-framework p_consensus values as observations and estimated posterior neighbourhood probabilities with a Beta-binomial approximation.",
  "",
  "## Probability interpretation",
  "",
  "`p_consensus_observed` is an observed support proportion from the thesis framework, not a posterior probability. Posterior quantities in this report are produced only after modelling representation evidence with priors and reliability weights.",
  "",
  "## Model",
  "",
  paste0("- Prior: Beta(", prior_alpha, ", ", prior_beta, ")."),
  paste0("- Representation evidence threshold for support calls: ", evidence_threshold, "."),
  "- Weighted support equals the sum of representation reliability weights for representations with observed support calls.",
  "",
  "## Outputs",
  "",
  paste0("- `", posterior_output, "`"),
  paste0("- `", summary_output, "`"),
  "",
  "## Manual inspection required",
  "",
  "- Inspect pairs with wide posterior intervals before treating them as stable tumour-neighbourhood relations.",
  "- Inspect pairs where high observed p_consensus values do not translate into high posterior probabilities because the contributing representations were downweighted."
)
write_markdown(report_lines, report_output)
cat("Posterior neighbourhood probabilities completed: ", posterior_output, "\n", sep = "")
