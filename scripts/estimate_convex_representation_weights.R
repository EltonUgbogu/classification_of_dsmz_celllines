#!/usr/bin/env Rscript

# The script estimates non-negative simplex-normalised representation weights.
# It uses a deterministic pilot objective when no external convex solver is used.

parse_arguments <- function() {
  argument_values <- commandArgs(trailingOnly = TRUE)
  parsed_arguments <- list()
  argument_index <- 1L
  while (argument_index <= length(argument_values)) {
    argument_name <- argument_values[[argument_index]]
    if (!startsWith(argument_name, "--")) {
      stop("Unexpected argument: ", argument_name)
    }
    option_name <- sub("^--", "", argument_name)
    if (argument_index == length(argument_values)) {
      stop("Missing value for argument: ", argument_name)
    }
    parsed_arguments[[option_name]] <- argument_values[[argument_index + 1L]]
    argument_index <- argument_index + 2L
  }
  parsed_arguments
}

require_argument <- function(arguments, name) {
  value <- arguments[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required argument --", name)
  }
  value
}

numeric_argument <- function(arguments, name, default_value) {
  value <- arguments[[name]]
  if (is.null(value) || !nzchar(value)) {
    return(default_value)
  }
  numeric_value <- suppressWarnings(as.numeric(value))
  if (!is.finite(numeric_value)) {
    stop("Argument --", name, " must be numeric")
  }
  numeric_value
}

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, quote = "", comment.char = "")
}

write_tsv <- function(data_frame, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data_frame, path, sep = "\t", row.names = FALSE, quote = FALSE,
              na = "")
}

as_numeric_safe <- function(values) {
  numeric_values <- suppressWarnings(as.numeric(values))
  numeric_values[!is.finite(numeric_values)] <- 0
  numeric_values
}

arguments <- parse_arguments()
representation_metrics_path <- require_argument(arguments, "representation_metrics")
output_path <- require_argument(arguments, "output")

stability_reward <- numeric_argument(arguments, "stability-reward", 0.40)
agreement_reward <- numeric_argument(arguments, "agreement-reward", 0.35)
redundancy_penalty <- numeric_argument(arguments, "redundancy-penalty", 0.15)
batch_penalty <- numeric_argument(arguments, "batch-penalty", 0.05)
instability_penalty <- numeric_argument(arguments, "instability-penalty", 0.05)

representation_metrics <- read_tsv(representation_metrics_path)
required_columns <- c(
  "representation_id",
  "neighbourhood_stability_score",
  "cancer_type_agreement_score",
  "representation_redundancy_index",
  "batch_association_score",
  "instability_score"
)
missing_columns <- setdiff(required_columns, names(representation_metrics))
if (length(missing_columns)) {
  stop("Representation metrics file is missing required columns: ", paste(missing_columns, collapse = ", "))
}

objective_raw <- (
  stability_reward * as_numeric_safe(representation_metrics$neighbourhood_stability_score) +
    agreement_reward * as_numeric_safe(representation_metrics$cancer_type_agreement_score) -
    redundancy_penalty * as_numeric_safe(representation_metrics$representation_redundancy_index) -
    batch_penalty * as_numeric_safe(representation_metrics$batch_association_score) -
    instability_penalty * as_numeric_safe(representation_metrics$instability_score)
)

objective_shifted <- pmax(objective_raw, 0)
if (sum(objective_shifted) <= 0) {
  representation_weight <- rep(1 / nrow(representation_metrics), nrow(representation_metrics))
  optimisation_status <- "uniform_fallback_all_objectives_nonpositive"
} else {
  representation_weight <- objective_shifted / sum(objective_shifted)
  optimisation_status <- "deterministic_simplex_normalised_pilot"
}

weights <- data.frame(
  representation_id = representation_metrics$representation_id,
  representation_weight = representation_weight,
  objective_raw = objective_raw,
  objective_nonnegative = objective_shifted,
  neighbourhood_stability_score = as_numeric_safe(representation_metrics$neighbourhood_stability_score),
  cancer_type_agreement_score = as_numeric_safe(representation_metrics$cancer_type_agreement_score),
  representation_redundancy_index = as_numeric_safe(representation_metrics$representation_redundancy_index),
  batch_association_score = as_numeric_safe(representation_metrics$batch_association_score),
  instability_score = as_numeric_safe(representation_metrics$instability_score),
  stability_reward = stability_reward,
  cancer_type_agreement_reward = agreement_reward,
  redundancy_penalty = redundancy_penalty,
  batch_penalty = batch_penalty,
  instability_penalty = instability_penalty,
  optimisation_status = optimisation_status,
  simplex_constraint = "representation_weight >= 0; sum(representation_weight) = 1",
  scientific_status = "pilot_approximation_not_scientifically_validated",
  stringsAsFactors = FALSE
)

write_tsv(weights, output_path)

