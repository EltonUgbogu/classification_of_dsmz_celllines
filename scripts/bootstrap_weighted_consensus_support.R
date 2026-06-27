#!/usr/bin/env Rscript

# The script assigns pilot bootstrap-support intervals to weighted consensus
# edges. The interval is a deterministic threshold-margin approximation until
# a resampling implementation is added.

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

arguments <- parse_arguments()
weighted_edges_path <- require_argument(arguments, "weighted-edges")
output_path <- require_argument(arguments, "output")
retained_min <- numeric_argument(arguments, "retained-min", 0.70)
rejected_max <- numeric_argument(arguments, "rejected-max", 0.40)
interval_half_width <- numeric_argument(arguments, "interval-half-width", 0.08)
boundary_margin_band <- numeric_argument(arguments, "boundary-margin-band", 0.05)
weak_weighted_consensus <- numeric_argument(arguments, "weak-weighted-consensus", 0.55)

weighted_edges <- read_tsv(weighted_edges_path)
if (!"weighted_consensus_score" %in% names(weighted_edges)) {
  stop("Weighted edge table must contain weighted_consensus_score")
}

weighted_score <- suppressWarnings(as.numeric(weighted_edges$weighted_consensus_score))
weighted_score[!is.finite(weighted_score)] <- 0

bootstrap_support_lower <- pmax(0, weighted_score - interval_half_width)
bootstrap_support_upper <- pmin(1, weighted_score + interval_half_width)
threshold_overlap_flag <- (
  bootstrap_support_lower < retained_min &
    bootstrap_support_upper >= rejected_max
)

edge_resolution_status <- ifelse(
  bootstrap_support_lower >= retained_min,
  "retained",
  ifelse(bootstrap_support_upper < rejected_max, "rejected", "boundary")
)

boundary_reason <- ifelse(
  edge_resolution_status == "boundary" & threshold_overlap_flag,
  "threshold_ci_overlap",
  ifelse(
    edge_resolution_status == "boundary" & weighted_score < weak_weighted_consensus,
    "weak_weighted_consensus",
    ifelse(
      abs(weighted_score - retained_min) <= boundary_margin_band |
        abs(weighted_score - rejected_max) <= boundary_margin_band,
      "low_rank_margin",
      ""
    )
  )
)

weighted_edges$bootstrap_support_lower <- bootstrap_support_lower
weighted_edges$bootstrap_support_upper <- bootstrap_support_upper
weighted_edges$edge_resolution_status <- edge_resolution_status
weighted_edges$boundary_status <- ifelse(edge_resolution_status == "boundary", "boundary_edge", "not_boundary")
weighted_edges$boundary_reason <- boundary_reason
weighted_edges$empirical_threshold <- retained_min
weighted_edges$threshold_overlap_flag <- threshold_overlap_flag
weighted_edges$bootstrap_method <- "pilot_threshold_margin_band"
weighted_edges$scientific_status <- "pilot_approximation_not_scientifically_validated"

write_tsv(weighted_edges, output_path)

