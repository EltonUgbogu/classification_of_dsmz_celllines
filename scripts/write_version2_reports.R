#!/usr/bin/env Rscript

# The script writes Version 2 method-development reports from generated tables.
# It writes documentation only in docs/version2/ and the completion flag only
# in results/version2/.

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

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, quote = "", comment.char = "")
}

write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
}

arguments <- parse_arguments()
representation_metrics <- read_tsv(require_argument(arguments, "representation_metrics"))
weights <- read_tsv(require_argument(arguments, "weights"))
weighted_edges <- read_tsv(require_argument(arguments, "weighted-edges"))
boundary_summary <- read_tsv(require_argument(arguments, "boundary-summary"))
graph_summary <- read_tsv(require_argument(arguments, "graph-summary"))
edge_comparison <- read_tsv(require_argument(arguments, "edge-comparison"))
ranking_comparison <- read_tsv(require_argument(arguments, "ranking-comparison"))
node_status_comparison <- read_tsv(require_argument(arguments, "node-status-comparison"))
version1_inputs_path <- require_argument(arguments, "version1-inputs")
documentation_reorganisation_path <- require_argument(arguments, "documentation-reorganisation")
methodology_out <- require_argument(arguments, "methodology-out")
weighting_report_out <- require_argument(arguments, "weighting-report-out")
checklist_out <- require_argument(arguments, "checklist-out")
status_out <- require_argument(arguments, "status-out")
complete_flag <- require_argument(arguments, "complete-flag")

methodology_lines <- c(
  "# Version 2 methodological summary",
  "",
  "Version 2 is a development workflow controlled by `Snakefile.v2`.",
  "",
  "## Scope",
  "",
  "- Convex representation weighting is implemented as a deterministic simplex-normalised pilot approximation.",
  "- Weighted consensus evidence is derived from protected Version 1 graph and ranking outputs.",
  "- Boundary cases are classified explicitly instead of being forced into binary decisions.",
  "- The probability-dependence graph is a development graph layer and not a validated biological network.",
  "- Version 1 versus Version 2 comparison tables are written for inspection.",
  "",
  "## Namespace",
  "",
  "- Version 2 results folder: `results/version2/`",
  "- Version 2 documentation folder: `docs/version2/`",
  "",
  "## Scientific status",
  "",
  "The implementation is suitable for development inspection. It should not be described as scientifically validated without external validation and review."
)
write_lines(methodology_lines, methodology_out)

weighting_lines <- c(
  "# Representation weighting report",
  "",
  "The Version 2 weighting layer estimates non-negative representation weights constrained to sum to one.",
  "",
  "## Objective",
  "",
  "The pilot objective rewards neighbourhood stability and cancer-type agreement, and penalises representation redundancy, batch association, and instability.",
  "",
  "## Weight table",
  "",
  paste0("- Representations: ", nrow(weights)),
  paste0("- Sum of representation weights: ", round(sum(as.numeric(weights$representation_weight), na.rm = TRUE), 8)),
  "",
  "## Representation metrics",
  "",
  paste0("- Metric rows: ", nrow(representation_metrics)),
  "",
  "## Scientific status",
  "",
  "No external convex solver is required for the current pilot. The method is documented as a deterministic simplex-normalised approximation."
)
write_lines(weighting_lines, weighting_report_out)

checklist_lines <- c(
  "# User inspection checklist",
  "",
  "- Confirm that Version 1 inputs listed in `docs/version2/version1_inputs_used_by_version2.md` are the intended protected thesis outputs.",
  "- Confirm that `results/version2/representation_weights.tsv` has non-negative weights summing to one.",
  "- Inspect `results/version2/weighted_consensus_edges.tsv` before interpreting retained, boundary, or rejected edges.",
  "- Inspect `results/version2/boundary_edges.tsv` and `results/version2/boundary_rankings.tsv` before any thesis-facing claims.",
  "- Confirm that the probability-dependence graph is described as a development graph layer, not a final biological interaction network.",
  "- Compare Version 1 and Version 2 outputs using the three `v1_v2_*_comparison.tsv` tables.",
  "- Do not rename `results/version2/` or `docs/version2/` to a method-specific name until a stable scientific name is approved."
)
write_lines(checklist_lines, checklist_out)

status_lines <- c(
  "# Version 2 implementation status",
  "",
  "## Controllers and namespaces",
  "",
  "- Version 1 workflow controller: `Snakefile`",
  "- Version 2 workflow controller: `Snakefile.v2`",
  "- Version 1 documentation folder: `docs/version1/`",
  "- Version 2 documentation folder: `docs/version2/`",
  "- Version 2 results folder: `results/version2/`",
  "- Version 1 outputs modified: none by the Version 2 workflow",
  "",
  "## Generated table counts",
  "",
  paste0("- Representation metric rows: ", nrow(representation_metrics)),
  paste0("- Representation weights rows: ", nrow(weights)),
  paste0("- Weighted consensus edge rows: ", nrow(weighted_edges)),
  paste0("- Boundary summary rows: ", nrow(boundary_summary)),
  paste0("- Graph summary rows: ", nrow(graph_summary)),
  paste0("- Edge comparison rows: ", nrow(edge_comparison)),
  paste0("- Ranking comparison rows: ", nrow(ranking_comparison)),
  paste0("- Node comparison rows: ", nrow(node_status_comparison)),
  "",
  "## Required inspection documents",
  "",
  paste0("- Version 1 input inventory: `", version1_inputs_path, "`"),
  paste0("- Documentation reorganisation report: `", documentation_reorganisation_path, "`"),
  "",
  "## Scientific assumptions",
  "",
  "- Representation weighting is a pilot simplex-normalised approximation.",
  "- Bootstrap support intervals are pilot threshold-margin bands unless replaced by resampling.",
  "- Boundary classes are development labels requiring user inspection.",
  "- The probability-dependence graph is not a causal graph and not a validated biological network.",
  "",
  "## Items requiring user inspection",
  "",
  "- Whether the Version 1 inputs selected for Version 2 are the intended protected thesis outputs.",
  "- Whether group-level ranking evidence should be mapped to profile-level graph nodes in a later implementation.",
  "- Whether a true convex solver should replace the pilot weighting approximation.",
  "- Whether bootstrap resampling should replace the threshold-margin support intervals."
)
write_lines(status_lines, status_out)

dir.create(dirname(complete_flag), recursive = TRUE, showWarnings = FALSE)
write_lines(c("version2_complete"), complete_flag)

