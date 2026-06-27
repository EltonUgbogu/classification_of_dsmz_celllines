command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- command_arguments[grepl("^--file=", command_arguments)]
current_script_path <- if (length(file_argument) > 0) sub("^--file=", "", file_argument[[1]]) else file.path("research_framework", "optimisation", "scripts", "08_write_optimisation_summary.R")
source(file.path(dirname(normalizePath(current_script_path, mustWork = TRUE)), "optimisation_utils.R"))
load_required_packages(c("data.table"))

arguments <- parse_key_value_arguments()
config <- load_optimisation_config(arguments$config)

output_root <- config_value(config, c("outputs", "optimisation_results_root"), file.path("research_framework", "optimisation", "results"))
docs_root <- config_value(config, c("outputs", "optimisation_docs_root"), file.path("research_framework", "optimisation", "docs"))
summary_report <- file.path(docs_root, "07_scientific_interpretation_notes.md")
portability_output <- file.path(output_root, "input_audit", "public_portability_audit.tsv")

read_optional <- function(path) {
  table <- read_table_if_exists(path)
  if (is.null(table)) data.frame() else table
}

audit_table <- read_optional(file.path(output_root, "input_audit", "optimisation_input_audit.tsv"))
redundancy_table <- read_optional(file.path(output_root, "representation_diagnostics", "representation_redundancy.tsv"))
stability_table <- read_optional(file.path(output_root, "representation_diagnostics", "representation_stability.tsv"))
weights_table <- read_optional(file.path(output_root, "representation_weights", "representation_weights.tsv"))
neighbourhood_summary <- read_optional(file.path(output_root, "probabilistic_neighbourhoods", "neighbourhood_probability_summary.tsv"))
edge_table <- read_optional(file.path(output_root, "probabilistic_graphs", "probabilistic_cellline_edges.tsv"))
resolved_edges <- read_optional(file.path(output_root, "probabilistic_graphs", "probabilistic_resolved_edges.tsv"))
isolates <- read_optional(file.path(output_root, "probabilistic_graphs", "probabilistic_isolates.tsv"))
bridge_anchors <- read_optional(file.path(output_root, "probabilistic_graphs", "probabilistic_bridge_anchors.tsv"))
ranking_summary <- read_optional(file.path(output_root, "ranking_uncertainty", "ranking_uncertainty_summary.tsv"))

scan_public_files <- function() {
  scan_roots <- c(
    "README.md",
    file.path("research_framework", "optimisation", "README.md"),
    file.path("research_framework", "optimisation", "config"),
    file.path("research_framework", "optimisation", "workflow"),
    file.path("research_framework", "optimisation", "scripts"),
    file.path("research_framework", "optimisation", "docs")
  )
  public_files <- unique(unlist(lapply(scan_roots, function(path) {
    if (file.exists(path) && !dir.exists(path)) {
      return(path)
    }
    if (dir.exists(path)) {
      return(list.files(path, recursive = TRUE, full.names = TRUE))
    }
    character(0)
  })))
  path_roots <- c("Users", "home", "work", "scratch", "mnt", "private")
  machine_path_pattern <- paste0("(/(", paste(path_roots, collapse = "|"), ")/|[A-Za-z]:\\\\\\\\)")
  rows <- lapply(public_files, function(path) {
    lines <- readLines(path, warn = FALSE)
    hits <- grep(machine_path_pattern, lines, value = TRUE)
    data.frame(
      file = path,
      n_machine_specific_path_hits = length(hits),
      example_hit = ifelse(length(hits) > 0, hits[[1]], ""),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

portability_table <- scan_public_files()
write_tsv_table(portability_table, portability_output)
if (any(portability_table$n_machine_specific_path_hits > 0)) {
  stop("Portability audit detected local or machine-specific absolute paths in public-facing files.", call. = FALSE)
}

high_weight <- if (nrow(weights_table) > 0) weights_table[weights_table$weight_class == "high_weight", , drop = FALSE] else data.frame()
high_edges <- if (nrow(edge_table) > 0) edge_table[edge_table$edge_probability_class == "high_probability", , drop = FALSE] else data.frame()
uncertain_edges <- if (nrow(edge_table) > 0) edge_table[edge_table$edge_probability_class == "uncertain_probability", , drop = FALSE] else data.frame()
probabilistic_isolates <- if (nrow(isolates) > 0) isolates[as.character(isolates$probabilistic_isolate_status) == "TRUE", , drop = FALSE] else data.frame()
bridge_hits <- if (nrow(bridge_anchors) > 0) bridge_anchors[as.character(bridge_anchors$probabilistic_bridge_anchor_status) == "TRUE", , drop = FALSE] else data.frame()

report_lines <- c(
  "# Scientific interpretation notes",
  "",
  "## Scope statement",
  "",
  "The optimisation framework is a research layer above the thesis framework. It does not replace the thesis framework.",
  "",
  "## Probability interpretation",
  "",
  "`p_consensus` is an observed support proportion from the thesis framework and must not be interpreted as a posterior probability. Posterior probabilities are generated only by the optimisation model.",
  "",
  "## Representation diagnostics",
  "",
  paste0("- Representations scored for redundancy: ", nrow(redundancy_table), "."),
  paste0("- Representations scored for stability: ", nrow(stability_table), "."),
  paste0("- High-weight representations: ", nrow(high_weight), "."),
  "",
  "## Threshold validation",
  "",
  "Probabilistic graph thresholds were selected by expected false edge-rate control. The threshold sweep table reports the expected false edge rate for each candidate threshold.",
  "",
  "## Probabilistic neighbourhoods",
  "",
  paste0("- Neighbourhood probability classes reported: ", paste(unique(neighbourhood_summary$probabilistic_neighbourhood_class), collapse = ", "), "."),
  "",
  "## Probabilistic graph model",
  "",
  paste0("- High-probability graph edges: ", nrow(high_edges), "."),
  paste0("- Uncertain graph edges: ", nrow(uncertain_edges), "."),
  "The current implementation is PDG-inspired and focuses on probabilistic neighbourhood variables and graph-edge variables.",
  "",
  "## Isolates and bridge-like anchors",
  "",
  paste0("- Probabilistic isolates: ", nrow(probabilistic_isolates), "."),
  paste0("- Probabilistic bridge-like anchors: ", nrow(bridge_hits), "."),
  "",
  "## Ranking uncertainty",
  "",
  if (nrow(ranking_summary) == 0) "- Ranking uncertainty was not estimated." else paste0("- ", ranking_summary$ranking_uncertainty_class, ": ", ranking_summary$n_cell_lines, " cell line(s)."),
  "",
  "## Questions the user must be prepared to answer",
  "",
  "### Is p_consensus a probability?",
  "",
  "No. It is an observed thesis-framework support proportion. Posterior probabilities are estimated only in the optimisation layer.",
  "",
  "### Why was a probabilistic graph model not used as the main thesis framework?",
  "",
  "The thesis framework was designed as a deterministic and consensus-based empirical baseline. The probabilistic graph model is a research extension layered above those outputs.",
  "",
  "### How would the deterministic thesis framework be extended probabilistically?",
  "",
  "The deterministic outputs can be treated as observations, representation reliability can be estimated, and posterior neighbourhood and graph-edge probabilities can be inferred from weighted evidence.",
  "",
  "### How are thresholds selected without being arbitrary?",
  "",
  "The probabilistic graph threshold is selected by controlling the expected false edge rate across retained edges.",
  "",
  "### How does the probabilistic model handle disagreement between feature-distance representations?",
  "",
  "Disagreement is represented through representation-specific evidence, redundancy diagnostics, reliability weights, and posterior uncertainty intervals.",
  "",
  "### What is the difference between representation support and posterior support?",
  "",
  "Representation support is an observed call or value from a thesis representation. Posterior support is inferred after modelling representation evidence with priors and reliability weights.",
  "",
  "### What does a probabilistic isolate mean?",
  "",
  "A probabilistic isolate is a cell line with no retained probabilistic graph edge at the selected threshold and an estimated probability of being disconnected under the model.",
  "",
  "### What does a probabilistic bridge-like anchor mean?",
  "",
  "A probabilistic bridge-like anchor is a node with high bridge score in the retained probabilistic graph topology. It is a graph uncertainty result and requires manual scientific inspection.",
  "",
  "## Scientific limitations",
  "",
  "- The current graph model is PDG-inspired, not a full formal probabilistic dependency graph solver.",
  "- Representation reliability weights depend on available thesis outputs and metadata completeness.",
  "- The null model currently uses within-cancer tumour-label permutation; other null models should be compared before publication-level claims.",
  "- Ranking uncertainty is limited by the structure of available thesis ranking tables.",
  "- A full PDG implementation would require explicit variable domains, conditional dependency claims, consistency objectives, and formal uncertainty propagation across all graph variables.",
  "",
  "## Public-release portability",
  "",
  "- The public README contains no detected local or HPC-specific absolute paths.",
  "- The optimisation README contains no detected local or HPC-specific absolute paths.",
  "- All config paths are repository-relative.",
  "- The optimisation Snakefile uses repository-relative paths.",
  "- Scripts accept paths from config or command-line arguments.",
  "- The workflow is intended to run from a fresh clone using repository-relative paths.",
  "- The thesis framework was not modified by this research layer.",
  "",
  "## Machine-readable portability audit",
  "",
  paste0("- `", portability_output, "`")
)
write_markdown(report_lines, summary_report)
cat("Optimisation scientific summary completed: ", summary_report, "\n", sep = "")
