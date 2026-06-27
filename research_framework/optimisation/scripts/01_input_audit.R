command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- command_arguments[grepl("^--file=", command_arguments)]
current_script_path <- if (length(file_argument) > 0) sub("^--file=", "", file_argument[[1]]) else file.path("research_framework", "optimisation", "scripts", "01_input_audit.R")
source(file.path(dirname(normalizePath(current_script_path, mustWork = TRUE)), "optimisation_utils.R"))
load_required_packages(c("data.table", "yaml"))

arguments <- parse_key_value_arguments()
config <- load_optimisation_config(arguments$config)

results_root <- config_value(config, c("inputs", "thesis_results_root"), "results")
output_root <- config_value(config, c("outputs", "optimisation_results_root"), file.path("research_framework", "optimisation", "results"))
docs_root <- config_value(config, c("outputs", "optimisation_docs_root"), file.path("research_framework", "optimisation", "docs"))

audit_output <- file.path(output_root, "input_audit", "optimisation_input_audit.tsv")
report_output <- file.path(docs_root, "01_input_audit_report.md")

expected_inputs <- data.frame(
  input_key = c(
    "p_consensus_direction_summary",
    "p_consensus_cellline_direction_summary",
    "p_consensus_best_cell_lines_ranked",
    "final_consensus_neighbourhoods",
    "representation_graph_edges",
    "representation_node_annotations",
    "resolved_graph_metadata",
    "isolate_tables",
    "bridge_like_anchor_tables",
    "marker_feature_panel",
    "tumour_to_cellline_rankings",
    "cellline_to_tumour_rankings",
    "sample_metadata"
  ),
  requirement_level = c(
    "required",
    "required",
    "optional",
    "required",
    "required",
    "optional",
    "optional",
    "optional",
    "optional",
    "optional",
    "optional",
    "optional",
    "optional"
  ),
  expected_pattern = c(
    "^p_consensus_direction_summary\\.tsv$",
    "^p_consensus_cellline_direction_summary\\.tsv$",
    "^p_consensus_best_cell_lines_ranked\\.tsv$",
    "^Final_consensus_tumour_neighbourhoods_.*\\.tsv$",
    "^cell_line_similarity_graph_edges_.*\\.tsv$",
    "^cell_line_similarity_graph_node_annotations_.*\\.tsv$",
    "metadata_with_components\\.tsv$|component_cell_line_summary\\.tsv$|pan_cancer_components\\.tsv$",
    "isolate.*\\.(tsv|csv)$",
    "bridge_like_anchor.*\\.tsv$",
    "ranked_marker_source_panel.*\\.(tsv|md)$",
    "tumour_to_cellline.*rank.*\\.tsv$|tumour_to_cell_line.*rank.*\\.tsv$",
    "cellline_to_tumour.*rank.*\\.tsv$|cell_line_to_tumour.*rank.*\\.tsv$",
    "metadata.*\\.tsv$|joint_metadata\\.tsv$"
  ),
  analysis_enabled_if_present = c(
    "representation weighting; model-selection context",
    "representation weighting; model-selection context",
    "ranking uncertainty",
    "representation diagnostics; posterior neighbourhood model",
    "probabilistic graph comparison",
    "graph metadata comparison",
    "component/isolate comparison",
    "probabilistic isolate comparison",
    "bridge-like anchor comparison",
    "feature-panel context",
    "tumour-to-cell-line ranking uncertainty",
    "cell-line-to-tumour ranking uncertainty",
    "source and cancer-type metadata checks"
  ),
  stringsAsFactors = FALSE
)

audit_rows <- lapply(seq_len(nrow(expected_inputs)), function(row_index) {
  pattern <- expected_inputs$expected_pattern[[row_index]]
  files_found <- discover_files(results_root, pattern)
  status <- if (length(files_found) > 0) {
    "present"
  } else if (expected_inputs$requirement_level[[row_index]] == "required") {
    "missing_required"
  } else {
    "missing_optional"
  }
  blocks_probabilistic_implementation <- status == "missing_required"
  analyses_skipped <- if (length(files_found) > 0) {
    ""
  } else if (blocks_probabilistic_implementation) {
    expected_inputs$analysis_enabled_if_present[[row_index]]
  } else {
    expected_inputs$analysis_enabled_if_present[[row_index]]
  }
  data.frame(
    input_key = expected_inputs$input_key[[row_index]],
    requirement_level = expected_inputs$requirement_level[[row_index]],
    expected_pattern = pattern,
    n_files_found = length(files_found),
    status = status,
    blocks_probabilistic_implementation = blocks_probabilistic_implementation,
    analyses_can_proceed = ifelse(length(files_found) > 0, expected_inputs$analysis_enabled_if_present[[row_index]], ""),
    analyses_skipped = ifelse(length(files_found) > 0, "", analyses_skipped),
    example_files = collapse_examples(files_found),
    stringsAsFactors = FALSE
  )
})

audit_table <- do.call(rbind, audit_rows)
write_tsv_table(audit_table, audit_output)

required_missing <- audit_table[audit_table$status == "missing_required", , drop = FALSE]
optional_missing <- audit_table[audit_table$status == "missing_optional", , drop = FALSE]
present_inputs <- audit_table[audit_table$status == "present", , drop = FALSE]

report_lines <- c(
  "# Input audit report",
  "",
  "## Scope",
  "",
  "The input audit located thesis-framework outputs available to the optimisation layer. The audit did not modify thesis-framework outputs.",
  "",
  "## Thesis outputs found",
  "",
  if (nrow(present_inputs) == 0) "- No expected thesis outputs were found." else paste0("- `", present_inputs$input_key, "`: ", present_inputs$n_files_found, " file(s). Example: `", present_inputs$example_files, "`."),
  "",
  "## Expected outputs missing",
  "",
  if (nrow(required_missing) == 0 && nrow(optional_missing) == 0) "- No expected inputs were missing." else c(
    if (nrow(required_missing) > 0) paste0("- Required missing `", required_missing$input_key, "` blocks: ", required_missing$analysis_enabled_if_present, ".") else "- No required inputs were missing.",
    if (nrow(optional_missing) > 0) paste0("- Optional missing `", optional_missing$input_key, "` skips: ", optional_missing$analysis_enabled_if_present, ".") else "- No optional inputs were missing."
  ),
  "",
  "## Blocking assessment",
  "",
  if (nrow(required_missing) == 0) "No missing required input currently blocks the probabilistic implementation." else "At least one required input is missing. The blocked analyses are listed in the audit table.",
  "",
  "## Analyses that can proceed",
  "",
  paste0("- ", present_inputs$analyses_can_proceed),
  "",
  "## Files or rules to inspect if inputs are missing",
  "",
  "- Check thesis workflow rules that write final p-consensus tables, representation-specific graph edges, and tumour-neighbourhood final consensus files.",
  "- Check whether optional ranking, isolate, bridge-like anchor, or marker-panel outputs were generated for the current repository state.",
  "",
  "## Machine-readable output",
  "",
  paste0("- `", audit_output, "`")
)

write_markdown(report_lines, report_output)
cat("Input audit completed: ", audit_output, "\n", sep = "")
