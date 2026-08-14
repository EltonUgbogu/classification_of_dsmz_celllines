#!/usr/bin/env Rscript

# Validate prepared graph-defined DESeq2 inputs before any DESeq2 model fitting.
# Unit of analysis: prepared count sample and resolved biological cell-line
# identity. The validator checks interface contracts only and writes numerical
# count metrics separately from PASS/FAIL checks.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

option_list <- list(
  make_option("--cohort", type = "character"),
  make_option("--counts", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--isolate_list", type = "character"),
  make_option("--anchor_list", type = "character"),
  make_option("--anchor_components", type = "character"),
  make_option("--components_list", type = "character"),
  make_option("--node_stats", type = "character"),
  make_option("--anchor_audit", type = "character"),
  make_option("--provenance", type = "character"),
  make_option("--sample_id_col", type = "character", default = "sample_id"),
  make_option("--cell_line_col", type = "character", default = "cell_line"),
  make_option("--component_col", type = "character", default = "component"),
  make_option("--expected_samples", type = "character", default = ""),
  make_option("--cell_line_verification_mode", type = "character", default = "source_population_declaration"),
  make_option("--cell_line_verification_column", type = "character", default = ""),
  make_option("--accepted_cell_line_values", type = "character", default = "Cell Line"),
  make_option("--output", type = "character"),
  make_option("--metrics_output", type = "character")
)

opts <- parse_args(OptionParser(option_list = option_list))

checks <- data.frame(check = character(), status = character(), details = character(), stringsAsFactors = FALSE)
metrics <- data.frame(metric = character(), value = character(), stringsAsFactors = FALSE)

add_check <- function(name, ok, details) {
  checks <<- rbind(
    checks,
    data.frame(
      check = name,
      status = if (isTRUE(ok)) "PASS" else "FAIL",
      details = as.character(details %||% ""),
      stringsAsFactors = FALSE
    )
  )
}

add_metric <- function(name, value) {
  metrics <<- rbind(
    metrics,
    data.frame(metric = name, value = as.character(value), stringsAsFactors = FALSE)
  )
}

canonicalise_cell_line_identifier <- function(value) {
  value <- toupper(trimws(as.character(value)))
  gsub("[^A-Z0-9]+", "", value)
}

truthy <- function(value) {
  toupper(trimws(as.character(value))) %in% c("TRUE", "T", "1", "YES", "Y")
}

read_required_table <- function(path, label) {
  if (is.null(path) || !file.exists(path)) {
    add_check(paste0(label, "_exists"), FALSE, path %||% "")
    return(NULL)
  }
  table <- tryCatch(
    fread(path, sep = "\t", data.table = FALSE, check.names = FALSE),
    error = function(error) error
  )
  add_check(paste0(label, "_readable"), !inherits(table, "error"),
            if (inherits(table, "error")) table$message else path)
  if (inherits(table, "error")) NULL else table
}

read_text_list <- function(path, sep = "\n") {
  if (!file.exists(path)) return(character(0))
  text <- paste(readLines(path, warn = FALSE), collapse = sep)
  values <- unique(trimws(unlist(strsplit(text, "[,\n\r]+"))))
  values[nzchar(values)]
}

require_columns <- function(table, columns, check_name) {
  missing <- setdiff(columns, names(table))
  add_check(check_name, length(missing) == 0L,
            if (length(missing)) paste(missing, collapse = ",") else paste(columns, collapse = ","))
  length(missing) == 0L
}

set_detail <- function(left, right, left_name, right_name) {
  paste0(
    "missing_from_", left_name, "=",
    paste(setdiff(right, left), collapse = ","),
    "; missing_from_", right_name, "=",
    paste(setdiff(left, right), collapse = ",")
  )
}

counts <- read_required_table(opts$counts, "counts")
metadata <- read_required_table(opts$metadata, "metadata")
node_stats <- read_required_table(opts$node_stats, "node_stats")
anchor_audit <- read_required_table(opts$anchor_audit, "anchor_audit")
provenance <- read_required_table(opts$provenance, "preparation_provenance")

count_sample_ids <- character()
if (is.data.frame(counts)) {
  add_check("counts_has_one_gene_identifier_column", ncol(counts) >= 2L && names(counts)[1L] == "gene_id",
            paste(names(counts)[seq_len(min(ncol(counts), 4L))], collapse = ","))
  if (ncol(counts) >= 2L) {
    count_sample_ids <- names(counts)[-1L]
    add_check("count_sample_identifiers_present",
              all(!is.na(count_sample_ids)) && all(nzchar(trimws(count_sample_ids))),
              paste(length(count_sample_ids), "count sample columns"))
    add_check("count_sample_identifiers_unique",
              anyDuplicated(count_sample_ids) == 0L,
              paste("duplicates", paste(unique(count_sample_ids[duplicated(count_sample_ids)]), collapse = ",")))
    gene_ids <- as.character(counts[[1L]])
    add_check("gene_identifiers_present",
              all(!is.na(gene_ids)) && all(nzchar(trimws(gene_ids))),
              paste("genes", length(gene_ids)))
    add_check("gene_identifiers_unique",
              anyDuplicated(gene_ids) == 0L,
              paste("duplicates", sum(duplicated(gene_ids))))
    count_matrix <- as.matrix(counts[, -1L, drop = FALSE])
    suppressWarnings(storage.mode(count_matrix) <- "numeric")
    add_check("counts_numeric", !any(is.na(count_matrix)), "all count columns coercible to numeric")
    add_check("counts_finite", all(is.finite(count_matrix)), "finite count values")
    add_check("counts_non_negative", all(count_matrix >= 0), "non-negative count values")
    add_check("counts_integer_valued", all(abs(count_matrix - round(count_matrix)) <= 1e-8),
              "integer-valued count values")
    add_metric("minimum_count", suppressWarnings(min(count_matrix, na.rm = TRUE)))
    add_metric("maximum_count", suppressWarnings(max(count_matrix, na.rm = TRUE)))
    add_metric("number_of_genes", nrow(counts))
    add_metric("number_of_samples", length(count_sample_ids))
  }
}

metadata_sample_ids <- character()
if (is.data.frame(metadata)) {
  if (require_columns(metadata, c(opts$sample_id_col, opts$cell_line_col, opts$component_col, "is_isolate"),
                      "metadata_has_required_columns")) {
    metadata_sample_ids <- trimws(as.character(metadata[[opts$sample_id_col]]))
    add_check("metadata_sample_identifiers_present",
              all(!is.na(metadata_sample_ids)) && all(nzchar(metadata_sample_ids)),
              paste(length(metadata_sample_ids), "metadata samples"))
    add_check("metadata_sample_identifiers_unique",
              anyDuplicated(metadata_sample_ids) == 0L,
              paste("duplicates", paste(unique(metadata_sample_ids[duplicated(metadata_sample_ids)]), collapse = ",")))
    if (length(count_sample_ids) > 0L && anyDuplicated(count_sample_ids) == 0L &&
        anyDuplicated(metadata_sample_ids) == 0L) {
      add_check("count_metadata_sample_sets_equal",
                setequal(count_sample_ids, metadata_sample_ids),
                set_detail(count_sample_ids, metadata_sample_ids, "counts", "metadata"))
    }
    if (nzchar(opts$expected_samples)) {
      expected_samples <- suppressWarnings(as.integer(opts$expected_samples))
      add_check("counts_expected_sample_count",
                !is.na(expected_samples) && length(count_sample_ids) == expected_samples,
                paste("observed", length(count_sample_ids), "expected", expected_samples))
    }
  }
}

node_map <- NULL
if (is.data.frame(node_stats) &&
    require_columns(node_stats, c("cell_line", "component", "is_isolate"), "node_stats_has_required_columns")) {
  graph_cell_lines <- trimws(as.character(node_stats$cell_line))
  node_map <- unique(data.frame(
    cell_line = graph_cell_lines,
    component = trimws(as.character(node_stats$component)),
    is_isolate = truthy(node_stats$is_isolate),
    stringsAsFactors = FALSE
  ))
  conflicts <- aggregate(
    paste(component, is_isolate, sep = "::") ~ cell_line,
    node_map,
    function(values) length(unique(values))
  )
  names(conflicts) <- c("cell_line", "n_assignments")
  add_check("graph_mapping_unique_component_and_isolate_status",
            all(conflicts$n_assignments == 1L),
            paste(conflicts$cell_line[conflicts$n_assignments > 1L], collapse = ","))
}

if (is.data.frame(metadata) && !is.null(node_map) &&
    all(c(opts$cell_line_col, opts$component_col, "is_isolate") %in% names(metadata))) {
  metadata_cell_lines <- trimws(as.character(metadata[[opts$cell_line_col]]))
  graph_cell_lines <- trimws(as.character(node_map$cell_line))
  metadata_canonical <- canonicalise_cell_line_identifier(metadata_cell_lines)
  graph_canonical <- canonicalise_cell_line_identifier(graph_cell_lines)
  add_check("graph_cell_lines_missing_from_staged_metadata",
            length(setdiff(unique(graph_canonical), unique(metadata_canonical))) == 0L,
            paste(setdiff(unique(graph_canonical), unique(metadata_canonical)), collapse = ","))
  add_check("staged_metadata_cell_lines_missing_from_graph",
            length(setdiff(unique(metadata_canonical), unique(graph_canonical))) == 0L,
            paste(setdiff(unique(metadata_canonical), unique(graph_canonical)), collapse = ","))

  metadata_map <- unique(data.frame(
    cell_line = metadata_cell_lines,
    component = trimws(as.character(metadata[[opts$component_col]])),
    is_isolate = truthy(metadata$is_isolate),
    stringsAsFactors = FALSE
  ))
  merged <- merge(metadata_map, node_map, by = "cell_line", all.x = TRUE, suffixes = c("_metadata", "_graph"))
  add_check("metadata_cell_lines_map_to_current_graph",
            !any(is.na(merged$component_graph)),
            paste(unique(merged$cell_line[is.na(merged$component_graph)]), collapse = ","))
  add_check("metadata_components_match_current_graph",
            all(merged$component_metadata == merged$component_graph, na.rm = TRUE),
            paste(unique(merged$cell_line[merged$component_metadata != merged$component_graph]), collapse = ","))
  add_check("metadata_isolate_flags_match_current_graph",
            all(merged$is_isolate_metadata == merged$is_isolate_graph, na.rm = TRUE),
            paste(unique(merged$cell_line[merged$is_isolate_metadata != merged$is_isolate_graph]), collapse = ","))

  component_cell_lines <- unique(
    node_map[
      !node_map$is_isolate & !node_map$component %in% c("", "NA", "-1"),
      c("component", "cell_line")
    ]
  )
  if (nrow(component_cell_lines) == 0L) {
    eligible_components <- character(0)
  } else {
    component_counts <- aggregate(
      cell_line ~ component,
      component_cell_lines,
      function(values) length(unique(values))
    )
    names(component_counts) <- c("component", "n_cell_lines")
    eligible_components <- sort(component_counts$component[component_counts$n_cell_lines >= 2L])
  }
  staged_components <- sort(read_text_list(opts$components_list))
  add_check("component_eligibility_uses_unique_biological_cell_lines",
            setequal(staged_components, eligible_components),
            set_detail(staged_components, eligible_components, "staged_components", "graph_eligible_components"))
}

if (is.data.frame(provenance) && require_columns(provenance, c("field", "value"), "provenance_has_required_columns")) {
  provenance_values <- setNames(as.character(provenance$value), as.character(provenance$field))
  add_check("raw_count_source_kind_declared",
            identical(provenance_values[["source_count_kind"]], "raw_gene_level_counts"),
            provenance_values[["source_count_kind"]] %||% "")
  add_check("cell_line_only_provenance_established",
            identical(provenance_values[["sample_population"]], "cell_line_only"),
            provenance_values[["sample_population"]] %||% "")
  add_check("value_transformation_subset_and_reorder_only",
            identical(provenance_values[["value_transformation"]], "subset_and_reorder_only"),
            provenance_values[["value_transformation"]] %||% "")
  add_check("normalisation_not_applied",
            identical(provenance_values[["normalisation_applied"]], "FALSE"),
            provenance_values[["normalisation_applied"]] %||% "")
  add_check("aggregation_not_applied",
            identical(provenance_values[["aggregation_applied"]], "FALSE"),
            provenance_values[["aggregation_applied"]] %||% "")
  add_check("value_preservation_check_passed",
            identical(provenance_values[["value_preservation_check"]], "prepared_values_identical_to_selected_source_values"),
            provenance_values[["value_preservation_check"]] %||% "")
}

if (opts$cell_line_verification_mode == "metadata_column" && is.data.frame(metadata)) {
  verification_column <- opts$cell_line_verification_column
  accepted_values <- trimws(unlist(strsplit(opts$accepted_cell_line_values, ",")))
  if (require_columns(metadata, verification_column, "metadata_has_cell_line_verification_column")) {
    observed <- trimws(as.character(metadata[[verification_column]]))
    add_check("metadata_cell_line_verification_values_accepted",
              all(observed %in% accepted_values),
              paste(unique(observed[!(observed %in% accepted_values)]), collapse = ","))
  }
}

staged_isolates <- read_text_list(opts$isolate_list, sep = ",")
staged_anchors <- read_text_list(opts$anchor_list, sep = ",")
has_graph_derived_contrasts <- length(staged_isolates) > 0L || length(staged_anchors) > 0L
add_check("graph_derived_contrasts_nonempty",
          has_graph_derived_contrasts,
          if (has_graph_derived_contrasts) {
            paste("n_isolates", length(staged_isolates), "n_anchors", length(staged_anchors))
          } else {
            "No isolate or anchor contrasts were derived for this cancer type. This is a methodological state requiring review."
          })

if (is.data.frame(metadata) && opts$cell_line_col %in% names(metadata)) {
  prepared_cell_lines <- unique(trimws(as.character(metadata[[opts$cell_line_col]])))
  add_check("isolate_list_members_present_in_prepared_metadata",
            all(staged_isolates %in% prepared_cell_lines),
            paste(setdiff(staged_isolates, prepared_cell_lines), collapse = ","))
  add_check("anchor_list_members_present_in_prepared_metadata",
            all(staged_anchors %in% prepared_cell_lines),
            paste(setdiff(staged_anchors, prepared_cell_lines), collapse = ","))
}

if (length(staged_anchors) > 0L && is.data.frame(anchor_audit) && !is.null(node_map)) {
  anchor_components <- read_required_table(opts$anchor_components, "anchor_components")
  if (is.data.frame(anchor_components) &&
      require_columns(anchor_components, c("anchor", "component"), "anchor_components_has_required_columns") &&
      require_columns(anchor_audit, c("node_id", "component_id", "anchor_selected"), "anchor_audit_has_required_columns")) {
    anchor_components$anchor <- trimws(as.character(anchor_components$anchor))
    anchor_components$component <- trimws(as.character(anchor_components$component))
    selected_audit <- anchor_audit[truthy(anchor_audit$anchor_selected), , drop = FALSE]
    selected_audit$node_id <- trimws(as.character(selected_audit$node_id))
    selected_audit$component_id <- trimws(as.character(selected_audit$component_id))
    add_check("anchor_list_matches_current_anchor_audit",
              setequal(staged_anchors, selected_audit$node_id),
              set_detail(staged_anchors, selected_audit$node_id, "staged_anchor_list", "selected_anchor_audit"))
    merged_anchor <- merge(anchor_components, selected_audit, by.x = "anchor", by.y = "node_id", all.x = TRUE)
    add_check("anchor_component_membership_matches_graph_statistics",
              all(merged_anchor$component == merged_anchor$component_id, na.rm = TRUE),
              paste(merged_anchor$anchor[merged_anchor$component != merged_anchor$component_id], collapse = ","))
  }
}

fwrite(checks, opts$output, sep = "\t")
fwrite(metrics, opts$metrics_output, sep = "\t")

failures <- checks[checks$status == "FAIL", , drop = FALSE]
if (nrow(failures) > 0L) {
  apply(
    failures,
    1,
    function(row) {
      message("[VALIDATION FAIL] ", row[["check"]], ": ", row[["details"]])
    }
  )
  message(sprintf("[VALIDATION FAIL] %d failed checks", nrow(failures)))
  quit(save = "no", status = 1L)
}

message(sprintf("[DESeq2 validation] PASS cohort=%s checks=%d", opts$cohort, nrow(checks)))
