#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }
    key <- sub("^--", "", key)
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- x[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

opts <- parse_args(args)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x) || identical(x, "")) y else x

required <- c(
  "cohort", "counts", "metadata", "isolate_list", "anchor_list",
  "anchor_components", "components_list", "node_stats", "anchor_audit",
  "sample_id_col", "cell_line_col", "component_col", "output"
)
missing_args <- setdiff(required, names(opts))
if (length(missing_args)) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "), call. = FALSE)
}

cohort <- opts$cohort
sample_id_col <- opts$sample_id_col
cell_line_col <- opts$cell_line_col
component_col <- opts$component_col
expected_samples <- suppressWarnings(as.integer(opts$expected_samples %||% NA_character_))

checks <- data.frame(
  cohort = character(),
  check = character(),
  status = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

add_check <- function(name, ok, detail) {
  checks[nrow(checks) + 1L, ] <<- list(
    cohort = cohort,
    check = name,
    status = if (isTRUE(ok)) "PASS" else "FAIL",
    detail = as.character(detail)
  )
}

normalise_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[!is.na(x) & nzchar(x)]
}

set_detail <- function(a, b, label_a = "staged", label_b = "current") {
  miss <- setdiff(a, b)
  extra <- setdiff(b, a)
  paste0(
    label_a, "_only=", if (length(miss)) paste(miss, collapse = ",") else "none",
    "; ", label_b, "_only=", if (length(extra)) paste(extra, collapse = ",") else "none"
  )
}

read_csv_list <- function(path) {
  txt <- if (file.exists(path)) paste(readLines(path, warn = FALSE), collapse = ",") else ""
  normalise_id(unlist(strsplit(txt, ",", fixed = TRUE), use.names = FALSE))
}

read_text_list <- function(path) {
  normalise_id(readLines(path, warn = FALSE))
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

file_paths <- c(
  counts = opts$counts,
  metadata = opts$metadata,
  isolate_list = opts$isolate_list,
  anchor_list = opts$anchor_list,
  anchor_components = opts$anchor_components,
  components_list = opts$components_list,
  node_stats = opts$node_stats,
  anchor_audit = opts$anchor_audit
)
for (nm in names(file_paths)) {
  add_check(paste0(nm, "_exists"), file.exists(file_paths[[nm]]), file_paths[[nm]])
}

counts_df <- NULL
metadata <- NULL
node_stats <- NULL
anchor_audit <- NULL

if (file.exists(opts$counts)) {
  counts_df <- tryCatch(
    read.delim(opts$counts, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) e
  )
  add_check("counts_readable", !inherits(counts_df, "error"), if (inherits(counts_df, "error")) counts_df$message else opts$counts)
}

if (file.exists(opts$metadata)) {
  metadata <- tryCatch(
    read.delim(opts$metadata, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) e
  )
  add_check("metadata_readable", !inherits(metadata, "error"), if (inherits(metadata, "error")) metadata$message else opts$metadata)
}

if (file.exists(opts$node_stats)) {
  node_stats <- tryCatch(
    read.delim(opts$node_stats, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) e
  )
  add_check("node_stats_readable", !inherits(node_stats, "error"), if (inherits(node_stats, "error")) node_stats$message else opts$node_stats)
}

if (file.exists(opts$anchor_audit)) {
  anchor_audit <- tryCatch(
    read.delim(opts$anchor_audit, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) e
  )
  add_check("anchor_audit_readable", !inherits(anchor_audit, "error"), if (inherits(anchor_audit, "error")) anchor_audit$message else opts$anchor_audit)
}

sample_ids <- character()
if (is.data.frame(counts_df)) {
  add_check("counts_has_gene_identifier_column", ncol(counts_df) >= 2L, paste("columns", ncol(counts_df)))
  genes <- normalise_id(counts_df[[1L]])
  add_check("gene_identifiers_present", length(genes) == nrow(counts_df), paste("genes", length(genes), "rows", nrow(counts_df)))
  add_check("gene_identifiers_unique", !anyDuplicated(genes), paste("duplicates", sum(duplicated(genes))))
  sample_ids <- colnames(counts_df)[-1L]
  numeric_counts <- suppressWarnings(as.matrix(data.frame(lapply(counts_df[-1L], as.numeric), check.names = FALSE)))
  na_count <- sum(is.na(numeric_counts))
  finite_count <- sum(is.finite(numeric_counts))
  total_count <- length(numeric_counts)
  negative_count <- sum(numeric_counts < 0, na.rm = TRUE)
  fractional_count <- sum(abs(numeric_counts - round(numeric_counts)) > 1e-8, na.rm = TRUE)
  rng <- range(numeric_counts, na.rm = TRUE)
  if (!all(is.finite(rng))) {
    rng <- c(NA_real_, NA_real_)
  }
  add_check("counts_orientation_genes_x_samples", length(sample_ids) > 0L && nrow(counts_df) > length(sample_ids),
            paste("genes", nrow(counts_df), "samples", length(sample_ids)))
  if (!is.na(expected_samples)) {
    add_check("counts_expected_sample_count", length(sample_ids) == expected_samples,
              paste("observed", length(sample_ids), "expected", expected_samples))
  }
  add_check("counts_no_missing_values", na_count == 0L, paste("NA_values", na_count))
  add_check("counts_no_nonfinite_values", finite_count == total_count, paste("finite", finite_count, "total", total_count))
  add_check("counts_no_negative_values", negative_count == 0L, paste("negative_values", negative_count))
  add_check("counts_integer_like", fractional_count == 0L, paste("fractional_values", fractional_count))
  add_check("counts_value_range", TRUE, paste(rng, collapse = " to "))
}

if (is.data.frame(metadata)) {
  metadata_cols <- colnames(metadata)
  add_check("metadata_has_sample_id_col", sample_id_col %in% metadata_cols, sample_id_col)
  add_check("metadata_has_cell_line_col", cell_line_col %in% metadata_cols, cell_line_col)
  add_check("metadata_has_component_col", component_col %in% metadata_cols, component_col)
  add_check("metadata_has_is_isolate_col", "is_isolate" %in% metadata_cols, "is_isolate")
  if (sample_id_col %in% metadata_cols && length(sample_ids)) {
    meta_samples <- normalise_id(metadata[[sample_id_col]])
    add_check("sample_ids_match_metadata", setequal(sample_ids, meta_samples),
              set_detail(sample_ids, meta_samples, "counts", "metadata"))
  }
}

if (is.data.frame(metadata) && is.data.frame(node_stats) &&
    all(c(cell_line_col, component_col, "is_isolate") %in% colnames(metadata)) &&
    all(c("cell_line", "component", "is_isolate") %in% colnames(node_stats))) {
  meta_map <- metadata[, c(cell_line_col, component_col, "is_isolate"), drop = FALSE]
  colnames(meta_map) <- c("cell_line", "component", "is_isolate")
  meta_map$cell_line <- as.character(meta_map$cell_line)
  meta_map$component <- as.character(meta_map$component)
  meta_map$is_isolate_norm <- truthy(meta_map$is_isolate)
  node_map <- node_stats[, c("cell_line", "component", "is_isolate"), drop = FALSE]
  node_map$cell_line <- as.character(node_map$cell_line)
  node_map$component <- as.character(node_map$component)
  node_map$is_isolate_norm <- truthy(node_map$is_isolate)
  merged <- merge(meta_map, node_map, by = "cell_line", all.x = TRUE, suffixes = c("_metadata", "_node"))
  unmapped <- merged$cell_line[is.na(merged$component_node)]
  component_mismatch <- merged$cell_line[!is.na(merged$component_node) & merged$component_metadata != merged$component_node]
  isolate_mismatch <- merged$cell_line[!is.na(merged$is_isolate_norm_node) & merged$is_isolate_norm_metadata != merged$is_isolate_norm_node]
  add_check("metadata_cell_lines_map_to_current_graph", length(unmapped) == 0L,
            if (length(unmapped)) paste(unique(unmapped), collapse = ",") else "all metadata cell lines found")
  add_check("metadata_components_match_current_graph", length(component_mismatch) == 0L,
            if (length(component_mismatch)) paste(unique(component_mismatch), collapse = ",") else "metadata components current")
  add_check("metadata_isolate_flags_match_current_graph", length(isolate_mismatch) == 0L,
            if (length(isolate_mismatch)) paste(unique(isolate_mismatch), collapse = ",") else "metadata isolate flags current")
}

if (file.exists(opts$isolate_list) && is.data.frame(node_stats) &&
    all(c("cell_line", "is_isolate") %in% colnames(node_stats))) {
  staged_isolates <- sort(unique(read_csv_list(opts$isolate_list)))
  current_isolates <- sort(unique(normalise_id(node_stats$cell_line[truthy(node_stats$is_isolate)])))
  add_check("isolate_list_matches_current_graph", setequal(staged_isolates, current_isolates),
            set_detail(staged_isolates, current_isolates, "staged", "current"))
}

if (file.exists(opts$components_list) && is.data.frame(node_stats) && component_col %in% colnames(node_stats)) {
  staged_components <- sort(unique(read_text_list(opts$components_list)))
  current_components_raw <- suppressWarnings(as.integer(as.character(node_stats[[component_col]])))
  tab <- table(current_components_raw[!is.na(current_components_raw) & current_components_raw >= 0L])
  current_components <- sort(names(tab[tab >= 2L]))
  add_check("components_list_matches_current_graph", setequal(staged_components, current_components),
            set_detail(staged_components, current_components, "staged", "current"))
}

if (file.exists(opts$anchor_list) && file.exists(opts$anchor_components) && is.data.frame(anchor_audit)) {
  staged_anchors <- sort(unique(read_csv_list(opts$anchor_list)))
  staged_anchor_comp <- tryCatch(
    read.delim(opts$anchor_components, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) e
  )
  add_check("anchor_components_readable", !inherits(staged_anchor_comp, "error"),
            if (inherits(staged_anchor_comp, "error")) staged_anchor_comp$message else opts$anchor_components)
  selected <- rep(FALSE, nrow(anchor_audit))
  if ("anchor_selected" %in% colnames(anchor_audit)) {
    selected <- truthy(anchor_audit$anchor_selected)
  } else {
    if ("most_connected_selected" %in% colnames(anchor_audit)) {
      selected <- selected | truthy(anchor_audit$most_connected_selected)
    }
    if ("canonical_bridge_selected" %in% colnames(anchor_audit)) {
      selected <- selected | truthy(anchor_audit$canonical_bridge_selected)
    }
  }
  if (all(c("node_id", "component_id") %in% colnames(anchor_audit))) {
    current_anchor_comp <- unique(data.frame(
      anchor = as.character(anchor_audit$node_id[selected]),
      component = as.character(anchor_audit$component_id[selected]),
      stringsAsFactors = FALSE
    ))
    current_anchors <- sort(unique(normalise_id(current_anchor_comp$anchor)))
    add_check("anchor_list_matches_current_graph", setequal(staged_anchors, current_anchors),
              set_detail(staged_anchors, current_anchors, "staged", "current"))
    if (is.data.frame(staged_anchor_comp) && all(c("anchor", "component") %in% colnames(staged_anchor_comp))) {
      staged_pairs <- sort(unique(paste(as.character(staged_anchor_comp$anchor), as.character(staged_anchor_comp$component), sep = "::")))
      current_pairs <- sort(unique(paste(current_anchor_comp$anchor, current_anchor_comp$component, sep = "::")))
      add_check("anchor_components_match_current_graph", setequal(staged_pairs, current_pairs),
                set_detail(staged_pairs, current_pairs, "staged", "current"))
    } else if (is.data.frame(staged_anchor_comp)) {
      add_check("anchor_components_has_required_columns", FALSE, paste(colnames(staged_anchor_comp), collapse = ","))
    }
  } else {
    add_check("anchor_audit_has_required_columns", FALSE, paste(colnames(anchor_audit), collapse = ","))
  }
}

dir.create(dirname(opts$output), showWarnings = FALSE, recursive = TRUE)
write.table(checks, opts$output, sep = "\t", quote = FALSE, row.names = FALSE)

failures <- checks[checks$status != "PASS", , drop = FALSE]
if (nrow(failures)) {
  message("Staged DESeq2 validation completed with ", nrow(failures), " failed check(s): ", opts$output)
} else {
  message("Staged DESeq2 validation passed: ", opts$output)
}
