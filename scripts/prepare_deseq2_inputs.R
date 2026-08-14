#!/usr/bin/env Rscript

# Prepare graph-defined, cancer-type-specific DESeq2 inputs.
# Unit of analysis: source count profile mapped to one resolved-graph biological
# cell line. Transformation: source counts + source metadata + analytical graph
# node statistics -> prepared count matrix, component-annotated metadata, sample
# mapping, and raw-count provenance.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--profile", type = "character"),
  make_option("--dsmz_counts", type = "character"),
  make_option("--dsmz_meta", type = "character"),
  make_option("--node_stats", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--source_count_kind", type = "character", default = "raw_gene_level_counts"),
  make_option("--source_sample_id_col", type = "character", default = "sample_id"),
  make_option("--source_cell_line_col", type = "character", default = "cell_line"),
  make_option("--sample_population", type = "character", default = "cell_line_only"),
  make_option("--cell_line_verification_mode", type = "character", default = "source_population_declaration"),
  make_option("--cell_line_verification_column", type = "character", default = ""),
  make_option("--accepted_cell_line_values", type = "character", default = "Cell Line"),
  make_option("--declared_population", type = "character", default = "cell_line_only")
)

opt <- parse_args(OptionParser(option_list = option_list))

required_options <- c("profile", "dsmz_counts", "dsmz_meta", "node_stats", "outdir")
for (option_name in required_options) {
  value <- opt[[option_name]]
  if (is.null(value) || is.na(value) || !nzchar(trimws(value))) {
    stop(sprintf("[DESeq2 preparation] --%s is required", option_name), call. = FALSE)
  }
}
if (!identical(opt$source_count_kind, "raw_gene_level_counts")) {
  stop("[DESeq2 preparation] source_count_kind must be raw_gene_level_counts for DESeq2 Wald inference",
       call. = FALSE)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

canonicalise_cell_line_identifier <- function(value) {
  value <- toupper(trimws(as.character(value)))
  gsub("[^A-Z0-9]+", "", value)
}

require_columns <- function(table, columns, table_name) {
  missing <- setdiff(columns, names(table))
  if (length(missing) > 0L) {
    stop(sprintf(
      "[DESeq2 preparation] %s is missing required column(s): %s",
      table_name,
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
}

read_count_source <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    object <- readRDS(path)
    if (is.matrix(object)) {
      table <- as.data.frame(object, check.names = FALSE)
      table <- data.frame(gene_id = rownames(object), table, check.names = FALSE)
    } else if (is.data.frame(object)) {
      table <- as.data.frame(object, check.names = FALSE)
      if (!("gene_id" %in% names(table)) && !is.null(rownames(table))) {
        table <- data.frame(gene_id = rownames(table), table, check.names = FALSE)
      }
    } else {
      stop("[DESeq2 preparation] Count RDS must contain a matrix or data.frame", call. = FALSE)
    }
  } else {
    table <- fread(path, sep = "\t", data.table = FALSE, check.names = FALSE)
  }
  if (!("gene_id" %in% names(table))) {
    names(table)[1L] <- "gene_id"
  }
  table
}

sha256_file <- function(path) {
  commands <- list(c("sha256sum", path), c("shasum", "-a", "256", path))
  for (command in commands) {
    result <- tryCatch(
      system2(command[[1L]], command[-1L], stdout = TRUE, stderr = FALSE),
      error = function(error) character(0)
    )
    if (length(result) > 0L && nzchar(result[[1L]])) {
      return(strsplit(result[[1L]], "[[:space:]]+")[[1L]][[1L]])
    }
  }
  ""
}

counts_source <- read_count_source(opt$dsmz_counts)
metadata_source <- fread(opt$dsmz_meta, sep = "\t", data.table = FALSE, check.names = FALSE)
node_stats <- fread(opt$node_stats, sep = "\t", data.table = FALSE, check.names = FALSE)

require_columns(counts_source, "gene_id", "source count matrix")
require_columns(metadata_source, c(opt$source_sample_id_col, opt$source_cell_line_col), "source metadata")
require_columns(node_stats, c("cell_line", "component", "is_isolate"), "resolved graph node statistics")

gene_ids <- as.character(counts_source$gene_id)
if (any(is.na(gene_ids)) || any(!nzchar(trimws(gene_ids)))) {
  stop("[DESeq2 preparation] Source count matrix contains missing or empty gene identifiers", call. = FALSE)
}
if (anyDuplicated(gene_ids)) {
  stop("[DESeq2 preparation] Source count matrix contains duplicate gene identifiers", call. = FALSE)
}
count_sample_ids <- names(counts_source)[names(counts_source) != "gene_id"]
if (length(count_sample_ids) == 0L) {
  stop("[DESeq2 preparation] Source count matrix contains no sample columns", call. = FALSE)
}
if (any(is.na(count_sample_ids)) || any(!nzchar(trimws(count_sample_ids)))) {
  stop("[DESeq2 preparation] Count sample identifiers are missing or empty", call. = FALSE)
}
if (anyDuplicated(count_sample_ids)) {
  stop("[DESeq2 preparation] Count sample identifiers are not unique", call. = FALSE)
}

count_values <- as.matrix(counts_source[, count_sample_ids, drop = FALSE])
suppressWarnings(storage.mode(count_values) <- "numeric")
if (any(is.na(count_values))) {
  stop("[DESeq2 preparation] Source count matrix contains NA values", call. = FALSE)
}
if (any(!is.finite(count_values))) {
  stop("[DESeq2 preparation] Source count matrix contains non-finite values", call. = FALSE)
}
if (any(count_values < 0)) {
  stop("[DESeq2 preparation] Source count matrix contains negative values", call. = FALSE)
}
if (any(abs(count_values - round(count_values)) > 1e-8)) {
  stop("[DESeq2 preparation] Source count matrix contains non-integer values", call. = FALSE)
}

metadata_source[[opt$source_sample_id_col]] <- trimws(as.character(metadata_source[[opt$source_sample_id_col]]))
metadata_source[[opt$source_cell_line_col]] <- trimws(as.character(metadata_source[[opt$source_cell_line_col]]))
if (any(is.na(metadata_source[[opt$source_sample_id_col]])) ||
    any(!nzchar(metadata_source[[opt$source_sample_id_col]]))) {
  stop("[DESeq2 preparation] Source metadata sample identifiers are missing or empty", call. = FALSE)
}
if (any(is.na(metadata_source[[opt$source_cell_line_col]])) ||
    any(!nzchar(metadata_source[[opt$source_cell_line_col]]))) {
  stop("[DESeq2 preparation] Source metadata biological cell-line identifiers are missing or empty", call. = FALSE)
}
count_bearing_metadata <- metadata_source[metadata_source[[opt$source_sample_id_col]] %in% count_sample_ids, , drop = FALSE]
if (anyDuplicated(count_bearing_metadata[[opt$source_sample_id_col]])) {
  duplicated_ids <- unique(count_bearing_metadata[[opt$source_sample_id_col]][duplicated(count_bearing_metadata[[opt$source_sample_id_col]])])
  stop("[DESeq2 preparation] Source metadata has duplicate count-bearing sample identifiers: ",
       paste(duplicated_ids, collapse = ", "), call. = FALSE)
}

graph_cell_lines <- sort(unique(trimws(as.character(node_stats$cell_line))))
graph_cell_lines <- graph_cell_lines[nzchar(graph_cell_lines)]
if (length(graph_cell_lines) == 0L) {
  stop("[DESeq2 preparation] Resolved graph node statistics contain no cell-line identifiers", call. = FALSE)
}
graph_canonical <- canonicalise_cell_line_identifier(graph_cell_lines)
if (anyDuplicated(graph_canonical)) {
  collision_table <- split(graph_cell_lines, graph_canonical)
  collisions <- collision_table[vapply(collision_table, function(ids) length(unique(ids)) > 1L, logical(1))]
  stop("[DESeq2 preparation] Canonical graph cell-line identifiers collide: ",
       paste(vapply(collisions, paste, character(1), collapse = "|"), collapse = "; "),
       call. = FALSE)
}
graph_id_by_canonical <- setNames(graph_cell_lines, graph_canonical)

metadata_source$source_sample_id_for_matching <- metadata_source[[opt$source_sample_id_col]]
metadata_source$source_cell_line_for_matching <- metadata_source[[opt$source_cell_line_col]]
metadata_source$source_cell_line_canonical <- canonicalise_cell_line_identifier(
  metadata_source$source_cell_line_for_matching
)

mapping_rows <- list()
for (target_cell_line in graph_cell_lines) {
  target_canonical <- canonicalise_cell_line_identifier(target_cell_line)
  candidate_samples <- character(0)
  candidate_sources <- character(0)

  exact_count_matches <- count_sample_ids[count_sample_ids == target_cell_line]
  if (length(exact_count_matches) > 0L) {
    candidate_samples <- c(candidate_samples, exact_count_matches)
    candidate_sources <- c(candidate_sources, rep("exact_count_column", length(exact_count_matches)))
  }

  exact_metadata_matches <- metadata_source$source_sample_id_for_matching[
    metadata_source$source_sample_id_for_matching == target_cell_line &
      metadata_source$source_sample_id_for_matching %in% count_sample_ids
  ]
  if (length(exact_metadata_matches) > 0L) {
    candidate_samples <- c(candidate_samples, exact_metadata_matches)
    candidate_sources <- c(candidate_sources, rep("exact_metadata_sample", length(exact_metadata_matches)))
  }

  metadata_cell_line_matches <- metadata_source$source_sample_id_for_matching[
    metadata_source$source_cell_line_canonical == target_canonical &
      metadata_source$source_sample_id_for_matching %in% count_sample_ids
  ]
  if (length(metadata_cell_line_matches) > 0L) {
    candidate_samples <- c(candidate_samples, metadata_cell_line_matches)
    candidate_sources <- c(candidate_sources, rep("metadata_cell_line_identifier", length(metadata_cell_line_matches)))
  }

  if (length(candidate_samples) > 0L) {
    candidate_table <- unique(data.frame(
      count_sample_id = candidate_samples,
      biological_cell_line = target_cell_line,
      mapping_evidence = candidate_sources,
      stringsAsFactors = FALSE
    ))
    mapping_rows[[target_cell_line]] <- candidate_table
  }
}

if (length(mapping_rows) == 0L) {
  stop("[DESeq2 preparation] No resolved graph cell lines mapped to count-bearing source samples", call. = FALSE)
}
sample_mapping <- unique(rbindlist(mapping_rows, use.names = TRUE))
same_pair <- unique(sample_mapping[, c("count_sample_id", "biological_cell_line")])
conflicts <- split(same_pair$biological_cell_line, same_pair$count_sample_id)
conflicts <- conflicts[vapply(conflicts, function(values) length(unique(values)) > 1L, logical(1))]
if (length(conflicts) > 0L) {
  details <- paste(
    sprintf("%s -> %s", names(conflicts), vapply(conflicts, paste, character(1), collapse = ",")),
    collapse = "; "
  )
  stop("[DESeq2 preparation] Count sample maps to multiple biological cell-line targets: ",
       details, call. = FALSE)
}
sample_mapping <- same_pair[order(same_pair$biological_cell_line, same_pair$count_sample_id), , drop = FALSE]
names(sample_mapping) <- c("sample_id", "cell_line")

unmatched_graph_cell_lines <- setdiff(graph_cell_lines, unique(sample_mapping$cell_line))
if (length(unmatched_graph_cell_lines) > 0L) {
  stop("[DESeq2 preparation] Resolved graph cell lines lacking count-bearing samples: ",
       paste(unmatched_graph_cell_lines, collapse = ", "), call. = FALSE)
}

metadata_join <- metadata_source[match(sample_mapping$sample_id, metadata_source[[opt$source_sample_id_col]]), , drop = FALSE]
if (nrow(metadata_join) != nrow(sample_mapping) ||
    any(is.na(metadata_join[[opt$source_sample_id_col]]))) {
  stop("[DESeq2 preparation] Metadata join failed to return exactly one row per prepared sample", call. = FALSE)
}
if (anyDuplicated(sample_mapping$sample_id)) {
  stop("[DESeq2 preparation] Prepared sample identifiers are duplicated after metadata join", call. = FALSE)
}

verification_mode <- opt$cell_line_verification_mode
accepted_values <- trimws(unlist(strsplit(opt$accepted_cell_line_values, ",")))
accepted_values <- accepted_values[nzchar(accepted_values)]
if (verification_mode == "metadata_column") {
  verification_column <- opt$cell_line_verification_column
  require_columns(metadata_join, verification_column, "source metadata for cell-line verification")
  observed_values <- trimws(as.character(metadata_join[[verification_column]]))
  rejected <- observed_values[!(observed_values %in% accepted_values)]
  if (length(rejected) > 0L) {
    stop("[DESeq2 preparation] Non-cell-line samples selected by metadata_column verification: ",
         paste(unique(rejected), collapse = ", "), call. = FALSE)
  }
  sample_population <- "cell_line_only"
} else if (verification_mode == "source_population_declaration") {
  if (opt$declared_population != "cell_line_only") {
    stop("[DESeq2 preparation] Source population declaration is not cell_line_only", call. = FALSE)
  }
  sample_population <- "cell_line_only"
} else {
  stop("[DESeq2 preparation] Unknown cell-line verification mode: ", verification_mode, call. = FALSE)
}

node_map <- unique(node_stats[, c("cell_line", "component", "is_isolate", "degree", "betweenness"), drop = FALSE])
if (anyDuplicated(node_map$cell_line)) {
  stop("[DESeq2 preparation] Graph node statistics have duplicate cell-line rows", call. = FALSE)
}
source_graph_annotation_columns <- c("component", "is_isolate", "degree", "betweenness")
metadata_join_extra <- metadata_join[
  ,
  setdiff(names(metadata_join), c(names(sample_mapping), source_graph_annotation_columns)),
  drop = FALSE
]
prepared_metadata <- cbind(sample_mapping, metadata_join_extra, stringsAsFactors = FALSE)
prepared_metadata <- merge(
  prepared_metadata,
  node_map,
  by = "cell_line",
  all.x = TRUE,
  sort = FALSE
)
prepared_metadata <- prepared_metadata[match(sample_mapping$sample_id, prepared_metadata$sample_id), , drop = FALSE]
if (any(is.na(prepared_metadata$component)) || any(is.na(prepared_metadata$is_isolate))) {
  missing_components <- prepared_metadata$cell_line[is.na(prepared_metadata$component) | is.na(prepared_metadata$is_isolate)]
  stop("[DESeq2 preparation] Component/isolate graph annotations missing after join: ",
       paste(unique(missing_components), collapse = ", "), call. = FALSE)
}

prepared_cell_lines <- sort(unique(prepared_metadata$cell_line))
if (!setequal(prepared_cell_lines, graph_cell_lines)) {
  stop("[DESeq2 preparation] Prepared biological cell-line set does not exactly equal resolved graph cell-line set; missing_from_prepared=",
       paste(setdiff(graph_cell_lines, prepared_cell_lines), collapse = ","),
       "; extra_in_prepared=", paste(setdiff(prepared_cell_lines, graph_cell_lines), collapse = ","),
       call. = FALSE)
}

prepared_counts <- counts_source[, c("gene_id", sample_mapping$sample_id), drop = FALSE]
source_slice <- counts_source[, c("gene_id", sample_mapping$sample_id), drop = FALSE]
if (!identical(prepared_counts, source_slice)) {
  stop("[DESeq2 preparation] Prepared count values do not exactly equal selected source count values", call. = FALSE)
}

fwrite(prepared_counts, file.path(opt$outdir, "counts.tsv"), sep = "\t")
fwrite(prepared_metadata, file.path(opt$outdir, "metadata.tsv"), sep = "\t")
fwrite(prepared_metadata, file.path(opt$outdir, "metadata_with_components.tsv"), sep = "\t")
fwrite(sample_mapping, file.path(opt$outdir, "sample_mapping.tsv"), sep = "\t")

provenance <- data.frame(
  field = c(
    "source_counts_path",
    "source_metadata_path",
    "source_count_kind",
    "sample_population",
    "cell_line_verification_mode",
    "cell_line_verification_column",
    "accepted_cell_line_values",
    "value_transformation",
    "normalisation_applied",
    "aggregation_applied",
    "value_preservation_check",
    "source_counts_sha256",
    "source_metadata_sha256"
  ),
  value = c(
    opt$dsmz_counts,
    opt$dsmz_meta,
    opt$source_count_kind,
    sample_population,
    verification_mode,
    opt$cell_line_verification_column,
    paste(accepted_values, collapse = ";"),
    "subset_and_reorder_only",
    "FALSE",
    "FALSE",
    "prepared_values_identical_to_selected_source_values",
    sha256_file(opt$dsmz_counts),
    sha256_file(opt$dsmz_meta)
  ),
  stringsAsFactors = FALSE
)
fwrite(provenance, file.path(opt$outdir, "input_provenance.tsv"), sep = "\t")

message(sprintf(
  "[DESeq2 preparation] PASS profile=%s samples=%d biological_cell_lines=%d genes=%d",
  opt$profile,
  nrow(sample_mapping),
  length(prepared_cell_lines),
  nrow(prepared_counts)
))
