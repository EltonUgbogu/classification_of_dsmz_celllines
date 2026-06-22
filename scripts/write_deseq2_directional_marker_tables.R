#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1L
    } else {
      out[[name]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

require_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required argument --", name, call. = FALSE)
  }
  value
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path, call. = FALSE)
  }
}

path_under <- function(root, maybe_relative) {
  if (is.na(maybe_relative) || !nzchar(maybe_relative)) {
    stop("Encountered an empty manifest path", call. = FALSE)
  }
  if (grepl("^/", maybe_relative)) {
    maybe_relative
  } else {
    file.path(root, maybe_relative)
  }
}

confirm_orientation <- function(deseq2_script) {
  stop_if_missing(deseq2_script, "DESeq2 source script")
  source_text <- paste(readLines(deseq2_script, warn = FALSE), collapse = "\n")

  isolate_pattern <- paste0(
    "results\\s*\\(\\s*dds_tmp\\s*,\\s*",
    "contrast\\s*=\\s*c\\s*\\(\\s*\"grp_tmp\"\\s*,\\s*iso\\s*,\\s*\"REST\"\\s*\\)"
  )
  anchor_pattern <- paste0(
    "results\\s*\\(\\s*dds_k\\s*,\\s*",
    "contrast\\s*=\\s*c\\s*\\(\\s*\"grp_tmp\"\\s*,\\s*anc\\s*,\\s*\"OUTSIDE_COMP\"\\s*\\)"
  )
  isolate_comment <- "Positive LFC indicates higher expression in the isolate"

  ok <- grepl(isolate_pattern, source_text, perl = TRUE) &&
    grepl(anchor_pattern, source_text, perl = TRUE) &&
    grepl(isolate_comment, source_text, fixed = TRUE)

  if (!ok) {
    stop(
      "Could not confirm DESeq2 contrast orientation in ", deseq2_script,
      ". Directional marker tables were not written.",
      call. = FALSE
    )
  }
  TRUE
}

parse_contrast <- function(contrast_id) {
  if (grepl("^isolate_.+_vs_rest$", contrast_id)) {
    focal <- sub("^isolate_", "", contrast_id)
    focal <- sub("_vs_rest$", "", focal)
    return(list(
      focal_cell_line = focal,
      focal_type = "isolate",
      reference_group = "rest",
      component = NA_character_
    ))
  }

  anchor_match <- regexec("^anchor_(.+)_vs_outside_component_(.+)$", contrast_id)
  parts <- regmatches(contrast_id, anchor_match)[[1L]]
  if (length(parts) == 3L) {
    return(list(
      focal_cell_line = parts[[2L]],
      focal_type = "anchor",
      reference_group = "outside_component",
      component = parts[[3L]]
    ))
  }

  stop("Could not parse contrast_id: ", contrast_id, call. = FALSE)
}

build_marker_table <- function(results_table, marker_genes, marker_rank, metadata) {
  if (!("gene_id" %in% names(results_table))) {
    stop("DESeq2 table is missing required gene_id column", call. = FALSE)
  }
  if (!("log2FoldChange" %in% names(results_table))) {
    stop("DESeq2 table is missing required log2FoldChange column", call. = FALSE)
  }

  missing_genes <- setdiff(marker_genes, results_table$gene_id)
  if (length(missing_genes) > 0L) {
    stop(
      "Marker genes were not found in the DESeq2 table: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  marker_dt <- data.table(gene_id = marker_genes, marker_rank = marker_rank)
  retained <- merge(marker_dt, results_table, by = "gene_id", all.x = TRUE, sort = FALSE)
  setorder(retained, marker_rank)

  retained[, `:=`(
    cohort = metadata$cohort,
    contrast_id = metadata$contrast_id,
    focal_cell_line = metadata$focal_cell_line,
    focal_type = metadata$focal_type,
    reference_group = metadata$reference_group,
    component = metadata$component
  )]

  meta_cols <- c(
    "cohort", "contrast_id", "focal_cell_line", "focal_type",
    "reference_group", "component", "direction", "marker_rank"
  )
  setcolorder(retained, c(intersect(meta_cols, names(retained)),
                          setdiff(names(retained), meta_cols)))
  retained
}

sort_top_markers <- function(dt) {
  if (nrow(dt) == 0L) {
    return(dt)
  }

  gene_col <- if ("gene_id" %in% names(dt)) "gene_id" else names(dt)[[1L]]
  sorted <- copy(dt)
  if ("padj" %in% names(sorted)) {
    sorted[, padj_sort := suppressWarnings(as.numeric(padj))]
  } else {
    sorted[, padj_sort := Inf]
  }
  sorted[is.na(padj_sort), padj_sort := Inf]
  sorted[, abs_lfc_sort := abs(suppressWarnings(as.numeric(log2FoldChange)))]
  sorted[is.na(abs_lfc_sort), abs_lfc_sort := -Inf]
  sorted[, gene_sort := as.character(get(gene_col))]
  sorted <- sorted[order(padj_sort, -abs_lfc_sort, marker_rank, gene_sort)]
  sorted[, c("padj_sort", "abs_lfc_sort", "gene_sort") := NULL]
  sorted
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

cohort <- require_arg(args, "cohort")
deseq2_dir <- normalizePath(require_arg(args, "deseq2_dir"), mustWork = TRUE)
deseq2_script <- normalizePath(require_arg(args, "deseq2_script"), mustWork = TRUE)
outdir <- args[["outdir"]]
if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- file.path(deseq2_dir, "directional_markers")
}

confirm_orientation(deseq2_script)

manifest_path <- file.path(deseq2_dir, "markers", "marker_sets_manifest.tsv")
stop_if_missing(manifest_path, "Marker manifest")

manifest <- fread(manifest_path, sep = "\t", data.table = TRUE)
required_manifest_cols <- c("contrast", "marker_file", "table_file", "n_markers")
missing_manifest_cols <- setdiff(required_manifest_cols, names(manifest))
if (length(missing_manifest_cols) > 0L) {
  stop(
    "Marker manifest is missing required columns: ",
    paste(missing_manifest_cols, collapse = ", "),
    call. = FALSE
  )
}

contrast_dir <- file.path(outdir, "contrast_level")
cohort_dir <- file.path(outdir, "cohort_level")
dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)

up_tables <- list()
down_tables <- list()
count_rows <- list()
zero_lfc_tables <- list()

for (i in seq_len(nrow(manifest))) {
  contrast_id <- manifest$contrast[[i]]
  marker_path <- path_under(deseq2_dir, manifest$marker_file[[i]])
  table_path <- path_under(deseq2_dir, manifest$table_file[[i]])
  stop_if_missing(marker_path, paste0("Marker file for ", contrast_id))
  stop_if_missing(table_path, paste0("DESeq2 table for ", contrast_id))

  marker_genes <- readLines(marker_path, warn = FALSE)
  marker_genes <- marker_genes[nzchar(marker_genes)]
  marker_rank <- seq_along(marker_genes)

  contrast_meta <- parse_contrast(contrast_id)
  contrast_meta$cohort <- cohort
  contrast_meta$contrast_id <- contrast_id

  results_table <- fread(table_path, sep = "\t", data.table = TRUE)
  retained <- build_marker_table(results_table, marker_genes, marker_rank, contrast_meta)

  up <- copy(retained[!is.na(log2FoldChange) & log2FoldChange > 0])
  down <- copy(retained[!is.na(log2FoldChange) & log2FoldChange < 0])
  zero <- copy(retained[!is.na(log2FoldChange) & log2FoldChange == 0])

  up[, direction := "upregulated"]
  down[, direction := "downregulated"]
  if (nrow(zero) > 0L) {
    zero[, direction := "non_directional"]
    zero_lfc_tables[[contrast_id]] <- zero
  }

  meta_cols <- c(
    "cohort", "contrast_id", "focal_cell_line", "focal_type",
    "reference_group", "component", "direction", "marker_rank"
  )
  setcolorder(up, c(intersect(meta_cols, names(up)), setdiff(names(up), meta_cols)))
  setcolorder(down, c(intersect(meta_cols, names(down)), setdiff(names(down), meta_cols)))

  fwrite(
    up,
    file.path(contrast_dir, paste0(contrast_id, "_upregulated.tsv")),
    sep = "\t"
  )
  fwrite(
    down,
    file.path(contrast_dir, paste0(contrast_id, "_downregulated.tsv")),
    sep = "\t"
  )

  up_tables[[contrast_id]] <- up
  down_tables[[contrast_id]] <- down
  count_rows[[contrast_id]] <- data.table(
    cohort = cohort,
    contrast_id = contrast_id,
    focal_cell_line = contrast_meta$focal_cell_line,
    focal_type = contrast_meta$focal_type,
    reference_group = contrast_meta$reference_group,
    n_upregulated = nrow(up),
    n_downregulated = nrow(down),
    n_total_directional = nrow(up) + nrow(down)
  )
}

all_up <- rbindlist(up_tables, use.names = TRUE, fill = TRUE)
all_down <- rbindlist(down_tables, use.names = TRUE, fill = TRUE)
counts <- rbindlist(count_rows, use.names = TRUE, fill = TRUE)

fwrite(all_up, file.path(cohort_dir, "all_upregulated_markers.tsv"), sep = "\t")
fwrite(all_down, file.path(cohort_dir, "all_downregulated_markers.tsv"), sep = "\t")
fwrite(counts, file.path(cohort_dir, "directional_marker_counts.tsv"), sep = "\t")
fwrite(sort_top_markers(all_up), file.path(cohort_dir, "top_upregulated_markers.tsv"), sep = "\t")
fwrite(sort_top_markers(all_down), file.path(cohort_dir, "top_downregulated_markers.tsv"), sep = "\t")

if (length(zero_lfc_tables) > 0L) {
  zero_lfc <- rbindlist(zero_lfc_tables, use.names = TRUE, fill = TRUE)
  fwrite(zero_lfc, file.path(cohort_dir, "non_directional_markers.tsv"), sep = "\t")
}

orientation_text <- c(
  "Positive log2FoldChange denotes higher expression in the focal isolate/anchor relative to the reference group.",
  "Negative log2FoldChange denotes lower expression in the focal isolate/anchor relative to the reference group.",
  paste0("Orientation verified from: ", deseq2_script)
)
writeLines(orientation_text, file.path(outdir, "directional_marker_orientation.txt"))

message(sprintf(
  "[OK] Wrote directional marker tables for %s: %d contrasts",
  cohort,
  nrow(manifest)
))
