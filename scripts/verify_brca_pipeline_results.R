#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list(root = "results/unsupervised/brca")
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args)) stop("Missing value for argument: ", key)
    out[[gsub("-", "_", substring(key, 3), fixed = TRUE)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

read_table_safe <- function(path) {
  tryCatch(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) NULL)
}

check_file <- function(path, label, errors) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    errors <- c(errors, paste0(label, " missing or empty: ", path))
    cat("[FAIL] ", label, ": ", path, "\n", sep = "")
  } else {
    cat("[OK] ", label, ": ", path, "\n", sep = "")
  }
  errors
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
root <- args$root
tn_root <- file.path(root, "tumour_neighbourhoods")
all_root <- file.path(tn_root, "final_consensus_all")

directions <- c(
  "Variance_euc", "Variance_corr", "MAD_euc", "MAD_corr",
  "MeanAbsDev_euc", "MeanAbsDev_corr", "Entropy_euc", "Entropy_corr",
  "PCA_euc", "PCA_corr", "Spearman_euc", "Spearman_corr",
  "MX_euc", "MX_corr", "kTotal_euc", "kTotal_corr",
  "HVG_euc", "HVG_corr", "pam50_euc", "pam50_corr"
)

cat("=== BRCA pipeline result verification ===\n")
errors <- character()

for (direction in directions) {
  cons <- file.path(tn_root, direction, "final_consensus",
                    paste0("Final_consensus_tumour_neighbourhoods_", direction, ".tsv"))
  edges <- file.path(tn_root, direction, "final_consensus",
                     paste0("cell_line_similarity_graph_edges_", direction, ".tsv"))
  errors <- check_file(cons, paste0(direction, " consensus"), errors)
  errors <- check_file(edges, paste0(direction, " similarity edges"), errors)

  tbl <- read_table_safe(cons)
  if (!is.null(tbl)) {
    required <- c("tumour_id", "p_consensus")
    missing <- setdiff(required, colnames(tbl))
    if (length(missing)) {
      errors <- c(errors, paste0(direction, " consensus missing columns: ", paste(missing, collapse = ", ")))
    }
    cat(sprintf("[INFO] %s consensus rows: %d\n", direction, nrow(tbl)))
  }
}

final_files <- c(
  "ranked cell lines" = file.path(all_root, "p_consensus_best_cell_lines_ranked.tsv"),
  "resolved neighbours" = file.path(all_root, "resolved_dsmz_neighbours.tsv"),
  "node stats" = file.path(all_root, "patient_referenced_aggregated_cell_line_similarity_graph_node_stats.tsv"),
  "model selection" = file.path(root, "validation", "model_selection_summary.tsv"),
  "permutation validation" = file.path(root, "validation", "neighbourhood_permutation_summary.tsv"),
  "random baseline" = file.path(root, "validation", "random_baseline_summary.tsv"),
  "silhouette report" = file.path(root, "validation", "silhouette_report.tsv"),
  "DESeq2 component done" = file.path(root, "deseq2", "component_vs_rest", ".done")
)

for (nm in names(final_files)) {
  errors <- check_file(final_files[[nm]], nm, errors)
}

ranked <- read_table_safe(final_files[["ranked cell lines"]])
if (!is.null(ranked)) {
  cat(sprintf("[INFO] Ranked table rows: %d\n", nrow(ranked)))
  expected_any <- c("cell_line", "composite_score", "frac_ge_thr", "median_p", "max_p")
  missing <- setdiff(expected_any, colnames(ranked))
  if (length(missing)) {
    cat("[WARN] Ranked table missing expected columns: ", paste(missing, collapse = ", "), "\n", sep = "")
  }
}

cat("\n=== Summary ===\n")
cat(sprintf("Errors: %d\n", length(errors)))
if (length(errors)) {
  cat(paste0(" - ", errors, collapse = "\n"), "\n")
  quit(status = 1)
}
quit(status = 0)
