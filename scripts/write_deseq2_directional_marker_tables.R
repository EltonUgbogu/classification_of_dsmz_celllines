#!/usr/bin/env Rscript

# Collate retained canonical DESeq2 markers into direction-specific enrichment
# tables. Unit of recurrence here is retained contrasts within this enrichment
# grouping; feature-selection recurrence is handled separately by the pan-cancer
# feature module.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--cohort", type = "character"),
  make_option("--deseq2_dir", type = "character"),
  make_option("--deseq2_script", type = "character", default = ""),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

for (name in c("cohort", "deseq2_dir", "outdir")) {
  value <- opt[[name]]
  if (is.null(value) || is.na(value) || !nzchar(trimws(value))) {
    stop(sprintf("[Directional markers] --%s is required", name), call. = FALSE)
  }
}

manifest_path <- file.path(opt$deseq2_dir, "markers", "contrast_level_marker_manifest.tsv")
if (!file.exists(manifest_path)) {
  stop("[Directional markers] Missing canonical manifest: ", manifest_path, call. = FALSE)
}

manifest <- fread(manifest_path, sep = "\t", data.table = FALSE, check.names = FALSE)
required_manifest <- c(
  "contrast_id", "marker_evidence_stratum", "marker_table_path",
  "marker_gene_list_path", "n_markers_after_cap"
)
missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest) > 0L) {
  stop("[Directional markers] Canonical manifest missing column(s): ",
       paste(missing_manifest, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(manifest$contrast_id)) {
  stop("[Directional markers] Duplicate contrast_id values in canonical manifest", call. = FALSE)
}

required_marker_columns <- c(
  "gene_id", "baseMean", "wald_statistic", "p_value", "adjusted_p_value",
  "log2_fold_change_unshrunk", "log2_fold_change_shrunken",
  "log2_fold_change_posterior_sd", "absolute_shrunken_log2_fold_change",
  "effect_direction", "contrast_marker_rank"
)

resolve_path <- function(value) {
  path <- as.character(value)
  if (grepl("^/", path)) return(path)
  file.path(opt$deseq2_dir, path)
}

load_gene_list <- function(path) {
  if (!file.exists(path)) stop("[Directional markers] Missing marker gene list: ", path, call. = FALSE)
  values <- trimws(readLines(path, warn = FALSE))
  values <- values[nzchar(values)]
  sub("\\..*$", "", values)
}

records <- list()
for (row_index in seq_len(nrow(manifest))) {
  row <- manifest[row_index, , drop = FALSE]
  contrast_id <- as.character(row$contrast_id)
  marker_table_path <- resolve_path(row$marker_table_path)
  marker_gene_list_path <- resolve_path(row$marker_gene_list_path)
  if (!file.exists(marker_table_path)) {
    stop("[Directional markers] Missing retained marker table for ", contrast_id, ": ",
         marker_table_path, call. = FALSE)
  }
  marker_table <- fread(marker_table_path, sep = "\t", data.table = FALSE, check.names = FALSE)
  missing_marker <- setdiff(required_marker_columns, names(marker_table))
  if (length(missing_marker) > 0L) {
    stop("[Directional markers] Retained marker table for ", contrast_id,
         " missing canonical column(s): ", paste(missing_marker, collapse = ", "),
         call. = FALSE)
  }
  marker_table$gene_id <- sub("\\..*$", "", as.character(marker_table$gene_id))
  gene_list <- load_gene_list(marker_gene_list_path)
  if (nrow(marker_table) != as.integer(row$n_markers_after_cap) ||
      length(gene_list) != as.integer(row$n_markers_after_cap) ||
      !setequal(marker_table$gene_id, gene_list)) {
    stop("[Directional markers] Manifest/table/gene-list mismatch for ", contrast_id,
         "; manifest=", row$n_markers_after_cap,
         " table=", nrow(marker_table),
         " gene_list=", length(gene_list), call. = FALSE)
  }
  if (nrow(marker_table) == 0L) next
  marker_table$cohort <- opt$cohort
  marker_table$contrast_id <- contrast_id
  marker_table$marker_evidence_stratum <- as.character(row$marker_evidence_stratum)
  marker_table$contrast_type <- as.character(row$contrast_type)
  marker_table$focal_profile_id <- as.character(row$focal_profile_id)
  marker_table$focal_component_id <- as.character(row$focal_component_id)
  records[[contrast_id]] <- marker_table
}

cohort_dir <- file.path(opt$outdir, "cohort_level")
dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)

output_columns <- c(
  "cohort", "gene_id", "contrast_id", "contrast_type", "marker_evidence_stratum",
  "focal_profile_id", "focal_component_id", "baseMean", "wald_statistic",
  "p_value", "adjusted_p_value", "log2_fold_change_unshrunk",
  "log2_fold_change_shrunken", "log2_fold_change_posterior_sd",
  "absolute_shrunken_log2_fold_change", "effect_direction", "contrast_marker_rank",
  "direction"
)
empty_output <- data.frame(matrix(ncol = length(output_columns), nrow = 0))
names(empty_output) <- output_columns

if (length(records) == 0L) {
  retained <- empty_output
} else {
  retained <- rbindlist(records, fill = TRUE)
  retained$log2_fold_change_shrunken <- as.numeric(retained$log2_fold_change_shrunken)
  retained$direction <- ifelse(
    retained$log2_fold_change_shrunken > 0,
    "UP",
    ifelse(retained$log2_fold_change_shrunken < 0, "DOWN", "MIXED")
  )
  if (any(retained$direction == "MIXED")) {
    stop("[Directional markers] Retained marker has zero shrunken LFC", call. = FALSE)
  }
}

gene_direction <- if (nrow(retained) > 0L) {
  aggregate(direction ~ gene_id, retained, function(values) {
    values <- unique(values)
    if (all(values == "UP")) return("UP")
    if (all(values == "DOWN")) return("DOWN")
    "MIXED"
  })
} else {
  data.frame(gene_id = character(), direction = character())
}
mixed_genes <- gene_direction$gene_id[gene_direction$direction == "MIXED"]

all_up <- retained[retained$direction == "UP", output_columns, drop = FALSE]
all_down <- retained[retained$direction == "DOWN", output_columns, drop = FALSE]
all_mixed <- retained[retained$gene_id %in% mixed_genes, output_columns, drop = FALSE]

rank_rows <- function(table) {
  if (nrow(table) == 0L) return(table)
  table[order(table$adjusted_p_value, -table$absolute_shrunken_log2_fold_change, table$contrast_marker_rank, table$gene_id), , drop = FALSE]
}
top_n <- 200L
top_up <- head(rank_rows(all_up), top_n)
top_down <- head(rank_rows(all_down), top_n)
top_mixed <- head(rank_rows(all_mixed), top_n)

write_table <- function(table, filename) {
  fwrite(table[, output_columns, drop = FALSE], file.path(cohort_dir, filename), sep = "\t")
}

write_table(all_up, "all_upregulated_markers.tsv")
write_table(all_down, "all_downregulated_markers.tsv")
write_table(all_mixed, "all_mixed_direction_markers.tsv")
write_table(top_up, "top_upregulated_markers.tsv")
write_table(top_down, "top_downregulated_markers.tsv")
write_table(top_mixed, "top_mixed_direction_markers.tsv")

counts <- data.frame(
  cohort = opt$cohort,
  n_retained_marker_rows = nrow(retained),
  n_up_rows = nrow(all_up),
  n_down_rows = nrow(all_down),
  n_mixed_rows = nrow(all_mixed),
  n_mixed_genes = length(mixed_genes),
  stringsAsFactors = FALSE
)
fwrite(counts, file.path(cohort_dir, "directional_marker_counts.tsv"), sep = "\t")
writeLines(c(
  "apeglm-shrunken log2 fold change > 0: UP relative to the focal isolate or anchor contrast.",
  "apeglm-shrunken log2 fold change < 0: DOWN relative to the focal isolate or anchor contrast.",
  "MIXED indicates a gene with retained UP and DOWN evidence within the enrichment grouping."
), file.path(opt$outdir, "directional_marker_orientation.txt"))

message(sprintf(
  "[Directional markers] PASS cohort=%s retained_rows=%d mixed_genes=%d",
  opt$cohort,
  nrow(retained),
  length(mixed_genes)
))
