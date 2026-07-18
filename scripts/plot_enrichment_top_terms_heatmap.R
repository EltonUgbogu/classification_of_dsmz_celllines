#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option("--input", type = "character"),
  make_option("--query-manifest", type = "character"),
  make_option("--display-mapping", type = "character"),
  make_option("--excluded-term-mapping", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--top-terms-per-query", type = "integer"),
  make_option("--maximum-terms", type = "integer"),
  make_option("--score-cap", type = "numeric"),
  make_option("--term-wrap-width", type = "integer"),
  make_option("--figure-width", type = "numeric"),
  make_option("--minimum-figure-height", type = "numeric"),
  make_option("--row-height-increment", type = "numeric"),
  make_option("--correction-method", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (nm in c("input", "query-manifest", "display-mapping", "excluded-term-mapping",
             "outdir", "top-terms-per-query", "maximum-terms", "score-cap",
             "term-wrap-width", "figure-width", "minimum-figure-height",
             "row-height-increment", "correction-method")) {
  if (is.null(opt[[nm]]) || (is.character(opt[[nm]]) && !nzchar(opt[[nm]]))) {
    stop("--", nm, " is required", call. = FALSE)
  }
}

input_file <- normalizePath(opt$input, mustWork = TRUE)
manifest_file <- normalizePath(opt$`query-manifest`, mustWork = TRUE)
display_mapping_file <- normalizePath(opt$`display-mapping`, mustWork = TRUE)
excluded_mapping_file <- normalizePath(opt$`excluded-term-mapping`, mustWork = FALSE)
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

wrap_label <- function(x, width) {
  vapply(strwrap(as.character(x), width = width, simplify = FALSE), paste, character(1), collapse = "\n")
}

terms <- fread(input_file, sep = "\t")
manifest <- fread(manifest_file, sep = "\t")
mapping <- fread(display_mapping_file, sep = "\t")

required_terms <- c("query_id", "analysis_mode", "source", "term_id", "term_name",
                    "adjusted_enrichment_p_value", "significant", "intersection_size",
                    "enrichment_ratio")
missing_terms <- setdiff(required_terms, names(terms))
if (length(missing_terms)) stop("Significant-term input missing columns: ", paste(missing_terms, collapse = ","), call. = FALSE)
if (anyDuplicated(manifest$query_id)) stop("Manifest query IDs are not unique", call. = FALSE)

required_mapping <- c("query_family", "attributed_cancer_type", "marker_evidence_stratum",
                      "evidence_class", "direction_consistency_class", "display_group_id",
                      "display_label", "display_order", "include_in_primary_heatmap")
missing_mapping <- setdiff(required_mapping, names(mapping))
if (length(missing_mapping)) stop("Display mapping missing columns: ", paste(missing_mapping, collapse = ","), call. = FALSE)

mapping <- mapping[tolower(as.character(include_in_primary_heatmap)) %in% c("true", "t", "1", "yes")]
if (!nrow(mapping)) stop("No primary heatmap display mappings are enabled", call. = FALSE)

metadata_cols <- c("query_family", "attributed_cancer_type", "marker_evidence_stratum",
                   "evidence_class", "direction_consistency_class")
for (col in metadata_cols) {
  if (!col %in% names(manifest)) manifest[, (col) := ""]
  manifest[[col]][is.na(manifest[[col]])] <- ""
  mapping[[col]][is.na(mapping[[col]])] <- ""
}

selected_queries <- list()
for (i in seq_len(nrow(mapping))) {
  map_row <- mapping[i]
  candidates <- manifest[query_execution_status == "runnable"]
  for (col in metadata_cols) {
    value <- trimws(as.character(map_row[[col]]))
    if (nzchar(value)) candidates <- candidates[as.character(get(col)) == value]
  }
  if (nrow(candidates) == 0L) next
  if (nrow(candidates) != 1L) {
    stop("Display mapping does not resolve to exactly one query for display_group_id=",
         map_row$display_group_id, " matched=", paste(candidates$query_id, collapse = ","), call. = FALSE)
  }
  selected_queries[[length(selected_queries) + 1L]] <- data.table(
    query_id = candidates$query_id,
    display_group_id = map_row$display_group_id,
    display_label = map_row$display_label,
    display_order = as.numeric(map_row$display_order)
  )
}
selected_queries <- rbindlist(selected_queries, fill = TRUE)
if (!nrow(selected_queries)) stop("No display mappings resolved to current runnable queries", call. = FALSE)
if (anyDuplicated(selected_queries$display_group_id)) {
  stop("display_group_id maps to multiple query IDs in the primary heatmap", call. = FALSE)
}
if (anyDuplicated(selected_queries$query_id)) {
  stop("A primary heatmap query maps to multiple display groups", call. = FALSE)
}

terms <- terms[analysis_mode == "primary" & significant == TRUE & query_id %in% selected_queries$query_id]
excluded <- data.table()
if (file.exists(excluded_mapping_file) && file.info(excluded_mapping_file)$size > 0) {
  excluded_map <- fread(excluded_mapping_file, sep = "\t")
  if (!all(c("source", "term_id") %in% names(excluded_map))) {
    stop("Excluded-term mapping must contain source and term_id", call. = FALSE)
  }
  excluded <- merge(terms, unique(excluded_map[, .(source, term_id)]), by = c("source", "term_id"))
  terms <- terms[!paste(source, term_id, sep = "\r") %in% paste(excluded_map$source, excluded_map$term_id, sep = "\r")]
}
if (!nrow(terms)) stop("No significant primary terms remain after exact exclusions", call. = FALSE)

# Per-query nomination order is deterministic: adjusted p-value ascending,
# enrichment ratio descending, intersection size descending, source, then term ID.
setorder(terms, query_id, adjusted_enrichment_p_value, -enrichment_ratio, -intersection_size, source, term_id)
terms[, nomination_rank_within_query := seq_len(.N), by = query_id]
nominees <- unique(terms[nomination_rank_within_query <= opt$`top-terms-per-query`, .(source, term_id)])
nominee_counts <- terms[
  nomination_rank_within_query <= opt$`top-terms-per-query`,
  .(n_primary_queries_nominating = uniqueN(query_id)),
  by = .(source, term_id)
]

candidate_stats <- merge(
  nominees,
  terms[, .(
    n_primary_queries_significant = uniqueN(query_id),
    minimum_adjusted_enrichment_p_value = min(adjusted_enrichment_p_value, na.rm = TRUE),
    representative_term_name = term_name[which.min(adjusted_enrichment_p_value)]
  ), by = .(source, term_id)],
  by = c("source", "term_id"),
  all.x = TRUE
)
candidate_stats <- merge(candidate_stats, nominee_counts, by = c("source", "term_id"), all.x = TRUE)
setorder(
  candidate_stats,
  -n_primary_queries_nominating,
  -n_primary_queries_significant,
  minimum_adjusted_enrichment_p_value,
  source,
  term_id
)
selected_terms <- candidate_stats[seq_len(min(.N, opt$`maximum-terms`))]

complete_selected <- merge(terms, selected_terms[, .(source, term_id)], by = c("source", "term_id"))
plot_grid <- CJ(
  query_id = selected_queries$query_id,
  term_key = paste(selected_terms$source, selected_terms$term_id, sep = "::"),
  unique = TRUE
)
selected_terms[, term_key := paste(source, term_id, sep = "::")]
complete_selected[, term_key := paste(source, term_id, sep = "::")]
plot_data <- merge(plot_grid, complete_selected, by = c("query_id", "term_key"), all.x = TRUE)
plot_data <- merge(plot_data, selected_queries, by = "query_id", all.x = TRUE)
plot_data <- merge(
  plot_data,
  selected_terms[, .(term_key, source_selected = source, term_id_selected = term_id,
                     representative_term_name, term_order = seq_len(.N))],
  by = "term_key",
  all.x = TRUE
)
plot_data[, plotting_score := -log10(adjusted_enrichment_p_value)]
plot_data[is.infinite(plotting_score), plotting_score := opt$`score-cap`]
plot_data[, plotting_score_capped := pmin(plotting_score, opt$`score-cap`)]
plot_data[, display_label := factor(display_label, levels = selected_queries[order(display_order)]$display_label)]
plot_data[, term_label := wrap_label(representative_term_name, opt$`term-wrap-width`)]
plot_data[, term_label := factor(term_label, levels = rev(unique(plot_data[order(term_order)]$term_label)))]

legend_title <- if (identical(opt$`correction-method`, "g_SCS")) {
  "-log10(g:SCS-adjusted p-value)"
} else {
  paste0("-log10(", opt$`correction-method`, "-adjusted p-value)")
}

p <- ggplot(plot_data, aes(x = display_label, y = term_label, fill = plotting_score_capped)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#d9dde4", high = "#b23a48", na.value = "#f0f1f3", limits = c(0, opt$`score-cap`), name = legend_title) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

height <- max(opt$`minimum-figure-height`, length(unique(plot_data$term_key)) * opt$`row-height-increment` + 2)
pdf_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap.pdf")
png_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap.png")
matrix_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_matrix.tsv")
selected_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_selected_terms.tsv")
excluded_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_excluded_terms.tsv")
provenance_path <- file.path(outdir, "Fig_enrichment_top_terms_heatmap_provenance.tsv")
ggsave(pdf_path, p, width = opt$`figure-width`, height = height, units = "in")
ggsave(png_path, p, width = opt$`figure-width`, height = height, units = "in", dpi = 300)

fwrite(plot_data, matrix_path, sep = "\t", na = "NA")
fwrite(selected_terms, selected_path, sep = "\t", na = "NA")
fwrite(excluded, excluded_path, sep = "\t", na = "NA")
provenance <- data.table(
  figure_name = "Fig_enrichment_top_terms_heatmap.pdf",
  figure_type = "primary_enrichment_heatmap",
  script = "scripts/plot_enrichment_top_terms_heatmap.R",
  input_files = paste(normalizePath(c(input_file, manifest_file, display_mapping_file, excluded_mapping_file), mustWork = FALSE), collapse = ";"),
  output_files = paste(normalizePath(c(pdf_path, png_path, matrix_path, selected_path, excluded_path), mustWork = FALSE), collapse = ";"),
  term_selection = "all_significant_terms_exact_exclusions_per_query_nomination_union_then_complete_matrix",
  one_display_group_per_query = TRUE,
  correction_method = opt$`correction-method`,
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
fwrite(provenance, provenance_path, sep = "\t")
message("[Primary heatmap] Each column corresponds to one actual enrichment query.")
