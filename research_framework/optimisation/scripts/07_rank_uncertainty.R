command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- command_arguments[grepl("^--file=", command_arguments)]
current_script_path <- if (length(file_argument) > 0) sub("^--file=", "", file_argument[[1]]) else file.path("research_framework", "optimisation", "scripts", "07_rank_uncertainty.R")
source(file.path(dirname(normalizePath(current_script_path, mustWork = TRUE)), "optimisation_utils.R"))
load_required_packages(c("data.table"))

arguments <- parse_key_value_arguments()
config <- load_optimisation_config(arguments$config)

output_root <- config_value(config, c("outputs", "optimisation_results_root"), file.path("research_framework", "optimisation", "results"))
docs_root <- config_value(config, c("outputs", "optimisation_docs_root"), file.path("research_framework", "optimisation", "docs"))
thesis_results_root <- config_value(config, c("inputs", "thesis_results_root"), "results")

summary_output <- file.path(output_root, "ranking_uncertainty", "ranking_uncertainty_summary.tsv")
confidence_output <- file.path(output_root, "ranking_uncertainty", "cellline_rank_confidence.tsv")
report_output <- file.path(docs_root, "06_ranking_uncertainty_report.md")
posterior_path <- file.path(output_root, "probabilistic_neighbourhoods", "posterior_neighbourhood_probabilities.tsv.gz")
components_path <- file.path(output_root, "probabilistic_graphs", "probabilistic_components.tsv")

ranking_paths <- discover_files(thesis_results_root, "cellline_to_tumour.*rank.*\\.tsv$|cell_line_to_tumour.*rank.*\\.tsv$")
if (length(ranking_paths) == 0) {
  warning("No cell-line-to-tumour ranking tables were found. The ranking uncertainty table will contain no thesis ranking rows.")
}

ranking_tables <- lapply(ranking_paths, function(path) {
  table <- read_table_if_exists(path)
  if (is.null(table) || !"cell_line" %in% names(table) || !"score" %in% names(table)) {
    return(NULL)
  }
  tumour_lineage_column <- c("tumour_lineage", "tum_lineage", "cancer_type")[c("tumour_lineage", "tum_lineage", "cancer_type") %in% names(table)][1]
  cell_lineage_column <- c("cell_lineage", "known_cancer_type")[c("cell_lineage", "known_cancer_type") %in% names(table)][1]
  if (is.na(tumour_lineage_column)) {
    table$tumour_lineage_for_uncertainty <- "unknown"
  } else {
    table$tumour_lineage_for_uncertainty <- as.character(table[[tumour_lineage_column]])
  }
  if (is.na(cell_lineage_column)) {
    table$known_cancer_type_for_uncertainty <- "unknown"
  } else {
    table$known_cancer_type_for_uncertainty <- as.character(table[[cell_lineage_column]])
  }
  table$source_file <- path
  table
})
ranking_data <- do.call(rbind, ranking_tables[!vapply(ranking_tables, is.null, logical(1))])

posterior_table <- if (file.exists(posterior_path)) data.table::fread(posterior_path, sep = "\t", data.table = FALSE, showProgress = FALSE) else NULL
posterior_by_cell <- if (!is.null(posterior_table) && nrow(posterior_table) > 0) {
  aggregate(posterior_mean_neighbourhood_probability ~ cell_line_id, data = posterior_table, FUN = mean)
} else {
  data.frame(cell_line_id = character(), posterior_mean_neighbourhood_probability = numeric())
}
names(posterior_by_cell)[names(posterior_by_cell) == "posterior_mean_neighbourhood_probability"] <- "posterior_neighbourhood_support"

components_table <- read_table_if_exists(components_path)
component_by_cell <- if (!is.null(components_table) && nrow(components_table) > 0) {
  components_table[, c("cell_line_id", "probabilistic_component_id")]
} else {
  data.frame(cell_line_id = character(), probabilistic_component_id = character())
}

if (is.null(ranking_data) || nrow(ranking_data) == 0) {
  confidence_table <- data.frame(
    cell_line_id = character(),
    known_cancer_type = character(),
    top_ranked_cancer_type = character(),
    top1_score = numeric(),
    top2_score = numeric(),
    rank_margin = numeric(),
    top1_confidence_class = character(),
    graph_ranking_agreement = character(),
    posterior_neighbourhood_support = numeric(),
    ranking_uncertainty_class = character(),
    notes = character(),
    stringsAsFactors = FALSE
  )
} else {
  ranking_data$score <- safe_numeric(ranking_data$score, default = NA_real_)
  lineage_scores <- aggregate(score ~ cell_line + known_cancer_type_for_uncertainty + tumour_lineage_for_uncertainty, data = ranking_data, FUN = mean)
  confidence_rows <- lapply(split(lineage_scores, lineage_scores$cell_line), function(cell_table) {
    cell_table <- cell_table[order(-cell_table$score), ]
    top1 <- cell_table[1, ]
    top2_score <- if (nrow(cell_table) >= 2) cell_table$score[[2]] else NA_real_
    rank_margin <- if (is.na(top2_score)) NA_real_ else top1$score - top2_score
    confidence_class <- if (is.na(rank_margin)) {
      "single_class_observed"
    } else if (rank_margin >= 0.10) {
      "high_top1_margin"
    } else if (rank_margin >= 0.03) {
      "moderate_top1_margin"
    } else {
      "low_top1_margin"
    }
    data.frame(
      cell_line_id = top1$cell_line,
      known_cancer_type = top1$known_cancer_type_for_uncertainty,
      top_ranked_cancer_type = top1$tumour_lineage_for_uncertainty,
      top1_score = top1$score,
      top2_score = top2_score,
      rank_margin = rank_margin,
      top1_confidence_class = confidence_class,
      stringsAsFactors = FALSE
    )
  })
  confidence_table <- do.call(rbind, confidence_rows)
  confidence_table <- merge(confidence_table, posterior_by_cell, by = "cell_line_id", all.x = TRUE)
  confidence_table <- merge(confidence_table, component_by_cell, by = "cell_line_id", all.x = TRUE)
  confidence_table$graph_ranking_agreement <- ifelse(
    is.na(confidence_table$probabilistic_component_id),
    "not_estimated",
    ifelse(confidence_table$known_cancer_type == confidence_table$top_ranked_cancer_type, "lineage_agrees_with_top_rank", "lineage_disagrees_with_top_rank")
  )
  confidence_table$ranking_uncertainty_class <- ifelse(
    confidence_table$top1_confidence_class == "high_top1_margin" & confidence_table$graph_ranking_agreement != "lineage_disagrees_with_top_rank",
    "low_uncertainty",
    ifelse(confidence_table$top1_confidence_class == "low_top1_margin" | confidence_table$graph_ranking_agreement == "lineage_disagrees_with_top_rank", "high_uncertainty", "moderate_uncertainty")
  )
  confidence_table$notes <- ifelse(is.na(confidence_table$posterior_neighbourhood_support), "posterior neighbourhood evidence unavailable", "")
  confidence_table <- confidence_table[, c(
    "cell_line_id", "known_cancer_type", "top_ranked_cancer_type", "top1_score", "top2_score",
    "rank_margin", "top1_confidence_class", "graph_ranking_agreement",
    "posterior_neighbourhood_support", "ranking_uncertainty_class", "notes"
  )]
}

summary_table <- as.data.frame(table(confidence_table$ranking_uncertainty_class), stringsAsFactors = FALSE)
if (nrow(summary_table) == 0) {
  summary_table <- data.frame(ranking_uncertainty_class = "not_estimated", n_cell_lines = 0L)
} else {
  names(summary_table) <- c("ranking_uncertainty_class", "n_cell_lines")
}

write_tsv_table(summary_table, summary_output)
write_tsv_table(confidence_table, confidence_output)

report_lines <- c(
  "# Ranking uncertainty report",
  "",
  "## Scope",
  "",
  "The ranking uncertainty step consumed available thesis ranking tables and combined rank margins with posterior neighbourhood and probabilistic graph context when available.",
  "",
  "## Assumptions",
  "",
  "- Rank margin was computed as the difference between the mean score for the top-ranked tumour lineage and the second-ranked tumour lineage for each cell line.",
  "- Graph-ranking agreement was recorded as lineage agreement where probabilistic component information was available.",
  "",
  "## Outputs",
  "",
  paste0("- `", summary_output, "`"),
  paste0("- `", confidence_output, "`"),
  "",
  "## Manual inspection required",
  "",
  "- Inspect high-uncertainty cell lines before making cancer-type matching claims.",
  "- Inspect cases where rank margin is low despite high posterior neighbourhood evidence."
)
write_markdown(report_lines, report_output)
cat("Ranking uncertainty completed: ", confidence_output, "\n", sep = "")

