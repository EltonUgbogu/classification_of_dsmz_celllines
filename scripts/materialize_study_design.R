#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(optparse)
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml"),
  make_option("--study-design", type = "character", default = "config/study_design.yaml"),
  make_option("--profile", type = "character", default = NULL),
  make_option("--out-question", type = "character"),
  make_option("--out-cohorts", type = "character"),
  make_option("--out-labels", type = "character"),
  make_option("--out-inference", type = "character"),
  make_option("--out-endpoints", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

config_dir <- dirname(normalizePath(opt$config))
lib_config_path <- file.path(config_dir, "..", "scripts", "lib_config.R")
source(lib_config_path)

raw_cfg <- yaml::read_yaml(opt$config)
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)
study <- yaml::read_yaml(opt$`study-design`)
cohort <- study$cohorts[[profile]]
if (is.null(cohort)) stop("Profile not found in study design: ", profile)

feature_methods <- cfg$feature_sets$methods %||% names(raw_cfg$defaults$feature_selection$method_topn)
distances <- cfg$feature_sets$distances %||% c("euc", "corr")
directions <- cfg$tumour_neighbourhoods$directions %||% as.character(outer(feature_methods, distances, paste, sep = "_"))
k_grid <- raw_cfg$defaults$clustering$k_grid %||% c(2, 3, 4, 5, 6, 7, 8)
consensus <- study$candidate_inference$consensus_metric
ranking_fields <- paste(study$candidate_inference$final_ranking$ranking_fields %||% character(), collapse = ";")

question_lines <- c(
  paste0("Study: ", study$study$title),
  paste0("Profile: ", profile),
  paste0("Primary question: ", study$study$primary_question),
  paste0("Primary endpoint: ", study$study$primary_endpoint),
  paste0("Scope note: ", study$study$scope_note)
)
writeLines(question_lines, opt$`out-question`)

cohort_tbl <- tibble(
  profile = profile,
  disease = cohort$disease %||% NA_character_,
  tumour_source = unlist(cohort$tumour_sources %||% list()),
  cell_line_panel = cohort$cell_line_panel %||% study$candidate_inference$cell_line_panel %||% NA_character_,
  current_analysis_input = cohort$current_analysis_input %||% NA_character_,
  current_preprocessing_output = cohort$current_preprocessing_output %||% NA_character_,
  source_evidence = cohort$source_evidence %||% NA_character_
)
write_tsv(cohort_tbl, opt$`out-cohorts`)

labels <- cohort$runtime_labels %||% list()
if (length(labels) == 0) {
  label_tbl <- tibble(
    profile = character(),
    label = character(),
    entity = character(),
    source_file = character(),
    status = character(),
    purpose = character()
  )
} else {
  label_tbl <- bind_rows(lapply(labels, function(x) {
    tibble(
      profile = profile,
      label = x$label %||% NA_character_,
      entity = x$entity %||% NA_character_,
      source_file = x$source_file %||% NA_character_,
      status = x$status %||% NA_character_,
      purpose = x$purpose %||% NA_character_
    )
  }))
}
write_tsv(label_tbl, opt$`out-labels`)

infer_tbl <- tibble(direction = directions) %>%
  mutate(
    feature_method = sub("_(euc|corr)$", "", direction),
    distance_metric = sub("^.*_", "", direction),
    clustering_algorithms = paste(study$candidate_inference$clustering_algorithms %||% c("hierarchical", "kmeans"), collapse = ";"),
    k_grid = paste(k_grid, collapse = ";"),
    consensus_metric = consensus$name %||% "p_consensus",
    strong_support_threshold = consensus$strong_support_threshold %||% 0.7,
    ranking_fields = ranking_fields
  ) %>%
  mutate(profile = profile, .before = 1)
write_tsv(infer_tbl, opt$`out-inference`)

endpoints <- bind_rows(
  tibble(endpoint_class = "primary", endpoint = study$study$primary_endpoint),
  tibble(endpoint_class = "secondary", endpoint = unlist(study$study$secondary_endpoints %||% list())),
  tibble(endpoint_class = "validation", endpoint = unlist(study$validation_framework$biological_validation_axes %||% list())),
  tibble(endpoint_class = "baseline", endpoint = unlist(study$validation_framework$baseline_analyses %||% list()))
) %>% mutate(profile = profile, .before = 1)
write_tsv(endpoints, opt$`out-endpoints`)
