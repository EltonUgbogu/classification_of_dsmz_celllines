#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

option_list <- list(
  make_option("--config", type = "character", default = "config/config.yaml"),
  make_option("--profile", type = "character", default = NULL),
  make_option("--out-tsv", type = "character"),
  make_option("--out-notes", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

config_dir <- dirname(normalizePath(opt$config))
pipe_root <- normalizePath(file.path(config_dir, ".."))
source(file.path(pipe_root, "scripts", "lib_config.R"))
profile <- opt$profile %||% Sys.getenv("SNAKEMAKE_PROFILE", "default")
cfg <- read_profiled_config(opt$config, profile)
abs_from_root <- function(p) if (grepl("^/", p)) p else file.path(pipe_root, p)
unsup_root <- abs_from_root(cfg$paths$unsup_root)

sil_files <- list.files(unsup_root, pattern = "silhouette", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (length(sil_files) == 0) {
  report_tbl <- tibble(
    cohort = profile,
    silhouette_artifacts_found = 0L,
    status = "no_persisted_silhouette_artifacts",
    note = "Silhouette is used in clustering code, but no standalone silhouette files were found under the profile output root."
  )
  notes <- c(
    sprintf("Profile: %s", profile),
    "No persisted silhouette artifacts were found under the current unsupervised output root.",
    "This does not mean silhouette was unused; the clustering code uses silhouette during model selection.",
    "It means the present profile outputs do not expose a dedicated silhouette report that can be summarised post hoc."
  )
} else {
  report_tbl <- tibble(
    cohort = profile,
    silhouette_artifacts_found = length(sil_files),
    status = "persisted_silhouette_artifacts_detected",
    file = sub(paste0("^", unsup_root, "/?"), "", sil_files)
  )
  notes <- c(
    sprintf("Profile: %s", profile),
    sprintf("Found %d persisted silhouette-related artifact(s) under the profile output root.", length(sil_files)),
    "This report is file-based: it enumerates persisted silhouette outputs when they exist.",
    "Interpretation should use the underlying artifact contents and associated clustering context."
  )
}
write_tsv(report_tbl, opt$`out-tsv`)
writeLines(notes, opt$`out-notes`)
