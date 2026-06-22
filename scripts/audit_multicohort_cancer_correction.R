#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- sub("^--", "", key)
    value <- TRUE
    if (i + 1L <= length(args) && !startsWith(args[[i + 1L]], "--")) {
      value <- args[[i + 1L]]
      i <- i + 1L
    }
    out[[name]] <- value
    i <- i + 1L
  }
  out
}

deep_merge <- function(base, override) {
  merged <- base
  for (nm in names(override)) {
    if (is.list(override[[nm]]) && is.list(merged[[nm]])) {
      merged[[nm]] <- deep_merge(merged[[nm]], override[[nm]])
    } else {
      merged[[nm]] <- override[[nm]]
    }
  }
  merged
}

status <- function(ok, warn = FALSE) {
  if (isTRUE(ok)) return("PASS")
  if (isTRUE(warn)) return("WARN")
  "FAIL"
}

line_count <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(0L)
  length(readLines(path, warn = FALSE))
}

safe_cols <- function(df) {
  paste(names(df), collapse = ", ")
}

method_id_from_kind <- function(kind) {
  kind <- sub("^ccp_hc_", "CCP_HC_", kind)
  kind <- sub("^ccp_kmeans_", "CCP_KM_", kind)
  kind
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_dir <- normalizePath(args[["project-dir"]] %||% getwd(), mustWork = TRUE)
profile <- args[["profile"]] %||% "multicohort_cancer"
config_file <- normalizePath(args[["config"]] %||% file.path(project_dir, "config", "config.yaml"), mustWork = TRUE)

cfg_full <- yaml::read_yaml(config_file)
cfg <- deep_merge(cfg_full$defaults %||% list(), cfg_full$profiles[[profile]] %||% list())

abs_path <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(path)
  if (grepl("^/", path)) return(path)
  file.path(project_dir, path)
}

directions <- cfg$tumour_neighbourhoods$directions %||% character()
unsup_root <- abs_path(cfg$paths$unsup_root %||% file.path("results", "unsupervised", profile))
tn_root <- file.path(unsup_root, "tumour_neighbourhoods")
cons_root <- file.path(unsup_root, "consensus")
log_root <- file.path(project_dir, "logs", profile)
validation_dir <- file.path(project_dir, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_tsv <- args[["output"]] %||% file.path(validation_dir, paste0("multicohort_cancer_correction_audit_", timestamp, ".tsv"))
out_txt <- sub("\\.tsv$", ".txt", out_tsv)

meta_file <- abs_path(cfg$paths$meta_tsv %||% cfg$paths$dsmz_meta_csv %||% "")
meta <- if (file.exists(meta_file)) {
  read.delim(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame()
}
sample_id_col <- cfg$deseq2$sample_id_col %||% NA_character_
staged_sample_id_col <- cfg$deseq2$staged_sample_id_col %||% NA_character_

marker_cfg <- cfg_full$defaults$marker_postprocessing %||% list()
pan_profiles <- marker_cfg$pan_cancer$disease_profiles %||% character()
include_heme <- isTRUE(marker_cfg$expression_matrix$include_heme)
mc_profiles <- cfg_full$multicohort_cancer$profiles %||% character()

summary_rows <- list()
add_row <- function(section, metric, expected, observed, stat, details = "") {
  summary_rows[[length(summary_rows) + 1L]] <<- data.frame(
    section = section,
    metric = metric,
    expected = as.character(expected),
    observed = as.character(observed),
    status = stat,
    details = as.character(details),
    stringsAsFactors = FALSE
  )
}

add_row("preflight", "metadata_sample_id_col", sample_id_col, sample_id_col %in% names(meta),
        status(sample_id_col %in% names(meta)),
        paste0("metadata=", meta_file, "; columns=", if (ncol(meta)) safe_cols(meta) else "<missing>"))
add_row("preflight", "staged_sample_id_col", staged_sample_id_col, staged_sample_id_col %in% names(meta),
        status(is.na(staged_sample_id_col) || staged_sample_id_col %in% names(meta), warn = TRUE),
        paste0("metadata=", meta_file))
add_row("preflight", "multicohort_profiles_no_heme", "BRCA,NBL,RBL only", paste(mc_profiles, collapse = ","),
        status(!"heme" %in% tolower(mc_profiles)))
add_row("preflight", "pan_cancer_disease_profiles_no_heme", "brca,nbl,rbl only", paste(pan_profiles, collapse = ","),
        status(!"heme" %in% tolower(pan_profiles)))
add_row("preflight", "pan_cancer_expression_include_heme", "false", include_heme,
        status(!include_heme))

profile_paths <- vapply(c("brca", "nbl", "rbl"), function(p) {
  pcfg <- deep_merge(cfg_full$defaults %||% list(), cfg_full$profiles[[p]] %||% list())
  abs_path(pcfg$paths$vst_joint_rds %||% "")
}, character(1))
for (p in names(profile_paths)) {
  add_row("preflight", paste0(p, "_vst_exists"), profile_paths[[p]], file.exists(profile_paths[[p]]),
          status(file.exists(profile_paths[[p]])))
}
heme_paths_planned <- grep("/heme/|heme_", profile_paths, value = TRUE, ignore.case = TRUE)
add_row("preflight", "planned_per_disease_paths_contain_heme", "none", length(heme_paths_planned),
        status(length(heme_paths_planned) == 0L),
        paste(heme_paths_planned, collapse = ";"))

hc_kinds <- c("ccp_hc_expr_cell_tumour", "ccp_hc_pca_cell_tumour")
km_kinds <- c("ccp_kmeans_expr_cell_tumour", "ccp_kmeans_pca_cell_tumour")
tn_jobs <- list()

for (direction in directions) {
  kinds <- hc_kinds
  families <- rep("hc", length(hc_kinds))
  if (grepl("_euc$", direction)) {
    kinds <- c(kinds, km_kinds)
    families <- c(families, rep("km", length(km_kinds)))
  }
  for (j in seq_along(kinds)) {
    kind <- kinds[[j]]
    method_id <- method_id_from_kind(kind)
    cluster <- file.path(cons_root, direction, kind, paste0(kind, "_clusters_optimal.rds"))
    csv <- file.path(tn_root, direction, method_id, paste0("Top_m_long_", method_id, ".csv"))
    rows <- max(0L, line_count(csv) - 1L)
    tn_jobs[[length(tn_jobs) + 1L]] <- data.frame(
      direction = direction,
      family = families[[j]],
      method = method_id,
      cluster_exists = file.exists(cluster),
      top_m_long_csv = csv,
      csv_exists = file.exists(csv),
      csv_rows = rows,
      output_dir_nonempty = dir.exists(dirname(csv)) &&
        length(list.files(dirname(csv), all.files = FALSE, no.. = TRUE)) > 0L,
      stringsAsFactors = FALSE
    )
  }
}

tn_df <- if (length(tn_jobs)) do.call(rbind, tn_jobs) else data.frame()
tn_expected <- nrow(tn_df)
tn_completed <- sum(tn_df$cluster_exists & tn_df$csv_exists & tn_df$csv_rows > 0)
tn_missing_csv <- sum(tn_df$cluster_exists & (!tn_df$csv_exists | tn_df$csv_rows == 0))
tn_missing_cluster <- sum(!tn_df$cluster_exists)
tn_nonempty_dirs <- sum(tn_df$output_dir_nonempty)

add_row("tumour_neighbourhoods", "expected_jobs", tn_expected, tn_expected, "INFO")
add_row("tumour_neighbourhoods", "completed_jobs", tn_expected, tn_completed,
        status(tn_completed == tn_expected))
add_row("tumour_neighbourhoods", "non_empty_output_directories", tn_expected, tn_nonempty_dirs,
        status(tn_nonempty_dirs == tn_expected))
add_row("tumour_neighbourhoods", "missing_or_empty_top_m_long_csv", 0, tn_missing_csv,
        status(tn_missing_csv == 0L))
add_row("tumour_neighbourhoods", "missing_upstream_cluster_rds", 0, tn_missing_cluster,
        status(tn_missing_cluster == 0L))

tn_job_path <- sub("\\.tsv$", "_tumour_neighbourhood_jobs.tsv", out_tsv)
write.table(tn_df, tn_job_path, sep = "\t", quote = FALSE, row.names = FALSE)

nh_logs <- list.files(log_root, pattern = "^tumour_nh_(hc|km)_.*\\.log$", full.names = TRUE)
log_text <- unlist(lapply(nh_logs, readLines, warn = FALSE), use.names = FALSE)
sample_errors <- grep(
  "sample_id_col 'sample_name' not found|Metadata sample ID column validation failed|No configured sample ID column found",
  log_text,
  value = TRUE
)
dsmz_errors <- grep(
  "raw_cell_line_col '.*' not found in metadata columns|DSMZ cell-line labels are missing|Non-DSMZ metadata rows",
  log_text,
  value = TRUE
)
failed_logs <- nh_logs[vapply(nh_logs, function(path) {
  txt <- tail(readLines(path, warn = FALSE), 40)
  any(grepl("Error|Execution halted|No Top_m_long|No KM-derived", txt))
}, logical(1))]
add_row("logs", "tumour_neighbourhood_sample_id_errors", 0, length(sample_errors),
        status(length(sample_errors) == 0L), paste(basename(nh_logs), collapse = ";"))
add_row("logs", "tumour_neighbourhood_dsmz_label_errors", 0, length(dsmz_errors),
        status(length(dsmz_errors) == 0L))
add_row("logs", "failed_tumour_neighbourhood_logs", 0, length(failed_logs),
        status(length(failed_logs) == 0L), paste(basename(failed_logs), collapse = ";"))

downstream <- c(
  final_consensus_all = file.path(tn_root, "final_consensus_all", "resolved_dsmz_neighbours.tsv"),
  pan_expr = abs_path(marker_cfg$expression_matrix$output_rds %||% "results/unsupervised/pan_cancer/inputs/pan_cancer_feature_expr.rds"),
  pan_cor = file.path(project_dir, "results", "unsupervised", "pan_cancer", "pan_cancer_cor.rds"),
  pan_graph = file.path(project_dir, "results", "unsupervised", "pan_cancer", "graph", "pan_cancer_graph_edges.tsv"),
  mapping_metrics = file.path(project_dir, "results", "unsupervised", "pan_cancer", "tumour_mapping", "metrics_summary.tsv"),
  pan_graph_pdf = file.path(project_dir, "results", "unsupervised", "pan_cancer", "figures", "Fig_pan_cancer_graph.pdf"),
  ecdf_pdf = file.path(project_dir, "results", "unsupervised", "pan_cancer", "figures", "ecdf_plots", "ecdf_model_prioritisation_combined.pdf")
)
for (nm in names(downstream)) {
  add_row("downstream", nm, downstream[[nm]], file.exists(downstream[[nm]]) && file.info(downstream[[nm]])$size > 0,
          status(file.exists(downstream[[nm]]) && file.info(downstream[[nm]])$size > 0))
}

audit_files <- list.files(file.path(tn_root, "validation"), pattern = "^dsmz_cell_line_label_audit_.*\\.tsv$", full.names = TRUE)
latest_audit <- if (length(audit_files)) audit_files[which.max(file.info(audit_files)$mtime)] else NA_character_
if (!is.na(latest_audit)) {
  audit <- read.delim(latest_audit, stringsAsFactors = FALSE, check.names = FALSE)
  dsmz_rows <- grepl("^NG[-_]", audit$sample_id) & toupper(audit$cohort) == "DSMZ" &
    grepl("cell", audit$sample_type, ignore.case = TRUE)
  blank_dsmz <- sum(dsmz_rows & (!nzchar(audit$derived_cell_line) | is.na(audit$derived_cell_line)))
  labelled_non_dsmz <- sum(!dsmz_rows & nzchar(audit$derived_cell_line))
  add_row("dsmz_label_audit", "latest_audit_file", "present", latest_audit, "INFO")
  add_row("dsmz_label_audit", "blank_dsmz_labels", 0, blank_dsmz, status(blank_dsmz == 0L))
  add_row("dsmz_label_audit", "labelled_non_dsmz_rows", 0, labelled_non_dsmz, status(labelled_non_dsmz == 0L))
  for (sample in c("NG-30919_Y_79_lib626626_10098_1", "NG-30919_RBL_14_lib628470_10098_1", "NG-30919_WERI_RB1_lib628472_10098_1")) {
    observed <- audit$derived_cell_line[match(sample, audit$sample_id)]
    expected <- sub("_lib.*$", "", sub("^NG[-_][^_]+_", "", sample))
    add_row("dsmz_label_audit", paste0("expected_label_", sample), expected, observed,
            status(length(observed) == 1L && !is.na(observed) && observed == expected))
  }
  for (cell_line in c("RBL_15", "RBL_20")) {
    n_ids <- sum(audit$derived_cell_line == cell_line, na.rm = TRUE)
    add_row("dsmz_label_audit", paste0(cell_line, "_library_count"), ">=2", n_ids,
            status(n_ids >= 2L))
  }
} else {
  add_row("dsmz_label_audit", "latest_audit_file", "present", "missing", "FAIL")
}

summary <- do.call(rbind, summary_rows)
write.table(summary, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

overall <- if (any(summary$status == "FAIL")) "FAIL" else if (any(summary$status == "WARN")) "WARN" else "PASS"
txt <- c(
  "MULTICOHORT CANCER CORRECTION AUDIT",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Project:", project_dir),
  paste("Profile:", profile),
  paste("Overall:", overall),
  "",
  paste("Summary TSV:", out_tsv),
  paste("Tumour-neighbourhood job TSV:", tn_job_path),
  "",
  capture.output(print(summary, row.names = FALSE))
)
writeLines(txt, out_txt)
cat(paste(txt, collapse = "\n"), "\n")
