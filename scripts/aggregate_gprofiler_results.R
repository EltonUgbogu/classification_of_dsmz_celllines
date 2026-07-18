#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--pipeline-root", type = "character"),
  make_option("--query-manifest", type = "character"),
  make_option("--gprofiler-dir", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--alpha", type = "numeric"),
  make_option("--run-iea-sensitivity", type = "logical")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (nm in c("pipeline-root", "query-manifest", "gprofiler-dir", "outdir", "alpha", "run-iea-sensitivity")) {
  if (is.null(opt[[nm]]) || (is.character(opt[[nm]]) && !nzchar(opt[[nm]]))) {
    stop("--", nm, " is required", call. = FALSE)
  }
}

pipeline_root <- normalizePath(opt$`pipeline-root`, mustWork = TRUE)
query_manifest_path <- normalizePath(opt$`query-manifest`, mustWork = TRUE)
gprofiler_dir <- normalizePath(opt$`gprofiler-dir`, mustWork = FALSE)
query_results_root <- file.path(gprofiler_dir, "query_results")
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rel_to_abs <- function(path) {
  value <- trimws(as.character(path))
  if (!nzchar(value)) return("")
  if (grepl("^/", value)) return(value)
  file.path(pipeline_root, value)
}

sha256_file <- function(path) tools::sha256sum(path)[[1]]

read_genes <- function(path) {
  values <- readLines(path, warn = FALSE)
  values <- trimws(sub("\t.*$", "", values))
  values <- values[nzchar(values) & values != "gene_id"]
  unique(values)
}

read_fingerprint <- function(path) {
  if (!file.exists(path)) stop("Missing query-mode fingerprint: ", path, call. = FALSE)
  x <- fread(path, sep = "\t")
  if (!all(c("field", "value") %in% names(x))) stop("Invalid fingerprint schema: ", path, call. = FALSE)
  stats::setNames(as.character(x$value), as.character(x$field))
}

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", na = "NA")
}

manifest <- fread(query_manifest_path, sep = "\t")
required_manifest <- c("query_id", "query_family", "gene_count", "gene_list_path",
                       "background_count", "background_path", "background_sha256",
                       "query_execution_status")
missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest)) {
  stop("Query manifest missing required columns: ", paste(missing_manifest, collapse = ","), call. = FALSE)
}
if (anyDuplicated(manifest$query_id)) stop("Duplicate query_id values in query manifest", call. = FALSE)
manifest[, gene_list_abs := vapply(gene_list_path, rel_to_abs, character(1))]
manifest[, background_abs := vapply(background_path, rel_to_abs, character(1))]

runnable <- manifest[query_execution_status == "runnable"]
required_modes <- c("primary")
if (isTRUE(opt$`run-iea-sensitivity`)) required_modes <- c(required_modes, "iea_sensitivity")

extra_dirs <- character()
for (mode in required_modes) {
  mode_dir <- file.path(query_results_root, mode)
  if (dir.exists(mode_dir)) {
    observed <- basename(list.dirs(mode_dir, full.names = TRUE, recursive = FALSE))
    extra <- setdiff(observed, runnable$query_id)
    if (length(extra)) extra_dirs <- c(extra_dirs, file.path(mode, extra))
  }
}
if (length(extra_dirs)) {
  stop("Found query result directories absent from current manifest: ", paste(extra_dirs, collapse = ","), call. = FALSE)
}

term_tables <- list()
run_rows <- list()
data_versions <- character()

for (i in seq_len(nrow(runnable))) {
  row <- runnable[i]
  genes <- read_genes(row$gene_list_abs)
  background <- read_genes(row$background_abs)
  if (length(genes) != as.integer(row$gene_count)) {
    stop("Manifest gene_count mismatch during aggregation: ", row$query_id, call. = FALSE)
  }
  if (length(background) != as.integer(row$background_count)) {
    stop("Manifest background_count mismatch during aggregation: ", row$query_id, call. = FALSE)
  }
  absent <- setdiff(genes, background)
  if (length(absent)) stop("Q is not a subset of B during aggregation for ", row$query_id, call. = FALSE)
  if (sha256_file(row$background_abs) != row$background_sha256) {
    stop("Background SHA256 mismatch during aggregation for ", row$query_id, call. = FALSE)
  }
  for (mode in required_modes) {
    qdir <- file.path(query_results_root, mode, row$query_id)
    full_path <- file.path(qdir, "gprofiler_full.tsv")
    sig_path <- file.path(qdir, "gprofiler_significant.tsv")
    fingerprint_path <- file.path(qdir, "query_mode_fingerprint.tsv")
    if (!file.exists(full_path)) stop("Missing runnable query-mode result: ", full_path, call. = FALSE)
    if (!file.exists(sig_path)) stop("Missing runnable query-mode significant result: ", sig_path, call. = FALSE)
    fp <- read_fingerprint(fingerprint_path)
    if (!identical(unname(fp["query_id"]), as.character(row$query_id))) {
      stop("Fingerprint query_id mismatch for ", row$query_id, " ", mode, call. = FALSE)
    }
    if (!identical(unname(fp["analysis_mode"]), mode)) {
      stop("Fingerprint analysis_mode mismatch for ", row$query_id, " ", mode, call. = FALSE)
    }
    if (!identical(unname(fp["query_file_sha256"]), sha256_file(row$gene_list_abs))) {
      stop("Query SHA256 mismatch for ", row$query_id, " ", mode, call. = FALSE)
    }
    if (!identical(unname(fp["background_file_sha256"]), row$background_sha256)) {
      stop("Background SHA256 fingerprint mismatch for ", row$query_id, " ", mode, call. = FALSE)
    }
    data_versions <- c(data_versions, unname(fp["gprofiler_annotation_data_version"]))
    terms <- fread(full_path, sep = "\t")
    if (nrow(terms)) {
      if (!all(c("query_id", "analysis_mode", "source", "term_id", "adjusted_enrichment_p_value", "significant") %in% names(terms))) {
        stop("Term result has non-canonical schema: ", full_path, call. = FALSE)
      }
      if (any(terms$query_id != row$query_id) || any(terms$analysis_mode != mode)) {
        stop("Term result query_id/analysis_mode mismatch: ", full_path, call. = FALSE)
      }
      manifest_meta <- row[, setdiff(names(row), c("gene_list_abs", "background_abs")), with = FALSE]
      for (nm in names(manifest_meta)) {
        if (!nm %in% names(terms)) terms[, (nm) := manifest_meta[[nm]][1]]
      }
      term_tables[[length(term_tables) + 1L]] <- terms
    }
    run_rows[[length(run_rows) + 1L]] <- data.table(
      query_id = row$query_id,
      analysis_mode = mode,
      result_dir = file.path("query_results", mode, row$query_id),
      full_tsv = file.path("query_results", mode, row$query_id, "gprofiler_full.tsv"),
      significant_tsv = file.path("query_results", mode, row$query_id, "gprofiler_significant.tsv"),
      fingerprint_tsv = file.path("query_results", mode, row$query_id, "query_mode_fingerprint.tsv"),
      query_gene_count = length(genes),
      background_gene_count = length(background),
      n_terms = if (exists("terms")) nrow(terms) else 0L,
      n_significant = if (exists("terms") && nrow(terms)) sum(terms$significant, na.rm = TRUE) else 0L,
      gprofiler_annotation_data_version = unname(fp["gprofiler_annotation_data_version"])
    )
    if (exists("terms")) rm(terms)
  }
}

data_versions <- unique(data_versions[nzchar(data_versions)])
if (length(data_versions) != 1L) {
  stop("Aggregated query modes do not share one g:Profiler annotation/data version: ",
       paste(data_versions, collapse = ";"), call. = FALSE)
}

all_terms <- if (length(term_tables)) rbindlist(term_tables, fill = TRUE) else data.table()
if (nrow(all_terms)) {
  key_cols <- c("query_id", "analysis_mode", "source", "term_id")
  if (anyDuplicated(all_terms[, ..key_cols])) {
    stop("Duplicate query-mode term outputs detected during aggregation", call. = FALSE)
  }
}
significant <- all_terms[significant == TRUE]
run_manifest <- if (length(run_rows)) rbindlist(run_rows, fill = TRUE) else data.table()

primary <- all_terms[analysis_mode == "primary"]
iea <- all_terms[analysis_mode == "iea_sensitivity"]
if (isTRUE(opt$`run-iea-sensitivity`)) {
  p <- primary[, .(
    query_id, source, term_id,
    term_name_primary = term_name,
    adjusted_enrichment_p_value_primary = adjusted_enrichment_p_value,
    intersection_genes_primary = intersection_genes
  )]
  e <- iea[, .(
    query_id, source, term_id,
    term_name_iea_excluded = term_name,
    adjusted_enrichment_p_value_iea_excluded = adjusted_enrichment_p_value,
    intersection_genes_iea_excluded = intersection_genes
  )]
  iea_join <- merge(p, e, by = c("query_id", "source", "term_id"), all = TRUE)
  iea_join[, present_primary := !is.na(adjusted_enrichment_p_value_primary)]
  iea_join[, present_iea_excluded := !is.na(adjusted_enrichment_p_value_iea_excluded)]
  iea_join[, significant_primary := present_primary & adjusted_enrichment_p_value_primary <= opt$alpha]
  iea_join[, significant_iea_excluded := present_iea_excluded & adjusted_enrichment_p_value_iea_excluded <= opt$alpha]
  iea_join[, changed_significance := significant_primary != significant_iea_excluded]
  iea_join[, presence_status := fifelse(
    present_primary & present_iea_excluded,
    "present_in_both",
    fifelse(present_primary, "primary_only", "iea_excluded_only")
  )]
  both <- iea_join[present_primary == TRUE & present_iea_excluded == TRUE]
  name_conflicts <- both[
    !is.na(term_name_primary) & !is.na(term_name_iea_excluded) &
      term_name_primary != term_name_iea_excluded
  ]
  if (nrow(name_conflicts)) {
    stop("Term-name conflicts in IEA sensitivity join", call. = FALSE)
  }
  meta_cols <- intersect(
    c("query_id", "query_family", "cancer_type", "attributed_cancer_type",
      "marker_evidence_stratum", "evidence_class", "contrast_id",
      "contrast_type", "direction", "direction_consistency_class",
      "background_id"),
    names(manifest)
  )
  iea_join <- merge(iea_join, unique(manifest[, ..meta_cols]), by = "query_id", all.x = TRUE)
} else {
  iea_join <- data.table()
}

run_summary <- data.table(
  metric = c(
    "manifest_queries_total",
    "runnable_queries_total",
    "required_analysis_modes",
    "primary_significant_terms_total",
    "primary_top_term_rows",
    "gprofiler_annotation_data_version",
    "aggregation_status"
  ),
  value = c(
    as.character(nrow(manifest)),
    as.character(nrow(runnable)),
    paste(required_modes, collapse = ","),
    as.character(nrow(significant[analysis_mode == "primary"])),
    "not_computed_by_aggregation_step",
    data_versions,
    "success"
  )
)

write_tsv(all_terms, file.path(outdir, "gprofiler_all_term_results.tsv"))
write_tsv(significant, file.path(outdir, "gprofiler_significant_term_results.tsv"))
write_tsv(run_manifest, file.path(outdir, "gprofiler_query_run_manifest.tsv"))
write_tsv(iea_join, file.path(outdir, "gprofiler_iea_sensitivity.tsv"))
write_tsv(run_summary, file.path(outdir, "gprofiler_run_summary.tsv"))
writeLines("success\n", file.path(outdir, ".success"))

message("[Aggregation] Manifest-driven g:Profiler aggregation completed.")
message("[IEA sensitivity] Full term join with explicit presence/significance status implemented.")
