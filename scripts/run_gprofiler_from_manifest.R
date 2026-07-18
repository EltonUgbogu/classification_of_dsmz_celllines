#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(gprofiler2)
})

option_list <- list(
  make_option("--pipeline-root", type = "character"),
  make_option("--query-manifest", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--organism", type = "character"),
  make_option("--sources", type = "character"),
  make_option("--alpha", type = "numeric"),
  make_option("--correction-method", type = "character"),
  make_option("--endpoint", type = "character"),
  make_option("--primary-exclude-iea", type = "logical"),
  make_option("--run-iea-sensitivity", type = "logical"),
  make_option("--reuse-matching-query-results", type = "logical"),
  make_option("--fail-on-query-error", type = "logical"),
  make_option("--as-short-link", type = "logical", default = FALSE),
  make_option("--query-delay-seconds", type = "numeric", default = 2),
  make_option("--retry-count", type = "integer", default = 1),
  make_option("--retry-delay-seconds", type = "numeric", default = 30),
  make_option("--batch-start", type = "integer", default = 1),
  make_option("--batch-end", type = "integer", default = NA)
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("pipeline-root", "query-manifest", "outdir", "organism",
              "sources", "alpha", "correction-method", "endpoint",
              "primary-exclude-iea", "run-iea-sensitivity",
              "reuse-matching-query-results", "fail-on-query-error")
for (nm in required) {
  if (is.null(opt[[nm]]) || (is.character(opt[[nm]]) && !nzchar(opt[[nm]]))) {
    stop("--", nm, " is required", call. = FALSE)
  }
}

pipeline_root <- normalizePath(opt$`pipeline-root`, mustWork = TRUE)
query_manifest_path <- normalizePath(opt$`query-manifest`, mustWork = TRUE)
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
query_results_root <- file.path(outdir, "query_results")
batch_root <- file.path(outdir, "batches")
dir.create(query_results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(batch_root, recursive = TRUE, showWarnings = FALSE)

sources <- strsplit(opt$sources, ",", fixed = TRUE)[[1]]
sources <- sources[nzchar(sources)]
if (!length(sources)) stop("--sources must contain at least one g:Profiler source", call. = FALSE)

if (!identical(opt$endpoint, "https://biit.cs.ut.ee/gprofiler/api") && nzchar(opt$endpoint)) {
  tryCatch(
    gprofiler2::set_base_url(opt$endpoint),
    error = function(e) stop("Failed to set g:Profiler endpoint: ", conditionMessage(e), call. = FALSE)
  )
}

is_true <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

rel_to_abs <- function(path) {
  value <- trimws(as.character(path))
  if (!nzchar(value)) return("")
  if (grepl("^/", value)) return(value)
  file.path(pipeline_root, value)
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

read_genes <- function(path) {
  if (!file.exists(path)) stop("Missing gene file: ", path, call. = FALSE)
  values <- readLines(path, warn = FALSE)
  values <- trimws(sub("\t.*$", "", values))
  values <- values[nzchar(values) & values != "gene_id"]
  unique(values)
}

read_ranked_genes <- function(path) {
  if (!file.exists(path)) stop("Missing ranked gene file: ", path, call. = FALSE)
  x <- fread(path, sep = "\t", data.table = FALSE)
  if (!"gene_id" %in% names(x)) names(x)[1] <- "gene_id"
  if ("rank_stat" %in% names(x)) {
    x$rank_stat <- suppressWarnings(as.numeric(x$rank_stat))
    x <- x[order(-x$rank_stat, x$gene_id), , drop = FALSE]
  }
  unique(trimws(as.character(x$gene_id)))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", na = "NA")
}

empty_terms <- function() {
  data.table(
    query_id = character(),
    analysis_mode = character(),
    source = character(),
    term_id = character(),
    term_name = character(),
    adjusted_enrichment_p_value = numeric(),
    significant = logical(),
    intersection_size = integer(),
    query_size = integer(),
    background_count = integer(),
    term_size = integer(),
    effective_domain_size = integer(),
    precision = numeric(),
    recall = numeric(),
    enrichment_ratio = numeric(),
    intersection_genes = character(),
    rank_within_query = integer(),
    rank_within_source = integer()
  )
}

normalise_result <- function(obj, row, analysis_mode) {
  if (is.null(obj) || is.null(obj$result) || nrow(obj$result) == 0L) {
    return(empty_terms())
  }
  x <- as.data.table(obj$result)
  needed <- c("source", "term_id", "term_name", "p_value", "intersection_size",
              "query_size", "term_size", "effective_domain_size", "precision",
              "recall", "intersection")
  for (nm in needed) {
    if (!nm %in% names(x)) x[, (nm) := NA]
  }
  x[, adjusted_enrichment_p_value := as.numeric(p_value)]
  x[, significant := !is.na(adjusted_enrichment_p_value) & adjusted_enrichment_p_value <= opt$alpha]
  x[, background_count := as.integer(row$background_count)]
  x[, enrichment_ratio := precision / (term_size / effective_domain_size)]
  x[, intersection_genes := as.character(intersection)]
  x[, query_id := row$query_id]
  x[, analysis_mode := analysis_mode]
  setorder(
    x,
    query_id,
    adjusted_enrichment_p_value,
    -enrichment_ratio,
    -intersection_size,
    term_id
  )
  x[, rank_within_query := seq_len(.N), by = .(query_id, analysis_mode)]
  setorder(
    x,
    query_id,
    analysis_mode,
    source,
    adjusted_enrichment_p_value,
    -enrichment_ratio,
    -intersection_size,
    term_id
  )
  x[, rank_within_source := seq_len(.N), by = .(query_id, analysis_mode, source)]
  cols <- names(empty_terms())
  x[, ..cols]
}

extract_data_version <- function(obj) {
  candidates <- character()
  if (!is.null(obj) && !is.null(obj$meta)) {
    meta <- obj$meta
    for (nm in c("data_version", "version", "release", "timestamp", "organism")) {
      if (!is.null(meta[[nm]]) && length(meta[[nm]]) > 0) {
        candidates <- c(candidates, paste(nm, paste(as.character(meta[[nm]]), collapse = ","), sep = "="))
      }
    }
  }
  value <- paste(unique(candidates[nzchar(candidates)]), collapse = ";")
  if (nzchar(value)) value else "unavailable_from_gprofiler_response"
}

get_current_gprofiler_data_version <- function() {
  ns <- asNamespace("gprofiler2")
  for (fn in c("get_version_info", "get_version")) {
    if (exists(fn, envir = ns, inherits = FALSE)) {
      value <- tryCatch(get(fn, envir = ns)(), error = function(e) NULL)
      if (!is.null(value)) return(paste(capture.output(str(value)), collapse = " "))
    }
  }
  ""
}

json_string <- function(x) {
  paste0("{", paste(sprintf('"%s":"%s"', names(x), gsub('"', '\\"', as.character(x))), collapse = ","), "}")
}

fingerprint_for <- function(row, analysis_mode, exclude_iea, data_version) {
  ranked_path <- rel_to_abs(row$ranked_genes_path)
  ranked_sha <- if (nzchar(ranked_path) && file.exists(ranked_path)) sha256_file(ranked_path) else ""
  c(
    query_id = row$query_id,
    analysis_mode = analysis_mode,
    query_file_sha256 = sha256_file(rel_to_abs(row$gene_list_path)),
    background_file_sha256 = sha256_file(rel_to_abs(row$background_path)),
    ranked_file_sha256 = ranked_sha,
    query_gene_count = as.character(row$gene_count),
    background_gene_count = as.character(row$background_count),
    ordered = as.character(is_true(row$ordered)),
    ranking_used_for_enrichment = as.character(is_true(row$ranking_used_for_enrichment)),
    organism = opt$organism,
    sources = paste(sources, collapse = ","),
    correction_method = opt$`correction-method`,
    alpha = as.character(opt$alpha),
    endpoint = opt$endpoint,
    iea_mode = ifelse(exclude_iea, "iea_excluded", "iea_included"),
    gprofiler2_package_version = as.character(utils::packageVersion("gprofiler2")),
    gprofiler_annotation_data_version = data_version
  )
}

write_fingerprint <- function(fingerprint, path) {
  write_tsv(data.table(field = names(fingerprint), value = unname(fingerprint)), path)
}

read_fingerprint <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- fread(path, sep = "\t")
  if (!all(c("field", "value") %in% names(x))) return(NULL)
  stats::setNames(as.character(x$value), as.character(x$field))
}

fingerprints_equal <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  identical(a[names(b)], b)
}

validate_manifest <- function(manifest) {
  required <- c("query_id", "query_family", "gene_count", "gene_list_path",
                "background_count", "background_path", "background_sha256",
                "ordered", "ranking_used_for_enrichment", "ranked_genes_path",
                "query_execution_status")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) stop("Query manifest missing required columns: ", paste(missing, collapse = ","), call. = FALSE)
  if (anyDuplicated(manifest$query_id)) {
    stop("Query manifest contains duplicate query_id values", call. = FALSE)
  }
  manifest[, gene_list_abs := vapply(gene_list_path, rel_to_abs, character(1))]
  manifest[, background_abs := vapply(background_path, rel_to_abs, character(1))]
  manifest[, ranked_abs := vapply(ranked_genes_path, rel_to_abs, character(1))]
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i]
    genes <- read_genes(row$gene_list_abs)
    background <- read_genes(row$background_abs)
    if (length(genes) != as.integer(row$gene_count)) {
      stop("Manifest gene_count mismatch for ", row$query_id, call. = FALSE)
    }
    if (length(background) != as.integer(row$background_count)) {
      stop("Manifest background_count mismatch for ", row$query_id, call. = FALSE)
    }
    absent <- setdiff(genes, background)
    if (length(absent)) {
      stop("Q is not a subset of B for ", row$query_id, ": ", paste(head(absent, 50), collapse = ","), call. = FALSE)
    }
    if (sha256_file(row$background_abs) != row$background_sha256) {
      stop("background SHA256 mismatch for ", row$query_id, call. = FALSE)
    }
    if (is_true(row$ordered) || is_true(row$ranking_used_for_enrichment)) {
      if (!nzchar(row$ranked_abs) || !file.exists(row$ranked_abs)) {
        stop("Ordered/ranked query lacks ranked_genes_path: ", row$query_id, call. = FALSE)
      }
    }
  }
  manifest
}

run_gost <- function(genes, background, ordered, exclude_iea) {
  gprofiler2::gost(
    query = genes,
    organism = opt$organism,
    ordered_query = ordered,
    multi_query = FALSE,
    significant = FALSE,
    exclude_iea = exclude_iea,
    evcodes = TRUE,
    user_threshold = opt$alpha,
    correction_method = opt$`correction-method`,
    domain_scope = "custom",
    custom_bg = background,
    sources = sources,
    as_short_link = isTRUE(opt$`as-short-link`)
  )
}

write_query_mode <- function(row, analysis_mode, exclude_iea, current_data_version) {
  qdir <- file.path(query_results_root, analysis_mode, row$query_id)
  fingerprint_path <- file.path(qdir, "query_mode_fingerprint.tsv")
  full_path <- file.path(qdir, "gprofiler_full.tsv")
  sig_path <- file.path(qdir, "gprofiler_significant.tsv")
  provenance_path <- file.path(qdir, "query_mode_provenance.tsv")
  expected_fingerprint <- fingerprint_for(row, analysis_mode, exclude_iea, current_data_version)
  existing_fingerprint <- read_fingerprint(fingerprint_path)
  can_consider_reuse <- isTRUE(opt$`reuse-matching-query-results`) &&
    nzchar(current_data_version) &&
    !identical(current_data_version, "unavailable_from_gprofiler_response")
  if (can_consider_reuse && file.exists(full_path) && file.exists(sig_path) &&
      fingerprints_equal(existing_fingerprint, expected_fingerprint)) {
    result <- fread(full_path, sep = "\t")
    return(list(status = "reused_matching_fingerprint", result = result, data_version = current_data_version, error = ""))
  }
  if (dir.exists(qdir)) unlink(qdir, recursive = TRUE, force = TRUE)
  tmpdir <- paste0(qdir, ".tmp.", Sys.getpid())
  if (dir.exists(tmpdir)) unlink(tmpdir, recursive = TRUE, force = TRUE)
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  genes <- if (is_true(row$ordered) && nzchar(row$ranked_abs)) {
    read_ranked_genes(row$ranked_abs)
  } else {
    read_genes(row$gene_list_abs)
  }
  background <- read_genes(row$background_abs)
  obj <- run_gost(genes, background, is_true(row$ordered), exclude_iea)
  saveRDS(obj, file.path(tmpdir, "gprofiler_raw.rds"))
  data_version <- extract_data_version(obj)
  result <- normalise_result(obj, row, analysis_mode)
  write_tsv(result, file.path(tmpdir, "gprofiler_full.tsv"))
  write_tsv(result[significant == TRUE], file.path(tmpdir, "gprofiler_significant.tsv"))
  fingerprint <- fingerprint_for(row, analysis_mode, exclude_iea, data_version)
  write_fingerprint(fingerprint, file.path(tmpdir, "query_mode_fingerprint.tsv"))
  provenance <- data.table(
    field = c("query_id", "analysis_mode", "endpoint", "gprofiler2_package_version",
              "gprofiler_annotation_data_version", "organism", "sources",
              "correction_method", "alpha", "exclude_iea", "custom_background_strategy",
              "run_timestamp"),
    value = c(row$query_id, analysis_mode, opt$endpoint,
              as.character(utils::packageVersion("gprofiler2")), data_version,
              opt$organism, paste(sources, collapse = ","), opt$`correction-method`,
              as.character(opt$alpha), as.character(exclude_iea),
              "query_specific_custom_analytical_background",
              format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  )
  write_tsv(provenance, file.path(tmpdir, "query_mode_provenance.tsv"))
  file.rename(tmpdir, qdir)
  list(status = "success", result = result, data_version = data_version, error = "")
}

manifest <- validate_manifest(fread(query_manifest_path, sep = "\t"))
runnable <- manifest[query_execution_status == "runnable"]
runnable[, runnable_index := seq_len(.N)]
batch_start <- max(1L, as.integer(opt$`batch-start`))
requested_batch_end <- opt$`batch-end`
if (is.na(requested_batch_end)) requested_batch_end <- nrow(runnable)
effective_batch_end <- min(as.integer(requested_batch_end), nrow(runnable))
todo <- runnable[runnable_index >= batch_start & runnable_index <= effective_batch_end]
batch_dir <- file.path(batch_root, sprintf("batch_%04d_%04d", batch_start, as.integer(requested_batch_end)))
dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)

current_data_version <- get_current_gprofiler_data_version()
if (!nzchar(current_data_version)) current_data_version <- "unavailable_from_gprofiler_response"

run_rows <- list()
failed_rows <- list()
required_modes <- c("primary")
if (isTRUE(opt$`run-iea-sensitivity`)) required_modes <- c(required_modes, "iea_sensitivity")

write_batch_records <- function() {
  run_manifest <- if (length(run_rows)) rbindlist(run_rows, fill = TRUE) else data.table()
  failed <- if (length(failed_rows)) rbindlist(failed_rows, fill = TRUE) else data.table()
  write_tsv(run_manifest, file.path(batch_dir, "query_mode_run_manifest.tsv"))
  write_tsv(failed, file.path(batch_dir, "failed_query_modes.tsv"))
}

write_batch_records()
for (i in seq_len(nrow(todo))) {
  row <- todo[i]
  message(sprintf("[g:Profiler execution] %d/%d %s", i, nrow(todo), row$query_id))
  for (mode in required_modes) {
    exclude_iea <- if (mode == "primary") isTRUE(opt$`primary-exclude-iea`) else TRUE
    started <- Sys.time()
    status <- "failed"
    result <- data.table()
    data_version <- ""
    error_message <- ""
    attempts <- max(1L, as.integer(opt$`retry-count`) + 1L)
    for (attempt in seq_len(attempts)) {
      run <- tryCatch(
        write_query_mode(row, mode, exclude_iea, current_data_version),
        error = function(e) list(status = "failed", result = data.table(), data_version = "", error = conditionMessage(e))
      )
      status <- run$status
      result <- run$result
      data_version <- run$data_version
      error_message <- run$error
      if (!identical(status, "failed")) break
      if (attempt < attempts && opt$`retry-delay-seconds` > 0) Sys.sleep(opt$`retry-delay-seconds`)
    }
    runtime <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 3)
    run_rows[[length(run_rows) + 1L]] <- data.table(
      query_id = row$query_id,
      analysis_mode = mode,
      status = status,
      n_terms = nrow(result),
      n_significant = if (nrow(result)) sum(result$significant, na.rm = TRUE) else 0L,
      attempts = attempts,
      runtime_seconds = runtime,
      gprofiler_annotation_data_version = data_version,
      result_dir = file.path("query_results", mode, row$query_id),
      error_message = error_message
    )
    if (identical(status, "failed")) {
      failed_rows[[length(failed_rows) + 1L]] <- data.table(
        query_id = row$query_id,
        analysis_mode = mode,
        status = status,
        error_message = error_message
      )
    }
    write_batch_records()
  }
  if (i < nrow(todo) && opt$`query-delay-seconds` > 0) Sys.sleep(opt$`query-delay-seconds`)
}

write_batch_records()
if (length(failed_rows) && isTRUE(opt$`fail-on-query-error`)) {
  stop("One or more runnable g:Profiler query modes failed; see ", file.path(batch_dir, "failed_query_modes.tsv"), call. = FALSE)
}

writeLines("success\n", file.path(batch_dir, ".success"))
message("[g:Profiler execution] Batch completed without failed required query modes.")
