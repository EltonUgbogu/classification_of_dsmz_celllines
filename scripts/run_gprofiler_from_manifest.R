#!/usr/bin/env Rscript

# Run g:Profiler enrichment for query sets listed in a query manifest.
#
# Each valid manifest row is submitted with its matched custom background. The
# manifest is a job table, not a gene list: gene IDs are read from the row-level
# gene_list_path/genes_path/genes_tsv file and backgrounds from the row-level
# background_path/background_tsv file.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(gprofiler2)
})

option_list <- list(
  make_option("--query-manifest", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--organism", type = "character", default = "hsapiens"),
  make_option("--sources", type = "character", default = "GO:BP,GO:MF,GO:CC,KEGG,REAC,WP,TF,HPA,CORUM"),
  make_option("--alpha", type = "numeric", default = 0.05),
  make_option("--correction-method", type = "character", default = "g_SCS"),
  make_option("--archive-url", type = "character", default = ""),
  make_option("--require-archive", type = "logical", default = TRUE),
  make_option("--run-iea-sensitivity", type = "logical", default = TRUE),
  make_option("--as-short-link", type = "logical", default = FALSE),
  make_option("--top-terms-per-source", type = "integer", default = 5),
  make_option("--query-delay-seconds", type = "numeric", default = 2),
  make_option("--retry-count", type = "integer", default = 1),
  make_option("--retry-delay-seconds", type = "numeric", default = 30),
  make_option("--resume", type = "logical", default = TRUE),
  make_option("--batch-start", type = "integer", default = 1),
  make_option("--batch-end", type = "integer", default = NA)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`query-manifest`) || is.null(opt$outdir)) {
  stop("--query-manifest and --outdir are required", call. = FALSE)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "primary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "iea_sensitivity"), recursive = TRUE, showWarnings = FALSE)

sources <- strsplit(opt$sources, ",", fixed = TRUE)[[1]]
sources <- sources[nzchar(sources)]

# Use a configured g:Profiler archive only when a fixed annotation snapshot is
# requested. When --archive-url is empty and --require-archive is FALSE, the
# runner uses the public endpoint.
if (nzchar(opt$`archive-url`)) {
  tryCatch(
    gprofiler2::set_base_url(opt$`archive-url`),
    error = function(e) {
      if (isTRUE(opt$`require-archive`)) {
        stop("Failed to set requested g:Profiler archive URL: ", conditionMessage(e), call. = FALSE)
      }
      warning("Failed to set requested g:Profiler archive URL: ", conditionMessage(e))
    }
  )
} else if (isTRUE(opt$`require-archive`)) {
  stop("require_archive is TRUE but --archive-url is empty", call. = FALSE)
}

as_flag <- function(x, default = FALSE) {
  if (length(x) == 0 || is.na(x)) return(default)
  value <- tolower(trimws(as.character(x)))
  if (!nzchar(value)) return(default)
  value %in% c("true", "t", "1", "yes", "y")
}

first_col <- function(dt, candidates, default = "") {
  for (candidate in candidates) {
    if (candidate %in% names(dt)) return(dt[[candidate]])
  }
  rep(default, nrow(dt))
}

coalesce_text <- function(...) {
  values <- list(...)
  out <- rep("", length(values[[1]]))
  for (value in values) {
    value <- as.character(value)
    value[is.na(value)] <- ""
    take <- !nzchar(out) & nzchar(value)
    out[take] <- value[take]
  }
  out
}

resolve_path <- function(path, manifest_dir) {
  path <- trimws(as.character(path))
  if (!nzchar(path)) return("")
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
  if (file.exists(path)) return(path)
  file.path(manifest_dir, path)
}

# Read one-column gene files. The first column is treated as gene_id when needed.
read_genes <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(character())
  x <- fread(path)
  if (!"gene_id" %in% names(x)) names(x)[1] <- "gene_id"
  genes <- trimws(as.character(x$gene_id))
  unique(genes[nzchar(genes) & !is.na(genes)])
}

# Ordered queries use ranked_genes.tsv sorted by decreasing rank_stat.
read_ranked <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(character())
  x <- fread(path)
  if (!"gene_id" %in% names(x)) names(x)[1] <- "gene_id"
  if ("rank_stat" %in% names(x)) setorder(x, -rank_stat)
  genes <- trimws(as.character(x$gene_id))
  unique(genes[nzchar(genes) & !is.na(genes)])
}

normalise_manifest <- function(path) {
  manifest <- fread(path)
  manifest_dir <- dirname(normalizePath(path, mustWork = FALSE))
  n <- nrow(manifest)
  if (n == 0) stop("query manifest has no rows: ", path, call. = FALSE)

  query_id <- coalesce_text(first_col(manifest, c("query_id", "query_name")), sprintf("query_%04d", seq_len(n)))
  gene_list_path <- coalesce_text(first_col(manifest, c("gene_list_path", "genes_path", "genes_tsv")))
  background_path <- coalesce_text(first_col(manifest, c("background_path", "background_tsv")))
  ranked_genes_path <- coalesce_text(first_col(manifest, c("ranked_genes_path", "ranked_genes_tsv")))

  gene_list_path <- vapply(gene_list_path, resolve_path, character(1), manifest_dir = manifest_dir)
  background_path <- vapply(background_path, resolve_path, character(1), manifest_dir = manifest_dir)
  ranked_genes_path <- vapply(ranked_genes_path, resolve_path, character(1), manifest_dir = manifest_dir)

  normalised <- data.table(
    query_id = query_id,
    query_name = coalesce_text(first_col(manifest, c("query_name", "query_id")), query_id),
    category = coalesce_text(first_col(manifest, c("category", "query_family"))),
    query_family = coalesce_text(first_col(manifest, c("query_family", "category"))),
    cohort = coalesce_text(first_col(manifest, c("cohort", "owner_profile", "disease", "profile"))),
    profile = coalesce_text(first_col(manifest, c("profile", "owner_profile", "cohort", "disease"))),
    disease = coalesce_text(first_col(manifest, c("disease", "owner_profile", "cohort", "profile"))),
    group_id = coalesce_text(first_col(manifest, c("group_id", "group", "category"))),
    contrast = coalesce_text(first_col(manifest, c("contrast", "source_contrast", "contrast_id"))),
    direction = coalesce_text(first_col(manifest, c("direction"))),
    source_marker_table = coalesce_text(first_col(manifest, c("source_marker_table_path", "source_table", "marker_file"))),
    gene_list_path = gene_list_path,
    background_path = background_path,
    ranked_genes_path = ranked_genes_path,
    rank_source = coalesce_text(first_col(manifest, c("rank_source"))),
    gene_count = suppressWarnings(as.integer(first_col(manifest, c("gene_count"), NA_character_))),
    background_count = suppressWarnings(as.integer(first_col(manifest, c("background_count"), NA_character_))),
    skip = vapply(first_col(manifest, c("skip"), "FALSE"), as_flag, logical(1), default = FALSE),
    ordered = vapply(first_col(manifest, c("ordered"), "FALSE"), as_flag, logical(1), default = FALSE),
    custom_background_available = vapply(first_col(manifest, c("custom_background_available"), "TRUE"), as_flag, logical(1), default = TRUE),
    background_strategy = coalesce_text(first_col(manifest, c("background_strategy"))),
    notes = coalesce_text(first_col(manifest, c("notes", "interpretation_hint"))),
    original_manifest_row = seq_len(n)
  )
  normalised[!nzchar(query_family), query_family := category]
  normalised[!nzchar(category), category := query_family]

  manifest_abs <- normalizePath(path, mustWork = FALSE)
  normalised[, skip_reason := ""]
  normalised[skip == TRUE, skip_reason := "skip_true_in_manifest"]
  normalised[!nzchar(gene_list_path), skip_reason := "missing_gene_list_path"]
  normalised[!nzchar(background_path), skip_reason := "missing_background_path"]
  normalised[nzchar(gene_list_path) & normalizePath(gene_list_path, mustWork = FALSE) == manifest_abs,
             skip_reason := "gene_list_path_points_to_query_manifest"]
  normalised[nzchar(background_path) & normalizePath(background_path, mustWork = FALSE) == manifest_abs,
             skip_reason := "background_path_points_to_query_manifest"]
  normalised[nzchar(gene_list_path) & !file.exists(gene_list_path), skip_reason := "gene_list_path_missing"]
  normalised[nzchar(background_path) & !file.exists(background_path), skip_reason := "background_path_missing"]
  normalised[nzchar(gene_list_path) & file.exists(gene_list_path) & file.info(gene_list_path)$size == 0,
             skip_reason := "empty_gene_list_file"]
  normalised[nzchar(background_path) & file.exists(background_path) & file.info(background_path)$size == 0,
             skip_reason := "empty_background_file"]
  normalised[is.na(gene_count), gene_count := vapply(gene_list_path, function(x) length(read_genes(x)), integer(1))]
  normalised[is.na(background_count), background_count := vapply(background_path, function(x) length(read_genes(x)), integer(1))]
  normalised[gene_count < 1, skip_reason := "zero_query_genes"]
  normalised[background_count < 1, skip_reason := "zero_background_genes"]
  normalised[!nzchar(skip_reason), skip := FALSE]
  normalised[nzchar(skip_reason), skip := TRUE]
  normalised
}

# g:Profiler is run with domain_scope = "custom" so adjusted p-values are
# calculated against the query-specific custom background: either the eligible
# marker-selection background or a contrast-specific DESeq2-tested background.
run_gost <- function(genes, background, ordered = FALSE, exclude_iea = FALSE) {
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
    as_short_link = opt$`as-short-link`
  )
}

# Stable empty schema keeps downstream aggregation readable when a query returns
# no terms.
empty_result <- function() {
  data.table(
    source = character(), term_id = character(), term_name = character(),
    p_value = numeric(), term_size = integer(), query_size = integer(),
    intersection_size = integer(), effective_domain_size = integer(),
    precision = numeric(), recall = numeric(), intersection = character()
  )
}

normalise_result <- function(gost_obj) {
  if (is.null(gost_obj) || is.null(gost_obj$result) || nrow(gost_obj$result) == 0) {
    return(empty_result())
  }
  x <- as.data.table(gost_obj$result)
  needed <- names(empty_result())
  for (nm in needed) if (!nm %in% names(x)) x[, (nm) := NA]
  x[, ..needed]
}

# Execute one manifest row and write the raw g:Profiler object, full term-level
# enrichment table, and significant-term table for the selected enrichment mode.
write_one <- function(row, mode, exclude_iea) {
  query_id <- row$query_id
  qdir <- file.path(opt$outdir, mode, query_id)
  dir.create(qdir, recursive = TRUE, showWarnings = FALSE)
  genes <- if (isTRUE(row$ordered) && nzchar(row$ranked_genes_path)) {
    read_ranked(row$ranked_genes_path)
  } else {
    read_genes(row$gene_list_path)
  }
  background <- read_genes(row$background_path)
  obj <- run_gost(genes, background, ordered = isTRUE(row$ordered), exclude_iea = exclude_iea)
  saveRDS(obj, file.path(qdir, "gprofiler_raw.rds"))
  res <- normalise_result(obj)
  res[, `:=`(
    query_id = query_id,
    query_name = row$query_name,
    category = row$category,
    query_family = row$query_family,
    cohort = row$cohort,
    contrast = row$contrast,
    profile = row$profile,
    disease = row$disease,
    group_id = row$group_id,
    direction = row$direction,
    source_marker_table = row$source_marker_table,
    gene_list_path = row$gene_list_path,
    background_path = row$background_path,
    ordered = row$ordered,
    rank_source = row$rank_source,
    gene_count = row$gene_count,
    background_count = row$background_count,
    iea_mode = ifelse(exclude_iea, "iea_excluded", "primary")
  )]
  res[, enrichment_ratio := precision / (term_size / effective_domain_size)]
  res[, leading_intersection_genes := vapply(strsplit(as.character(intersection), ","), function(z) {
    paste(head(z[nzchar(z)], 20), collapse = ",")
  }, character(1))]
  metadata_cols <- c("query_id", "query_name", "category", "query_family",
                     "cohort", "contrast", "profile", "disease", "group_id",
                     "direction", "source_marker_table", "gene_list_path",
                     "background_path", "ordered", "rank_source", "gene_count",
                     "background_count", "iea_mode")
  setcolorder(res, c(metadata_cols, setdiff(names(res), metadata_cols)))
  fwrite(res, file.path(qdir, "gprofiler_full.tsv"), sep = "\t")
  fwrite(res[p_value <= opt$alpha], file.path(qdir, "gprofiler_sig.tsv"), sep = "\t")
  res
}

corpus_schema <- c("query_id", "query_name", "category", "query_family",
                   "cohort", "contrast", "profile", "disease", "group_id",
                   "direction", "source_marker_table", "gene_list_path",
                   "background_path", "ordered", "rank_source", "iea_mode",
                   "n_terms", "n_significant", "gene_count", "background_count",
                   "full_tsv", "sig_tsv", "raw_rds", "status", "attempts",
                   "runtime_seconds", "error_message")

top_schema <- c("query_id", "query_name", "category", "query_family",
                "cohort", "contrast", "profile", "disease", "group_id",
                "direction", "source_marker_table", "gene_list_path",
                "background_path", "ordered", "rank_source", "iea_mode", "source",
                "term_id", "term_name", "p_value", "intersection_size",
                "query_size", "gene_count", "background_count",
                "term_size", "effective_domain_size",
                "precision", "recall", "enrichment_ratio", "intersection",
                "top_intersection_genes", "rank_within_query", "rank_within_source")

sens_schema <- c("query_id", "query_family", "category", "cohort", "contrast",
                 "direction", "source_marker_table", "source", "term_id",
                 "term_name", "p_primary", "p_iea_excluded",
                 "significant_primary", "significant_iea_excluded",
                 "delta_log10_p", "changed_significance",
                 "intersection_primary", "intersection_iea_excluded")

failed_schema <- c("query_id", "query_name", "category", "iea_mode", "status",
                   "attempts", "runtime_seconds", "error_message",
                   "gene_list_path", "background_path", "gene_count",
                   "background_count")

coerce_schema <- function(dt, schema) {
  if (nrow(dt) == 0 && length(names(dt)) == 0) {
    return(as.data.table(setNames(rep(list(logical()), length(schema)), schema)))
  }
  for (nm in schema) if (!nm %in% names(dt)) dt[, (nm) := NA]
  dt[, ..schema]
}

result_paths <- function(row, mode) {
  qdir <- file.path(opt$outdir, mode, row$query_id)
  list(
    full_tsv = file.path(qdir, "gprofiler_full.tsv"),
    sig_tsv = file.path(qdir, "gprofiler_sig.tsv"),
    raw_rds = file.path(qdir, "gprofiler_raw.rds")
  )
}

read_existing_result <- function(row, mode) {
  paths <- result_paths(row, mode)
  if (!isTRUE(opt$resume) || !file.exists(paths$full_tsv)) return(NULL)
  fread(paths$full_tsv)
}

corpus_entry <- function(row, mode, result, status, attempts, runtime_seconds, error_message) {
  paths <- result_paths(row, mode)
  data.table(
    query_id = row$query_id,
    query_name = row$query_name,
    category = row$category,
    query_family = row$query_family,
    cohort = row$cohort,
    contrast = row$contrast,
    profile = row$profile,
    disease = row$disease,
    group_id = row$group_id,
    direction = row$direction,
    source_marker_table = row$source_marker_table,
    gene_list_path = row$gene_list_path,
    background_path = row$background_path,
    ordered = row$ordered,
    rank_source = row$rank_source,
    iea_mode = ifelse(mode == "iea_sensitivity", "iea_excluded", "primary"),
    n_terms = nrow(result),
    n_significant = if (nrow(result) > 0 && "p_value" %in% names(result)) {
      sum(result$p_value <= opt$alpha, na.rm = TRUE)
    } else {
      0L
    },
    gene_count = row$gene_count,
    background_count = row$background_count,
    full_tsv = paths$full_tsv,
    sig_tsv = paths$sig_tsv,
    raw_rds = paths$raw_rds,
    status = status,
    attempts = attempts,
    runtime_seconds = runtime_seconds,
    error_message = error_message
  )
}

failed_entry <- function(row, mode, status, attempts, runtime_seconds, error_message) {
  data.table(
    query_id = row$query_id,
    query_name = row$query_name,
    category = row$category,
    iea_mode = ifelse(mode == "iea_sensitivity", "iea_excluded", "primary"),
    status = status,
    attempts = attempts,
    runtime_seconds = runtime_seconds,
    error_message = error_message,
    gene_list_path = row$gene_list_path,
    background_path = row$background_path,
    gene_count = row$gene_count,
    background_count = row$background_count
  )
}

run_with_retries <- function(row, mode, exclude_iea) {
  existing <- read_existing_result(row, mode)
  if (!is.null(existing)) {
    message(sprintf("  %s: reusing existing result", mode))
    return(list(result = existing, status = "resumed", attempts = 0L,
                runtime_seconds = 0, error_message = ""))
  }

  max_attempts <- max(1L, as.integer(opt$`retry-count`) + 1L)
  last_error <- ""
  total_start <- Sys.time()
  for (attempt in seq_len(max_attempts)) {
    attempt_start <- Sys.time()
    result <- tryCatch(
      write_one(row, mode, exclude_iea = exclude_iea),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(result)) {
      runtime <- round(as.numeric(difftime(Sys.time(), total_start, units = "secs")), 3)
      return(list(result = result, status = "success", attempts = attempt,
                  runtime_seconds = runtime, error_message = ""))
    }
    warning(sprintf("%s failed for %s on attempt %d/%d: %s",
                    mode, row$query_id, attempt, max_attempts, last_error))
    if (attempt < max_attempts && opt$`retry-delay-seconds` > 0) {
      Sys.sleep(opt$`retry-delay-seconds`)
    }
  }
  runtime <- round(as.numeric(difftime(Sys.time(), total_start, units = "secs")), 3)
  list(result = data.table(), status = "failed", attempts = max_attempts,
       runtime_seconds = runtime, error_message = last_error)
}

append_top_rows <- function(primary) {
  if (nrow(primary) == 0) return(NULL)
  top <- copy(primary[p_value <= opt$alpha])
  if (nrow(top) == 0) return(NULL)
  setorder(top, query_id, source, p_value)
  top[, rank_within_source := seq_len(.N), by = .(query_id, source, iea_mode)]
  top <- top[rank_within_source <= opt$`top-terms-per-source`]
  top[, top_intersection_genes := vapply(strsplit(as.character(intersection), ","), function(z) {
    paste(head(z[nzchar(z)], 20), collapse = ",")
  }, character(1))]
  top[, rank_within_query := frank(p_value, ties.method = "first"), by = .(query_id, iea_mode)]
  top
}

append_sensitivity_rows <- function(primary, iea) {
  if (nrow(primary) == 0 && nrow(iea) == 0) return(NULL)
  p <- primary[, .(query_id, query_family, direction, source, term_id, term_name,
                   category, cohort, contrast, source_marker_table,
                   p_primary = p_value, significant_primary = p_value <= opt$alpha,
                   intersection_primary = intersection)]
  e <- iea[, .(query_id, source, term_id,
               p_iea_excluded = p_value,
               significant_iea_excluded = p_value <= opt$alpha,
               intersection_iea_excluded = intersection)]
  joined <- merge(p, e, by = c("query_id", "source", "term_id"), all = TRUE)
  joined[, p_primary := as.numeric(p_primary)]
  joined[, p_iea_excluded := as.numeric(p_iea_excluded)]
  joined[, delta_log10_p := -log10(p_iea_excluded) - (-log10(p_primary))]
  joined[, changed_significance := significant_primary != significant_iea_excluded]
  joined[, .(
    query_id, query_family, category, cohort, contrast, direction,
    source_marker_table, source, term_id, term_name,
    p_primary, p_iea_excluded, significant_primary,
    significant_iea_excluded, delta_log10_p, changed_significance,
    intersection_primary, intersection_iea_excluded
  )]
}

manifest <- normalise_manifest(opt$`query-manifest`)
todo_all <- manifest[skip == FALSE]
todo_all[, runnable_index := .I]
batch_start <- max(1L, as.integer(opt$`batch-start`))
batch_end <- opt$`batch-end`
if (is.na(batch_end)) batch_end <- nrow(todo_all)
batch_end <- min(as.integer(batch_end), nrow(todo_all))
todo <- todo_all[runnable_index >= batch_start & runnable_index <= batch_end]
skipped <- manifest[skip == TRUE]

corpus_rows <- list()
top_rows <- list()
sens_rows <- list()
failed_rows <- list()

write_outputs <- function() {
  corpus <- if (length(corpus_rows)) rbindlist(corpus_rows, fill = TRUE) else data.table()
  tops <- if (length(top_rows)) rbindlist(top_rows, fill = TRUE) else data.table()
  sens <- if (length(sens_rows)) rbindlist(sens_rows, fill = TRUE) else data.table()
  failed <- if (length(failed_rows)) rbindlist(failed_rows, fill = TRUE) else data.table()

  corpus <- coerce_schema(corpus, corpus_schema)
  tops <- coerce_schema(tops, top_schema)
  sens <- coerce_schema(sens, sens_schema)
  failed <- coerce_schema(failed, failed_schema)

  fwrite(corpus, file.path(opt$outdir, "corpus_manifest.tsv"), sep = "\t")
  fwrite(tops, file.path(opt$outdir, "top_terms.tsv"), sep = "\t")
  fwrite(sens, file.path(opt$outdir, "iea_sensitivity_summary.tsv"), sep = "\t")
  fwrite(skipped, file.path(opt$outdir, "skipped_queries.tsv"), sep = "\t")
  fwrite(failed, file.path(opt$outdir, "failed_queries.tsv"), sep = "\t")

  run_summary <- data.table(
    metric = c(
      "manifest_queries",
      "manifest_runnable_queries",
      "batch_start",
      "batch_end",
      "submitted_queries",
      "skipped_queries",
      "failed_query_modes",
      "resumed_query_modes",
      "primary_queries_with_significant_terms",
      "primary_significant_terms",
      "iea_sensitivity_enabled",
      "alpha",
      "correction_method",
      "query_delay_seconds",
      "retry_count",
      "retry_delay_seconds",
      "resume_enabled",
      "batch_status"
    ),
    value = c(
      as.character(nrow(manifest)),
      as.character(nrow(todo_all)),
      as.character(batch_start),
      as.character(batch_end),
      as.character(nrow(todo)),
      as.character(nrow(skipped)),
      as.character(nrow(failed)),
      as.character(nrow(corpus[status == "resumed"])),
      as.character(uniqueN(corpus[iea_mode == "primary" & n_significant > 0, query_id])),
      as.character(nrow(tops[iea_mode == "primary"])),
      as.character(isTRUE(opt$`run-iea-sensitivity`)),
      as.character(opt$alpha),
      opt$`correction-method`,
      as.character(opt$`query-delay-seconds`),
      as.character(opt$`retry-count`),
      as.character(opt$`retry-delay-seconds`),
      as.character(isTRUE(opt$resume)),
      if (nrow(failed) > 0) "completed_with_failed_query_modes" else "completed"
    )
  )
  fwrite(run_summary, file.path(opt$outdir, "run_summary.tsv"), sep = "\t")

  version_info <- data.table(
    item = c("gprofiler2_version", "organism", "sources", "archive_url",
             "correction_method", "alpha", "as_short_link"),
    value = c(as.character(utils::packageVersion("gprofiler2")), opt$organism,
              paste(sources, collapse = ","), opt$`archive-url`,
              opt$`correction-method`, as.character(opt$alpha),
              as.character(opt$`as-short-link`))
  )
  fwrite(version_info, file.path(opt$outdir, "gprofiler_version.tsv"), sep = "\t")
}

write_outputs()

# Submit each runnable query serially. Completed query files are reused when
# resume is enabled, and batch-level outputs are refreshed after every row.
for (i in seq_len(nrow(todo))) {
  row <- todo[i]
  message(sprintf("[%d/%d runnable=%d] %s", i, nrow(todo), row$runnable_index, row$query_id))

  primary_run <- run_with_retries(row, "primary", exclude_iea = FALSE)
  primary <- primary_run$result
  corpus_rows[[length(corpus_rows) + 1]] <- corpus_entry(
    row, "primary", primary, primary_run$status, primary_run$attempts,
    primary_run$runtime_seconds, primary_run$error_message
  )
  if (identical(primary_run$status, "failed")) {
    failed_rows[[length(failed_rows) + 1]] <- failed_entry(
      row, "primary", primary_run$status, primary_run$attempts,
      primary_run$runtime_seconds, primary_run$error_message
    )
  }

  top <- append_top_rows(primary)
  if (!is.null(top)) top_rows[[length(top_rows) + 1]] <- top

  if (isTRUE(opt$`run-iea-sensitivity`)) {
    if (identical(primary_run$status, "failed")) {
      iea_run <- list(result = data.table(), status = "not_run_primary_failed",
                      attempts = 0L, runtime_seconds = 0,
                      error_message = primary_run$error_message)
    } else {
      iea_run <- run_with_retries(row, "iea_sensitivity", exclude_iea = TRUE)
    }
    iea <- iea_run$result
    corpus_rows[[length(corpus_rows) + 1]] <- corpus_entry(
      row, "iea_sensitivity", iea, iea_run$status, iea_run$attempts,
      iea_run$runtime_seconds, iea_run$error_message
    )
    if (identical(iea_run$status, "failed")) {
      failed_rows[[length(failed_rows) + 1]] <- failed_entry(
        row, "iea_sensitivity", iea_run$status, iea_run$attempts,
        iea_run$runtime_seconds, iea_run$error_message
      )
    }
    sens <- append_sensitivity_rows(primary, iea)
    if (!is.null(sens)) sens_rows[[length(sens_rows) + 1]] <- sens
  }

  write_outputs()
  if (i < nrow(todo) && opt$`query-delay-seconds` > 0) {
    Sys.sleep(opt$`query-delay-seconds`)
  }
}

write_outputs()

writeLines(c(
  "g:Profiler enrichment report",
  "",
  "Primary interpretation files:",
  "  top_terms.tsv: significant terms, capped per source and query for review.",
  "  corpus_manifest.tsv: one row per submitted query and mode, with result counts and file paths.",
  "  iea_sensitivity_summary.tsv: compares primary results with electronic annotations excluded.",
  "  run_summary.tsv: batch-level counts and statistical settings.",
  "  primary/<query_id>/gprofiler_full.tsv: all returned terms for one query.",
  "  primary/<query_id>/gprofiler_sig.tsv: significant terms for one query.",
  "",
  "Column notes:",
  "  p_value is g:Profiler's adjusted p-value under the selected g:SCS multiple-testing correction.",
  "  enrichment_ratio = precision / (term_size / effective_domain_size).",
  "  precision = intersection_size / query_size.",
  "  recall = intersection_size / term_size.",
  "  leading_intersection_genes/top_intersection_genes list up to 20 query genes driving a term.",
  "",
  "Scientific notes:",
  "  query_manifest.tsv is treated as a functional-enrichment job manifest, not as a gene list.",
  "  Each runnable query uses its row-level gene_list_path/genes_path/genes_tsv file.",
  "  Each runnable query uses its row-level background_path/background_tsv file as the query-specific custom background.",
  "  category, cohort, contrast, direction, and source_marker_table metadata are preserved in tabular outputs.",
  "  ordered=TRUE uses ranked_genes.tsv and ordered_query=TRUE in g:Profiler.",
  "  Primary enrichment mode includes IEA annotations; IEA-sensitivity mode reruns the query with exclude_iea = TRUE."
), file.path(opt$outdir, "README_gprofiler_results.txt"))
