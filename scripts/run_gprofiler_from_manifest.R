#!/usr/bin/env Rscript

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
  make_option("--top-terms-per-source", type = "integer", default = 5)
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

read_genes <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(character())
  x <- fread(path)
  if (!"gene_id" %in% names(x)) names(x)[1] <- "gene_id"
  unique(x$gene_id[nzchar(x$gene_id)])
}

read_ranked <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(character())
  x <- fread(path)
  if (!"gene_id" %in% names(x)) names(x)[1] <- "gene_id"
  if ("rank_stat" %in% names(x)) setorder(x, -rank_stat)
  unique(x$gene_id[nzchar(x$gene_id)])
}

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

write_one <- function(row, mode, exclude_iea) {
  query_id <- row$query_id
  qdir <- file.path(opt$outdir, mode, query_id)
  dir.create(qdir, recursive = TRUE, showWarnings = FALSE)
  genes <- if (isTRUE(row$ordered)) read_ranked(row$ranked_genes_tsv) else read_genes(row$genes_tsv)
  background <- read_genes(row$background_tsv)
  obj <- run_gost(genes, background, ordered = isTRUE(row$ordered), exclude_iea = exclude_iea)
  saveRDS(obj, file.path(qdir, "gprofiler_raw.rds"))
  res <- normalise_result(obj)
  res[, `:=`(
    query_id = query_id,
    query_family = row$query_family,
    profile = row$profile,
    disease = row$disease,
    group_id = row$group_id,
    direction = row$direction,
    ordered = row$ordered,
    rank_source = row$rank_source,
    iea_mode = ifelse(exclude_iea, "iea_excluded", "primary")
  )]
  setcolorder(res, c("query_id", "query_family", "profile", "disease", "group_id",
                     "direction", "ordered", "rank_source", "iea_mode",
                     setdiff(names(res), c("query_id", "query_family", "profile",
                                           "disease", "group_id", "direction",
                                           "ordered", "rank_source", "iea_mode"))))
  fwrite(res, file.path(qdir, "gprofiler_full.tsv"), sep = "\t")
  fwrite(res[p_value <= opt$alpha], file.path(qdir, "gprofiler_sig.tsv"), sep = "\t")
  res
}

manifest <- fread(opt$`query-manifest`)
todo <- manifest[skip == FALSE]

corpus_rows <- list()
top_rows <- list()
sens_rows <- list()

for (i in seq_len(nrow(todo))) {
  row <- todo[i]
  message(sprintf("[%d/%d] %s", i, nrow(todo), row$query_id))
  primary <- tryCatch(
    write_one(row, "primary", exclude_iea = FALSE),
    error = function(e) {
      warning("g:Profiler failed for ", row$query_id, ": ", conditionMessage(e))
      data.table()
    }
  )
  corpus_rows[[length(corpus_rows) + 1]] <- data.table(
    query_id = row$query_id,
    query_family = row$query_family,
    profile = row$profile,
    disease = row$disease,
    group_id = row$group_id,
    direction = row$direction,
    ordered = row$ordered,
    rank_source = row$rank_source,
    iea_mode = "primary",
    n_terms = nrow(primary),
    n_significant = if (nrow(primary) > 0) sum(primary$p_value <= opt$alpha, na.rm = TRUE) else 0,
    full_tsv = file.path(opt$outdir, "primary", row$query_id, "gprofiler_full.tsv"),
    sig_tsv = file.path(opt$outdir, "primary", row$query_id, "gprofiler_sig.tsv"),
    raw_rds = file.path(opt$outdir, "primary", row$query_id, "gprofiler_raw.rds")
  )
  if (nrow(primary) > 0) {
    top <- copy(primary[p_value <= opt$alpha])
    if (nrow(top) > 0) {
      setorder(top, query_id, source, p_value)
      top[, rank_within_source := seq_len(.N), by = .(query_id, source, iea_mode)]
      top <- top[rank_within_source <= opt$`top-terms-per-source`]
      top[, enrichment_ratio := precision / (term_size / effective_domain_size)]
      top[, top_intersection_genes := vapply(strsplit(as.character(intersection), ","), function(z) {
        paste(head(z[nzchar(z)], 20), collapse = ",")
      }, character(1))]
      top[, rank_within_query := frank(p_value, ties.method = "first"), by = .(query_id, iea_mode)]
      top_rows[[length(top_rows) + 1]] <- top
    }
  }

  if (isTRUE(opt$`run-iea-sensitivity`)) {
    iea <- tryCatch(
      write_one(row, "iea_sensitivity", exclude_iea = TRUE),
      error = function(e) {
        warning("IEA sensitivity failed for ", row$query_id, ": ", conditionMessage(e))
        data.table()
      }
    )
    corpus_rows[[length(corpus_rows) + 1]] <- data.table(
      query_id = row$query_id,
      query_family = row$query_family,
      profile = row$profile,
      disease = row$disease,
      group_id = row$group_id,
      direction = row$direction,
      ordered = row$ordered,
      rank_source = row$rank_source,
      iea_mode = "iea_excluded",
      n_terms = nrow(iea),
      n_significant = if (nrow(iea) > 0) sum(iea$p_value <= opt$alpha, na.rm = TRUE) else 0,
      full_tsv = file.path(opt$outdir, "iea_sensitivity", row$query_id, "gprofiler_full.tsv"),
      sig_tsv = file.path(opt$outdir, "iea_sensitivity", row$query_id, "gprofiler_sig.tsv"),
      raw_rds = file.path(opt$outdir, "iea_sensitivity", row$query_id, "gprofiler_raw.rds")
    )
    if (nrow(primary) > 0 || nrow(iea) > 0) {
      p <- primary[, .(query_id, query_family, direction, source, term_id, term_name,
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
      sens_rows[[length(sens_rows) + 1]] <- joined[, .(
        query_id, query_family, direction, source, term_id, term_name,
        p_primary, p_iea_excluded, significant_primary,
        significant_iea_excluded, delta_log10_p, changed_significance,
        intersection_primary, intersection_iea_excluded
      )]
    }
  }
}

corpus <- if (length(corpus_rows)) rbindlist(corpus_rows, fill = TRUE) else data.table()
tops <- if (length(top_rows)) rbindlist(top_rows, fill = TRUE) else data.table()
sens <- if (length(sens_rows)) rbindlist(sens_rows, fill = TRUE) else data.table()

top_schema <- c("query_id", "query_family", "profile", "disease", "group_id",
                "direction", "ordered", "rank_source", "iea_mode", "source",
                "term_id", "term_name", "p_value", "intersection_size",
                "query_size", "term_size", "effective_domain_size",
                "precision", "recall", "enrichment_ratio", "intersection",
                "top_intersection_genes", "rank_within_query", "rank_within_source")
for (nm in top_schema) if (!nm %in% names(tops)) tops[, (nm) := NA]
tops <- tops[, ..top_schema]

sens_schema <- c("query_id", "query_family", "direction", "source", "term_id",
                 "term_name", "p_primary", "p_iea_excluded",
                 "significant_primary", "significant_iea_excluded",
                 "delta_log10_p", "changed_significance",
                 "intersection_primary", "intersection_iea_excluded")
for (nm in sens_schema) if (!nm %in% names(sens)) sens[, (nm) := NA]
sens <- sens[, ..sens_schema]

fwrite(corpus, file.path(opt$outdir, "corpus_manifest.tsv"), sep = "\t")
fwrite(tops, file.path(opt$outdir, "top_terms.tsv"), sep = "\t")
fwrite(sens, file.path(opt$outdir, "iea_sensitivity_summary.tsv"), sep = "\t")

version_info <- data.table(
  item = c("gprofiler2_version", "organism", "sources", "archive_url",
           "correction_method", "alpha", "as_short_link"),
  value = c(as.character(utils::packageVersion("gprofiler2")), opt$organism,
            paste(sources, collapse = ","), opt$`archive-url`,
            opt$`correction-method`, as.character(opt$alpha),
            as.character(opt$`as-short-link`))
)
fwrite(version_info, file.path(opt$outdir, "gprofiler_version.tsv"), sep = "\t")
