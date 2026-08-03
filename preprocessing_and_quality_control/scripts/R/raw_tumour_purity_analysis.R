options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(tidyestimate)
})

for (path in unname(snakemake@output)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

## Shared figure module: palette, typography, gridless panel frame, threshold
## line and the atomic PDF device used by the purity figures further down.
##
## Sourced here, before this script's own helpers, deliberately. Both files
## define `sha256_file()`; loading the module first means the definition below
## wins, so the provenance digests keep the exact implementation they have
## always used. The module contributes presentation only.
source(snakemake@input[["shared_figure_module"]])

log_path <- snakemake@output[["execution_log"]]
log_connection <- file(log_path, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_connection)
}, add = TRUE)

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

save_rds_atomic <- function(object, path) {
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Failed to atomically promote RDS output: ", path)
  }
  invisible(path)
}

sha256_file <- function(path) {
  if (is.null(path) || !length(path) || !nzchar(path) || !file.exists(path)) {
    return(NA_character_)
  }
  result <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path, ": ", paste(result, collapse = " "))
  }
  strsplit(result[[1]], "[[:space:]]+")[[1]][1]
}

read_count_matrix <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
  object <- readRDS(path)
  if (methods::is(object, "SummarizedExperiment")) {
    object <- SummarizedExperiment::assay(object)
  } else if (is.list(object) && !is.null(object$counts)) {
    object <- object$counts
  }
  if (!is.matrix(object) && !is.data.frame(object)) {
    stop(
      label,
      " must be a matrix, data.frame, SummarizedExperiment, or list(counts=); class=",
      paste(class(object), collapse = ",")
    )
  }
  counts <- as.matrix(object)
  if (!is.numeric(counts)) stop(label, " must be numeric")
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop(label, " must have row and column names")
  }
  storage.mode(counts) <- "double"
  counts
}

validate_raw_counts <- function(counts, label, chunk_rows = 1000L) {
  failures <- character(0)
  if (anyDuplicated(rownames(counts))) failures <- c(failures, "duplicated gene identifiers")
  if (anyDuplicated(colnames(counts))) failures <- c(failures, "duplicated sample identifiers")
  for (start in seq.int(1L, nrow(counts), by = chunk_rows)) {
    rows <- start:min(nrow(counts), start + chunk_rows - 1L)
    block <- counts[rows, , drop = FALSE]
    if (anyNA(block)) failures <- c(failures, "missing values")
    if (any(!is.finite(block))) failures <- c(failures, "non-finite values")
    if (any(block < 0)) failures <- c(failures, "negative values")
    if (any(abs(block - round(block)) > 1e-8)) {
      failures <- c(failures, "non-integer-like values")
    }
    if (length(failures)) break
  }
  if (length(failures)) {
    stop(label, " failed raw-count validation: ", paste(unique(failures), collapse = "; "))
  }
  invisible(TRUE)
}

read_metadata <- function(path, sample_ids, default_cohort) {
  if (is.null(path) || !length(path) || !nzchar(path)) {
    return(data.frame(
      sample_id = sample_ids,
      cohort = rep(default_cohort, length(sample_ids)),
      stringsAsFactors = FALSE
    ))
  }
  if (!file.exists(path)) stop("Tumour metadata not found: ", path)
  metadata <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) {
    utils::read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    )
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  sample_column <- intersect(
    c("sample_id", "sample", "aliquot_id", "run"),
    colnames(metadata)
  )
  if (!length(sample_column)) {
    stop("Tumour metadata lacks a recognised sample identifier column: ", path)
  }
  metadata_ids <- as.character(metadata[[sample_column[[1]]]])
  if (anyNA(metadata_ids) || any(!nzchar(metadata_ids)) || anyDuplicated(metadata_ids)) {
    stop("Tumour metadata has missing, empty, or duplicated sample identifiers")
  }
  index <- match(sample_ids, metadata_ids)
  if (anyNA(index)) {
    stop(
      "Tumour metadata lacks eligible samples: ",
      paste(sample_ids[is.na(index)], collapse = ", ")
    )
  }
  metadata <- metadata[index, , drop = FALSE]
  if (!identical(as.character(metadata[[sample_column[[1]]]]), sample_ids)) {
    stop("Tumour metadata could not be aligned exactly to eligible raw-count columns")
  }
  cohort_column <- intersect(c("cohort", "source", "project"), colnames(metadata))
  cohort <- if (length(cohort_column)) {
    as.character(metadata[[cohort_column[[1]]]])
  } else {
    rep(default_cohort, length(sample_ids))
  }
  cohort[is.na(cohort) | !nzchar(cohort)] <- default_cohort
  data.frame(
    sample_id = sample_ids,
    cohort = cohort,
    stringsAsFactors = FALSE
  )
}

map_raw_counts_to_hgnc <- function(counts, map_path) {
  if (!file.exists(map_path)) stop("Gene map not found: ", map_path)
  mapping <- utils::read.delim(
    map_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  ensembl_column <- intersect(
    c("Ensembl_ID", "ensembl_id", "gene_id", "ensembl_gene_id"),
    colnames(mapping)
  )
  symbol_column <- intersect(
    c("HGNC_Symbol", "symbol", "gene_symbol", "hgnc_symbol"),
    colnames(mapping)
  )
  if (!length(ensembl_column) || !length(symbol_column)) {
    stop("Gene map must contain recognised Ensembl and HGNC-symbol columns")
  }

  ensembl <- sub("\\..*$", "", trimws(as.character(mapping[[ensembl_column[[1]]]])))
  symbols <- toupper(trimws(as.character(mapping[[symbol_column[[1]]]])))
  valid <- !is.na(ensembl) & nzchar(ensembl) & !is.na(symbols) & nzchar(symbols)
  symbol_by_ensembl <- setNames(symbols[valid], ensembl[valid])

  count_ensembl <- sub("\\..*$", "", rownames(counts))
  mapped_symbols <- unname(symbol_by_ensembl[count_ensembl])
  keep <- !is.na(mapped_symbols) & nzchar(mapped_symbols)
  if (sum(keep) < 1000L) {
    stop("Only ", sum(keep), " raw-count genes map to HGNC symbols")
  }
  mapped <- counts[keep, , drop = FALSE]
  rownames(mapped) <- mapped_symbols[keep]
  if (anyDuplicated(rownames(mapped))) {
    mapped <- rowsum(mapped, group = rownames(mapped), reorder = TRUE)
  }
  storage.mode(mapped) <- "double"
  validate_raw_counts(mapped, "HGNC-aggregated tumour counts")
  mapped
}

rank_normalise_for_tidyestimate <- function(expression) {
  normalised <- vapply(
    seq_len(ncol(expression)),
    function(column) {
      rank(expression[, column], ties.method = "average", na.last = "keep") *
        (10000 / nrow(expression))
    },
    numeric(nrow(expression))
  )
  if (!is.matrix(normalised)) {
    normalised <- matrix(normalised, nrow = nrow(expression))
  }
  rownames(normalised) <- rownames(expression)
  colnames(normalised) <- colnames(expression)
  storage.mode(normalised) <- "double"
  if (anyNA(normalised) || any(!is.finite(normalised))) {
    stop("Tumour-only ESTIMATE rank normalisation returned invalid values")
  }
  normalised
}

standardise_scores <- function(scores, expected_samples, force_sample_order = FALSE) {
  scores <- as.data.frame(scores, check.names = FALSE)
  explicit_sample <- intersect(c("sample_id", "sample", "Sample"), colnames(scores))
  if (length(explicit_sample)) {
    rownames(scores) <- as.character(scores[[explicit_sample[[1]]]])
    scores[[explicit_sample[[1]]]] <- NULL
  }
  lower_names <- tolower(colnames(scores))
  names <- colnames(scores)
  names[grepl("strom", lower_names)] <- "StromalScore"
  names[grepl("immu", lower_names)] <- "ImmuneScore"
  names[grepl("estimate", lower_names)] <- "ESTIMATEScore"
  names[grepl("purity", lower_names)] <- "tumourPurity"
  colnames(scores) <- names

  default_rownames <- identical(rownames(scores), as.character(seq_len(nrow(scores))))
  if (nrow(scores) == length(expected_samples) &&
      (force_sample_order || is.null(rownames(scores)) || default_rownames)) {
    rownames(scores) <- expected_samples
  }
  required <- c("StromalScore", "ImmuneScore", "ESTIMATEScore")
  missing <- setdiff(required, colnames(scores))
  if (length(missing)) {
    stop("ESTIMATE result lacks score columns: ", paste(missing, collapse = ", "))
  }
  if (!"tumourPurity" %in% colnames(scores)) {
    scores$tumourPurity <- cos(
      0.6049872018 + 0.0001467884 * as.numeric(scores$ESTIMATEScore)
    )
  }
  scores$tumourPurity <- pmin(1, pmax(0, as.numeric(scores$tumourPurity)))
  if (!setequal(rownames(scores), expected_samples)) {
    stop("Purity-score sample IDs do not exactly match eligible raw tumour columns")
  }
  scores <- scores[
    expected_samples,
    c(required, "tumourPurity"),
    drop = FALSE
  ]
  for (column in colnames(scores)) scores[[column]] <- as.numeric(scores[[column]])
  if (anyNA(scores) || any(!is.finite(as.matrix(scores)))) {
    stop("Purity scores contain missing or non-finite values")
  }
  scores
}

compute_estimate_scores <- function(scoring_expression) {
  input <- data.frame(
    hgnc_symbol = rownames(scoring_expression),
    scoring_expression,
    check.names = FALSE
  )
  exports <- getNamespaceExports("tidyestimate")
  score_function <- if ("estimate_score" %in% exports) {
    tidyestimate::estimate_score
  } else if ("estimate_scores" %in% exports) {
    tidyestimate::estimate_scores
  } else {
    stop("Installed tidyestimate exposes neither estimate_score nor estimate_scores")
  }
  standardise_scores(
    score_function(input, is_affymetrix = FALSE),
    colnames(scoring_expression),
    force_sample_order = TRUE
  )
}

optional_input <- function(name) {
  if (!name %in% names(snakemake@input)) return("")
  value <- snakemake@input[[name]]
  if (!length(value)) return("")
  as.character(value[[1]])
}

cohort <- as.character(snakemake@params[["cohort"]])
threshold <- as.numeric(snakemake@params[["threshold"]])
expected_input_samples <- as.integer(snakemake@params[["expected_input_samples"]])
expected_eligible_samples <- as.integer(snakemake@params[["expected_eligible_samples"]])
expected_retained_samples <- as.integer(snakemake@params[["expected_retained_samples"]])
primary_sample_type <- as.character(snakemake@params[["primary_sample_type"]])
score_policy <- as.character(snakemake@params[["score_policy"]])
rule_name <- as.character(snakemake@params[["rule_name"]])
config_path <- as.character(snakemake@params[["config_path"]])

raw_path <- as.character(snakemake@input[["raw_tumour_counts"]])
gene_map_path <- as.character(snakemake@input[["gene_map"]])
metadata_path <- optional_input("tumour_metadata")
reference_scores_path <- optional_input("reference_scores")
reference_script_path <- optional_input("reference_script")
reference_log_path <- optional_input("reference_log")

cat("[INFO] Tumour-purity analysis from raw tumour patient counts\n")
cat("[INFO] cohort=", cohort, "\n", sep = "")
cat("[INFO] rule=", rule_name, "\n", sep = "")
cat("[INFO] raw tumour input=", raw_path, "\n", sep = "")
cat("[INFO] DSMZ influence=false\n")
cat("[INFO] transformation=per-sample average-rank scaling (10000 * rank / genes)\n")
cat("[INFO] purity threshold=", threshold, "\n", sep = "")
cat("[INFO] score policy=", score_policy, "\n", sep = "")

raw_counts <- read_count_matrix(raw_path, paste(cohort, "raw tumour counts"))
validate_raw_counts(raw_counts, paste(cohort, "raw tumour counts"))
if (expected_input_samples > 0L && ncol(raw_counts) != expected_input_samples) {
  stop(
    "Expected ", expected_input_samples, " raw tumour profiles; observed ",
    ncol(raw_counts)
  )
}

eligible_ids <- colnames(raw_counts)
if (nzchar(primary_sample_type)) {
  sample_type <- ifelse(
    grepl("^TCGA-[^-]+-[^-]+-[0-9]{2}", eligible_ids),
    substr(eligible_ids, 14L, 15L),
    NA_character_
  )
  eligible_ids <- eligible_ids[
    !is.na(sample_type) & sample_type == primary_sample_type
  ]
}
if (!length(eligible_ids)) stop("No eligible tumour samples remain after validation")
if (expected_eligible_samples > 0L && length(eligible_ids) != expected_eligible_samples) {
  stop(
    "Expected ", expected_eligible_samples, " eligible tumours; observed ",
    length(eligible_ids)
  )
}

eligible_raw <- raw_counts[, eligible_ids, drop = FALSE]
eligible_metadata <- read_metadata(metadata_path, eligible_ids, cohort)
eligible_manifest <- data.frame(
  raw_order = match(eligible_ids, colnames(raw_counts)),
  eligible_order = seq_along(eligible_ids),
  sample_id = eligible_ids,
  cohort = eligible_metadata$cohort,
  stringsAsFactors = FALSE
)
write_tsv(eligible_manifest, snakemake@output[["eligible_samples"]])

cat(
  "[INFO] raw dimensions=", nrow(raw_counts), "x", ncol(raw_counts),
  "; eligible dimensions=", nrow(eligible_raw), "x", ncol(eligible_raw), "\n",
  sep = ""
)

hgnc_counts <- map_raw_counts_to_hgnc(eligible_raw, gene_map_path)
common_symbols <- unique(toupper(as.character(tidyestimate::common_genes$hgnc_symbol)))
estimate_genes <- intersect(rownames(hgnc_counts), common_symbols)
if (length(estimate_genes) < 50L) {
  stop("Only ", length(estimate_genes), " ESTIMATE common genes are present")
}
estimate_counts <- hgnc_counts[estimate_genes, , drop = FALSE]
scoring_expression <- rank_normalise_for_tidyestimate(estimate_counts)
if (!identical(colnames(scoring_expression), eligible_ids)) {
  stop("Tumour-only scoring expression columns differ from eligible raw tumours")
}
save_rds_atomic(
  scoring_expression,
  snakemake@output[["purity_scoring_expression"]]
)
cat(
  "[INFO] purity scoring dimensions=", nrow(scoring_expression), "x",
  ncol(scoring_expression), "\n", sep = ""
)

if (identical(score_policy, "compute")) {
  scores <- compute_estimate_scores(scoring_expression)
  score_source <- snakemake@output[["purity_scoring_expression"]]
  score_source_sha256 <- sha256_file(score_source)
} else if (identical(score_policy, "reuse_verified")) {
  if (!nzchar(reference_scores_path) || !file.exists(reference_scores_path)) {
    stop("Verified purity-score table is missing: ", reference_scores_path)
  }
  reference_scores <- utils::read.csv(
    reference_scores_path,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  scores <- standardise_scores(
    reference_scores,
    eligible_ids,
    force_sample_order = FALSE
  )
  score_source <- reference_scores_path
  score_source_sha256 <- sha256_file(reference_scores_path)
} else {
  stop("Unsupported purity score policy: ", score_policy)
}

retained_ids <- eligible_ids[scores[eligible_ids, "tumourPurity"] >= threshold]
excluded_ids <- eligible_ids[scores[eligible_ids, "tumourPurity"] < threshold]
if (!length(retained_ids)) stop("Purity threshold excludes every eligible tumour")
if (expected_retained_samples > 0L &&
    length(retained_ids) != expected_retained_samples) {
  stop(
    "Expected ", expected_retained_samples, " retained tumours; observed ",
    length(retained_ids), " at threshold ", threshold
  )
}
if (length(retained_ids) + length(excluded_ids) != length(eligible_ids) ||
    length(intersect(retained_ids, excluded_ids))) {
  stop("Retained/excluded purity partition does not reconcile")
}

score_table <- data.frame(
  sample_id = eligible_ids,
  scores[eligible_ids, , drop = FALSE],
  purity_threshold = threshold,
  retained = eligible_ids %in% retained_ids,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)
utils::write.csv(
  score_table,
  snakemake@output[["purity_scores"]],
  row.names = FALSE
)

cohort_by_sample <- setNames(eligible_metadata$cohort, eligible_metadata$sample_id)
manifest_score_columns <- c(
  "StromalScore", "ImmuneScore", "ESTIMATEScore", "tumourPurity"
)
retained_manifest <- data.frame(
  retained_order = seq_along(retained_ids),
  raw_order = match(retained_ids, colnames(raw_counts)),
  sample_id = retained_ids,
  cohort = unname(cohort_by_sample[retained_ids]),
  scores[retained_ids, manifest_score_columns, drop = FALSE],
  purity_threshold = threshold,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)
excluded_manifest <- data.frame(
  raw_order = match(excluded_ids, colnames(raw_counts)),
  sample_id = excluded_ids,
  cohort = unname(cohort_by_sample[excluded_ids]),
  scores[excluded_ids, manifest_score_columns, drop = FALSE],
  purity_threshold = threshold,
  exclusion_reason = "tumour_purity_below_threshold",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)
write_tsv(retained_manifest, snakemake@output[["retained_samples"]])
write_tsv(excluded_manifest, snakemake@output[["excluded_samples"]])

retained_raw <- raw_counts[, retained_ids, drop = FALSE]
if (!identical(colnames(retained_raw), retained_manifest$sample_id)) {
  stop("Retained raw tumour columns do not exactly equal the retained manifest")
}
if (!identical(
  colnames(retained_raw),
  score_table$sample_id[
    score_table$tumourPurity >= threshold & !is.na(score_table$tumourPurity)
  ]
)) {
  stop("Retained raw tumour columns do not exactly equal samples above threshold")
}
save_rds_atomic(
  retained_raw,
  snakemake@output[["retained_raw_tumour_counts"]]
)

## ---------------------------------------------------------------------------
## Presentation settings for the purity figures.
##
## Appearance and figure composition only: the purity scores, the threshold and
## the retained/excluded partition are exactly the values computed above.
## ---------------------------------------------------------------------------

## The palette, the typography, the gridless panel frame, the threshold line and
## the atomic PDF device all come from the shared figure module loaded at the
## top of this script, so the purity figures of every cohort match the
## batch-correction figures of every cohort. Nothing here defines a colour.

retained_flag <- score_table$retained
purity_values <- score_table$tumourPurity
retained_total <- sum(retained_flag)
excluded_total <- sum(!retained_flag)
retained_percent <- 100 * retained_total / length(retained_flag)

## Composition: the purity distribution carries the decision, so it spans the
## top row; the three score relationships sit beneath it, with a shared legend.
figure_open(snakemake@output[["purity_diagnostics"]], width = 12, height = 8.6)
graphics::layout(
  matrix(c(1, 1, 1, 2, 3, 4, 5, 5, 5), nrow = 3L, byrow = TRUE),
  heights = c(1.05, 1, 0.14)
)
graphics::par(oma = c(0, 0, 3.2, 0))

figure_reset_par()
graphics::par(mar = c(4.3, 4.8, 3.0, 1.6))
histogram <- graphics::hist(purity_values, breaks = 30, plot = FALSE)
bin_count <- length(histogram$counts)
figure_panel_frame(
  xlim = range(histogram$breaks),
  ylim = c(0, max(histogram$counts) * 1.16),
  xlab = "Tumour purity (ESTIMATE)",
  ylab = "Tumour samples",
  main = "Tumour-purity distribution",
  subtitle = sprintf("Dashed line marks the retention threshold of %.2f", threshold)
)
graphics::rect(
  histogram$breaks[seq_len(bin_count)],
  0,
  histogram$breaks[seq_len(bin_count) + 1L],
  histogram$counts,
  col = ifelse(histogram$mids >= threshold, FIGURE_RETAINED, FIGURE_EXCLUDED),
  border = "#ffffff",
  lwd = 0.8
)
figure_threshold_line(threshold, "v")
graphics::text(
  x = threshold,
  y = max(histogram$counts) * 1.12,
  labels = sprintf(
    "  %s retained (%.1f%%)  |  %s excluded",
    figure_count(retained_total), retained_percent, figure_count(excluded_total)
  ),
  adj = c(0, 0.5),
  col = FIGURE_INK,
  cex = 1.0
)

for (score_name in c("StromalScore", "ImmuneScore", "ESTIMATEScore")) {
  figure_reset_par()
  graphics::par(mar = c(4.3, 4.8, 3.0, 1.6))
  score_values <- score_table[[score_name]]
  figure_panel_frame(
    xlim = range(score_values, finite = TRUE),
    ylim = range(purity_values, finite = TRUE),
    xlab = score_name,
    ylab = "Tumour purity",
    main = paste("Purity versus", score_name),
    subtitle = if (identical(score_name, "ESTIMATEScore")) {
      "Purity is a deterministic function of this score"
    } else {
      NULL
    }
  )
  graphics::points(
    score_values[retained_flag],
    purity_values[retained_flag],
    pch = 19, cex = 0.6,
    col = grDevices::adjustcolor(FIGURE_RETAINED, alpha.f = 0.45)
  )
  graphics::points(
    score_values[!retained_flag],
    purity_values[!retained_flag],
    pch = 21, cex = 0.9, lwd = 0.5,
    bg = FIGURE_EXCLUDED, col = "#ffffff"
  )
  figure_threshold_line(threshold, "h")
}

graphics::par(mar = c(0, 0, 0, 0))
graphics::plot.new()
graphics::legend(
  "center",
  legend = c(
    sprintf(FIGURE_LEGEND_RETAINED, threshold, figure_count(retained_total)),
    sprintf(FIGURE_LEGEND_EXCLUDED, threshold, figure_count(excluded_total)),
    sprintf(FIGURE_LEGEND_THRESHOLD, threshold)
  ),
  pch = c(21, 21, NA),
  lty = c(NA, NA, 2),
  lwd = c(NA, NA, 1.8),
  pt.bg = c(FIGURE_RETAINED, FIGURE_EXCLUDED, NA),
  col = c("#ffffff", "#ffffff", FIGURE_INK),
  pt.lwd = 0.6,
  pt.cex = 1.6,
  cex = 1.05,
  horiz = TRUE,
  bty = "n",
  text.col = FIGURE_INK,
  x.intersp = 0.9
)
graphics::mtext(
  sprintf("%s tumour-purity filtering", cohort),
  side = 3, line = 1.1, outer = TRUE, adj = 0.02,
  font = 2, cex = 1.35, col = FIGURE_INK
)
graphics::mtext(
  sprintf(
    "ESTIMATE purity for %s eligible primary tumours  |  threshold %.2f",
    figure_count(length(purity_values)), threshold
  ),
  side = 3, line = -0.3, outer = TRUE, adj = 0.02,
  cex = 0.95, col = FIGURE_INK_MUTED
)
figure_close()

cohort_levels <- sort(unique(eligible_metadata$cohort))
before_counts <- table(factor(eligible_metadata$cohort, levels = cohort_levels))
after_counts <- table(factor(
  unname(cohort_by_sample[retained_ids]),
  levels = cohort_levels
))
retention_counts <- rbind(
  Before = as.integer(before_counts),
  After = as.integer(after_counts)
)
colnames(retention_counts) <- cohort_levels
## Same counts as before; drawn with direct value labels and a per-cohort
## retention rate so the figure is readable without reading the axis.
bar_midpoints <- graphics::barplot(retention_counts, beside = TRUE, plot = FALSE)
bar_values <- as.vector(retention_counts)
bar_half_width <- 0.42
cohort_midpoints <- colMeans(bar_midpoints)
cohort_retention <- 100 * retention_counts["After", ] / pmax(retention_counts["Before", ], 1L)

figure_open(
  snakemake@output[["sample_retention"]],
  width = max(6.4, 2.6 + 1.9 * length(cohort_levels)),
  height = 6
)
figure_reset_par()
graphics::par(mar = c(5.2, 5.0, 4.6, 1.6))
figure_panel_frame(
  xlim = c(min(bar_midpoints) - 0.8, max(bar_midpoints) + 0.8),
  ylim = c(0, max(retention_counts) * 1.16),
  xlab = "",
  ylab = "Tumour samples",
  main = sprintf("%s retention at purity >= %.2f", cohort, threshold),
  subtitle = sprintf(
    "%s eligible primary tumours -> %s retained (%.1f%%); %s excluded",
    figure_count(length(eligible_ids)),
    figure_count(length(retained_ids)),
    100 * length(retained_ids) / length(eligible_ids),
    figure_count(length(eligible_ids) - length(retained_ids))
  ),
  x_axis = FALSE
)
graphics::rect(
  as.vector(bar_midpoints) - bar_half_width,
  0,
  as.vector(bar_midpoints) + bar_half_width,
  bar_values,
  col = rep(c(FIGURE_STAGE_BEFORE, FIGURE_STAGE_AFTER), times = length(cohort_levels)),
  border = "#ffffff",
  lwd = 1.2
)
graphics::text(
  as.vector(bar_midpoints),
  bar_values,
  labels = figure_count(bar_values),
  pos = 3,
  offset = 0.4,
  col = FIGURE_INK,
  cex = 1.0
)
graphics::mtext(
  cohort_levels,
  side = 1, line = 0.9, at = cohort_midpoints,
  col = FIGURE_INK, cex = 1.05, font = 2
)
graphics::mtext(
  sprintf("%.1f%% retained", cohort_retention),
  side = 1, line = 2.2, at = cohort_midpoints,
  col = FIGURE_INK_MUTED, cex = 0.9
)
graphics::legend(
  "topright",
  legend = c(FIGURE_LEGEND_STAGE_BEFORE, FIGURE_LEGEND_STAGE_AFTER),
  fill = c(FIGURE_STAGE_BEFORE, FIGURE_STAGE_AFTER),
  border = "#ffffff",
  bty = "n",
  cex = 1.0,
  text.col = FIGURE_INK
)
figure_close()

validation <- data.frame(
  check = c(
    "raw_count_class",
    "raw_dimensions",
    "raw_unique_gene_ids",
    "raw_unique_sample_ids",
    "raw_missing_values",
    "raw_nonfinite_values",
    "raw_nonnegative_integer_values",
    "score_table_equals_eligible_tumours",
    "scoring_object_tumour_only",
    "retained_excluded_disjoint",
    "retained_plus_excluded_equals_eligible",
    "retained_raw_equals_manifest",
    "retained_raw_equals_threshold_decision",
    "configured_threshold_applied"
  ),
  observed = c(
    paste(class(raw_counts), collapse = ","),
    paste(dim(raw_counts), collapse = "x"),
    sum(duplicated(rownames(raw_counts))),
    sum(duplicated(colnames(raw_counts))),
    sum(is.na(raw_counts)),
    sum(!is.finite(raw_counts)),
    all(raw_counts >= 0 & abs(raw_counts - round(raw_counts)) <= 1e-8),
    identical(score_table$sample_id, eligible_ids),
    identical(colnames(scoring_expression), eligible_ids),
    length(intersect(retained_ids, excluded_ids)) == 0L,
    setequal(c(retained_ids, excluded_ids), eligible_ids),
    identical(colnames(retained_raw), retained_manifest$sample_id),
    identical(
      colnames(retained_raw),
      score_table$sample_id[score_table$tumourPurity >= threshold]
    ),
    all(score_table$purity_threshold == threshold)
  ),
  expected = c(
    "matrix,array",
    paste(nrow(raw_counts), expected_input_samples, sep = "x"),
    0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
  ),
  status = "PASS",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(validation, snakemake@output[["validation_report"]])

provenance <- data.frame(
  field = c(
    "generated_at",
    "cohort",
    "rule",
    "configuration",
    "processing_order",
    "raw_tumour_count_input",
    "raw_tumour_count_sha256",
    "raw_tumour_dimensions",
    "tumour_metadata_input",
    "tumour_metadata_sha256",
    "gene_map_input",
    "gene_map_sha256",
    "eligible_tumour_samples",
    "purity_scoring_expression",
    "purity_scoring_expression_sha256",
    "purity_scoring_dimensions",
    "purity_transformation",
    "dsmz_influenced_transformation",
    "purity_method",
    "purity_score_policy",
    "purity_score_source",
    "purity_score_source_sha256",
    "reference_score_script",
    "reference_score_script_sha256",
    "reference_score_log",
    "reference_score_log_sha256",
    "purity_threshold",
    "tumours_retained",
    "tumours_excluded",
    "retained_raw_tumour_counts",
    "retained_raw_tumour_counts_sha256",
    "tidyestimate_version",
    "absolute_purity_calibration_limitation"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    cohort,
    rule_name,
    config_path,
    paste(
      "raw tumour patient counts -> tumour-only validation/gene mapping ->",
      "tumour-only rank-normalised scoring expression -> ESTIMATE purity ->",
      "purity threshold -> retained raw tumour patient counts"
    ),
    raw_path,
    sha256_file(raw_path),
    paste(dim(raw_counts), collapse = "x"),
    metadata_path,
    sha256_file(metadata_path),
    gene_map_path,
    sha256_file(gene_map_path),
    length(eligible_ids),
    snakemake@output[["purity_scoring_expression"]],
    sha256_file(snakemake@output[["purity_scoring_expression"]]),
    paste(dim(scoring_expression), collapse = "x"),
    paste(
      "HGNC-mapped raw tumour counts filtered to tidyestimate common genes;",
      "per-sample average ranks scaled as 10000 * rank / number of genes;",
      "no DSMZ samples"
    ),
    "false",
    "ESTIMATE via tidyestimate with is_affymetrix=FALSE",
    score_policy,
    score_source,
    score_source_sha256,
    reference_script_path,
    sha256_file(reference_script_path),
    reference_log_path,
    sha256_file(reference_log_path),
    threshold,
    length(retained_ids),
    length(excluded_ids),
    snakemake@output[["retained_raw_tumour_counts"]],
    sha256_file(snakemake@output[["retained_raw_tumour_counts"]]),
    as.character(utils::packageVersion("tidyestimate")),
    paste(
      "The cosine conversion from ESTIMATEScore to absolute tumour purity was",
      "developed for Affymetrix data; its use with RNA-seq and a 0.70 threshold",
      "remains a methodological assumption."
    )
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(provenance, snakemake@output[["provenance"]])

writeLines(
  c(
    "PASS",
    paste0("cohort=", cohort),
    paste0("eligible_tumours=", length(eligible_ids)),
    paste0("retained_tumours=", length(retained_ids)),
    paste0("excluded_tumours=", length(excluded_ids)),
    paste0("threshold=", threshold),
    paste0(
      "retained_raw_sha256=",
      sha256_file(snakemake@output[["retained_raw_tumour_counts"]])
    )
  ),
  snakemake@output[["validation_ok"]]
)

cat(
  "[SUCCESS] ", cohort, " tumour-purity stage completed: ",
  length(eligible_ids), " eligible -> ", length(retained_ids),
  " retained + ", length(excluded_ids), " excluded\n",
  sep = ""
)
