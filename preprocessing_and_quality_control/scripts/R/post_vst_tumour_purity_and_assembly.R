options(stringsAsFactors = FALSE)

for (path in unname(snakemake@output)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

log_path <- snakemake@output[["execution_log"]]
log_connection <- file(log_path, open = "wt")
sink(log_connection, split = TRUE)
on.exit({
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
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Failed to atomically promote RDS output: ", path)
  }
  invisible(path)
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path, ": ", paste(output, collapse = " "))
  }
  strsplit(output[[1]], "[[:space:]]+")[[1]][1]
}

read_numeric_matrix <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
  object <- readRDS(path)
  if (inherits(object, "SummarizedExperiment")) {
    object <- SummarizedExperiment::assay(object)
  } else if (is.list(object) && !is.null(object$counts)) {
    object <- object$counts
  }
  if (!is.matrix(object) && !is.data.frame(object)) {
    stop(label, " must be matrix/data.frame-like; class=", paste(class(object), collapse = ","))
  }
  matrix <- as.matrix(object)
  if (!is.numeric(matrix)) stop(label, " must be numeric")
  if (is.null(rownames(matrix)) || is.null(colnames(matrix))) {
    stop(label, " must have row and column names")
  }
  matrix
}

validate_expression_matrix <- function(
    matrix,
    label,
    min_genes,
    min_samples = 2L,
    expected_samples = 0L) {
  failures <- character(0)
  if (nrow(matrix) < min_genes) {
    failures <- c(failures, sprintf("%d genes is below minimum %d", nrow(matrix), min_genes))
  }
  if (ncol(matrix) < min_samples) {
    failures <- c(failures, sprintf("%d samples is below minimum %d", ncol(matrix), min_samples))
  }
  if (expected_samples > 0L && ncol(matrix) != expected_samples) {
    failures <- c(failures, sprintf("%d samples does not equal expected %d", ncol(matrix), expected_samples))
  }
  if (anyDuplicated(rownames(matrix))) failures <- c(failures, "duplicated gene identifiers")
  if (anyDuplicated(colnames(matrix))) failures <- c(failures, "duplicated sample identifiers")
  if (anyNA(matrix)) failures <- c(failures, "missing values")
  if (any(!is.finite(matrix))) failures <- c(failures, "non-finite values")
  if (length(failures)) {
    stop(label, " failed expression-matrix validation: ", paste(failures, collapse = "; "))
  }
  invisible(TRUE)
}

read_coldata <- function(path) {
  table <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  sample_candidates <- intersect(
    c("sample_id", "sample", "Sample", "run", "aliquot_id"),
    colnames(table)
  )
  if (length(sample_candidates)) {
    sample_ids <- as.character(table[[sample_candidates[[1]]]])
  } else if (ncol(table) > 0L &&
             (colnames(table)[[1]] == "" || grepl("^Unnamed", colnames(table)[[1]]))) {
    sample_ids <- as.character(table[[1]])
    table <- table[, -1, drop = FALSE]
  } else {
    stop("coldata must contain an explicit sample_id/sample column: ", path)
  }
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("coldata contains missing, empty, or duplicated sample identifiers")
  }
  table$sample_id <- sample_ids
  table
}

map_ensembl_to_hgnc <- function(expression, map_path) {
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
  ensembl <- sub("\\..*$", "", as.character(mapping[[ensembl_column[[1]]]]))
  symbols <- toupper(trimws(as.character(mapping[[symbol_column[[1]]]])))
  valid_mapping <- !is.na(ensembl) & nzchar(ensembl) & !is.na(symbols) & nzchar(symbols)
  symbol_by_ensembl <- setNames(symbols[valid_mapping], ensembl[valid_mapping])

  expression_ensembl <- sub("\\..*$", "", rownames(expression))
  mapped_symbols <- unname(symbol_by_ensembl[expression_ensembl])
  keep <- !is.na(mapped_symbols) & nzchar(mapped_symbols)
  if (sum(keep) < 1000L) {
    stop("Only ", sum(keep), " expression genes mapped to HGNC symbols")
  }
  mapped <- expression[keep, , drop = FALSE]
  rownames(mapped) <- mapped_symbols[keep]
  if (anyDuplicated(rownames(mapped))) {
    mapped <- rowsum(mapped, group = rownames(mapped), reorder = TRUE)
  }
  storage.mode(mapped) <- "double"
  mapped
}

standardise_scores <- function(scores, expected_samples, force_sample_order = FALSE) {
  scores <- as.data.frame(scores, check.names = FALSE)
  if ("sample" %in% colnames(scores)) {
    rownames(scores) <- as.character(scores$sample)
    scores$sample <- NULL
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
  if (length(missing)) stop("ESTIMATE result lacks score columns: ", paste(missing, collapse = ", "))
  if (!"tumourPurity" %in% colnames(scores)) {
    scores$tumourPurity <- cos(0.6049872018 + 0.0001467884 * scores$ESTIMATEScore)
  }
  scores$tumourPurity <- pmin(1, pmax(0, scores$tumourPurity))
  if (!setequal(rownames(scores), expected_samples)) {
    stop("ESTIMATE score sample IDs do not exactly match tumour VST columns")
  }
  scores <- scores[expected_samples, c(required, "tumourPurity"), drop = FALSE]
  if (anyNA(scores) || any(!is.finite(as.matrix(scores)))) {
    stop("ESTIMATE returned missing or non-finite scores")
  }
  scores
}

cohort <- as.character(snakemake@params[["cohort"]])
threshold <- as.numeric(snakemake@params[["threshold"]])
expected_tumours_before <- as.integer(snakemake@params[["expected_tumours_before"]])
expected_tumours_retained <- as.integer(snakemake@params[["expected_tumours_retained"]])
expected_dsmz <- as.integer(snakemake@params[["expected_dsmz"]])
min_genes <- as.integer(snakemake@params[["min_genes"]])
score_policy <- as.character(snakemake@params[["score_policy"]])

cat("[INFO] Post-VST tumour-purity workflow\n")
cat("[INFO] cohort=", cohort, "\n", sep = "")
cat("[INFO] method=ESTIMATE/tidyestimate\n")
cat("[INFO] score policy=", score_policy, "\n", sep = "")
cat("[INFO] purity threshold=", threshold, "\n", sep = "")
cat("[INFO] authoritative order: raw counts -> ComBat-seq -> VST -> purity -> final subset\n")

joint <- read_numeric_matrix(snakemake@input[["joint_vst_post_bc"]], "joint post-BC VST")
tumour <- read_numeric_matrix(snakemake@input[["tumour_vst_post_bc"]], "tumour post-BC VST")
dsmz <- read_numeric_matrix(snakemake@input[["dsmz_vst_post_bc"]], "DSMZ post-BC VST")
coldata <- read_coldata(snakemake@input[["coldata"]])

validate_expression_matrix(
  joint, "joint post-BC VST", min_genes,
  min_samples = expected_tumours_before + expected_dsmz,
  expected_samples = expected_tumours_before + expected_dsmz
)
validate_expression_matrix(
  tumour, "tumour post-BC VST", min_genes,
  min_samples = expected_tumours_before,
  expected_samples = expected_tumours_before
)
validate_expression_matrix(
  dsmz, "DSMZ post-BC VST", min_genes,
  min_samples = expected_dsmz,
  expected_samples = expected_dsmz
)

if (!identical(rownames(joint), rownames(tumour)) ||
    !identical(rownames(joint), rownames(dsmz))) {
  stop("Joint, tumour, and DSMZ VST matrices do not have identical genes in order")
}
if (!setequal(colnames(joint), c(colnames(tumour), colnames(dsmz)))) {
  stop("Joint VST samples do not exactly equal tumour plus DSMZ samples")
}
if (!setequal(coldata$sample_id, colnames(joint))) {
  stop("Batch coldata samples do not exactly equal joint VST samples")
}
coldata <- coldata[match(colnames(joint), coldata$sample_id), , drop = FALSE]

if (identical(score_policy, "compute_post_vst")) {
  tumour_hgnc <- map_ensembl_to_hgnc(
    tumour,
    snakemake@input[["gene_map"]]
  )
  common_symbols <- unique(toupper(as.character(tidyestimate::common_genes$hgnc_symbol)))
  estimate_genes <- intersect(rownames(tumour_hgnc), common_symbols)
  if (length(estimate_genes) < 50L) {
    stop("Only ", length(estimate_genes), " ESTIMATE common genes are present")
  }
  cat("[INFO] HGNC genes for ESTIMATE=", length(estimate_genes), "\n", sep = "")
  estimate_input <- data.frame(
    hgnc_symbol = estimate_genes,
    tumour_hgnc[estimate_genes, , drop = FALSE],
    check.names = FALSE
  )
  score_exports <- getNamespaceExports("tidyestimate")
  score_function <- if ("estimate_score" %in% score_exports) {
    tidyestimate::estimate_score
  } else if ("estimate_scores" %in% score_exports) {
    tidyestimate::estimate_scores
  } else {
    stop("Installed tidyestimate exposes neither estimate_score nor estimate_scores")
  }
  scores <- standardise_scores(
    score_function(estimate_input, is_affymetrix = FALSE),
    colnames(tumour),
    force_sample_order = TRUE
  )
  score_source_path <- snakemake@input[["tumour_vst_post_bc"]]
  score_source_sha256 <- sha256_file(score_source_path)
  purity_method_description <- "ESTIMATE/tidyestimate computed on post-ComBat-seq tumour VST"
} else if (identical(score_policy, "apply_verified_scores")) {
  reference_path <- as.character(snakemake@input[["reference_scores"]])
  if (!nzchar(reference_path) || !file.exists(reference_path)) {
    stop("Verified purity-score input is missing: ", reference_path)
  }
  reference_scores <- utils::read.csv(
    reference_path,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  scores <- standardise_scores(
    reference_scores,
    colnames(tumour),
    force_sample_order = FALSE
  )
  score_source_path <- reference_path
  score_source_sha256 <- sha256_file(reference_path)
  purity_method_description <- paste(
    "Verified ESTIMATE/tidyestimate scores applied after post-ComBat-seq VST;",
    "scores were not recomputed from transformed values"
  )
} else {
  stop("Unsupported purity score policy: ", score_policy)
}

retained_ids <- rownames(scores)[scores$tumourPurity >= threshold]
excluded_ids <- rownames(scores)[scores$tumourPurity < threshold]
if (!length(retained_ids)) stop("Purity threshold excluded every tumour sample")
if (expected_tumours_retained > 0L && length(retained_ids) != expected_tumours_retained) {
  stop(
    "Purity retained ", length(retained_ids), " tumour samples; expected ",
    expected_tumours_retained, " at threshold ", threshold
  )
}
if (length(retained_ids) + length(excluded_ids) != ncol(tumour)) {
  stop("Retained/excluded purity accounting does not reconcile")
}

score_table <- data.frame(
  sample_id = rownames(scores),
  scores,
  purity_threshold = threshold,
  retained = rownames(scores) %in% retained_ids,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)
utils::write.csv(score_table, snakemake@output[["purity_scores"]], row.names = FALSE)

tumour_coldata <- coldata[match(colnames(tumour), coldata$sample_id), , drop = FALSE]
cohort_values <- if ("cohort" %in% colnames(tumour_coldata)) {
  as.character(tumour_coldata$cohort)
} else {
  rep(cohort, nrow(tumour_coldata))
}
cohort_values[is.na(cohort_values) | !nzchar(cohort_values)] <- cohort
cohort_by_sample <- setNames(cohort_values, tumour_coldata$sample_id)

manifest_columns <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "tumourPurity")
retained_manifest <- data.frame(
  retained_order = seq_along(retained_ids),
  sample_id = retained_ids,
  cohort = unname(cohort_by_sample[retained_ids]),
  scores[retained_ids, manifest_columns, drop = FALSE],
  purity_threshold = threshold,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)
excluded_manifest <- data.frame(
  raw_order = match(excluded_ids, colnames(tumour)),
  sample_id = excluded_ids,
  cohort = unname(cohort_by_sample[excluded_ids]),
  scores[excluded_ids, manifest_columns, drop = FALSE],
  purity_threshold = threshold,
  exclusion_reason = "tumour_purity_below_threshold",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)
write_tsv(retained_manifest, snakemake@output[["retained_samples"]])
write_tsv(excluded_manifest, snakemake@output[["excluded_samples"]])

final_ids <- colnames(joint)[
  colnames(joint) %in% retained_ids | colnames(joint) %in% colnames(dsmz)
]
final_tumour <- tumour[, retained_ids, drop = FALSE]
final_joint <- joint[, final_ids, drop = FALSE]
final_coldata <- coldata[match(final_ids, coldata$sample_id), , drop = FALSE]

validate_expression_matrix(
  final_tumour, "final post-purity tumour VST", min_genes,
  min_samples = length(retained_ids),
  expected_samples = length(retained_ids)
)
validate_expression_matrix(
  final_joint, "final post-purity joint VST", min_genes,
  min_samples = length(retained_ids) + expected_dsmz,
  expected_samples = length(retained_ids) + expected_dsmz
)
if (!identical(final_coldata$sample_id, colnames(final_joint))) {
  stop("Final metadata order does not exactly match final joint VST columns")
}
if (!setequal(intersect(colnames(final_joint), colnames(tumour)), retained_ids)) {
  stop("Final tumour columns do not exactly match retained purity samples")
}
if (!setequal(intersect(colnames(final_joint), colnames(dsmz)), colnames(dsmz))) {
  stop("Final joint VST does not retain every DSMZ profile")
}

save_rds_atomic(final_tumour, snakemake@output[["tumour_vst_post_purity"]])
save_rds_atomic(final_joint, snakemake@output[["joint_vst_post_purity"]])
write_tsv(final_coldata, snakemake@output[["coldata_post_purity"]])

grDevices::pdf(snakemake@output[["purity_diagnostics"]], width = 10, height = 8)
graphics::par(mfrow = c(2, 2))
graphics::hist(
  scores$tumourPurity,
  breaks = 30,
  main = paste(cohort, "tumour-purity distribution"),
  xlab = "Tumour purity",
  col = "#5B8FF9",
  border = "white"
)
graphics::abline(v = threshold, col = "#D62728", lwd = 2, lty = 2)
for (score_name in c("StromalScore", "ImmuneScore", "ESTIMATEScore")) {
  graphics::plot(
    scores[[score_name]],
    scores$tumourPurity,
    pch = 19,
    cex = 0.6,
    col = grDevices::adjustcolor("#2F4B7C", alpha.f = 0.55),
    xlab = score_name,
    ylab = "Tumour purity",
    main = paste("Purity versus", score_name)
  )
  graphics::abline(h = threshold, col = "#D62728", lwd = 2, lty = 2)
}
grDevices::dev.off()

before_counts <- table(factor(cohort_by_sample, levels = sort(unique(cohort_by_sample))))
after_counts <- table(factor(cohort_by_sample[retained_ids], levels = names(before_counts)))
retention <- rbind(Before = as.integer(before_counts), After = as.integer(after_counts))
colnames(retention) <- names(before_counts)
grDevices::pdf(snakemake@output[["sample_retention"]], width = 9, height = 6)
graphics::barplot(
  retention,
  beside = TRUE,
  col = c("#BDBDBD", "#2CA25F"),
  ylab = "Tumour samples",
  main = sprintf("%s retention at purity >= %.2f (%d -> %d)", cohort, threshold, ncol(tumour), length(retained_ids)),
  legend.text = rownames(retention),
  args.legend = list(x = "topright", bty = "n")
)
grDevices::dev.off()

validation <- data.frame(
  check = c(
    "object_class",
    "minimum_gene_count",
    "final_dimensions",
    "duplicated_gene_ids",
    "duplicated_sample_ids",
    "missing_values",
    "nonfinite_values",
    "metadata_exact_order",
    "retained_purity_exact_set",
    "dsmz_exact_set"
  ),
  observed = c(
    paste(class(final_joint), collapse = ","),
    nrow(final_joint),
    paste(dim(final_joint), collapse = "x"),
    sum(duplicated(rownames(final_joint))),
    sum(duplicated(colnames(final_joint))),
    sum(is.na(final_joint)),
    sum(!is.finite(final_joint)),
    identical(final_coldata$sample_id, colnames(final_joint)),
    setequal(intersect(colnames(final_joint), colnames(tumour)), retained_ids),
    setequal(intersect(colnames(final_joint), colnames(dsmz)), colnames(dsmz))
  ),
  expected = c(
    "matrix,array",
    paste0(">=", min_genes),
    paste(nrow(final_joint), length(retained_ids) + expected_dsmz, sep = "x"),
    0,
    0,
    0,
    0,
    TRUE,
    TRUE,
    TRUE
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
    "processing_order",
    "batch_correction_method",
    "normalisation_method",
    "purity_method",
    "purity_score_policy",
    "purity_score_source",
    "purity_score_source_sha256",
    "purity_threshold",
    "tumours_before",
    "tumours_retained",
    "tumours_excluded",
    "dsmz_profiles",
    "genes",
    "joint_vst_post_bc_input",
    "joint_vst_post_bc_sha256",
    "tumour_vst_post_bc_input",
    "tumour_vst_post_bc_sha256",
    "dsmz_vst_post_bc_input",
    "dsmz_vst_post_bc_sha256",
    "purity_scores_output",
    "retained_manifest_output",
    "excluded_manifest_output",
    "joint_vst_post_purity_output",
    "joint_vst_post_purity_sha256",
    "tidyestimate_version"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    cohort,
    if (identical(score_policy, "compute_post_vst")) {
      "raw counts -> validation -> ComBat-seq -> DESeq2 VST -> ESTIMATE purity estimation/filtering -> post-purity joint VST"
    } else {
      "raw counts -> validation -> ComBat-seq -> DESeq2 VST -> verified ESTIMATE purity filtering -> post-purity joint VST"
    },
    "sva::ComBat_seq on raw joint counts",
    "DESeq2::vst after ComBat-seq",
    purity_method_description,
    score_policy,
    score_source_path,
    score_source_sha256,
    threshold,
    ncol(tumour),
    length(retained_ids),
    length(excluded_ids),
    ncol(dsmz),
    nrow(final_joint),
    snakemake@input[["joint_vst_post_bc"]],
    sha256_file(snakemake@input[["joint_vst_post_bc"]]),
    snakemake@input[["tumour_vst_post_bc"]],
    sha256_file(snakemake@input[["tumour_vst_post_bc"]]),
    snakemake@input[["dsmz_vst_post_bc"]],
    sha256_file(snakemake@input[["dsmz_vst_post_bc"]]),
    snakemake@output[["purity_scores"]],
    snakemake@output[["retained_samples"]],
    snakemake@output[["excluded_samples"]],
    snakemake@output[["joint_vst_post_purity"]],
    sha256_file(snakemake@output[["joint_vst_post_purity"]]),
    as.character(utils::packageVersion("tidyestimate"))
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(provenance, snakemake@output[["provenance"]])
writeLines(
  sprintf(
    "PASS\t%s\t%d genes\t%d tumours\t%d DSMZ\tsha256=%s",
    cohort,
    nrow(final_joint),
    length(retained_ids),
    ncol(dsmz),
    sha256_file(snakemake@output[["joint_vst_post_purity"]])
  ),
  snakemake@output[["validation_ok"]]
)

cat(sprintf(
  "[SUCCESS] %s: %d/%d tumours retained; %d DSMZ; final=%d x %d\n",
  cohort,
  length(retained_ids),
  ncol(tumour),
  ncol(dsmz),
  nrow(final_joint),
  ncol(final_joint)
))
