options(stringsAsFactors = FALSE)

for (path in unname(snakemake@output)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

log_path <- snakemake@output[["validation_log"]]
log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
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

read_numeric_matrix <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("%s not found: %s", label, path))
  obj <- readRDS(path)
  if (methods::is(obj, "SummarizedExperiment")) {
    obj <- SummarizedExperiment::assay(obj)
  }
  if (!is.matrix(obj) && !is.data.frame(obj)) {
    stop(sprintf("%s must be matrix/data.frame-like; class=%s", label, paste(class(obj), collapse = ",")))
  }
  mat <- as.matrix(obj)
  if (!is.numeric(mat) || is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(sprintf("%s must be a named numeric matrix", label))
  }
  mat
}

matrix_is_finite <- function(mat, chunk_rows = 1000L) {
  for (start in seq.int(1L, nrow(mat), by = chunk_rows)) {
    idx <- start:min(nrow(mat), start + chunk_rows - 1L)
    block <- mat[idx, , drop = FALSE]
    if (anyNA(block) || any(!is.finite(block))) return(FALSE)
  }
  TRUE
}

checks <- list()
add_check <- function(check, status, observed, expected, detail = "") {
  allowed <- c("PASS", "FAIL", "WARN", "INFO")
  if (!status %in% allowed) stop("Unknown validation status: ", status)
  checks[[length(checks) + 1L]] <<- data.frame(
    check = as.character(check),
    status = status,
    observed = as.character(observed),
    expected = as.character(expected),
    detail = as.character(detail),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cat(sprintf(
    "[%s] %s | observed=%s | expected=%s%s\n",
    status,
    check,
    observed,
    expected,
    if (nzchar(detail)) paste0(" | ", detail) else ""
  ))
}

pass_fail <- function(condition) if (isTRUE(condition)) "PASS" else "FAIL"

expected_tumours <- as.integer(snakemake@params[["expected_tumours"]])
expected_dsmz <- as.integer(snakemake@params[["expected_dsmz"]])
expected_genes <- as.integer(snakemake@params[["expected_shared_genes"]])
correlation_threshold <- as.numeric(snakemake@params[["correlation_threshold"]])
distribution_sample_size <- as.integer(snakemake@params[["distribution_sample_size"]])
expected_total <- expected_tumours + expected_dsmz

cat("[INFO] Validating regenerated BRCA outputs against the protected active matrix\n")
cat("[INFO] Validation is read-only with respect to all protected data/brca objects.\n")

active <- read_numeric_matrix(snakemake@input[["active_reference_rds"]], "active BRCA matrix")
regenerated <- read_numeric_matrix(snakemake@input[["joint_vst_post_bc_rds"]], "regenerated joint VST")
tumour <- read_numeric_matrix(snakemake@input[["tumour_vst_post_bc_rds"]], "regenerated tumour VST")
dsmz <- read_numeric_matrix(snakemake@input[["dsmz_vst_post_bc_rds"]], "regenerated DSMZ VST")

add_check(
  "regenerated joint dimensions",
  pass_fail(identical(dim(regenerated), c(expected_genes, expected_total))),
  paste(dim(regenerated), collapse = "x"),
  paste(expected_genes, expected_total, sep = "x")
)
add_check(
  "regenerated tumour dimensions",
  pass_fail(identical(dim(tumour), c(expected_genes, expected_tumours))),
  paste(dim(tumour), collapse = "x"),
  paste(expected_genes, expected_tumours, sep = "x")
)
add_check(
  "regenerated DSMZ dimensions",
  pass_fail(identical(dim(dsmz), c(expected_genes, expected_dsmz))),
  paste(dim(dsmz), collapse = "x"),
  paste(expected_genes, expected_dsmz, sep = "x")
)
add_check(
  "active reference dimensions",
  pass_fail(identical(dim(active), c(expected_genes, expected_total))),
  paste(dim(active), collapse = "x"),
  paste(expected_genes, expected_total, sep = "x")
)

add_check(
  "regenerated joint duplicated gene names",
  pass_fail(!anyDuplicated(rownames(regenerated))),
  sum(duplicated(rownames(regenerated))),
  0
)
add_check(
  "regenerated joint duplicated sample names",
  pass_fail(!anyDuplicated(colnames(regenerated))),
  sum(duplicated(colnames(regenerated))),
  0
)
add_check(
  "regenerated tumour duplicated names",
  pass_fail(!anyDuplicated(rownames(tumour)) && !anyDuplicated(colnames(tumour))),
  sum(duplicated(rownames(tumour))) + sum(duplicated(colnames(tumour))),
  0
)
add_check(
  "regenerated DSMZ duplicated names",
  pass_fail(!anyDuplicated(rownames(dsmz)) && !anyDuplicated(colnames(dsmz))),
  sum(duplicated(rownames(dsmz))) + sum(duplicated(colnames(dsmz))),
  0
)
joint_finite <- matrix_is_finite(regenerated)
tumour_finite <- matrix_is_finite(tumour)
dsmz_finite <- matrix_is_finite(dsmz)
add_check("regenerated joint finite values", pass_fail(joint_finite), joint_finite, TRUE)
add_check("regenerated tumour finite values", pass_fail(tumour_finite), tumour_finite, TRUE)
add_check("regenerated DSMZ finite values", pass_fail(dsmz_finite), dsmz_finite, TRUE)

active_tcga <- colnames(active)[grepl("^TCGA-", colnames(active))]
active_dsmz <- colnames(active)[grepl("^NG-", colnames(active))]
regenerated_tcga <- colnames(regenerated)[grepl("^TCGA-", colnames(regenerated))]
regenerated_dsmz <- colnames(regenerated)[grepl("^NG-", colnames(regenerated))]

joint_sample_set_identical <- identical(sort(colnames(regenerated)), sort(colnames(active)))
joint_sample_order_identical <- identical(colnames(regenerated), colnames(active))
gene_set_identical <- identical(sort(rownames(regenerated)), sort(rownames(active)))
gene_order_identical <- identical(rownames(regenerated), rownames(active))

add_check(
  "joint sample-set identity with active matrix",
  pass_fail(joint_sample_set_identical),
  joint_sample_set_identical,
  TRUE
)
add_check(
  "joint sample order identity with active matrix",
  if (joint_sample_order_identical) "PASS" else "WARN",
  joint_sample_order_identical,
  TRUE,
  "Order differences do not change sample-set identity but are reported for drop-in comparison"
)
add_check(
  "gene-set identity with active matrix",
  pass_fail(gene_set_identical),
  gene_set_identical,
  TRUE
)
add_check(
  "gene order identity with active matrix",
  if (gene_order_identical) "PASS" else "WARN",
  gene_order_identical,
  TRUE
)
add_check(
  "TCGA sample-set identity",
  pass_fail(identical(sort(regenerated_tcga), sort(active_tcga))),
  identical(sort(regenerated_tcga), sort(active_tcga)),
  TRUE
)
add_check(
  "DSMZ sample-set identity",
  pass_fail(identical(sort(regenerated_dsmz), sort(active_dsmz))),
  identical(sort(regenerated_dsmz), sort(active_dsmz)),
  TRUE
)
add_check(
  "tumour split identity",
  pass_fail(identical(sort(colnames(tumour)), sort(active_tcga))),
  identical(sort(colnames(tumour)), sort(active_tcga)),
  TRUE
)
add_check(
  "DSMZ split identity",
  pass_fail(identical(sort(colnames(dsmz)), sort(active_dsmz))),
  identical(sort(colnames(dsmz)), sort(active_dsmz)),
  TRUE
)

coldata <- utils::read.delim(
  snakemake@input[["coldata_tsv"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
retained_manifest <- utils::read.delim(
  snakemake@input[["retained_manifest"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
add_check(
  "coldata sample-set identity",
  pass_fail("sample_id" %in% colnames(coldata) &&
              identical(sort(coldata$sample_id), sort(colnames(regenerated)))),
  if ("sample_id" %in% colnames(coldata)) nrow(coldata) else "missing sample_id",
  expected_total
)
add_check(
  "retained manifest count",
  pass_fail(nrow(retained_manifest) == expected_tumours),
  nrow(retained_manifest),
  expected_tumours
)

required_files <- c(
  provenance = snakemake@input[["provenance_tsv"]],
  batch_log = snakemake@input[["batch_log"]],
  tumour_filtering_log = snakemake@input[["purity_log"]]
)
for (name in names(required_files)) {
  size <- if (file.exists(required_files[[name]])) file.info(required_files[[name]])$size else 0
  add_check(
    paste(name, "exists and is non-empty"),
    pass_fail(is.finite(size) && size > 0),
    size,
    ">0 bytes"
  )
}

common_genes <- intersect(rownames(active), rownames(regenerated))
common_samples <- intersect(colnames(active), colnames(regenerated))
active_gene_index <- match(common_genes, rownames(active))
regenerated_gene_index <- match(common_genes, rownames(regenerated))
active_sample_index <- match(common_samples, colnames(active))
regenerated_sample_index <- match(common_samples, colnames(regenerated))

sample_correlations <- vapply(seq_along(common_samples), function(i) {
  stats::cor(
    active[active_gene_index, active_sample_index[[i]]],
    regenerated[regenerated_gene_index, regenerated_sample_index[[i]]],
    method = "pearson"
  )
}, numeric(1))
correlation_table <- data.frame(
  sample_id = common_samples,
  sample_class = ifelse(grepl("^TCGA-", common_samples), "tumour", "cell_line"),
  common_genes = length(common_genes),
  pearson_correlation = sample_correlations,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(correlation_table, snakemake@output[["correlations_tsv"]])

median_correlation <- stats::median(sample_correlations, na.rm = TRUE)
minimum_correlation <- min(sample_correlations, na.rm = TRUE)
correlation_status <- if (
  is.finite(median_correlation) && median_correlation >= correlation_threshold
) "PASS" else "WARN"
add_check(
  "median per-sample Pearson correlation with active matrix",
  correlation_status,
  sprintf("%.6f", median_correlation),
  sprintf(">=%0.2f", correlation_threshold),
  sprintf("minimum=%0.6f across %d matched samples", minimum_correlation, length(common_samples))
)

sample_values <- function(mat, n) {
  total <- length(mat)
  index <- unique(as.integer(seq.int(1L, total, length.out = min(n, total))))
  mat[index]
}
probabilities <- c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
active_values <- sample_values(active, distribution_sample_size)
regenerated_values <- sample_values(regenerated, distribution_sample_size)
active_quantiles <- stats::quantile(active_values, probabilities, names = FALSE, na.rm = TRUE)
regenerated_quantiles <- stats::quantile(regenerated_values, probabilities, names = FALSE, na.rm = TRUE)
distribution <- data.frame(
  probability = probabilities,
  active_value = active_quantiles,
  regenerated_value = regenerated_quantiles,
  regenerated_minus_active = regenerated_quantiles - active_quantiles,
  sampled_values_per_matrix = c(length(active_values), rep(NA_integer_, length(probabilities) - 1L)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(distribution, snakemake@output[["distribution_tsv"]])
add_check(
  "sampled median VST shift",
  "INFO",
  sprintf("%.6f", regenerated_quantiles[probabilities == 0.5] - active_quantiles[probabilities == 0.5]),
  "reported, not silently thresholded",
  sprintf("systematic sample size per matrix=%d", length(active_values))
)

# The current count RDS files post-date the active matrix, and the historical
# producing chain is unresolved. Technical similarity cannot establish
# historical identity or authorise promotion.
ready_to_replace <- FALSE
add_check(
  "ready to replace protected active matrix",
  "WARN",
  ready_to_replace,
  FALSE,
  paste(
    "Regenerated outputs remain a reproducible current variant;",
    "historical source identity is unresolved and explicit promotion approval is required"
  )
)

report <- do.call(rbind, checks)
write_tsv(report, snakemake@output[["report_tsv"]])

failed <- report$check[report$status == "FAIL"]
if (length(failed)) {
  stop(sprintf("Regenerated BRCA output validation failed: %s", paste(failed, collapse = "; ")))
}

writeLines(
  c(
    "BRCA regenerated-output structural validation passed",
    sprintf("validated_at=%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    sprintf("joint_dimensions=%dx%d", nrow(regenerated), ncol(regenerated)),
    sprintf("sample_set_identical=%s", joint_sample_set_identical),
    sprintf("gene_set_identical=%s", gene_set_identical),
    sprintf("median_sample_pearson=%.8f", median_correlation),
    "ready_to_replace_protected_active_matrix=false",
    "promotion_requires_explicit_approval=true"
  ),
  con = snakemake@output[["validation_ok"]]
)
cat("[SUCCESS] Structural validation passed; regenerated outputs were not promoted\n")
