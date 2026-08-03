options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

for (path in unname(snakemake@output)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

log_path <- snakemake@output[["validation_log"]]
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

sha256_file <- function(path) {
  result <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path, ": ", paste(result, collapse = " "))
  }
  strsplit(result[[1]], "[[:space:]]+")[[1]][1]
}

read_matrix <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
  object <- readRDS(path)
  if (methods::is(object, "SummarizedExperiment")) {
    object <- SummarizedExperiment::assay(object)
  } else if (is.list(object) && !is.null(object$counts)) {
    object <- object$counts
  }
  if (!is.matrix(object) && !is.data.frame(object)) {
    stop(label, " is not matrix/data.frame-like; class=", paste(class(object), collapse = ","))
  }
  matrix <- as.matrix(object)
  if (!is.numeric(matrix)) stop(label, " is not numeric")
  if (is.null(rownames(matrix)) || is.null(colnames(matrix))) {
    stop(label, " lacks row or column names")
  }
  matrix
}

read_table <- function(path, label, csv = FALSE) {
  if (!file.exists(path)) stop(label, " not found: ", path)
  table <- if (csv) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    utils::read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    )
  }
  if (!"sample_id" %in% colnames(table)) {
    stop(label, " must contain sample_id")
  }
  table$sample_id <- as.character(table$sample_id)
  if (anyNA(table$sample_id) || any(!nzchar(table$sample_id)) ||
      anyDuplicated(table$sample_id)) {
    stop(label, " has missing, empty, or duplicated sample IDs")
  }
  table
}

validate_expression <- function(matrix, label) {
  failures <- character(0)
  if (anyDuplicated(rownames(matrix))) failures <- c(failures, "duplicated row names")
  if (anyDuplicated(colnames(matrix))) failures <- c(failures, "duplicated column names")
  if (anyNA(matrix)) failures <- c(failures, "missing values")
  if (any(!is.finite(matrix))) failures <- c(failures, "non-finite values")
  if (length(failures)) {
    stop(label, " failed validation: ", paste(failures, collapse = "; "))
  }
  invisible(TRUE)
}

validate_counts <- function(matrix, label, chunk_rows = 1000L) {
  validate_expression(matrix, label)
  for (start in seq.int(1L, nrow(matrix), by = chunk_rows)) {
    rows <- start:min(nrow(matrix), start + chunk_rows - 1L)
    block <- matrix[rows, , drop = FALSE]
    if (any(block < 0)) stop(label, " contains negative values")
    if (any(abs(block - round(block)) > 1e-8)) {
      stop(label, " contains non-integer-like values")
    }
  }
  invisible(TRUE)
}

add_check <- local({
  checks <- list()
  function(name = NULL, observed = NULL, expected = NULL, pass = NULL, finish = FALSE) {
    if (finish) return(do.call(rbind, checks))
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name,
      observed = as.character(observed),
      expected = as.character(expected),
      status = if (isTRUE(pass)) "PASS" else "FAIL",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    invisible(pass)
  }
})

cohort <- as.character(snakemake@params[["cohort"]])
expected_tumours <- as.integer(snakemake@params[["expected_tumours"]])
expected_dsmz <- as.integer(snakemake@params[["expected_dsmz"]])
expected_collapsed_dsmz <- as.integer(snakemake@params[["expected_collapsed_dsmz"]])

retained_manifest <- read_table(
  snakemake@input[["retained_samples"]],
  "retained-sample manifest"
)
excluded_manifest <- read_table(
  snakemake@input[["excluded_samples"]],
  "excluded-sample manifest"
)
score_table <- read_table(
  snakemake@input[["purity_scores"]],
  "purity-score table",
  csv = TRUE
)
if (!all(c("tumourPurity", "purity_threshold") %in% colnames(score_table))) {
  stop("Purity-score table lacks tumourPurity or purity_threshold")
}

retained_raw <- read_matrix(
  snakemake@input[["retained_raw_tumour_counts"]],
  "retained raw tumour counts"
)
dsmz_counts <- read_matrix(
  snakemake@input[["dsmz_counts"]],
  "harmonised raw DSMZ counts"
)
merged_raw <- read_matrix(
  snakemake@input[["merged_counts_raw"]],
  "merged filtered raw counts"
)
corrected_counts <- read_matrix(
  snakemake@input[["merged_counts_batch_corrected"]],
  "ComBat-seq corrected counts"
)
joint_pre <- read_matrix(
  snakemake@input[["joint_vst_pre_bc"]],
  "joint VST before batch correction"
)
joint_post <- read_matrix(
  snakemake@input[["joint_vst_post_bc"]],
  "joint VST after batch correction"
)
tumour_post <- read_matrix(
  snakemake@input[["tumour_vst_post_bc"]],
  "tumour VST after batch correction"
)
dsmz_post <- read_matrix(
  snakemake@input[["dsmz_vst_post_bc"]],
  "DSMZ VST after batch correction"
)
coldata <- read_table(snakemake@input[["coldata"]], "joint sample metadata")

validate_counts(retained_raw, "retained raw tumour counts")
validate_counts(dsmz_counts, "harmonised raw DSMZ counts")
validate_counts(merged_raw, "merged filtered raw counts")
validate_counts(corrected_counts, "ComBat-seq corrected counts")
validate_expression(joint_pre, "joint VST before batch correction")
validate_expression(joint_post, "joint VST after batch correction")
validate_expression(tumour_post, "tumour VST after batch correction")
validate_expression(dsmz_post, "DSMZ VST after batch correction")

retained_ids <- retained_manifest$sample_id
excluded_ids <- excluded_manifest$sample_id
threshold <- unique(as.numeric(score_table$purity_threshold))
if (length(threshold) != 1L || !is.finite(threshold)) {
  stop("Purity-score table does not contain exactly one finite threshold")
}
threshold_ids <- score_table$sample_id[
  as.numeric(score_table$tumourPurity) >= threshold
]
dsmz_ids <- colnames(dsmz_counts)
joint_ids <- c(retained_ids, dsmz_ids)

downstream_sample_sets <- list(
  retained_raw = colnames(retained_raw),
  merged_raw = colnames(merged_raw),
  corrected_counts = colnames(corrected_counts),
  joint_pre = colnames(joint_pre),
  joint_post = colnames(joint_post),
  tumour_post = colnames(tumour_post),
  dsmz_post = colnames(dsmz_post),
  coldata = coldata$sample_id
)
excluded_downstream_hits <- unique(unlist(lapply(
  downstream_sample_sets,
  function(ids) intersect(excluded_ids, ids)
)))

add_check(
  "retained_manifest_count",
  length(retained_ids),
  expected_tumours,
  expected_tumours <= 0L || length(retained_ids) == expected_tumours
)
add_check(
  "dsmz_profile_count",
  length(dsmz_ids),
  expected_dsmz,
  expected_dsmz <= 0L || length(dsmz_ids) == expected_dsmz
)
add_check(
  "retained_raw_columns_equal_manifest",
  paste(colnames(retained_raw), collapse = ","),
  paste(retained_ids, collapse = ","),
  identical(colnames(retained_raw), retained_ids)
)
add_check(
  "retained_manifest_equals_threshold_decision",
  paste(retained_ids, collapse = ","),
  paste(threshold_ids, collapse = ","),
  identical(retained_ids, threshold_ids)
)
add_check(
  "retained_excluded_disjoint",
  length(intersect(retained_ids, excluded_ids)),
  0,
  length(intersect(retained_ids, excluded_ids)) == 0L
)
add_check(
  "score_partition_complete",
  length(c(retained_ids, excluded_ids)),
  nrow(score_table),
  setequal(c(retained_ids, excluded_ids), score_table$sample_id)
)
add_check(
  "excluded_tumours_absent_downstream",
  paste(excluded_downstream_hits, collapse = ","),
  "",
  length(excluded_downstream_hits) == 0L
)
add_check(
  "dsmz_samples_unique",
  sum(duplicated(dsmz_ids)),
  0,
  !anyDuplicated(dsmz_ids)
)
add_check(
  "merged_sample_order",
  paste(colnames(merged_raw), collapse = ","),
  paste(joint_ids, collapse = ","),
  identical(colnames(merged_raw), joint_ids)
)
add_check(
  "corrected_sample_order",
  paste(colnames(corrected_counts), collapse = ","),
  paste(joint_ids, collapse = ","),
  identical(colnames(corrected_counts), joint_ids)
)
add_check(
  "joint_pre_sample_order",
  paste(colnames(joint_pre), collapse = ","),
  paste(joint_ids, collapse = ","),
  identical(colnames(joint_pre), joint_ids)
)
add_check(
  "joint_post_sample_order",
  paste(colnames(joint_post), collapse = ","),
  paste(joint_ids, collapse = ","),
  identical(colnames(joint_post), joint_ids)
)
add_check(
  "coldata_sample_order",
  paste(coldata$sample_id, collapse = ","),
  paste(joint_ids, collapse = ","),
  identical(coldata$sample_id, joint_ids)
)
add_check(
  "gene_order_merged_corrected",
  identical(rownames(merged_raw), rownames(corrected_counts)),
  TRUE,
  identical(rownames(merged_raw), rownames(corrected_counts))
)
add_check(
  "gene_order_merged_joint_pre",
  identical(rownames(merged_raw), rownames(joint_pre)),
  TRUE,
  identical(rownames(merged_raw), rownames(joint_pre))
)
add_check(
  "gene_order_merged_joint_post",
  identical(rownames(merged_raw), rownames(joint_post)),
  TRUE,
  identical(rownames(merged_raw), rownames(joint_post))
)
add_check(
  "tumour_split_exact",
  paste(colnames(tumour_post), collapse = ","),
  paste(retained_ids, collapse = ","),
  identical(colnames(tumour_post), retained_ids)
)
add_check(
  "dsmz_split_exact",
  paste(colnames(dsmz_post), collapse = ","),
  paste(dsmz_ids, collapse = ","),
  identical(colnames(dsmz_post), dsmz_ids)
)
add_check(
  "joint_post_class",
  paste(class(joint_post), collapse = ","),
  "matrix,array",
  is.matrix(joint_post)
)
add_check(
  "joint_post_dimensions",
  paste(dim(joint_post), collapse = "x"),
  paste(nrow(merged_raw), length(joint_ids), sep = "x"),
  identical(dim(joint_post), c(nrow(merged_raw), length(joint_ids)))
)
add_check(
  "joint_post_unique_rows",
  sum(duplicated(rownames(joint_post))),
  0,
  !anyDuplicated(rownames(joint_post))
)
add_check(
  "joint_post_unique_columns",
  sum(duplicated(colnames(joint_post))),
  0,
  !anyDuplicated(colnames(joint_post))
)
add_check(
  "joint_post_missing_values",
  sum(is.na(joint_post)),
  0,
  !anyNA(joint_post)
)
add_check(
  "joint_post_nonfinite_values",
  sum(!is.finite(joint_post)),
  0,
  all(is.finite(joint_post))
)

collapsed_path <- if ("collapsed_dsmz" %in% names(snakemake@input)) {
  value <- snakemake@input[["collapsed_dsmz"]]
  if (length(value)) as.character(value[[1]]) else ""
} else {
  ""
}
if (nzchar(collapsed_path)) {
  collapsed <- read_matrix(collapsed_path, "collapsed DSMZ VST")
  validate_expression(collapsed, "collapsed DSMZ VST")
  add_check(
    "collapsed_dsmz_profile_count",
    ncol(collapsed),
    expected_collapsed_dsmz,
    expected_collapsed_dsmz <= 0L || ncol(collapsed) == expected_collapsed_dsmz
  )
  add_check(
    "collapsed_dsmz_gene_order",
    identical(rownames(collapsed), rownames(dsmz_post)),
    TRUE,
    identical(rownames(collapsed), rownames(dsmz_post))
  )
}

figure_names <- grep("^figure_", names(snakemake@input), value = TRUE)
for (figure_name in figure_names) {
  figure_path <- as.character(snakemake@input[[figure_name]])
  figure_size <- if (file.exists(figure_path)) file.info(figure_path)$size else 0
  add_check(
    paste0("nonempty_", figure_name),
    figure_size,
    ">1000",
    is.finite(figure_size) && figure_size > 1000
  )
}

validation <- add_check(finish = TRUE)
write_tsv(validation, snakemake@output[["validation_report"]])

final_manifest <- data.frame(
  cohort = cohort,
  repository_relative_path = as.character(snakemake@params[["final_relative_path"]]),
  absolute_path = normalizePath(
    snakemake@input[["joint_vst_post_bc"]],
    mustWork = TRUE
  ),
  object_class = paste(class(joint_post), collapse = ","),
  genes = nrow(joint_post),
  columns = ncol(joint_post),
  tumour_samples = length(retained_ids),
  dsmz_samples = length(dsmz_ids),
  duplicated_rows = sum(duplicated(rownames(joint_post))),
  duplicated_columns = sum(duplicated(colnames(joint_post))),
  missing_values = sum(is.na(joint_post)),
  nonfinite_values = sum(!is.finite(joint_post)),
  metadata_exact_order = identical(coldata$sample_id, colnames(joint_post)),
  purity_manifest_exact_order = identical(retained_ids, colnames(tumour_post)),
  sha256 = sha256_file(snakemake@input[["joint_vst_post_bc"]]),
  meaning = paste(
    "VST of the joint matrix formed from purity-retained raw tumour counts",
    "and raw DSMZ counts after ComBat-seq correction."
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(final_manifest, snakemake@output[["final_matrix_manifest"]])

failed <- validation$check[validation$status != "PASS"]
if (length(failed)) {
  stop(
    cohort,
    " final purity-filtered preprocessing validation failed: ",
    paste(failed, collapse = "; ")
  )
}

writeLines(
  c(
    "PASS",
    paste0("cohort=", cohort),
    paste0("genes=", nrow(joint_post)),
    paste0("tumour_samples=", length(retained_ids)),
    paste0("dsmz_samples=", length(dsmz_ids)),
    paste0("sha256=", sha256_file(snakemake@input[["joint_vst_post_bc"]]))
  ),
  snakemake@output[["validation_ok"]]
)
cat("[SUCCESS] ", cohort, " final output validation passed\n", sep = "")
