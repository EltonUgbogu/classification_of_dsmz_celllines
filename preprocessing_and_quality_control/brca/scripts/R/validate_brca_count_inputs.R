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

checks <- list()
add_check <- function(check, observed, expected, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = as.character(check),
    observed = as.character(observed),
    expected = as.character(expected),
    status = if (isTRUE(pass)) "PASS" else "FAIL",
    detail = as.character(detail),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cat(sprintf(
    "[%s] %s | observed=%s | expected=%s%s\n",
    if (isTRUE(pass)) "PASS" else "FAIL",
    check,
    observed,
    expected,
    if (nzchar(detail)) paste0(" | ", detail) else ""
  ))
}

read_count_matrix <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("%s does not exist: %s", label, path))
  obj <- readRDS(path)
  if (methods::is(obj, "SummarizedExperiment")) {
    obj <- SummarizedExperiment::assay(obj)
  } else if (is.list(obj) && !is.null(obj$counts)) {
    obj <- obj$counts
  }
  if (!is.matrix(obj) && !is.data.frame(obj)) {
    stop(sprintf("%s must be matrix/data.frame-like; class=%s", label, paste(class(obj), collapse = ",")))
  }
  mat <- as.matrix(obj)
  if (!is.numeric(mat)) stop(sprintf("%s is not numeric", label))
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(sprintf("%s requires row names and column names", label))
  }
  mat
}

scan_values <- function(mat, chunk_rows = 1000L) {
  has_na <- FALSE
  has_nonfinite <- FALSE
  has_negative <- FALSE
  has_noninteger <- FALSE
  observed_min <- Inf
  observed_max <- -Inf
  starts <- seq.int(1L, nrow(mat), by = chunk_rows)
  for (start in starts) {
    idx <- start:min(nrow(mat), start + chunk_rows - 1L)
    block <- mat[idx, , drop = FALSE]
    has_na <- has_na || anyNA(block)
    has_nonfinite <- has_nonfinite || any(!is.finite(block))
    finite <- block[is.finite(block)]
    if (length(finite)) {
      observed_min <- min(observed_min, min(finite))
      observed_max <- max(observed_max, max(finite))
      has_negative <- has_negative || any(finite < 0)
      has_noninteger <- has_noninteger || any(abs(finite - round(finite)) > 1e-8)
    }
  }
  list(
    has_na = has_na,
    has_nonfinite = has_nonfinite,
    has_negative = has_negative,
    has_noninteger = has_noninteger,
    min = observed_min,
    max = observed_max
  )
}

strip_ensembl_version <- function(ids) {
  sub("\\..*$", "", ids)
}

tcga_sample_type <- function(ids) {
  ifelse(grepl("^TCGA-[^-]+-[^-]+-[0-9]{2}", ids), substr(ids, 14L, 15L), NA_character_)
}

cat("[INFO] Validating BRCA count-level inputs\n")
cat("[INFO] This rule performs validation only; it does not run purity scoring or batch correction.\n")

tumour <- read_count_matrix(snakemake@input[["tumour_rds"]], "TCGA-BRCA count matrix")
dsmz <- read_count_matrix(snakemake@input[["dsmz_rds"]], "DSMZ BRCA count matrix")
meta <- utils::read.csv(
  snakemake@input[["dsmz_meta_csv"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expected_raw <- as.integer(snakemake@params[["expected_raw_tcga"]])
expected_primary <- as.integer(snakemake@params[["expected_primary"]])
expected_dsmz <- as.integer(snakemake@params[["expected_dsmz"]])
expected_shared <- as.integer(snakemake@params[["expected_shared_genes"]])
primary_code <- as.character(snakemake@params[["primary_code"]])
dsmz_sample_col <- as.character(snakemake@params[["dsmz_sample_col"]])
dsmz_code_col <- as.character(snakemake@params[["dsmz_code_col"]])
dsmz_code <- as.character(snakemake@params[["dsmz_code"]])
dsmz_cell_line_col <- as.character(snakemake@params[["dsmz_cell_line_col"]])

add_check("TCGA object class", paste(class(tumour), collapse = ","), "numeric matrix", is.matrix(tumour))
add_check("DSMZ object class", paste(class(dsmz), collapse = ","), "numeric matrix", is.matrix(dsmz))

tumour_values <- scan_values(tumour)
dsmz_values <- scan_values(dsmz)
add_check("TCGA has no NA values", tumour_values$has_na, FALSE, !tumour_values$has_na)
add_check("TCGA has only finite values", tumour_values$has_nonfinite, FALSE, !tumour_values$has_nonfinite)
add_check("TCGA has no negative values", tumour_values$has_negative, FALSE, !tumour_values$has_negative)
add_check("TCGA is integer-like", tumour_values$has_noninteger, FALSE, !tumour_values$has_noninteger,
          sprintf("range=%s..%s", tumour_values$min, tumour_values$max))
add_check("DSMZ has no NA values", dsmz_values$has_na, FALSE, !dsmz_values$has_na)
add_check("DSMZ has only finite values", dsmz_values$has_nonfinite, FALSE, !dsmz_values$has_nonfinite)
add_check("DSMZ has no negative values", dsmz_values$has_negative, FALSE, !dsmz_values$has_negative)
add_check("DSMZ is integer-like", dsmz_values$has_noninteger, FALSE, !dsmz_values$has_noninteger,
          sprintf("range=%s..%s", dsmz_values$min, dsmz_values$max))

add_check("TCGA duplicated row names", sum(duplicated(rownames(tumour))), 0, !anyDuplicated(rownames(tumour)))
add_check("TCGA duplicated column names", sum(duplicated(colnames(tumour))), 0, !anyDuplicated(colnames(tumour)))
add_check("DSMZ duplicated row names", sum(duplicated(rownames(dsmz))), 0, !anyDuplicated(rownames(dsmz)))
add_check("DSMZ duplicated column names", sum(duplicated(colnames(dsmz))), 0, !anyDuplicated(colnames(dsmz)))

add_check("raw TCGA sample count", ncol(tumour), expected_raw, ncol(tumour) == expected_raw)
types <- tcga_sample_type(colnames(tumour))
type_table <- table(types, useNA = "ifany")
cat("[INFO] TCGA sample-type counts:\n")
print(type_table)
primary_ids <- colnames(tumour)[!is.na(types) & types == primary_code]
add_check("raw TCGA primary-tumour count", length(primary_ids), expected_primary, length(primary_ids) == expected_primary,
          sprintf("primary code=%s", primary_code))

required_meta <- c(dsmz_sample_col, dsmz_code_col, dsmz_cell_line_col)
missing_meta <- setdiff(required_meta, colnames(meta))
add_check("DSMZ metadata required columns", paste(missing_meta, collapse = ","), "none missing", length(missing_meta) == 0)
if (!length(missing_meta)) {
  match_counts <- vapply(colnames(dsmz), function(id) sum(meta[[dsmz_sample_col]] == id, na.rm = TRUE), integer(1))
  add_check(
    "DSMZ count columns match exactly one metadata row",
    sum(match_counts != 1L),
    0,
    all(match_counts == 1L)
  )
  matched_meta <- meta[match(colnames(dsmz), meta[[dsmz_sample_col]]), , drop = FALSE]
  code_ok <- !is.na(matched_meta[[dsmz_code_col]]) & matched_meta[[dsmz_code_col]] == dsmz_code
  add_check("DSMZ profiles have BRCA cancer code", sum(!code_ok), 0, all(code_ok))
  cell_lines <- matched_meta[[dsmz_cell_line_col]]
  valid_cell_lines <- !is.na(cell_lines) & nzchar(trimws(cell_lines))
  add_check("DSMZ profiles have cell-line labels", sum(!valid_cell_lines), 0, all(valid_cell_lines))
  add_check(
    "unique DSMZ cell-line profiles",
    length(unique(cell_lines[valid_cell_lines])),
    expected_dsmz,
    length(unique(cell_lines[valid_cell_lines])) == expected_dsmz
  )
}

tumour_genes <- strip_ensembl_version(rownames(tumour))
dsmz_genes <- strip_ensembl_version(rownames(dsmz))
add_check(
  "TCGA duplicate occurrences after Ensembl-version harmonisation",
  sum(duplicated(tumour_genes)),
  44,
  sum(duplicated(tumour_genes)) == 44,
  "The 44 colliding _PAR_Y rows are validated as zero before aggregation"
)
add_check(
  "DSMZ duplicated genes after Ensembl-version stripping",
  sum(duplicated(dsmz_genes)),
  0,
  !anyDuplicated(dsmz_genes)
)
shared_genes <- intersect(unique(tumour_genes), unique(dsmz_genes))
add_check("shared harmonised genes", length(shared_genes), expected_shared, length(shared_genes) == expected_shared)
tcga_duplicate_rows <- which(duplicated(tumour_genes))
tcga_duplicate_nonzero <- if (length(tcga_duplicate_rows)) {
  any(tumour[tcga_duplicate_rows, , drop = FALSE] != 0)
} else {
  FALSE
}
add_check(
  "TCGA harmonisation duplicate rows are zero",
  tcga_duplicate_nonzero,
  FALSE,
  !tcga_duplicate_nonzero,
  "Aggregation therefore preserves all non-zero counts"
)
add_check("TCGA-only harmonised genes", length(setdiff(unique(tumour_genes), unique(dsmz_genes))), 0,
          length(setdiff(unique(tumour_genes), unique(dsmz_genes))) == 0)

report <- do.call(rbind, checks)
utils::write.table(
  report,
  file = snakemake@output[["report_tsv"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

failed <- report$check[report$status == "FAIL"]
if (length(failed)) {
  stop(sprintf("BRCA count-input validation failed: %s", paste(failed, collapse = "; ")))
}

writeLines(
  c(
    "BRCA count-level input validation passed",
    sprintf("validated_at=%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    sprintf("tcga_dimensions=%dx%d", nrow(tumour), ncol(tumour)),
    sprintf("tcga_primary_tumours=%d", length(primary_ids)),
    sprintf("dsmz_dimensions=%dx%d", nrow(dsmz), ncol(dsmz)),
    sprintf("shared_harmonised_genes=%d", length(shared_genes))
  ),
  con = snakemake@output[["validation_ok"]]
)
cat("[SUCCESS] BRCA count-level inputs passed all declared checks\n")
