#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list(
    joint_rds = "data/brca/brca_vst_joint.rds",
    dsmz_rds = "data/brca/dsmz_brca_counts.rds",
    tcga_rds = "data/brca/tcga_brca_counts.rds",
    pam50_annotation_csv = "data/brca/annotations/TCGA_BRCA_project_subtype.csv",
    pam50_expr_rds = "results/unsupervised/brca/tumour_neighbourhoods_input/expr_pam50.rds",
    pam50_map_rds = "results/unsupervised/brca/tumour_neighbourhoods_input/cell_line_to_original_sample_id_pam50.rds"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args)) stop("Missing value for argument: ", key)
    name <- gsub("-", "_", substring(key, 3), fixed = TRUE)
    out[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

status <- new.env(parent = emptyenv())
status$errors <- character()
status$warnings <- character()

fail <- function(...) {
  msg <- paste0(...)
  status$errors <- c(status$errors, msg)
  cat("[FAIL] ", msg, "\n", sep = "")
}

warn <- function(...) {
  msg <- paste0(...)
  status$warnings <- c(status$warnings, msg)
  cat("[WARN] ", msg, "\n", sep = "")
}

pass <- function(...) {
  cat("[OK] ", paste0(...), "\n", sep = "")
}

read_rds <- function(path, label) {
  if (!file.exists(path)) {
    fail(label, " missing: ", path)
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) {
    fail(label, " unreadable RDS: ", path, " (", conditionMessage(e), ")")
    NULL
  })
}

as_expr_matrix <- function(obj, label) {
  if (is.null(obj)) return(NULL)
  if (is.matrix(obj)) return(obj)
  if (inherits(obj, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      fail(label, " is a Matrix object but the Matrix package is unavailable")
      return(NULL)
    }
    return(as.matrix(obj))
  }
  if (inherits(obj, "SummarizedExperiment") || inherits(obj, "DESeqTransform")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      fail(label, " is a SummarizedExperiment-like object but SummarizedExperiment is unavailable")
      return(NULL)
    }
    return(SummarizedExperiment::assay(obj))
  }
  fail(label, " has unsupported class: ", paste(class(obj), collapse = ", "))
  NULL
}

validate_matrix <- function(path, label, expect_ng = FALSE, expect_tcga = FALSE) {
  obj <- read_rds(path, label)
  mat <- as_expr_matrix(obj, label)
  if (is.null(mat)) return(NULL)

  pass(label, " readable as ", paste(class(obj), collapse = ", "))
  cat(sprintf("[INFO] %s dimensions: %d genes x %d samples\n", label, nrow(mat), ncol(mat)))

  if (!is.numeric(mat)) fail(label, " must be numeric")
  if (is.null(rownames(mat))) fail(label, " must have gene rownames")
  if (is.null(colnames(mat))) fail(label, " must have sample colnames")
  if (!is.null(rownames(mat)) && anyDuplicated(rownames(mat))) fail(label, " has duplicated gene rownames")
  if (!is.null(colnames(mat)) && anyDuplicated(colnames(mat))) fail(label, " has duplicated sample colnames")
  if (nrow(mat) < 1000) warn(label, " has fewer than 1000 gene rows")
  if (nrow(mat) <= ncol(mat)) warn(label, " orientation is unusual; expected genes in rows and samples in columns")

  if (is.numeric(mat)) {
    n_na <- sum(is.na(mat))
    n_inf <- sum(is.infinite(mat))
    if (n_na > 0) fail(label, " contains NA values: ", n_na)
    if (n_inf > 0) fail(label, " contains infinite values: ", n_inf)
  }

  cn <- colnames(mat)
  n_ng <- sum(grepl("^NG-", cn))
  n_tcga <- sum(grepl("^TCGA-", cn))
  n_other_tumour <- sum(grepl("^(GSE|SRP|TARGET-|SRR|GDC-)", cn))
  cat(sprintf("[INFO] %s sample IDs: NG=%d TCGA=%d other_tumour_like=%d\n",
              label, n_ng, n_tcga, n_other_tumour))

  if (expect_ng && n_ng == 0) fail(label, " has no DSMZ NG-* sample columns")
  if (expect_tcga && n_tcga == 0) fail(label, " has no TCGA-* tumour columns")
  if (identical(label, "BRCA joint VST") && n_tcga > 0 && n_other_tumour == 0) {
    warn(label, " uses TCGA-* tumour columns; split_joint_vst_by_sample_type.R must treat TCGA-* as tumour samples")
  }

  mat
}

validate_annotation <- function(path) {
  if (!file.exists(path)) {
    fail("PAM50 annotation CSV missing: ", path)
    return(NULL)
  }
  ann <- tryCatch(read.csv(path, check.names = FALSE), error = function(e) {
    fail("PAM50 annotation CSV unreadable: ", path, " (", conditionMessage(e), ")")
    NULL
  })
  if (is.null(ann)) return(NULL)
  pass("PAM50 annotation CSV readable")
  required <- c("sample_barcode", "PAM50_Subtype")
  missing <- setdiff(required, names(ann))
  if (length(missing)) fail("PAM50 annotation missing columns: ", paste(missing, collapse = ", "))
  if ("sample_barcode" %in% names(ann) && anyDuplicated(ann$sample_barcode)) {
    fail("PAM50 annotation has duplicated sample_barcode values")
  }
  cat(sprintf("[INFO] PAM50 annotation rows: %d\n", nrow(ann)))
  ann
}

validate_pam50_staged_inputs <- function(expr_path, map_path) {
  if (file.exists(expr_path) && file.exists(map_path)) {
    expr <- validate_matrix(expr_path, "Staged PAM50 expression", expect_ng = FALSE, expect_tcga = FALSE)
    mp <- read_rds(map_path, "Staged PAM50 mapping")
    if (!is.null(mp)) {
      if (is.null(names(mp))) fail("Staged PAM50 mapping must be a named vector/list")
      pass("Staged PAM50 mapping readable")
      cat(sprintf("[INFO] Staged PAM50 mapping entries: %d\n", length(mp)))
    }
    return(invisible(expr))
  }
  warn("BRCA config enables pam50_euc/pam50_corr, but staged PAM50 inputs are missing: ",
       expr_path, " and/or ", map_path)
  warn("Dry runs will fail unless these files are produced or PAM50 directions are disabled.")
  invisible(NULL)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

cat("=== BRCA RDS structure validation ===\n")
joint <- validate_matrix(args$joint_rds, "BRCA joint VST", expect_ng = TRUE, expect_tcga = TRUE)
dsmz <- validate_matrix(args$dsmz_rds, "DSMZ BRCA source", expect_ng = TRUE, expect_tcga = FALSE)
tcga <- validate_matrix(args$tcga_rds, "TCGA BRCA source", expect_ng = FALSE, expect_tcga = TRUE)
ann <- validate_annotation(args$pam50_annotation_csv)
validate_pam50_staged_inputs(args$pam50_expr_rds, args$pam50_map_rds)

if (!is.null(joint) && !is.null(dsmz)) {
  missing_ng <- setdiff(colnames(joint)[grepl("^NG-", colnames(joint))], colnames(dsmz))
  if (length(missing_ng)) fail("Joint VST has DSMZ columns not present in DSMZ source: ", paste(head(missing_ng, 10), collapse = ", "))
  common_genes <- length(intersect(rownames(joint), rownames(dsmz)))
  cat(sprintf("[INFO] Joint/DSMZ shared genes: %d\n", common_genes))
}

if (!is.null(joint) && !is.null(tcga)) {
  missing_tcga <- setdiff(colnames(joint)[grepl("^TCGA-", colnames(joint))], colnames(tcga))
  if (length(missing_tcga)) fail("Joint VST has TCGA columns not present in TCGA source: ", paste(head(missing_tcga, 10), collapse = ", "))
  common_genes <- length(intersect(rownames(joint), rownames(tcga)))
  cat(sprintf("[INFO] Joint/TCGA shared genes: %d\n", common_genes))
  if (common_genes == 0) {
    tcga_genes_unversioned <- sub("\\..*$", "", rownames(tcga))
    common_unversioned <- length(intersect(rownames(joint), tcga_genes_unversioned))
    cat(sprintf("[INFO] Joint/TCGA shared genes after Ensembl version stripping: %d\n", common_unversioned))
    if (common_unversioned == 0) {
      warn("Joint and TCGA source rownames do not overlap, even after Ensembl version stripping")
    }
  }
}

if (!is.null(joint) && !is.null(ann) && "sample_barcode" %in% names(ann)) {
  tcga_cols <- colnames(joint)[grepl("^TCGA-", colnames(joint))]
  annotated <- tcga_cols %in% ann$sample_barcode
  cat(sprintf("[INFO] Joint TCGA samples with exact PAM50 annotation: %d/%d\n",
              sum(annotated), length(tcga_cols)))
  if (sum(annotated) == 0) {
    warn("No exact match between joint TCGA column names and PAM50 sample_barcode values")
  }
}

cat("\n=== Summary ===\n")
cat(sprintf("Errors: %d\n", length(status$errors)))
cat(sprintf("Warnings: %d\n", length(status$warnings)))

if (length(status$errors)) {
  quit(status = 1)
}
quit(status = 0)
