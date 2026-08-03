#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

options <- list(
  make_option("--input", type = "character"),
  make_option("--output", type = "character"),
  make_option("--min-genes", type = "integer", default = 1000L),
  make_option("--min-samples", type = "integer", default = 3L)
)
args <- parse_args(OptionParser(option_list = options))
if (is.null(args$input) || is.null(args$output)) {
  stop("--input and --output are required")
}
if (!file.exists(args$input)) stop("Expression RDS not found: ", args$input)

object <- readRDS(args$input)
if (inherits(object, "SummarizedExperiment")) {
  object <- SummarizedExperiment::assay(object)
} else if (is.list(object) && !is.null(object$counts)) {
  object <- object$counts
}
if (!is.matrix(object) && !is.data.frame(object)) {
  stop("Expression object must be matrix/data.frame-like; class=", paste(class(object), collapse = ","))
}
matrix <- as.matrix(object)
if (!is.numeric(matrix)) stop("Expression matrix must be numeric")

failures <- character(0)
if (nrow(matrix) < args$`min-genes`) {
  failures <- c(failures, sprintf("%d genes < required %d", nrow(matrix), args$`min-genes`))
}
if (ncol(matrix) < args$`min-samples`) {
  failures <- c(failures, sprintf("%d samples < required %d", ncol(matrix), args$`min-samples`))
}
if (is.null(rownames(matrix)) || is.null(colnames(matrix))) {
  failures <- c(failures, "missing gene or sample identifiers")
} else {
  if (anyDuplicated(rownames(matrix))) failures <- c(failures, "duplicated gene identifiers")
  if (anyDuplicated(colnames(matrix))) failures <- c(failures, "duplicated sample identifiers")
}
if (anyNA(matrix)) failures <- c(failures, "missing values")
if (any(!is.finite(matrix))) failures <- c(failures, "non-finite values")
if (length(failures)) {
  stop("Active expression matrix is invalid: ", paste(failures, collapse = "; "))
}

sha <- strsplit(system2("sha256sum", args$input, stdout = TRUE), "[[:space:]]+")[[1]][1]
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
writeLines(
  sprintf(
    "PASS\tclass=%s\tdimensions=%dx%d\tsha256=%s\tinput=%s",
    paste(class(object), collapse = ","),
    nrow(matrix),
    ncol(matrix),
    sha,
    normalizePath(args$input)
  ),
  args$output
)
