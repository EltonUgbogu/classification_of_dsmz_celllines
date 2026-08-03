#!/usr/bin/env Rscript

parse_args <- function(args) {
  defaults <- list(
    counts_tsv = "results/unsupervised/nbl/deseq2_inputs/counts.tsv",
    metadata_csv = "data/dsmz/DSMZ_metadata.csv",
    gene_map_tsv = "data/nbl/count_data/nbl_ensembl_to_hgnc.tsv",
    output_rds = paste0(
      "preprocessing_and_quality_control/nbl/results/dsmz_input/",
      "DSMZ_nbl_raw_counts.rds"
    ),
    validation_tsv = paste0(
      "preprocessing_and_quality_control/nbl/results/dsmz_input/",
      "DSMZ_nbl_raw_counts.validation.tsv"
    )
  )
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args)) stop("Missing value for --", key)
    key <- gsub("-", "_", key)
    defaults[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  defaults
}

record <- function(check, status, details) {
  data.frame(
    check = check,
    status = status,
    details = as.character(details),
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c(args$counts_tsv, args$metadata_csv, args$gene_map_tsv)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required input(s): ", paste(missing, collapse = ", "))
if (file.exists(args$output_rds)) {
  stop("Refusing to overwrite existing output: ", args$output_rds)
}
if (file.exists(args$validation_tsv)) {
  stop("Refusing to overwrite existing validation report: ", args$validation_tsv)
}

counts <- utils::read.delim(
  args$counts_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(names(counts)[1], "gene_id")) {
  stop("Expected the first count-table column to be gene_id")
}
gene_id_with_version <- as.character(counts$gene_id)
if (anyNA(gene_id_with_version) || any(!nzchar(gene_id_with_version))) {
  stop("gene_id contains missing or empty values")
}
if (anyDuplicated(gene_id_with_version)) {
  stop("gene_id contains duplicates")
}

sample_ids <- setdiff(names(counts), "gene_id")
if (!length(sample_ids) || anyDuplicated(sample_ids)) {
  stop("Count table has no sample columns or has duplicated sample columns")
}
count_matrix <- as.matrix(counts[, sample_ids, drop = FALSE])
suppressWarnings(storage.mode(count_matrix) <- "numeric")
if (anyNA(count_matrix) || any(!is.finite(count_matrix))) {
  stop("Count payload contains missing or non-finite values")
}
if (any(count_matrix < 0)) stop("Count payload contains negative values")
if (any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
  stop("Count payload contains non-integer values")
}

metadata <- utils::read.csv(
  args$metadata_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_meta <- c("sample_name", "Disease")
if (!all(required_meta %in% names(metadata))) {
  stop("DSMZ metadata lacks required columns: ", paste(required_meta, collapse = ", "))
}
meta_idx <- match(sample_ids, metadata$sample_name)
if (anyNA(meta_idx)) {
  stop(
    "Count samples missing from DSMZ metadata: ",
    paste(sample_ids[is.na(meta_idx)], collapse = ", ")
  )
}
sample_diseases <- unique(metadata$Disease[meta_idx])
if (!identical(sample_diseases, "Neuroblastoma")) {
  stop(
    "Expected only Neuroblastoma DSMZ samples; observed: ",
    paste(sample_diseases, collapse = ", ")
  )
}

ensembl_id <- sub("\\..*$", "", gene_id_with_version)
n_duplicate_unversioned <- sum(duplicated(ensembl_id))
gene_map <- utils::read.delim(
  args$gene_map_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!all(c("Ensembl_ID", "HGNC_Symbol") %in% names(gene_map))) {
  stop("Gene map lacks Ensembl_ID and HGNC_Symbol")
}
if (anyDuplicated(gene_map$Ensembl_ID)) {
  stop("Gene map contains duplicated Ensembl_ID values")
}
gene_name <- gene_map$HGNC_Symbol[match(ensembl_id, gene_map$Ensembl_ID)]

output <- data.frame(
  Ensembl_ID = ensembl_id,
  gene_name = gene_name,
  Ensembl_ID_with_version = gene_id_with_version,
  count_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dir.create(dirname(args$output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(output, args$output_rds)

checks <- do.call(
  rbind,
  list(
    record("source_counts_path", "PASS", args$counts_tsv),
    record("source_is_integer_nonnegative", "PASS", TRUE),
    record("dimensions", "PASS", paste(nrow(count_matrix), ncol(count_matrix), sep = "x")),
    record("unique_gene_ids", "PASS", nrow(count_matrix)),
    record(
      "duplicate_unversioned_ensembl_ids",
      "PASS",
      paste0(
        n_duplicate_unversioned,
        " duplicate rows retained here; the declared batch loader aggregates ",
        "them with rowsum after version stripping"
      )
    ),
    record("unique_sample_ids", "PASS", length(sample_ids)),
    record("metadata_samples_matched", "PASS", sum(!is.na(meta_idx))),
    record("metadata_disease", "PASS", paste(sample_diseases, collapse = "|")),
    record("hgnc_symbols_mapped", "PASS", sum(!is.na(gene_name) & nzchar(gene_name))),
    record("output_object_class", "PASS", paste(class(output), collapse = ",")),
    record("output_path", "PASS", args$output_rds)
  )
)
utils::write.table(
  checks,
  args$validation_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Created validated NBL DSMZ raw-count staging object\n")
cat("Source:", args$counts_tsv, "\n")
cat("Output:", args$output_rds, "\n")
cat("Dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
cat("DSMZ disease:", paste(sample_diseases, collapse = "|"), "\n")
cat("Mapped HGNC symbols:", sum(!is.na(gene_name) & nzchar(gene_name)), "\n")
