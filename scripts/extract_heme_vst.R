#!/usr/bin/env Rscript
# =============================================================================
# extract_heme_vst.R
# =============================================================================
#
# Extract HEME samples from DSMZ count file and generate VST matrix
# Uses same preprocessing approach as BRCA/NBL/RBL
#
# USAGE:
# Rscript scripts/extract_heme_vst.R \
#   --counts <validated HEME-specific DSMZ raw-count RDS> \
#   --metadata data/dsmz/DSMZ_metadata.csv \
#   --output results/heme/inputs/vst_joint.rds
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(optparse)
})

option_list <- list(
  make_option(c("--counts"), type="character", default=NULL,
              help="Path to DSMZ count file (RDS)"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Path to DSMZ metadata CSV"),
  make_option(c("--output"), type="character", default=NULL,
              help="Output path for HEME VST matrix (RDS)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$counts) || is.null(opt$metadata) || is.null(opt$output)) {
  stop("--counts, --metadata, and --output are required")
}

# Create output directory
output_dir <- dirname(opt$output)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

cat("Loading DSMZ count file...\n")
DSMZ_RAW_COUNT_ERROR <- paste(
  "Configured DSMZ input is not a raw integer count matrix.",
  "Do not use transformed VST/log expression for batch-correction input."
)
stop_invalid_dsmz_counts <- function(detail) {
  stop(sprintf("%s %s", DSMZ_RAW_COUNT_ERROR, detail), call. = FALSE)
}
if (!file.exists(opt$counts)) {
  stop_invalid_dsmz_counts(sprintf("Configured file does not exist: %s.", opt$counts))
}
counts_df <- tryCatch(
  readRDS(opt$counts),
  error = function(e) stop_invalid_dsmz_counts(sprintf(
    "readRDS failed for %s: %s.",
    opt$counts,
    conditionMessage(e)
  ))
)
if (!(is.matrix(counts_df) || is.data.frame(counts_df))) {
  stop_invalid_dsmz_counts("Expected matrix/data.frame-like input.")
}
if (is.null(rownames(counts_df)) || is.null(colnames(counts_df))) {
  stop_invalid_dsmz_counts("Row and column names are required.")
}
if (anyDuplicated(rownames(counts_df)) || anyDuplicated(colnames(counts_df))) {
  stop_invalid_dsmz_counts("Duplicated row or column names were detected.")
}

# Check structure
if (!"Ensembl_ID" %in% colnames(counts_df)) {
  stop("Count file must have 'Ensembl_ID' column")
}
annotation <- c("Ensembl_ID", "gene_name", "Ensembl_ID_with_version")
sample_cols <- setdiff(colnames(counts_df), annotation)
if (!length(sample_cols) || !all(vapply(counts_df[sample_cols], is.numeric, logical(1)))) {
  stop_invalid_dsmz_counts("All sample columns must be numeric.")
}
count_payload <- as.matrix(counts_df[, sample_cols, drop = FALSE])
if (anyNA(count_payload) || any(!is.finite(count_payload))) {
  stop_invalid_dsmz_counts("NA or non-finite values were detected.")
}
if (any(count_payload < 0)) {
  stop_invalid_dsmz_counts("Negative values were detected.")
}
if (any(abs(count_payload - round(count_payload)) > 1e-8)) {
  stop_invalid_dsmz_counts("Non-integer values were detected.")
}

cat("  Count file: ", nrow(counts_df), " genes, ", ncol(counts_df), " columns\n", sep="")

cat("Loading metadata...\n")
meta <- fread(opt$metadata)

# Identify HEME samples
cat("Identifying HEME samples...\n")
heme_meta <- meta[
  grepl("hematologic|leukemia|lymphoma|LL-100", group2, ignore.case=TRUE) |
  grepl("hematologic|leukemia|lymphoma", Disease, ignore.case=TRUE)
]

cat("  Found ", nrow(heme_meta), " HEME samples in metadata\n", sep="")

# Get HEME sample names (match format in count file)
heme_sample_names <- heme_meta$sample_name

# Extract HEME columns from count file
count_cols <- colnames(counts_df)
heme_cols <- count_cols[count_cols %in% heme_sample_names]

if (length(heme_cols) == 0) {
  stop("No HEME samples found in count file. Check sample name matching.")
}

cat("  Found ", length(heme_cols), " HEME samples in count file\n", sep="")

# Extract count matrix
cat("Extracting count matrix...\n")
gene_col <- "Ensembl_ID"
# Extract columns (handle both data.table and data.frame)
if (inherits(counts_df, "data.table")) {
  count_matrix <- as.matrix(counts_df[, heme_cols, with=FALSE])
} else {
  count_matrix <- as.matrix(counts_df[, heme_cols])
}
rownames(count_matrix) <- counts_df[[gene_col]]

# Remove version numbers from gene IDs (to match other VST files)
rownames(count_matrix) <- sub("\\.[0-9]+$", "", rownames(count_matrix))

cat("  Count matrix: ", nrow(count_matrix), " genes x ", ncol(count_matrix), " samples\n", sep="")

# Create DESeqDataSet
cat("Creating DESeqDataSet...\n")
col_data <- data.frame(
  sample = colnames(count_matrix),
  row.names = colnames(count_matrix)
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ 1  # No design for VST only
)

# Run VST
cat("Running variance stabilizing transformation...\n")
vst <- vst(dds, blind=TRUE)
vst_matrix <- assay(vst)

cat("  VST matrix: ", nrow(vst_matrix), " genes x ", ncol(vst_matrix), " samples\n", sep="")

# Save
cat("Saving VST matrix...\n")
saveRDS(vst_matrix, file=opt$output)

cat("[OK] HEME VST matrix saved to:", opt$output, "\n")
cat("  Final dimensions: ", nrow(vst_matrix), " genes x ", ncol(vst_matrix), " samples\n", sep="")
