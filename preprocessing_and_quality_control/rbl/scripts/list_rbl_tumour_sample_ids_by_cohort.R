#!/usr/bin/env Rscript

# list_rbl_sample_ids_by_gse.R
#
# Reads tumour_vst_rbl_batch_corrected.rds, extracts sample IDs,
# derives cohort from the prefix before the first underscore,
# and writes one TSV per cohort:
#
# sample_ids_cohorts/
#   GSE111168_sample_ids.tsv
#   GSE196420_sample_ids.tsv
#   GSE268136_sample_ids.tsv
#
# Usage:
#   Rscript list_rbl_sample_ids_by_gse.R \
#       data/rbl/preprocessing_results/tumour_vst_rbl_batch_corrected.rds \
#       data/rbl/preprocessing_results/sample_ids_cohorts

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), "../../.."), mustWork = FALSE)
rbl_data_root <- Sys.getenv("RBL_DATA_ROOT", unset = file.path(repo_root, "data", "rbl"))

input_rds <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(rbl_data_root, "preprocessing_results", "tumour_vst_rbl_batch_corrected.rds")
}

outdir <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(rbl_data_root, "preprocessing_results", "sample_ids_cohorts")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)

extract_sample_ids <- function(x) {
  if (inherits(x, "SummarizedExperiment")) {
    return(colnames(x))
  }

  if (is.matrix(x) || is.data.frame(x)) {
    ids <- colnames(x)
    if (is.null(ids)) {
      stop("Object has no column names.")
    }
    return(ids)
  }

  if (is.list(x)) {
    for (nm in c("assay", "expr", "vst", "counts", "mat")) {
      if (!is.null(x[[nm]]) && !is.null(colnames(x[[nm]]))) {
        return(colnames(x[[nm]]))
      }
    }

    if (!is.null(x$sample_id)) {
      return(as.character(x$sample_id))
    }
  }

  stop("Unsupported object type: ", paste(class(x), collapse = ", "))
}

sample_ids <- extract_sample_ids(obj)
sample_ids <- as.character(sample_ids)

# Derive cohort from prefix before first underscore
cohort <- sub("_.*$", "", sample_ids)

df <- data.frame(
  cohort = cohort,
  sample_id = sample_ids,
  stringsAsFactors = FALSE
)

df <- df[order(df$cohort, df$sample_id), , drop = FALSE]

# Write combined table
write.table(
  df,
  file = file.path(outdir, "all_sample_ids_by_cohort.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Write counts
count_df <- as.data.frame(table(df$cohort), stringsAsFactors = FALSE)
colnames(count_df) <- c("cohort", "n_samples")
write.table(
  count_df,
  file = file.path(outdir, "sample_id_counts_by_cohort.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Write one file per cohort
for (coh in unique(df$cohort)) {
  subdf <- df[df$cohort == coh, "sample_id", drop = FALSE]
  outfile <- file.path(outdir, paste0(coh, "_sample_ids.tsv"))
  write.table(
    subdf,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

message("Done.")
message("Output directory: ", outdir)
message("")
message("Counts per cohort:")
print(count_df, row.names = FALSE)
