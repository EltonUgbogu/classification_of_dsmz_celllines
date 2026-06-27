#!/usr/bin/env Rscript

# list_nbl_sample_ids_by_cohort.R
#
# Reads tumour_vst_nbl_batch_corrected.rds, extracts sample IDs,
# derives cohort from the sample ID prefix, and writes one TSV per cohort.
#
# Expected cohort examples:
#   GSE100148
#   GSE189367
#   GSE49711
#   SRP409177
#   TARGET
#
# Usage:
#   Rscript list_nbl_sample_ids_by_cohort.R \
#       data/nbl/new_data/preprocessing_results/tumour_vst_nbl_batch_corrected.rds \
#       data/nbl/new_data/preprocessing_results/sample_id_exports

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), "../../.."), mustWork = FALSE)
nbl_data_root <- Sys.getenv("NBL_DATA_ROOT", unset = file.path(repo_root, "data", "nbl", "new_data"))

input_rds <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(nbl_data_root, "preprocessing_results", "tumour_vst_nbl_batch_corrected.rds")
}

outdir <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(nbl_data_root, "preprocessing_results", "sample_id_exports")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)

extract_sample_ids <- function(x) {
  if (inherits(x, "SummarizedExperiment")) {
    ids <- colnames(x)
    if (is.null(ids)) {
      stop("SummarizedExperiment has no column names.")
    }
    return(as.character(ids))
  }

  if (is.matrix(x) || is.data.frame(x)) {
    ids <- colnames(x)
    if (is.null(ids)) {
      stop("Object has no column names.")
    }
    return(as.character(ids))
  }

  if (is.list(x)) {
    for (nm in c("assay", "expr", "vst", "counts", "mat")) {
      if (!is.null(x[[nm]]) && !is.null(colnames(x[[nm]]))) {
        return(as.character(colnames(x[[nm]])))
      }
    }

    if (!is.null(x$sample_id)) {
      return(as.character(x$sample_id))
    }
  }

  stop("Unsupported object type: ", paste(class(x), collapse = ", "))
}

derive_cohort <- function(ids) {
  ids <- as.character(ids)

  cohort <- ifelse(
    grepl("^TARGET-", ids),
    "TARGET",
    sub("_.*$", "", ids)
  )

  cohort
}

sample_ids <- extract_sample_ids(obj)
cohort <- derive_cohort(sample_ids)

df <- data.frame(
  cohort = cohort,
  sample_id = sample_ids,
  stringsAsFactors = FALSE
)

df <- df[order(df$cohort, df$sample_id), , drop = FALSE]

write.table(
  df,
  file = file.path(outdir, "all_sample_ids_by_cohort.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

count_df <- as.data.frame(table(df$cohort), stringsAsFactors = FALSE)
colnames(count_df) <- c("cohort", "n_samples")
count_df <- count_df[order(count_df$cohort), , drop = FALSE]

write.table(
  count_df,
  file = file.path(outdir, "sample_id_counts_by_cohort.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

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
message("Input RDS: ", input_rds)
message("Output directory: ", outdir)
message("")
message("Counts per cohort:")
print(count_df, row.names = FALSE)
