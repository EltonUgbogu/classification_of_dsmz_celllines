#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
})

# ---------- CLI ARG / CONFIG ----------
args <- commandArgs(trailingOnly = TRUE)
in_rds <- if (length(args) >= 1) args[1] else "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds"
outdir <- if (length(args) >= 2) args[2] else "/home/chu25/TCGA/results/tumour_purity_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("[INFO] Reading:", in_rds, "\n")
obj <- readRDS(in_rds)

# ---------- Extract colData and sample IDs ----------
as_df <- function(x) {
  # Coerce S4 / list / factor to plain vectors for clean CSVs
  if (inherits(x, "Rle")) as.vector(x)
  else if (is.factor(x)) as.character(x)
  else if (is.list(x))  vapply(x, function(y) paste(y, collapse=";"), character(1))
  else x
}

if (inherits(obj, "SummarizedExperiment")) {
  cd <- as.data.frame(colData(obj))
  sample_id <- colnames(assay(obj)) # first assay
  project_hint <- NA_character_
} else if (is.list(obj)) {
  # If the RDS is a list of SEs (per project), row-bind their colData
  parts <- lapply(names(obj), function(p){
    se <- obj[[p]]
    stopifnot(inherits(se, "SummarizedExperiment"))
    cdf <- as.data.frame(colData(se))
    cdf$..sample_id_tmp <- colnames(assay(se))
    cdf$..source_project <- p
    cdf
  })
  cd <- do.call(rbind, parts)
  sample_id <- cd$..sample_id_tmp
  project_hint <- cd$..source_project
} else {
  stop("[ERROR] Object is neither SummarizedExperiment nor list of SEs.")
}

cd[] <- lapply(cd, as_df)

# ---------- Find candidate purity columns ----------
nm <- names(cd)
cand_idx <- grepl("purity", nm, ignore.case = TRUE) | nm %in% c("CPE", "purityEstimate", "Purity", "ABSOLUTE.Purity")
cands <- nm[cand_idx]

if (!length(cands)) {
  cat("[WARN] No columns with name containing 'purity'. Also checked {CPE, purityEstimate, Purity, ABSOLUTE.Purity}.\n")
} else {
  cat("[INFO] Candidate purity columns found:", paste(cands, collapse=", "), "\n")
}

# ---------- Choose the best purity column (most non-NA numeric) ----------
choose_best <- function(df, cols) {
  if (!length(cols)) return(NULL)
  stat <- lapply(cols, function(cl){
    x <- suppressWarnings(as.numeric(df[[cl]]))
    list(name = cl,
         non_na = sum(!is.na(x)),
         min = suppressWarnings(min(x, na.rm = TRUE)),
         max = suppressWarnings(max(x, na.rm = TRUE)),
         ex = x)
  })
  # pick column with the most non-NA values; tie-breaker: larger max
  non_na_counts <- sapply(stat, `[[`, "non_na")
  if (all(non_na_counts == 0)) return(NULL)
  best_i <- order(non_na_counts, sapply(stat, `[[`, "max"), decreasing = TRUE)[1]
  stat[[best_i]]
}

best <- choose_best(cd, cands)

# ---------- Build output table ----------
df <- data.frame(
  sample_id = sample_id,
  stringsAsFactors = FALSE
)

# Add a project_id if available
proj_col <- intersect(c("project_id","study","project","Project"), names(cd))
if (length(proj_col)) {
  df$project_id <- as.character(cd[[proj_col[1]]])
} else if (exists("project_hint")) {
  df$project_id <- project_hint
} else {
  df$project_id <- NA_character_
}

# Attach every candidate purity column (as character, for inspection)
for (cl in cands) {
  df[[cl]] <- as.character(cd[[cl]])
}

# Add the chosen purity as a clean numeric column
if (is.null(best)) {
  df$purity_selected <- NA_real_
  cat("[INFO] No usable numeric purity column found.\n")
} else {
  df$purity_selected <- best$ex
  cat(sprintf("[OK] Selected purity column: %s (non-NA: %d)\n",
              best$name, best$non_na))
  if (best$non_na > 0) {
    cat(sprintf("     Range (min..max): %s .. %s\n",
                format(best$min, digits = 4), format(best$max, digits = 4)))
  }
}

# ---------- Summaries & export ----------
summary_path <- file.path(outdir, "purity_presence_summary.txt")
sink(summary_path)
cat("=== Purity presence summary ===\n")
cat("RDS:", in_rds, "\n")
cat("Candidates:", if (length(cands)) paste(cands, collapse=", ") else "none", "\n")
if (!is.null(best)) {
  cat("Selected:", best$name, "\n")
  cat("Non-NA:", best$non_na, "\n")
  cat("Min..Max:", best$min, "..", best$max, "\n")
}
# quick table of non-NA by project
if ("project_id" %in% names(df)) {
  cat("\nNon-NA counts by project_id:\n")
  tab <- with(df, tapply(!is.na(purity_selected), project_id, sum))
  print(tab)
}
sink()
cat("[OK] Wrote summary:", summary_path, "\n")

csv_path <- file.path(outdir, "TCGA_purity_extract.csv")
write.csv(df, csv_path, row.names = FALSE)
cat("[OK] Wrote CSV:", csv_path, "\n")

# Optional: warn if purity values look out of typical range [0,1]
if (any(df$purity_selected < 0 | df$purity_selected > 1, na.rm = TRUE)) {
  cat("[WARN] Some purity values are outside [0,1]. Check the selected column and scaling.\n")
}

cat("[DONE] Purity check complete.\n")
