#!/usr/bin/env Rscript
# compute_representativeness_score.R
# Layer 3c: Compute single scalar representativeness score per cell line
# Derived from tumour-informed cell-line characterisation (Layer 3b)

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
})

option_list <- list(
  make_option(c("--cell_char_tsv"), type="character",
              help="cell_line_characterisation.tsv (from Layer 3b)"),
  make_option(c("--out_tsv"), type="character",
              help="Output TSV for ranked representativeness"),
  make_option(c("--score_mode"), type="character", default="pXspec",
              help="pXspec (default: p_consensus × specificity), p_only, or spec_only")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$cell_char_tsv) || !file.exists(opt$cell_char_tsv)) {
  stop("[ERROR] --cell_char_tsv is required and must exist: ", opt$cell_char_tsv)
}
if (is.null(opt$out_tsv) || !nzchar(opt$out_tsv)) {
  stop("[ERROR] --out_tsv is required")
}

cat("[INFO] Loading cell-line characterisation from:", opt$cell_char_tsv, "\n")
df <- read_tsv(opt$cell_char_tsv, show_col_types = FALSE)

req <- c("cell_line", "dominant_cluster", "p_consensus", "specificity")
miss <- setdiff(req, colnames(df))
if (length(miss) > 0) {
  stop("[ERROR] Missing required columns: ", paste(miss, collapse=", "))
}

cat("[INFO] Computing representativeness score (mode: ", opt$score_mode, ")...\n", sep="")

df <- df %>%
  mutate(
    p_consensus = as.numeric(p_consensus),
    specificity = as.numeric(specificity),
    rep_score = case_when(
      opt$score_mode == "p_only"    ~ p_consensus,
      opt$score_mode == "spec_only" ~ specificity,
      TRUE                          ~ p_consensus * specificity
    )
  ) %>%
  arrange(desc(rep_score), desc(p_consensus), desc(specificity))

cat("[INFO] Score range: [", min(df$rep_score, na.rm=TRUE), ", ", max(df$rep_score, na.rm=TRUE), "]\n", sep="")
cat("[INFO] Top 5 cell lines by representativeness:\n")
print(df %>% head(5) %>% select(cell_line, dominant_cluster, p_consensus, specificity, rep_score))

# Write ranked table
dir.create(dirname(opt$out_tsv), recursive = TRUE, showWarnings = FALSE)
write_tsv(df, opt$out_tsv)
cat("[SUCCESS] Wrote:", opt$out_tsv, "\n")

# Per-program top-k list (always write if output path suggests it)
out_top10 <- sub("ranked.tsv$", "top10_per_cluster.tsv", opt$out_tsv)
if (out_top10 != opt$out_tsv) {
  df_top10 <- df %>%
    group_by(dominant_cluster) %>%
    slice_max(rep_score, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  write_tsv(df_top10, out_top10)
  cat("[SUCCESS] Wrote top-10 per cluster:", out_top10, "\n")
} else {
  # If output path doesn't match pattern, write to same directory
  out_top10 <- file.path(dirname(opt$out_tsv), "cell_line_representativeness_top10_per_cluster.tsv")
  df_top10 <- df %>%
    group_by(dominant_cluster) %>%
    slice_max(rep_score, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  write_tsv(df_top10, out_top10)
  cat("[SUCCESS] Wrote top-10 per cluster:", out_top10, "\n")
}

