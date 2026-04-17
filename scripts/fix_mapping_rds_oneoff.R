#!/usr/bin/env Rscript
# One-off hotfix: rewrite cell_line_to_original_sample_id_*.rds so values are
# canonical cell line names (BC3, BJAB, ...) extracted from tech IDs, not "HEME".
# Run from pipeline dir: Rscript scripts/fix_mapping_rds_oneoff.R [path_to.rds]
# If no path given, fixes the kTotal mapping used in the doc.

infile <- if (length(commandArgs(trailingOnly = TRUE)) > 0L) {
  commandArgs(trailingOnly = TRUE)[1L]
} else {
  "results/unsupervised/heme/tumour_neighbourhoods_input/cell_line_to_original_sample_id_kTotal.rds"
}
outfile <- infile

if (!file.exists(infile)) {
  stop("File not found: ", infile)
}

m <- readRDS(infile)
keys <- names(m)

# Extract canonical cell line from key: NG-<digits>_<CELL>_lib...
cell <- sub("^NG-[0-9]+_", "", keys)
cell <- sub("_lib.*$", "", cell)
cell <- gsub("-", "_", cell)
cell <- gsub(" ", "_", cell)

# Only rewrite DSMZ keys (NG-*); leave others unchanged if any
dsmz <- grepl("^NG-", keys)
m2 <- m
m2[dsmz] <- cell[dsmz]

cat("[INFO] Rewriting mapping with extracted cell line names.\n")
cat("[INFO] First 10 (key -> value):\n")
n <- min(10L, length(m2))
for (i in seq_len(n)) {
  cat("  ", names(m2)[i], " -> ", unname(m2)[i], "\n", sep = "")
}

saveRDS(m2, outfile)
cat("[OK] Saved: ", outfile, "\n", sep = "")
cat("[INFO] Unique cell lines (DSMZ): ", length(unique(unname(m2[dsmz]))), "\n", sep = "")
