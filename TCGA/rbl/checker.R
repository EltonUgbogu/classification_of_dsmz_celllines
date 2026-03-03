#!/usr/bin/env Rscript

x <- readRDS("/home/chu25/data/rbl/rbl_raw_count_matrix.rds")
str(x)
# and check:
cat("class:", class(x), "\n")
cat("rownames length:", length(rownames(x)), "\n")
cat("colnames length:", length(colnames(x)), "\n")
