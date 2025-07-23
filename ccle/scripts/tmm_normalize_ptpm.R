#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
ptpm_file <- args[1]
output_file <- args[2]

library(NOISeq)

# Read TPM
expr <- read.table(ptpm_file, header=TRUE, row.names=1, sep="\t")

# TMM normalization using NOISeq::tmm
norm_expr <- tmm(expr)

# Write normalized expression
write.table(norm_expr, file=output_file, sep="\t", quote=FALSE)
