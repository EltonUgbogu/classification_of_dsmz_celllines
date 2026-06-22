#!/usr/bin/env Rscript
# mycn_gene_presence.R
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
counts_file <- if (length(args) >= 1) args[[1]] else file.path(
  repo_root, "results", "unsupervised", "nbl", "deseq2_inputs", "counts.tsv"
)
counts <- read.table(counts_file,
                     header=TRUE, sep="\t", check.names=FALSE)

# Check raw presence
any(grepl("ENSG00000134323", counts$gene_id))

# Print matching rows
counts[grep("ENSG00000134323", counts$gene_id), 1]
