#!/usr/bin/env Rscript
# mycn_gene_presence.R
counts <- read.table("/work/ugbogu/pipeline/results/unsupervised/nbl/deseq2_inputs/counts.tsv",
                     header=TRUE, sep="\t", check.names=FALSE)

# Check raw presence
any(grepl("ENSG00000134323", counts$gene_id))

# Print matching rows
counts[grep("ENSG00000134323", counts$gene_id), 1]
