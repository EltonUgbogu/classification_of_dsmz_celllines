# Assumes 'genes_gtf' is available from your previous R script's 'gene_map' chunk
library(rtracklayer)
library(dplyr)
library(readr)
gtf_path <- "/home/chu25/data/GTF/Homo_sapiens.GRCh38.107.gtf"
gtf <- rtracklayer::import(gtf_path)
genes_gtf <- gtf[gtf$type == "gene"]

gene_map_r <- as.data.frame(mcols(genes_gtf)) %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::rename(Ensembl_ID = gene_id, HGNC_Symbol = gene_name)

# Save the mapping file
write_tsv(gene_map_r, "/home/chu25/data/ensembl_to_hgnc_map.tsv")