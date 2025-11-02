# Normalization functions (VST per-dataset and ComBat_seq merge)
# Place in R/normalization.R

vst_per_dataset <- function(counts, design = ~1, blind = TRUE) {
  dds <- DESeqDataSetFromMatrix(round(counts), data.frame(row.names = colnames(counts)), design = design)
  vsd <- vst(dds, blind = blind)
  assay(vsd)
}

combat_seq_merge <- function(common_genes_counts, V_tcga_norm, V_dsmz_norm, seed = 42) {
  batch <- factor(c(rep("TCGA", ncol(V_tcga_norm)), rep("DSMZ", ncol(V_dsmz_norm))))
  set.seed(seed)
  all_counts_combat <- sva::ComBat_seq(counts = as.matrix(common_genes_counts),
                                      batch = batch, group = NULL,
                                      covar_mod = NULL, full_mod = TRUE,
                                      shrink = FALSE, shrink.disp = FALSE)
  # Convert back to VST
  coldata <- tibble(sample = colnames(all_counts_combat),
                    dataset = ifelse(colnames(all_counts_combat) %in% colnames(V_tcga_norm), "TCGA", "DSMZ")) %>%
    column_to_rownames("sample")
  dds_combat <- DESeqDataSetFromMatrix(round(all_counts_combat), coldata, design = ~ 1)
  vsd_combat <- vst(dds_combat, blind = TRUE)
  assay(vsd_combat)
}