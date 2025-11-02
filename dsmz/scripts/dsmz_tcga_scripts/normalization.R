# Normalization functions (VST per-dataset and ComBat_seq merge)
# Place in R/normalization.R

# counts_matrix - counts matrix 
# coldata - data frame of counts data
vst_per_dataset <- function(counts_matrix, coldata, design = ~1, blind = TRUE) {
  cat("[INFO] VST normalizing counts...\n")
  dds <- DESeqDataSetFromMatrix(round(counts_matrix), coldata, design = design)
  vsd <- vst(dds, blind = blind)
  assay(vsd)
}

# counts_matrix - counts matrix from combined counts matrix 
# batch_factor - batch factor from merge_counts function
#seed - random seed for ComBat_seq
batch_adjust <- function(counts_matrix, batch_factor, seed = 42) {
  cat("[INFO] Batch adjusting counts...\n")
  set.seed(seed)
  batch_corrected_counts <- sva::ComBat_seq(counts = as.matrix(counts_matrix),
                                      batch = batch_factor, group = NULL,
                                      covar_mod = NULL, full_mod = TRUE,
                                      shrink = FALSE, shrink.disp = FALSE)
  return(batch_corrected_counts)
}
