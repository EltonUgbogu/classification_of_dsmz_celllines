# Check for strand-specificity batch effects
library(readr)
library(dplyr)
library(pheatmap)

# Load metadata
meta <- read_csv("/home/chu25/data/nbl/final_nbl_sample_metadata.csv", show_col_types = FALSE)

# Count library types
print(table(meta$strand_info))

# Load VST matrix 
vst_mat <- readRDS("/home/chu25/TCGA/nbl/nbl_vst_matrix.rds")   # Fixed path

# Create annotation dataframe
annot <- meta %>%
  select(sample, strand_info) %>%
  column_to_rownames("sample")

# Sample 5000 genes
set.seed(123)
genes_plot <- sample(rownames(vst_mat), 5000)

# Generate heatmap with full path
pheatmap(vst_mat[genes_plot, ],
         annotation_col = annot,
         annotation_colors = list(strand_info = c(`reverse/-` = "#56B4E9", unstranded = "#E69F00")),
         show_colnames = FALSE,
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         filename = "/home/chu25/TCGA/nbl/results/check_strand_batch.pdf",
         width = 12, height = 10,
         main = "NBL Cohort: Strand Batch Effect Check (n=120)")

message("Strand batch check complete! Plot saved to /home/chu25/TCGA/nbl/results/check_strand_batch.pdf")