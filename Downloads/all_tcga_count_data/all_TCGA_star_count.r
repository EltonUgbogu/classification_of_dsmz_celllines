# Create personal library directory 
#dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
#.libPaths(Sys.getenv("R_LIBS_USER"))

# Set CRAN mirror explicitly 
#options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load necessary library
library(TCGAbiolinks)
library(SummarizedExperiment)

# Query all TCGA RNA-seq counts
query_all <- GDCquery(
  project = getGDCprojects()$project_id[grepl("TCGA", getGDCprojects()$project_id)],
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download data
#GDCdownload(query_all)

# Prepare data 
tcga_se <- GDCprepare(query_all)

# Save as .rds 
saveRDS(tcga_se, "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment.rds")






