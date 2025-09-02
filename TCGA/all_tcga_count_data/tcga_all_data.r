# ----------------------------------------
# Setup: Personal R Library & CRAN Mirror
# ----------------------------------------
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(Sys.getenv("R_LIBS_USER"))
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ----------------------------------------
# Load Required Libraries
# ----------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)

# ----------------------------------------
# Query: All TCGA Projects with RNA-seq (STAR - Counts)
# ----------------------------------------

# Get all GDC project IDs
gdc_projects <- getGDCprojects()$project_id

# Filter for TCGA-only projects
tcga_projects <- gdc_projects[grepl("^TCGA-", gdc_projects)]

# Define the query object
query_all <- GDCquery(
  project = tcga_projects,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# ----------------------------------------
# Download Data (can be very large)
# ----------------------------------------
# This may take hours and consume >100 GB
GDCdownload(query_all)

# ----------------------------------------
# Prepare SummarizedExperiment Object
# ----------------------------------------
# ⚠️ Warning: This may crash your session if RAM is insufficient
tcga_se <- GDCprepare(query_all)

# ----------------------------------------
# Save Result to .rds File
# ----------------------------------------
saveRDS(tcga_se, file = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment.rds")
