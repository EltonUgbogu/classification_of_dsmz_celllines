library(TCGAbiolinks)
library(SummarizedExperiment)

# Define a list of TCGA projects
projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD",
  "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC",
  "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
  "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ",
  "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)

# Initialize a list to store SummarizedExperiment objects
se_list <- list()

# Local data directory where GDCdownload stored the files
data_dir <- "~/TCGA/all_tcga_count_data/GDCdata"

# Loop over each project
for (project in projects) {
  message("Processing project: ", project)

  # Create GDC query
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  # Prepare data using GDCprepare and specify local directory
  se <- GDCprepare(query, directory = data_dir)

  # Store in the list
  se_list[[project]] <- se
}

# Save all SummarizedExperiment objects into one RDS file
saveRDS(se_list, file = "~/data/tcga/all_TCGA_SummarizedExperiment.rds")

message("Done saving all SE objects.")

