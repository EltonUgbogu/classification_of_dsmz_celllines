# Install and use TCGA data packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

# Query BRCA clinical data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Clinical",
                  file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, "patient")
