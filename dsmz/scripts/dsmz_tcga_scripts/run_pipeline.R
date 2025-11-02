.libPaths("/home/chu25/miniconda3/envs/tcga-r-env/lib/R/library")
Sys.setenv(R_LIBS_USER="/home/chu25/miniconda3/envs/tcga-r-env/lib/R/library")

suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(DESeq2)
  library(matrixStats)
  library(pheatmap)
  library(uwot)
  library(cluster)
  library(fpc)
  library(clusterCrit)
  library(limma)
  library(sva)
  library(dbscan)
  library(PMCMRplus)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dynamicTreeCut)  # for cutreeDynamic / cutreeHybrid
  library(WGCNA)           # for plotDendroAndColors
  library(paran)           # for Parallel Analysis for PCs
})
