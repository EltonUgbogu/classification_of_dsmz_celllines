suppressPackageStartupMessages({

  library(dplyr)
  library(tibble)
  library(data.table)
  library(SummarizedExperiment)
})

message("[INFO] Loading TARGET object ...")
target_obj <- readRDS("/home/chu25/data/nbl/TARGET-NBL/TARGET_NBL_STAR_Counts_data.rds")
#check coldata for strandedness
print(colData(target_obj))
