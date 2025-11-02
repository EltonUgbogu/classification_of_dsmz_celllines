#!/usr/bin/env Rscript
# install_pkgs.R — install all packages required by breast_cancer_subtypes.R

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus = max(1L, parallel::detectCores() - 1L)
)

message("==> Library paths BEFORE:\n", paste(" -", .libPaths(), collapse = "\n"))

# Prefer a user-writable library inside the conda env (safest on HPC)
conda_lib <- file.path(Sys.getenv("CONDA_PREFIX", unset = R.home("library")), "rlib")
if (!dir.exists(conda_lib)) dir.create(conda_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(conda_lib, .libPaths())))

message("==> Library paths AFTER:\n", paste(" -", .libPaths(), collapse = "\n"))

need <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

# --- Package sets -------------------------------------------------------------
cran_pkgs <- c(
  # Core
  "tidyverse","matrixStats","pheatmap","uwot","cluster",
  "fpc","clusterCrit","patchwork","ggpubr","rstatix","RColorBrewer",
  # sva CRAN deps
  "mgcv","nlme","survival"
)

bioc_pkgs <- c(
  "DESeq2","SummarizedExperiment","edgeR","limma",
  # sva + its Bioc dep
  "genefilter","sva"
)

# --- Install BiocManager if needed -------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager from CRAN...")
  install.packages("BiocManager")
}

# --- Install CRAN packages ----------------------------------------------------
cran_missing <- need(cran_pkgs)
if (length(cran_missing)) {
  message("Installing CRAN packages: ", paste(cran_missing, collapse = ", "))
  install.packages(cran_missing)
} else {
  message("All CRAN packages already present.")
}

# --- Install Bioconductor packages -------------------------------------------
bioc_missing <- need(bioc_pkgs)
if (length(bioc_missing)) {
  message("Installing Bioconductor packages: ", paste(bioc_missing, collapse = ", "))
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
} else {
  message("All Bioconductor packages already present.")
}

BiocManager::install(c("genefilter","sva"), ask=FALSE, update=FALSE, dependencies=TRUE)


# --- Verify installs & report versions ---------------------------------------
all_pkgs <- c(cran_pkgs, bioc_pkgs)
failed <- character()

cat("\n==> Package versions:\n")
for (p in all_pkgs) {
  ok <- suppressWarnings(requireNamespace(p, quietly = TRUE))
  if (!ok) {
    cat(sprintf(" - %-24s : NOT INSTALLED\n", p))
    failed <- c(failed, p)
  } else {
    v <- as.character(utils::packageVersion(p))
    cat(sprintf(" - %-24s : %s\n", p, v))
  }
}

# Extra sanity: can we attach the key ones?
attach_test <- c("DESeq2", "SummarizedExperiment", "limma", "edgeR", "sva", "tidyverse")
attach_fail <- character()
for (p in attach_test) {
  ok <- try(suppressPackageStartupMessages(library(p, character.only = TRUE)), silent = TRUE)
  if (inherits(ok, "try-error")) attach_fail <- c(attach_fail, p)
}
if (length(attach_fail)) {
  message("\n[ERROR] Unable to attach some core packages: ", paste(attach_fail, collapse = ", "))
  failed <- unique(c(failed, attach_fail))
}

cat("\n==> sessionInfo() snapshot:\n")
print(sessionInfo())

if (length(failed)) {
  message("\nInstallation FAILED for: ", paste(failed, collapse = ", "))
  quit(status = 1)
} else {
  message("\nAll required packages are installed and loadable. ✅")
}
