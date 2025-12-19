#!/usr/bin/env Rscript

# rebuild_joint_vst.R
# Build a true joint DSMZ+TCGA VST matrix from separate batch-corrected matrices.
# Ensures original sample IDs (NG-... for DSMZ, TCGA-... for TCGA) are preserved.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
})

cat("=== Rebuilding joint DSMZ+TCGA VST matrix for BRCA ===\n")

bcc_res <- "/home/chu25/dsmz/bcc_analysis/results"

tcga_vst_path <- file.path(bcc_res, "tcga_vst_brca_subtype_batch_corrected.rds")
dsmz_vst_path <- file.path(bcc_res, "dsmz_vst_bcc_batch_corrected.rds")
joint_out     <- file.path(bcc_res, "vst_dsmz_bcc_tcga_brca_batch_corrected.rds")

cat("[INFO] TCGA VST  :", tcga_vst_path, "\n")
cat("[INFO] DSMZ VST  :", dsmz_vst_path, "\n")
cat("[INFO] Joint out :", joint_out, "\n")

stopifnot(file.exists(tcga_vst_path))
stopifnot(file.exists(dsmz_vst_path))

tcga_vst <- readRDS(tcga_vst_path)
dsmz_vst <- readRDS(dsmz_vst_path)

cat("[INFO] TCGA VST dims:", nrow(tcga_vst), "genes x", ncol(tcga_vst), "samples\n")
cat("[INFO] DSMZ VST dims:", nrow(dsmz_vst), "genes x", ncol(dsmz_vst), "samples\n")

common_genes <- intersect(rownames(tcga_vst), rownames(dsmz_vst))
if (length(common_genes) < 10000) {
  stop("Too few common genes between DSMZ and TCGA VST matrices: ", length(common_genes))
}

tcga_vst  <- tcga_vst[common_genes, , drop = FALSE]
dsmz_vst  <- dsmz_vst[common_genes, , drop = FALSE]
joint_vst <- cbind(dsmz_vst, tcga_vst)

cat("[INFO] Joint VST dims:", nrow(joint_vst), "genes x", ncol(joint_vst), "samples\n")
cat("[INFO] DSMZ samples:", sum(grepl("^NG-", colnames(joint_vst))),
    " | TCGA samples:", sum(grepl("^TCGA-", colnames(joint_vst))), "\n")

saveRDS(joint_vst, joint_out)
cat("[OK] Saved joint VST to:\n  ", joint_out, "\n")

