#!/usr/bin/env Rscript

# rebuild_joint_vst.R
# Build a true joint DSMZ + tumour VST matrix from separate batch-corrected matrices.
# Uses config.yaml + --profile so it works for BRCA/NBL/RBL.
# Preserves original sample IDs (NG-... for DSMZ, TCGA-... or tumour IDs for tumours).

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

option_list <- list(
  make_option(c("--config"),  type="character", default="config/config.yaml",
              help="Path to config YAML"),
  make_option(c("--profile"), type="character", default=NA,
              help="Profile name (brca, nbl, rbl). Overrides config$profile")
)

opt <- parse_args(OptionParser(option_list = option_list))

cfg_all <- yaml::read_yaml(opt$config)

profile <- opt$profile
if (is.na(profile) || !nzchar(profile)) {
  profile <- cfg_all$profile
}
if (is.null(profile) || !nzchar(profile)) {
  stop("No profile given and config.yaml has no top-level `profile:`")
}
if (is.null(cfg_all$profiles[[profile]])) {
  stop("Unknown profile: ", profile)
}

# merge defaults + profile (simple deep merge)
deep_merge <- function(base, override) {
  out <- base
  for (nm in names(override)) {
    if (is.list(override[[nm]]) && is.list(out[[nm]])) {
      out[[nm]] <- deep_merge(out[[nm]], override[[nm]])
    } else {
      out[[nm]] <- override[[nm]]
    }
  }
  out
}

cfg <- deep_merge(cfg_all$defaults, cfg_all$profiles[[profile]])

cancer_type <- cfg$analysis$cancer_type
cell_vst_path   <- cfg$paths$cell_vst_rds
tumour_vst_path <- cfg$paths$tumour_vst_rds
joint_out       <- cfg$paths$vst_joint_rds

cat("=== Rebuilding joint VST matrix ===\n")
cat("[INFO] Profile   :", profile, "\n")
cat("[INFO] Cancer    :", cancer_type, "\n")
cat("[INFO] Tumour VST:", tumour_vst_path, "\n")
cat("[INFO] Cell  VST :", cell_vst_path, "\n")
cat("[INFO] Joint out :", joint_out, "\n")

stopifnot(file.exists(tumour_vst_path))
stopifnot(file.exists(cell_vst_path))

tumour_vst <- readRDS(tumour_vst_path)
cell_vst   <- readRDS(cell_vst_path)

cat("[INFO] Tumour VST dims:", nrow(tumour_vst), "genes x", ncol(tumour_vst), "samples\n")
cat("[INFO] Cell   VST dims:", nrow(cell_vst),   "genes x", ncol(cell_vst),   "samples\n")

common_genes <- intersect(rownames(tumour_vst), rownames(cell_vst))
cat("[INFO] Common genes:", length(common_genes), "\n")
if (length(common_genes) < 10000) {
  stop("Too few common genes between tumour and cell VST matrices: ", length(common_genes))
}

tumour_vst <- tumour_vst[common_genes, , drop = FALSE]
cell_vst   <- cell_vst[common_genes, , drop = FALSE]

joint_vst <- cbind(cell_vst, tumour_vst)

cat("[INFO] Joint VST dims:", nrow(joint_vst), "genes x", ncol(joint_vst), "samples\n")

n_dsmz <- sum(grepl("^NG-", colnames(joint_vst)))
n_tcga <- sum(grepl("^TCGA-", colnames(joint_vst)))
cat("[INFO] Detected DSMZ (NG-):", n_dsmz, "\n")
cat("[INFO] Detected TCGA (TCGA-):", n_tcga, "\n")

dir.create(dirname(joint_out), recursive = TRUE, showWarnings = FALSE)
saveRDS(joint_vst, joint_out)

cat("[OK] Saved joint VST to:\n  ", joint_out, "\n", sep = "")
