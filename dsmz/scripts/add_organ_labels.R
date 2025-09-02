#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)  # brings in readr, dplyr, purrr, tidyr, stringr
})

# -----------------------------------------------------------------------------#
# CLI args
# -----------------------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[[1]] else "/home/chu25/data/dsmz/DSMZ_metadata.csv"
out_csv <- if (length(args) >= 2) args[[2]] else sub("\\.csv$", "_with_organs.csv", in_csv)

message("[INFO] Input:  ", in_csv)
message("[INFO] Output: ", out_csv)

stopifnot(file.exists(in_csv))
meta <- readr::read_csv(in_csv, show_col_types = FALSE)

# -----------------------------------------------------------------------------#
# Helpers
# -----------------------------------------------------------------------------#

# Safe lower-case + trim; convert NA to ""
safe_lc <- function(x) {
  x <- tidyr::replace_na(x, "")
  x <- as.character(x)
  stringr::str_trim(stringr::str_to_lower(x))
}

# Treat these strings as "missing" (expects already lower-cased input)
is_missing_val <- function(x) {
  x %in% c("", "na", "n/a", "none", "unknown", "unspecified", "not reported", "nan")
}

# Derive high-level organ label from Description + Disease
# Returns one of: Breast, Adrenal/SNS, Uveal/Eye, Lymphoid, Myeloid, Hematologic,
# Brain/CNS, Lung, Colon/Rectum, ... (or Other/Unknown)
organ_from_fields <- function(desc, disease) {
  d <- safe_lc(desc)
  z <- safe_lc(disease)

  # --- Solid tumor heuristics by description or disease ----------------------
  if (stringr::str_detect(d, "\\bbreast\\b") || stringr::str_detect(z, "\\bbreast\\b"))
    return("Breast")
  if (stringr::str_detect(d, "neuroblastoma") || stringr::str_detect(z, "neuroblastoma"))
    return("Adrenal/SNS")
  if (stringr::str_detect(d, "retinoblastoma") || stringr::str_detect(z, "retinoblastoma"))
    return("Uveal/Eye")

  # --- Hematologic: disease-driven (overrides panels) ------------------------
  # Lymphoid patterns (incl. ALCL)
  if (stringr::str_detect(z, "primary effusion lymphoma")) return("Lymphoid")
  if (stringr::str_detect(z, "hairy cell leukemia"))       return("Lymphoid")
  if (stringr::str_detect(z, "t-?acute lymphoblastic|t-?lymphoblastic|\\bt[- ]?all\\b|\\bt-?ll\\b"))
    return("Lymphoid")
  if (stringr::str_detect(z, "b-?cell precursor|pre-?b-?acute lymphoblastic|\\bpre-?b-?all\\b"))
    return("Lymphoid")
  if (stringr::str_detect(z, "diffuse large b-?cell lymphoma|\\bdlbcl\\b|primary mediastinal b-?cell"))
    return("Lymphoid")
  if (stringr::str_detect(z, "mantle cell lymphoma"))      return("Lymphoid")
  if (stringr::str_detect(z, "hodgkin lymphoma"))          return("Lymphoid")
  if (stringr::str_detect(z, "multiple myeloma|plasma cell leukemia"))
    return("Lymphoid")
  if (stringr::str_detect(z, "\\bnatural\\s*killer\\b|\\bnk\\b"))
    return("Lymphoid")
  if (stringr::str_detect(z, "chronic myeloid leukemia.*lymphoid blast"))
    return("Lymphoid")
  if (stringr::str_detect(z, "anaplastic large cell lymphoma|\\balcl\\b"))
    return("Lymphoid")

  # Myeloid patterns
  if (stringr::str_detect(z, "^acute myeloid leukemia|\\baml\\b|myelocytic|monocytic|megakaryocytic|erythroid"))
    return("Myeloid")
  if (stringr::str_detect(z, "chronic myeloid leukemia.*myeloid blast"))
    return("Myeloid")

  # --- Broad heme panel with no disease detail -------------------------------
  if (stringr::str_detect(d, "^leukemia and lymphoma panel$") && is_missing_val(z))
    return("Hematologic")

  # --- Explicit organ/site hints (optional expansions) -----------------------
  if (stringr::str_detect(d, "\\blung\\b") || stringr::str_detect(z, "\\blung\\b"))
    return("Lung")
  if (stringr::str_detect(d, "\\bcolon|rectum\\b") || stringr::str_detect(z, "\\bcolon|rectum\\b"))
    return("Colon/Rectum")
  if (stringr::str_detect(d, "\\bbrain|cns\\b") || stringr::str_detect(z, "\\bbrain|cns\\b"))
    return("Brain/CNS")

  # --- Missing/unknown -------------------------------------------------------
  if (is_missing_val(d) && is_missing_val(z)) return("Unknown")

  # --- Fallback --------------------------------------------------------------
  "Other"
}

# Map normalized organ to candidate TCGA cohort codes
# Note: TCGA has DLBC (B-cell lymphoma) and LAML (AML). No T-cell lymphoma cohort.
tcga_from_organ <- function(org) {
  switch(org,
    "Adrenal"        = c("ACC"),
    "Adrenal/SNS"    = c("ACC","PCPG"),     # nearest adult cohorts for neuroblastoma/adr. axis
    "Bladder"        = c("BLCA"),
    "Brain/CNS"      = c("LGG","GBM"),
    "Breast"         = c("BRCA"),
    "Cervix"         = c("CESC"),
    "Liver/Biliary"  = c("LIHC","CHOL"),
    "Colon/Rectum"   = c("COAD","READ"),
    "Esophagus"      = c("ESCA"),
    "Head/Neck"      = c("HNSC"),
    "Kidney"         = c("KICH","KIRC","KIRP"),
    "Lung"           = c("LUAD","LUSC"),
    "Mesothelium"    = c("MESO"),
    "Ovary"          = c("OV"),
    "Pancreas"       = c("PAAD"),
    "Prostate"       = c("PRAD"),
    "Sarcoma"        = c("SARC"),
    "Skin"           = c("SKCM"),
    "Stomach"        = c("STAD"),
    "Testis"         = c("TGCT"),
    "Thyroid"        = c("THCA"),
    "Thymus"         = c("THYM"),
    "Uterus"         = c("UCEC","UCS"),
    "Uveal/Eye"      = c("UVM"),
    # Hematologic groupings
    "Lymphoid"       = c("DLBC"),
    "Myeloid"        = c("LAML"),           # keep to LAML; LCML is rarely present
    "Hematologic"    = c("DLBC","LAML"),
    # Fallback/unknown → no mapping
    character(0)
  )
}

# -----------------------------------------------------------------------------#
# Validate required columns and compute labels
# -----------------------------------------------------------------------------#
need_cols <- c("Description", "Disease")
missing_cols <- setdiff(need_cols, names(meta))
if (length(missing_cols)) {
  stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
}

meta2 <- meta %>%
  mutate(
    organ = purrr::pmap_chr(list(Description, Disease), organ_from_fields),
    tcga_codes = purrr::map_chr(organ, ~ paste(tcga_from_organ(.x), collapse = ";"))
  )

# -----------------------------------------------------------------------------#
# Write output & print a small summary
# -----------------------------------------------------------------------------#
readr::write_csv(meta2, out_csv)
message("[OK] Wrote: ", out_csv)

summary_tbl <- meta2 %>% count(organ, sort = TRUE)
message("[INFO] Organ label counts:")
print(summary_tbl, n = nrow(summary_tbl))
