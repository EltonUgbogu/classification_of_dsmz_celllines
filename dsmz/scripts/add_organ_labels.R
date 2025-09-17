#!/usr/bin/env Rscript
# Shebang line to run script with Rscript interpreter

suppressPackageStartupMessages({
# Suppress startup messages for cleaner output
  library(tidyverse)  # Load tidyverse: readr, dplyr, tidyr, stringr, purrr
})

# ----------------------------- CLI -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
# Get command-line arguments, excluding Rscript internals
in_csv <- if (length(args) >= 1) args[[1]] else "/home/chu25/data/dsmz/DSMZ_metadata.csv"
# Set input CSV path from args or use default
out_csv <- if (length(args) >= 2) args[[2]] else sub("\\.csv$", "_with_organs_tcga.csv", in_csv)
# Set output CSV path from args or derive from input
message("[INFO] Input:  ", in_csv)
# Print input file path for logging
message("[INFO] Output: ", out_csv)
# Print output file path for logging
stopifnot(file.exists(in_csv))
# Stop if input CSV file does not exist
meta <- readr::read_csv(in_csv, show_col_types = FALSE)
# Read input CSV into meta dataframe, hide column type messages

# ---------------------------- Helpers ----------------------------------------
safe_lc <- function(x) {
# Define function to safely lowercase and trim strings
  x <- tidyr::replace_na(x, "")
# Replace NA with empty string
  stringr::str_trim(stringr::str_to_lower(as.character(x)))
# Trim and convert to lowercase
}
is_missing_val <- function(x) x %in% c("", "na", "n/a", "none", "unknown", "unspecified", "not reported", "nan")
# Check if value indicates missing data

# Disease pattern helpers
is_aml <- function(x) str_detect(x, regex("^acute myeloid leukemia|\\baml\\b|myelocytic|monocytic|megakaryocytic|erythroid", TRUE))
# Detect AML patterns in strings
is_cml <- function(x) str_detect(x, regex("^chronic myeloid leukemia", TRUE))
# Detect CML patterns in strings
is_all <- function(x) str_detect(x, regex("acute lymphoblastic|lymphoblastic|\\ball\\b", TRUE))
# Detect ALL patterns in strings
is_dlbcl <- function(x) str_detect(x, regex("diffuse large b-?cell lymphoma|\\bdlbcl\\b|primary mediastinal b-?cell", TRUE))
# Detect DLBCL patterns in strings
is_hodg <- function(x) str_detect(x, regex("hodgkin", TRUE))
# Detect Hodgkin lymphoma patterns
is_mcl <- function(x) str_detect(x, regex("mantle cell lymphoma", TRUE))
# Detect mantle cell lymphoma patterns
is_mm <- function(x) str_detect(x, regex("myeloma|plasma cell leukemia", TRUE))
# Detect multiple myeloma patterns
is_nk <- function(x) str_detect(x, regex("natural\\s*killer|\\bnk\\b", TRUE))
# Detect NK cell malignancy patterns
is_pel <- function(x) str_detect(x, regex("primary effusion", TRUE))
# Detect primary effusion lymphoma patterns
is_hairy <- function(x) str_detect(x, regex("^hairy cell", TRUE))
# Detect hairy cell leukemia patterns
is_neuro <- function(x) str_detect(x, regex("^neuroblastoma", TRUE))
# Detect neuroblastoma patterns
is_rb <- function(x) str_detect(x, regex("^retinoblastoma", TRUE))
# Detect retinoblastoma patterns

# Valid TCGA project IDs (RNA-seq cohorts)
valid_tcga <- c(
# Define vector of valid TCGA project codes
  "ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
# First set of TCGA codes
  "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
# Second set of TCGA codes
  "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
# Third set of TCGA codes
  "UCEC","UCS","UVM","LAML"
# Final set of TCGA codes
)

# Diseases with no matching TCGA cohort -> force NONE
no_tcga_exact <- c(
# Define diseases with no TCGA match
  "Hairy cell leukemia",
# Hairy cell leukemia
  "Chronic myeloid leukemia, myeloid blast crisis",
# CML myeloid blast crisis
  "Chronic myeloid leukemia, lymphoid blast crisis",
# CML lymphoid blast crisis
  "B-cell precursor / pre-B-acute lymphoblastic leukemia",
# B-cell ALL
  "T-acute lymphoblastic leukemia / T-lymphoblastic lymphoma",
# T-cell ALL/lymphoma
  "Hodgkin lymphoma",
# Hodgkin lymphoma
  "Mantle cell lymphoma",
# Mantle cell lymphoma
  "Multiple myeloma, plasma cell leukemia",
# Multiple myeloma
  "Natural killer malignancy",
# NK cell malignancy
  "Primary effusion lymphoma",
# Primary effusion lymphoma
  "Retinoblastoma",
# Retinoblastoma
  "Neuroblastoma",
# Neuroblastoma
  "NONE"
# Default NONE value
)

# ---------------------- Organ: derive normalized label ------------------------
organ_from_fields <- function(desc, disease) {
# Define function to derive organ labels
  d <- safe_lc(desc); z <- safe_lc(disease)
# Lowercase and trim description and disease
  if (str_detect(d, "\\bbreast\\b") || str_detect(z, "\\bbreast\\b")) return("Breast")
# Assign Breast for breast-related terms
  if (is_neuro(z) || is_neuro(d)) return("Adrenal/SNS")
# Assign Adrenal/SNS for neuroblastoma
  if (is_rb(z) || is_rb(d)) return("Uveal/Eye")
# Assign Uveal/Eye for retinoblastoma
  if (is_dlbcl(z) || is_all(z) || is_hodg(z) || is_mcl(z) || is_mm(z) ||
# Check hematologic-lymphoid conditions
      is_nk(z) || is_pel(z) || is_hairy(z) ||
# Include NK, PEL, hairy cell
      str_detect(z, "primary effusion lymphoma") ||
# Include primary effusion lymphoma
      str_detect(d, "^leukemia and lymphoma panel$")) return("Hematologic-Lymphoid")
# Assign Hematologic-Lymphoid
  if (is_aml(z) || str_detect(z, "chronic myeloid leukemia.*myeloid blast") ||
# Check AML or CML myeloid blast
      is_cml(z)) return("Hematologic-Myeloid")
# Assign Hematologic-Myeloid
  if (str_detect(d, "\\blung\\b") || str_detect(z, "\\blung\\b")) return("Lung")
# Assign Lung for lung-related terms
  if (str_detect(d, "colon|rectum") || str_detect(z, "colon|rectum")) return("Colon/Rectum")
# Assign Colon/Rectum for related terms
  if (str_detect(d, "brain|cns") || str_detect(z, "brain|cns")) return("Brain/CNS")
# Assign Brain/CNS for related terms
  if (is_missing_val(d) && is_missing_val(z)) return("Unknown")
# Assign Unknown if both fields missing
  "Other"
# Default to Other
}

# -------------------- TCGA codes: derive/clean per row ------------------------
clean_codes_vec <- function(x) {
# Define function to clean TCGA codes
  if (length(x) == 0) return(character(0))
# Return empty vector if input is empty
  codes <- str_split(x, ";", simplify = FALSE)[[1]]
# Split codes by semicolon
  codes <- str_trim(codes)
# Trim whitespace from codes
  codes <- codes[codes != "" & codes %in% valid_tcga]
# Keep only valid TCGA codes
  unique(codes)
# Return unique codes
}

derive_tcga <- function(desc, disease, existing_codes = "") {
# Define function to derive TCGA codes
  d <- safe_lc(desc); z <- safe_lc(disease)
# Lowercase and trim description and disease
  if (disease %in% no_tcga_exact || is_cml(z) || is_all(z) || is_hodg(z) ||
# Check diseases with no TCGA match
      is_mcl(z) || is_mm(z) || is_nk(z) || is_pel(z) || is_hairy(z) ||
# Include MCL, MM, NK, PEL, hairy cell
      is_neuro(z) || is_rb(z)) {
# Include neuroblastoma, retinoblastoma
    return("NONE")
# Return NONE for no-match diseases
  }
  if (is_aml(z)) return("LAML")
# Assign LAML for AML
  if (is_dlbcl(z)) return("DLBC")
# Assign DLBC for DLBCL
  if (str_detect(d, "\\bbreast\\b") || str_detect(z, "\\bbreast\\b")) return("BRCA")
# Assign BRCA for breast-related terms
  kept <- clean_codes_vec(existing_codes)
# Clean existing TCGA codes
  if (length(kept) > 0) return(paste(kept, collapse = ";"))
# Return cleaned existing codes if present
  "NONE"
# Default to NONE
}

# -------------------------- Validate columns ----------------------------------
need_cols <- c("Description", "Disease")
# Define required column names
missing_cols <- setdiff(need_cols, names(meta))
# Find missing required columns
if (length(missing_cols)) stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
# Stop if required columns are missing
if (!"tcga_codes" %in% names(meta)) meta$tcga_codes <- ""
# Add empty tcga_codes column if missing

# ---------------------------- Process rows ------------------------------------
meta2 <- meta %>%
# Start dplyr pipeline with input dataframe
  mutate(
# Add new columns
    organ = pmap_chr(list(Description, Disease), organ_from_fields),
# Derive organ labels
    tcga_codes = pmap_chr(list(Description, Disease, tcga_codes), derive_tcga)
# Derive TCGA codes
  )

# ----------------------------- Write out --------------------------------------
readr::write_csv(meta2, out_csv)
# Write processed dataframe to output CSV
message("[OK] Wrote: ", out_csv)
# Log successful write to output file
message("[INFO] Organ counts:")
# Log header for organ counts
print(count(meta2, organ, sort = TRUE), n = 50)
# Print sorted organ counts, up to 50
message("[INFO] TCGA code counts (including NONE):")
# Log header for TCGA code counts
print(count(meta2, tcga_codes, sort = TRUE), n = 50)
# Print sorted TCGA code counts, up to 50
