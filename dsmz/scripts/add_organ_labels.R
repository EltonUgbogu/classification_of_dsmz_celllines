#!/usr/bin/env Rscript

# Load required packages silently
suppressPackageStartupMessages({
  library(tidyverse)  # brings in readr, dplyr, purrr, tidyr, stringr
})

# -----------------------------------------------------------------------------#
# CLI args                                      # Command line argument parsing
# -----------------------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)       # Get command line arguments
# Set input CSV path from first arg or default
in_csv  <- if (length(args) >= 1) args[[1]] else "/home/chu25/data/dsmz/DSMZ_metadata.csv"
# Set output CSV path from second arg or modify input filename  
out_csv <- if (length(args) >= 2) args[[2]] else sub("\\.csv$", "_with_organs.csv", in_csv)

message("[INFO] Input:  ", in_csv)             # Print input file path
message("[INFO] Output: ", out_csv)            # Print output file path

stopifnot(file.exists(in_csv))                 # Ensure input file exists
meta <- readr::read_csv(in_csv, show_col_types = FALSE)  # Read input CSV

# -----------------------------------------------------------------------------#
# Helpers                                       # Helper functions for processing
# -----------------------------------------------------------------------------#

# Safe lower-case + trim; convert NA to ""     # Handle missing values safely
safe_lc <- function(x) {
  x <- tidyr::replace_na(x, "")                # Replace NA with empty string
  x <- as.character(x)                         # Convert to character
  stringr::str_trim(stringr::str_to_lower(x)) # Trim whitespace and lowercase
}

# Treat these strings as "missing"             # Define what counts as missing
is_missing_val <- function(x) {
  # Check if value is in list of missing value indicators
  x %in% c("", "na", "n/a", "none", "unknown", "unspecified", "not reported", "nan")
}

# Derive high-level organ label from Description + Disease
# Returns one of: Breast, Adrenal/SNS, Uveal/Eye, Lymphoid, Myeloid, 
# Hematologic, Brain/CNS, Lung, Colon/Rectum, ... (or Other/Unknown)
organ_from_fields <- function(desc, disease) {
  d <- safe_lc(desc)                           # Clean description field
  z <- safe_lc(disease)                        # Clean disease field

  # --- Solid tumor heuristics by description or disease --------------------
  # Check for breast cancer patterns
  if (stringr::str_detect(d, "\\bbreast\\b") || stringr::str_detect(z, "\\bbreast\\b"))
    return("Breast")
  # Check for neuroblastoma (adrenal/sympathetic nervous system)
  if (stringr::str_detect(d, "neuroblastoma") || stringr::str_detect(z, "neuroblastoma"))
    return("Adrenal/SNS")
  # Check for retinoblastoma (eye cancer)
  if (stringr::str_detect(d, "retinoblastoma") || stringr::str_detect(z, "retinoblastoma"))
    return("Uveal/Eye")

  # --- Hematologic: disease-driven (overrides panels) ----------------------
  # Lymphoid patterns (including ALCL)
  # Primary effusion lymphoma
  if (stringr::str_detect(z, "primary effusion lymphoma")) return("Lymphoid")
  # Hairy cell leukemia
  if (stringr::str_detect(z, "hairy cell leukemia"))       return("Lymphoid")
  # T-cell acute lymphoblastic leukemia/lymphoma
  if (stringr::str_detect(z, "t-?acute lymphoblastic|t-?lymphoblastic|\\bt[- ]?all\\b|\\bt-?ll\\b"))
    return("Lymphoid")
  # B-cell precursor acute lymphoblastic leukemia
  if (stringr::str_detect(z, "b-?cell precursor|pre-?b-?acute lymphoblastic|\\bpre-?b-?all\\b"))
    return("Lymphoid")
  # Diffuse large B-cell lymphoma
  if (stringr::str_detect(z, "diffuse large b-?cell lymphoma|\\bdlbcl\\b|primary mediastinal b-?cell"))
    return("Lymphoid")
  # Mantle cell lymphoma
  if (stringr::str_detect(z, "mantle cell lymphoma"))      return("Lymphoid")
  # Hodgkin lymphoma
  if (stringr::str_detect(z, "hodgkin lymphoma"))          return("Lymphoid")
  # Multiple myeloma and plasma cell leukemia
  if (stringr::str_detect(z, "multiple myeloma|plasma cell leukemia"))
    return("Lymphoid")
  # Natural killer cell malignancies
  if (stringr::str_detect(z, "\\bnatural\\s*killer\\b|\\bnk\\b"))
    return("Lymphoid")
  # Chronic myeloid leukemia in lymphoid blast phase
  if (stringr::str_detect(z, "chronic myeloid leukemia.*lymphoid blast"))
    return("Lymphoid")
  # Anaplastic large cell lymphoma
  if (stringr::str_detect(z, "anaplastic large cell lymphoma|\\balcl\\b"))
    return("Lymphoid")

  # Myeloid patterns
  # Acute myeloid leukemia and related subtypes
  if (stringr::str_detect(z, "^acute myeloid leukemia|\\baml\\b|myelocytic|monocytic|megakaryocytic|erythroid"))
    return("Myeloid")
  # Chronic myeloid leukemia in myeloid blast phase
  if (stringr::str_detect(z, "chronic myeloid leukemia.*myeloid blast"))
    return("Myeloid")

  # --- Broad heme panel with no disease detail -----------------------------
  # Generic leukemia/lymphoma panel without specific disease
  if (stringr::str_detect(d, "^leukemia and lymphoma panel$") && is_missing_val(z))
    return("Hematologic")

  # --- Explicit organ/site hints (optional expansions) ---------------------
  # Lung cancer patterns
  if (stringr::str_detect(d, "\\blung\\b") || stringr::str_detect(z, "\\blung\\b"))
    return("Lung")
  # Colorectal cancer patterns
  if (stringr::str_detect(d, "\\bcolon|rectum\\b") || stringr::str_detect(z, "\\bcolon|rectum\\b"))
    return("Colon/Rectum")
  # Brain/central nervous system patterns
  if (stringr::str_detect(d, "\\bbrain|cns\\b") || stringr::str_detect(z, "\\bbrain|cns\\b"))
    return("Brain/CNS")

  # --- Missing/unknown -----------------------------------------------------
  # Both description and disease are missing/unknown
  if (is_missing_val(d) && is_missing_val(z)) return("Unknown")

  # --- Fallback ------------------------------------------------------------
  "Other"                                      # Default category
}

# -----------------------------------------------------------------------------#
# Validate required columns and compute labels # Main processing section
# -----------------------------------------------------------------------------#
need_cols <- c("Description", "Disease")       # Required column names
missing_cols <- setdiff(need_cols, names(meta)) # Find missing columns
if (length(missing_cols)) {                    # Check if any columns missing
  # Stop execution with error message listing missing columns
  stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Add organ classification column to the metadata
meta2 <- meta %>%
  mutate(
    # Apply organ classification function to each row
    organ = purrr::pmap_chr(list(Description, Disease), organ_from_fields)
  )

# -----------------------------------------------------------------------------#
# Write output & print a small summary         # Output generation and summary
# -----------------------------------------------------------------------------#
readr::write_csv(meta2, out_csv)               # Write results to CSV file
message("[OK] Wrote: ", out_csv)               # Confirm file written

# Create summary table of organ classifications
summary_tbl <- meta2 %>% count(organ, sort = TRUE)
message("[INFO] Organ label counts:")          # Print summary header
print(summary_tbl, n = nrow(summary_tbl))      # Print all rows of summary
