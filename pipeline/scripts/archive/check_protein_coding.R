library(rtracklayer)

gtf_path  <- "/work/ugbogu/pipeline/data/gtf/gencode.v44.basic.annotation.gtf"
# NBL gene list (ENSG gene IDs). For BRCA use e.g. .../brca/deseq2_markers/markers/brca_features.txt
#list_path <- "/work/ugbogu/pipeline/results/unsupervised/nbl/deseq2_markers/markers/unique_feature_set_recurrence_ge_2.txt"
list_path <- "/work/ugbogu/pipeline/results/unsupervised/nbl/deseq2_markers/markers/nbl_features.txt"
gtf <- import(gtf_path)
genes <- gtf[gtf$type == "gene"]

# DIAGNOSTIC: Check what gene_id actually looks like in GTF
cat("\n=== GTF gene_id DIAGNOSTICS ===\n")
gene_ids_gtf_raw <- mcols(genes)$gene_id
cat("First 20 GTF gene_id values:\n")
print(head(gene_ids_gtf_raw, 20))
cat("\nHow many GTF gene_id start with ENSG?\n")
print(table(grepl("^ENSG", as.character(gene_ids_gtf_raw))))
cat("\nHow many GTF gene_id contain a dot (version suffix)?\n")
print(table(grepl("\\.", as.character(gene_ids_gtf_raw))))
cat("\nGTF gene_id class:\n")
print(class(gene_ids_gtf_raw))
cat("\nGTF gene_id sample (first 5):\n")
print(head(as.character(gene_ids_gtf_raw), 5))

# See what columns are actually available
cols <- colnames(mcols(genes))
message("\nAvailable gene mcols: ", paste(cols, collapse = ", "))

# Pick the biotype column name automatically
biotype_col <- if ("gene_type" %in% cols) "gene_type" else if ("gene_biotype" %in% cols) "gene_biotype" else NA
if (is.na(biotype_col)) stop("Neither gene_type nor gene_biotype found in GTF attributes.")

# gene_name sometimes missing; fall back gracefully
name_col <- if ("gene_name" %in% cols) "gene_name" else NA

keep <- c("gene_id", name_col, biotype_col)
keep <- keep[!is.na(keep)]

gene_info <- as.data.frame(mcols(genes)[, keep, drop = FALSE])
colnames(gene_info)[colnames(gene_info) == biotype_col] <- "gene_biotype"
if (!("gene_name" %in% colnames(gene_info))) gene_info$gene_name <- NA_character_

# Read IDs
ids <- readLines(list_path, warn = FALSE)
ids <- unique(trimws(ids))
ids <- ids[ids != ""]
ids <- sub("\\..*$", "", ids)

# DIAGNOSTIC: Check what input IDs look like
cat("\n=== INPUT IDs DIAGNOSTICS ===\n")
cat("First 10 input IDs:\n")
print(head(ids, 10))
cat("\nHow many input IDs start with ENSG?\n")
print(table(grepl("^ENSG", ids)))
cat("\nHow many input IDs start with ENST?\n")
print(table(grepl("^ENST", ids)))
cat("\nInput ID lengths:\n")
print(table(nchar(ids)))
cat("\nFirst input ID raw bytes:\n")
print(charToRaw(ids[1]))

# Normalize both sides for robust matching
cat("\n=== NORMALIZED MATCHING ===\n")
gene_ids_gtf <- as.character(gene_info$gene_id)
gene_ids_gtf <- sub("\\..*$", "", gene_ids_gtf)   # strip versions if present
gene_ids_gtf <- trimws(gene_ids_gtf)

ids_normalized <- trimws(ids)

cat("Example GTF gene_id (normalized):\n")
print(head(gene_ids_gtf, 10))
cat("\nExample list ids (normalized):\n")
print(head(ids_normalized, 10))

# Try matching with normalized versions
matched_normalized <- intersect(ids_normalized, gene_ids_gtf)
missing_normalized <- setdiff(ids_normalized, gene_ids_gtf)

cat("\nMatched genes (normalized):", length(matched_normalized), "of", length(ids_normalized), "\n")
if (length(matched_normalized) > 0) {
  cat("First 10 matched IDs:\n")
  print(head(matched_normalized, 10))
}
if (length(missing_normalized) > 0 && length(missing_normalized) <= 20) {
  cat("\nMissing genes (not found in GTF gene_id):\n")
  print(missing_normalized)
} else if (length(missing_normalized) > 20) {
  cat("\nMissing genes (first 20 of", length(missing_normalized), "total):\n")
  print(head(missing_normalized, 20))
}

# Use normalized matching only (original gene_id has version suffixes and would give empty table)
tab_normalized <- table(character(0))  # default empty
if (length(matched_normalized) > 0) {
  gene_info_normalized <- gene_info
  gene_info_normalized$gene_id_normalized <- gene_ids_gtf
  matched_biotypes <- gene_info_normalized$gene_biotype[gene_info_normalized$gene_id_normalized %in% ids_normalized]
  tab_normalized <- table(matched_biotypes, useNA = "ifany")
  cat("\n=== BIOTYPE TABLE (normalized matching) ===\n")
  print(tab_normalized)
}

# Clean summary for reporting
pc <- if ("protein_coding" %in% names(tab_normalized)) tab_normalized[["protein_coding"]] else 0
total <- sum(tab_normalized)

cat("\n=== SUMMARY (for reporting) ===\n")
cat(sprintf("Input genes (unique): %d\n", length(ids_normalized)))
cat(sprintf("Matched in GTF: %d\n", length(matched_normalized)))
if (total > 0) {
  cat(sprintf("Protein-coding: %d (%.1f%% of matched)\n", pc, 100 * pc / total))
  cat(sprintf("Non–protein-coding: %d (%.1f%% of matched)\n", total - pc, 100 * (total - pc) / total))
} else {
  cat("No matched genes; cannot compute biotype percentages.\n")
}
