#!/usr/bin/env Rscript
#
# make_cell_cycle_ensg_from_gmt.R
#
# Extract Ensembl gene IDs associated with cell cycle GO terms from a GMT file.
# GMT (Gene Matrix Transposed) format: each line contains a gene set with ID, name, and genes.
#

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

# Input: GMT file containing GO:BP gene sets with Ensembl IDs
# Default path points to a g:Profiler GMT file for human GO Biological Process terms
gmt_file <- if (length(args) >= 1) args[1] else file.path(
  repo_root, "results", "unsupervised", "nbl", "gsea_ova",
  "gprofiler_hsapiens", "hsapiens.GO:BP.ENSG.gmt"
)

# Output: plain text file with one Ensembl gene ID per line
out_file <- if (length(args) >= 2) args[2] else "cell_cycle_ensg.txt"

# GO term IDs to extract genes from
# GO:0000278 = mitotic cell cycle
go_ids <- c("GO:0000278")

# Validate input file exists
if (!file.exists(gmt_file)) stop("File not found: ", gmt_file)

# Open file connection (handle both plain text and gzipped files)
con <- if (grepl("\\.gz$", gmt_file)) gzfile(gmt_file, "rt") else file(gmt_file, "rt")

# Read all lines from GMT file
lines <- readLines(con, warn = FALSE)
close(con)

# Find lines matching any of the specified GO IDs
# GMT format: GO_ID<TAB>GO_NAME<TAB>GENE1<TAB>GENE2<TAB>...
# Initialize boolean vector to track which lines to keep
keep <- rep(FALSE, length(lines))

# For each GO ID, mark lines that start with that ID (followed by tab)
for (id in go_ids) keep <- keep | grepl(paste0("^", id, "\t"), lines)

# Extract the matching lines
hits <- lines[keep]

# Error if no matching GO terms found
if (length(hits) == 0) stop("No matching GO IDs found in GMT: ", paste(go_ids, collapse = ", "))

# Extract gene IDs from matching GMT lines
genes <- unique(unlist(lapply(strsplit(hits, "\t", fixed = TRUE), function(x) {
  # GMT format: column 1 = GO ID, column 2 = GO name, columns 3+ = genes
  if (length(x) < 3) return(character(0))  # Skip malformed lines
  x[seq.int(3, length(x))]  # Extract all columns from position 3 onward (the genes)
})))

# Clean up gene list
genes <- genes[nzchar(genes)]  # Remove empty strings

# Filter to keep only valid Ensembl gene IDs
# Pattern matches: ENSG followed by digits, optionally with version (e.g., ENSG00000012345.2)
genes <- genes[grepl("^ENSG[0-9]+(\\.[0-9]+)?$", genes)]

# Remove version suffixes (e.g., ENSG00000012345.2 → ENSG00000012345)
genes <- sub("\\.[0-9]+$", "", genes)

# Get unique gene IDs and sort alphabetically
genes <- sort(unique(genes))

# Write output: one gene ID per line
writeLines(genes, con = out_file)

# Print summary statistics
cat("GMT:", gmt_file, "\n")
cat("GO IDs:", paste(go_ids, collapse = ", "), "\n")
cat("Wrote:", out_file, "\n")
cat("Genes:", length(genes), "\n")
