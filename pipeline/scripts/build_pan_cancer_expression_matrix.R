#!/usr/bin/env Rscript
# =============================================================================
# build_pan_cancer_expression_matrix.R
# =============================================================================
#
# Build combined pan-cancer expression matrix restricted to pan-cancer features
# Includes: BRCA, NBL, RBL cell lines + hematologic (HEME) as negative control
#
# USAGE:
# Rscript build_pan_cancer_expression_matrix.R \
#   --brca-vst results/brca/inputs/vst_joint.rds \
#   --nbl-vst results/nbl/inputs/vst_joint.rds \
#   --rbl-vst results/rbl/inputs/vst_joint.rds \
#   --heme-vst results/heme/inputs/vst_joint.rds \
#   --metadata data/dsmz/DSMZ_metadata.csv \
#   --genes results/unsupervised/pan_cancer/pan_cancer_features_clean.txt \
#   --output results/unsupervised/pan_cancer/pan_cancer_expr.rds
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--brca-vst"), type="character", default=NULL,
              help="Path to BRCA VST expression matrix (RDS)"),
  make_option(c("--nbl-vst"), type="character", default=NULL,
              help="Path to NBL VST expression matrix (RDS)"),
  make_option(c("--rbl-vst"), type="character", default=NULL,
              help="Path to RBL VST expression matrix (RDS)"),
  make_option(c("--heme-vst"), type="character", default=NULL,
              help="Path to HEME VST expression matrix (RDS)"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Path to DSMZ metadata CSV file"),
  make_option(c("--genes"), type="character", default=NULL,
              help="Path to pan-cancer feature gene list (one per line)"),
  make_option(c("--output"), type="character", default=NULL,
              help="Output path for combined expression matrix (RDS)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
required <- c("brca-vst", "nbl-vst", "rbl-vst", "metadata", "genes", "output")
missing <- setdiff(required, names(opt)[!sapply(opt, is.null)])
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse=", "))
}

cat("[1/5] Loading pan-cancer feature genes...\n")
genes <- fread(opt$genes, header=FALSE)$V1
cat("  Loaded", length(genes), "genes\n")

cat("[2/5] Loading expression matrices...\n")
expr_list <- list()

# Load BRCA
if (!is.null(opt$`brca-vst`) && file.exists(opt$`brca-vst`)) {
  cat("  Loading BRCA:", opt$`brca-vst`, "\n")
  brca <- readRDS(opt$`brca-vst`)
  if (is.matrix(brca)) {
    expr_list[["BRCA"]] <- brca
  } else if (is.list(brca) && "expr" %in% names(brca)) {
    expr_list[["BRCA"]] <- brca$expr
  } else {
    stop("BRCA VST file has unexpected structure")
  }
  cat("    BRCA: ", nrow(expr_list[["BRCA"]]), " genes, ", 
      ncol(expr_list[["BRCA"]]), " samples\n", sep="")
}

# Load NBL
if (!is.null(opt$`nbl-vst`) && file.exists(opt$`nbl-vst`)) {
  cat("  Loading NBL:", opt$`nbl-vst`, "\n")
  nbl <- readRDS(opt$`nbl-vst`)
  if (is.matrix(nbl)) {
    expr_list[["NBL"]] <- nbl
  } else if (is.list(nbl) && "expr" %in% names(nbl)) {
    expr_list[["NBL"]] <- nbl$expr
  } else {
    stop("NBL VST file has unexpected structure")
  }
  cat("    NBL: ", nrow(expr_list[["NBL"]]), " genes, ", 
      ncol(expr_list[["NBL"]]), " samples\n", sep="")
}

# Load RBL
if (!is.null(opt$`rbl-vst`) && file.exists(opt$`rbl-vst`)) {
  cat("  Loading RBL:", opt$`rbl-vst`, "\n")
  rbl <- readRDS(opt$`rbl-vst`)
  if (is.matrix(rbl)) {
    expr_list[["RBL"]] <- rbl
  } else if (is.list(rbl) && "expr" %in% names(rbl)) {
    expr_list[["RBL"]] <- rbl$expr
  } else {
    stop("RBL VST file has unexpected structure")
  }
  cat("    RBL: ", nrow(expr_list[["RBL"]]), " genes, ", 
      ncol(expr_list[["RBL"]]), " samples\n", sep="")
}

# Load HEME (optional)
if (!is.null(opt$`heme-vst`) && file.exists(opt$`heme-vst`)) {
  cat("  Loading HEME:", opt$`heme-vst`, "\n")
  heme <- readRDS(opt$`heme-vst`)
  if (is.matrix(heme)) {
    expr_list[["HEME"]] <- heme
  } else if (is.list(heme) && "expr" %in% names(heme)) {
    expr_list[["HEME"]] <- heme$expr
  } else {
    stop("HEME VST file has unexpected structure")
  }
  cat("    HEME: ", nrow(expr_list[["HEME"]]), " genes, ", 
      ncol(expr_list[["HEME"]]), " samples\n", sep="")
}

if (length(expr_list) == 0) {
  stop("No expression matrices loaded")
}

cat("[3/5] Filtering to pan-cancer features and merging...\n")
# Strip version numbers from gene IDs if present (ENSG00000000003.15 -> ENSG00000000003)
# VST matrices typically don't have version numbers
genes_no_version <- sub("\\.[0-9]+$", "", genes)

# Filter each matrix to pan-cancer features
expr_filtered <- lapply(expr_list, function(expr) {
  # Try matching with and without version numbers
  vst_genes <- rownames(expr)
  common_genes <- intersect(vst_genes, genes_no_version)
  
  if (length(common_genes) == 0) {
    # Try with version numbers if no match
    common_genes <- intersect(vst_genes, genes)
  }
  
  if (length(common_genes) == 0) {
    stop("No common genes found between expression matrix and feature list")
  }
  expr[common_genes, , drop=FALSE]
})

# Merge matrices (cbind, handling missing genes)
all_genes <- unique(unlist(lapply(expr_filtered, rownames)))
all_genes <- intersect(all_genes, genes_no_version)  # Only keep genes in feature list (no version)

cat("  Merging", length(all_genes), "genes across", length(expr_filtered), "diseases\n")

# Build combined matrix
combined_expr <- matrix(NA, nrow=length(all_genes), ncol=0)
rownames(combined_expr) <- all_genes

for (disease in names(expr_filtered)) {
  expr_sub <- expr_filtered[[disease]]
  common <- intersect(rownames(combined_expr), rownames(expr_sub))
  combined_expr <- cbind(combined_expr, expr_sub[common, , drop=FALSE])
}

# Fill missing values with 0 (shouldn't happen if all matrices have same genes)
combined_expr[is.na(combined_expr)] <- 0

cat("  Combined matrix: ", nrow(combined_expr), " genes, ", 
    ncol(combined_expr), " samples\n", sep="")

cat("[4/5] Loading and aligning metadata...\n")
meta <- fread(opt$metadata)

# Determine lineage from available metadata columns
if ("lineage" %in% colnames(meta)) {
  meta[, lineage := as.character(lineage)]
} else if ("cancer_type" %in% colnames(meta)) {
  meta[, lineage := fcase(
    cancer_type %like% "BRCA|Breast", "BRCA",
    cancer_type %like% "NBL|Neuroblastoma", "NBL",
    cancer_type %like% "RBL|Retinoblastoma", "RBL",
    cancer_type %like% "HEME|Hematologic", "HEME",
    default = NA_character_
  )]
} else if ("Disease" %in% colnames(meta)) {
  meta[, lineage := fcase(
    Disease %like% "Breast", "BRCA",
    Disease %like% "Neuroblastoma", "NBL",
    Disease %like% "Retinoblastoma", "RBL",
    Disease %like% "leukemia|lymphoma|Hematologic", "HEME",
    default = NA_character_
  )]
  if ("group2" %in% colnames(meta)) {
    meta[is.na(lineage), lineage := fcase(
      group2 %like% "Breast", "BRCA",
      group2 %like% "Neuroblastoma", "NBL",
      group2 %like% "Retinoblastoma", "RBL",
      group2 %like% "LL-100|leukemia|lymphoma|Hematologic", "HEME",
      default = NA_character_
    )]
  }
} else if ("group2" %in% colnames(meta)) {
  meta[, lineage := fcase(
    group2 %like% "Breast", "BRCA",
    group2 %like% "Neuroblastoma", "NBL",
    group2 %like% "Retinoblastoma", "RBL",
    group2 %like% "LL-100|leukemia|lymphoma|Hematologic", "HEME",
    default = NA_character_
  )]
}

# Determine sample_id column
sample_id_col <- NULL
for (col in c("sample_id", "sample", "sample_name", "Cell_Line", "DSMZ_Cell_line_norm")) {
  if (col %in% colnames(meta)) {
    sample_id_col <- col
    break
  }
}
if (is.null(sample_id_col)) {
  stop("Could not find sample ID column in metadata")
}

meta[, sample_id := get(sample_id_col)]
if ("sample_type" %in% colnames(meta)) {
  meta[, type := fcase(
    sample_type %like% "Cell Line|cell_line|cell line", "cell_line",
    sample_type %like% "Tumour|Tumor|tumour|tumor", "tumour",
    default = "cell_line"
  )]
} else {
  meta[, type := "cell_line"]  # Default for DSMZ-only metadata
}

# Build expression-derived metadata for all samples (including tumours)
expr_meta <- rbindlist(lapply(names(expr_list), function(disease) {
  data.table(sample_id = colnames(expr_list[[disease]]), lineage_expr = disease)
}), use.names = TRUE)

meta <- meta[, .(sample_id, lineage_meta = lineage, type_meta = type)]
meta <- unique(meta)

meta_all <- merge(expr_meta, meta, by = "sample_id", all.x = TRUE)
meta_all[, lineage := ifelse(!is.na(lineage_meta), lineage_meta, lineage_expr)]
meta_all[, type := ifelse(!is.na(type_meta), type_meta, "tumour")]

# Override type based on sample ID patterns: force tumour classification for
# tumour-like IDs (GSE*, SRR*, SRP*, TCGA*, TARGET*, GDC*) that don't match DSMZ cell line patterns
# This override applies even if metadata says "cell_line" to ensure consistency
is_tumour_pattern <- grepl("^GSE[0-9]+|^SRR[0-9]+|^SRP[0-9]+|^TCGA-|^TARGET-|^GDC-", meta_all$sample_id)
is_dsmz_cellline <- grepl("^NG-[0-9]+_.*_lib[0-9]+_|^NG-[0-9]+_.*_libLAB", meta_all$sample_id)
meta_all[is_tumour_pattern & !is_dsmz_cellline, type := "tumour"]

meta_all <- meta_all[, .(sample_id, lineage, type)]
meta_all <- unique(meta_all)

# Align with expression matrix and preserve column order
meta_all <- meta_all[match(colnames(combined_expr), sample_id)]
if (any(is.na(meta_all$sample_id))) {
  stop("Missing metadata for one or more expression samples")
}

combined_expr <- combined_expr[, meta_all$sample_id, drop = FALSE]
meta <- meta_all

cat("  Aligned metadata: ", nrow(meta), " samples\n", sep="")
cat("  Lineage distribution:\n")
print(table(meta$lineage, useNA = "ifany"))
cat("  Type distribution:\n")
print(table(meta$type, useNA = "ifany"))

cat("[5/5] Saving combined expression matrix...\n")
result <- list(
  expr = combined_expr,
  meta = meta,
  genes = all_genes
)

saveRDS(result, file = opt$output)
cat("[OK] Combined pan-cancer expression matrix saved to:", opt$output, "\n")
cat("  Final dimensions: ", nrow(combined_expr), " genes x ", 
    ncol(combined_expr), " samples\n", sep="")
