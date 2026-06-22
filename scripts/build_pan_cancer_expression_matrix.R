#!/usr/bin/env Rscript
# =============================================================================
# build_pan_cancer_expression_matrix.R
# =============================================================================
#
# Build a feature-restricted pan-cancer expression object by combining the
# BRCA, NBL, and RBL joint VST matrices (optional HEME) into a single list
# containing expr, meta, and genes. The gene set is supplied by
# pan_cancer_features_clean.txt and the output preserves that ordering.
# Metadata embedded in the input RDS files is used when available; an
# external --metadata CSV acts only as a fallback to fill missing lineage
# or type annotations.
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

strip_version <- function(x) {
  sub("\\.[0-9]+$", "", x)
}

ensure_sample_id <- function(dt) {
  if ("sample_id" %in% names(dt)) {
    dt[, sample_id := as.character(sample_id)]
    return(dt)
  }
  for (col in c("sample", "sample_name", "Cell_Line", "DSMZ_Cell_line_norm",
                "Sample.ID", "SampleID", "Barcode", "id")) {
    if (col %in% names(dt)) {
      dt[, sample_id := as.character(get(col))]
      return(dt)
    }
  }
  stop("Metadata is missing a recognizable sample identifier column")
}

lineage_from_values <- function(values) {
  out <- rep(NA_character_, length(values))
  if (is.null(values)) {
    return(out)
  }
  out[grepl("BRCA|Breast", values, ignore.case = TRUE)] <- "BRCA"
  out[grepl("NBL|Neuroblastoma", values, ignore.case = TRUE)] <- "NBL"
  out[grepl("RBL|Retinoblastoma", values, ignore.case = TRUE)] <- "RBL"
  out[grepl("HEME|Hema|LL-100|Leukemia|Lymphoma", values, ignore.case = TRUE)] <- "HEME"
  out
}

derive_lineage <- function(dt, default_lineage) {
  lineage <- if ("lineage" %in% names(dt)) as.character(dt$lineage) else rep(NA_character_, nrow(dt))
  if ("lineage_fbk" %in% names(dt)) {
    idx <- is.na(lineage) | lineage == ""
    lineage[idx] <- as.character(dt$lineage_fbk[idx])
  }
  for (col in c("cancer_type", "cancer_type_fbk", "Disease", "Disease_fbk",
                "group2", "group2_fbk", "organ", "organ_fbk")) {
    if (col %in% names(dt)) {
      guess <- lineage_from_values(dt[[col]])
      idx <- is.na(lineage) | lineage == ""
      lineage[idx] <- guess[idx]
    }
  }
  lineage[is.na(lineage) | lineage == ""] <- default_lineage
  toupper(lineage)
}

derive_type <- function(dt) {
  type <- if ("type" %in% names(dt)) as.character(dt$type) else rep(NA_character_, nrow(dt))
  if ("type_fbk" %in% names(dt)) {
    idx <- is.na(type) | type == ""
    type[idx] <- as.character(dt$type_fbk[idx])
  }
  for (col in c("sample_type", "sample_type_fbk")) {
    if (col %in% names(dt)) {
      idx <- is.na(type) | type == ""
      type[idx] <- as.character(dt[[col]][idx])
    }
  }
  type <- tolower(type)
  type[type %in% c("cell", "cells", "cell line", "cellline", "cell_line")] <- "cell_line"
  type[type %in% c("tumour", "tumor", "tumors", "tumours")] <- "tumour"
  type
}

prepare_metadata <- function(meta_dt, fallback_meta, default_lineage, source_label, sample_ids) {
  if (is.null(meta_dt)) {
    meta_dt <- data.table(sample_id = sample_ids)
  } else {
    meta_dt <- as.data.table(meta_dt)
    meta_dt <- ensure_sample_id(meta_dt)
    meta_dt <- unique(meta_dt, by = "sample_id")
    meta_dt <- meta_dt[sample_id %in% sample_ids]
  }
  missing_ids <- setdiff(sample_ids, meta_dt$sample_id)
  if (length(missing_ids) > 0) {
    meta_dt <- rbind(meta_dt, data.table(sample_id = missing_ids), fill = TRUE)
  }
  if (!is.null(fallback_meta)) {
    fbk <- copy(fallback_meta)
    fallback_cols <- setdiff(names(fallback_meta), "sample_id")
    setnames(fbk, fallback_cols, paste0(fallback_cols, "_fbk"))
    meta_dt <- merge(meta_dt, fbk, by = "sample_id", all.x = TRUE, sort = FALSE)
  }
  meta_dt <- meta_dt[match(sample_ids, meta_dt$sample_id)]
  meta_dt[, lineage := derive_lineage(.SD, default_lineage)]
  meta_dt[, type := derive_type(.SD)]
  tumour_pattern <- grepl("^(GSE|SRR|SRP|TCGA-|TARGET-|GDC-)", meta_dt$sample_id, ignore.case = TRUE)
  dsmz_pattern <- grepl("^NG-[0-9]+", meta_dt$sample_id, ignore.case = TRUE)
  meta_dt[(is.na(type) | type == "") & tumour_pattern & !dsmz_pattern, type := "tumour"]
  meta_dt[(is.na(type) | type == ""), type := "cell_line"]
  meta_dt[, source_profile := toupper(source_label)]
  drop_cols <- grep("_fbk$", names(meta_dt), value = TRUE)
  if (length(drop_cols) > 0) {
    meta_dt[, (drop_cols) := NULL]
  }
  meta_dt[, .(
    sample_id,
    lineage,
    cancer_type = lineage,
    type,
    sample_type = type,
    source_profile
  )]
}

load_expr_with_meta <- function(path, profile_label) {
  obj <- readRDS(path)
  expr <- NULL
  meta <- NULL
  if (is.matrix(obj)) {
    expr <- obj
  } else if (is.list(obj) && "expr" %in% names(obj)) {
    expr <- obj$expr
    if (!is.null(obj$meta)) {
      meta <- as.data.table(obj$meta)
    }
  } else {
    stop(sprintf("Unsupported RDS structure for %s: %s", profile_label, path))
  }
  expr <- as.matrix(expr)
  rownames(expr) <- strip_version(rownames(expr))
  list(expr = expr, meta = meta)
}

load_fallback_metadata <- function(path) {
  if (is.null(path) || path == "" || !file.exists(path)) {
    return(NULL)
  }
  dt <- fread(path)
  dt <- ensure_sample_id(dt)
  keep_cols <- intersect(c("sample_id", "lineage", "cancer_type", "Disease",
                           "group2", "organ", "sample_type", "type"), names(dt))
  unique(dt[, ..keep_cols])
}

option_list <- list(
  make_option(c("--brca-vst"), type = "character", default = NULL,
              help = "Path to BRCA joint VST RDS"),
  make_option(c("--nbl-vst"), type = "character", default = NULL,
              help = "Path to NBL joint VST RDS"),
  make_option(c("--rbl-vst"), type = "character", default = NULL,
              help = "Path to RBL joint VST RDS"),
  make_option(c("--heme-vst"), type = "character", default = NULL,
              help = "Optional path to HEME VST RDS"),
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Fallback metadata CSV/TSV"),
  make_option(c("--genes"), type = "character", default = NULL,
              help = "Feature gene list (pan_cancer_features_clean.txt)"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Output RDS path"),
  make_option(c("--output-metadata"), type = "character", default = NULL,
              help = "Optional TSV path for metadata export"),
  make_option(c("--sample-type-filter"), type = "character", default = "all",
              help = "Sample type to retain after metadata harmonisation: all, cell_line, or tumour")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("brca-vst", "nbl-vst", "rbl-vst", "genes", "output")
missing <- required[!nzchar(unlist(opt[required]))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

if (!file.exists(opt$genes)) {
  stop("Gene list not found: ", opt$genes)
}

gene_list <- fread(opt$genes, header = FALSE)$V1
gene_list <- gene_list[gene_list != ""]
gene_list <- unique(gene_list)
genes_clean <- strip_version(gene_list)
cat("[1/4] Loaded", length(genes_clean), "pan-cancer features\n")

sample_type_filter <- tolower(gsub("-", "_", opt$`sample-type-filter`))
if (!sample_type_filter %in% c("all", "cell_line", "tumour")) {
  stop("Unsupported --sample-type-filter: ", opt$`sample-type-filter`)
}
if (sample_type_filter != "all") {
  cat("[INFO] Retaining only samples with type = ", sample_type_filter, "\n", sep = "")
}

profile_inputs <- list(
  BRCA = opt$`brca-vst`,
  NBL  = opt$`nbl-vst`,
  RBL  = opt$`rbl-vst`
)
if (!is.null(opt$`heme-vst`) && opt$`heme-vst` != "") {
  profile_inputs$HEME <- opt$`heme-vst`
}

expr_list <- list()
meta_list <- list()
fallback_meta <- load_fallback_metadata(opt$metadata)
if (!is.null(fallback_meta)) {
  cat("[INFO] Loaded fallback metadata for", nrow(fallback_meta), "samples\n")
}

cat("[2/4] Loading expression matrices...\n")
for (profile in names(profile_inputs)) {
  path <- profile_inputs[[profile]]
  if (is.null(path) || path == "") {
    next
  }
  if (!file.exists(path)) {
    stop(sprintf("Expression file for %s not found: %s", profile, path))
  }
  cat(sprintf("  Loading %s: %s\n", profile, path))
  obj <- load_expr_with_meta(path, profile)
  expr_profile <- obj$expr
  sample_ids <- colnames(expr_profile)
  meta_processed <- prepare_metadata(
    meta_dt = obj$meta,
    fallback_meta = fallback_meta,
    default_lineage = profile,
    source_label = profile,
    sample_ids = sample_ids
  )
  if (sample_type_filter != "all") {
    keep <- meta_processed$type == sample_type_filter
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) {
      cat(sprintf("    %s: 0 samples match type = %s; skipping\n",
                  profile, sample_type_filter))
      next
    }
    keep_ids <- meta_processed$sample_id[keep]
    expr_profile <- expr_profile[, keep_ids, drop = FALSE]
    meta_processed <- meta_processed[keep]
  }
  expr_list[[profile]] <- expr_profile
  meta_list[[profile]] <- meta_processed
  cat(sprintf("    %s: %d genes x %d samples\n",
              profile, nrow(expr_profile), ncol(expr_profile)))
}

if (length(expr_list) == 0) {
  stop("No expression matrices loaded")
}

gene_sets <- lapply(expr_list, rownames)
common_genes <- Reduce(intersect, gene_sets)
common_genes <- intersect(common_genes, genes_clean)
if (length(common_genes) == 0) {
  stop("No genes present across all profiles and the feature list")
}
ordered_genes <- genes_clean[genes_clean %in% common_genes]

expr_list <- lapply(expr_list, function(mat) {
  mat[ordered_genes, , drop = FALSE]
})
combined_expr <- do.call(cbind, expr_list)

meta_combined <- rbindlist(meta_list, use.names = TRUE, fill = TRUE)
meta_combined <- meta_combined[match(colnames(combined_expr), sample_id)]
if (any(is.na(meta_combined$sample_id))) {
  stop("Metadata is missing annotations for one or more samples")
}

result <- list(
  expr = combined_expr,
  meta = meta_combined,
  genes = ordered_genes
)

saveRDS(result, file = opt$output)
cat(sprintf("[OK] Saved feature-restricted expression matrix: %s\n", opt$output))
cat(sprintf("      Dimensions: %d genes x %d samples\n",
            nrow(combined_expr), ncol(combined_expr)))

if (!is.null(opt$`output-metadata`) && opt$`output-metadata` != "") {
  fwrite(meta_combined, opt$`output-metadata`, sep = "\t")
  cat(sprintf("[INFO] Wrote metadata TSV: %s\n", opt$`output-metadata`))
}

cat("[SUMMARY] Sample counts by lineage:\n")
print(table(meta_combined$lineage))
cat("[SUMMARY] Sample counts by type:\n")
print(table(meta_combined$type))
