#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--significant-terms", type = "character"),
  make_option("--query-manifest", type = "character"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (nm in c("significant-terms", "query-manifest", "output")) {
  if (is.null(opt[[nm]]) || !nzchar(opt[[nm]])) stop("--", nm, " is required", call. = FALSE)
}

terms <- fread(normalizePath(opt$`significant-terms`, mustWork = TRUE), sep = "\t")
manifest <- fread(normalizePath(opt$`query-manifest`, mustWork = TRUE), sep = "\t")
required_manifest <- c("query_id", "query_family", "query_execution_status", "cancer_type", "marker_evidence_stratum")
missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest)) stop("Manifest missing columns: ", paste(missing_manifest, collapse = ","), call. = FALSE)
required_terms <- c("query_id", "analysis_mode", "source", "term_id", "term_name", "significant")
missing_terms <- setdiff(required_terms, names(terms))
if (length(missing_terms)) stop("Significant term table missing columns: ", paste(missing_terms, collapse = ","), call. = FALSE)

contrast_queries <- manifest[
  query_family == "contrast_level_marker_set" &
    query_execution_status == "runnable",
  .(query_id, cancer_type, marker_evidence_stratum)
]
if (!nrow(contrast_queries)) stop("No runnable contrast-level marker-set queries in manifest", call. = FALSE)

denominators <- contrast_queries[
  ,
  .(runnable_contrast_count = uniqueN(query_id)),
  by = .(cancer_type, marker_evidence_stratum)
]
sig <- merge(
  terms[analysis_mode == "primary" & significant == TRUE],
  contrast_queries,
  by = "query_id",
  all = FALSE
)
support <- sig[
  ,
  .(
    significant_contrast_count = uniqueN(query_id),
    term_name = term_name[1L]
  ),
  by = .(cancer_type, marker_evidence_stratum, source, term_id)
]
support <- merge(support, denominators, by = c("cancer_type", "marker_evidence_stratum"), all.x = TRUE)
support[is.na(significant_contrast_count), significant_contrast_count := 0L]
support[, term_support_fraction := significant_contrast_count / runnable_contrast_count]
setorder(support, cancer_type, marker_evidence_stratum, -term_support_fraction, source, term_id)

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
fwrite(support, opt$output, sep = "\t", na = "NA")
fwrite(
  data.table(
    artifact_name = basename(opt$output),
    artifact_type = "enrichment_contrast_term_support_table",
    numerator = "significant_contrast_count",
    denominator = "runnable_contrast_count",
    support_definition = "term_support_fraction = significant_contrast_count / runnable_contrast_count within cancer_type x marker_evidence_stratum"
  ),
  sub("\\.tsv$", "_provenance.tsv", opt$output),
  sep = "\t"
)
message("[Contrast interpretation] Contrast-level term support fractions written.")
