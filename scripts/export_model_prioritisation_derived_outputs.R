#!/usr/bin/env Rscript
# Export derived model-prioritisation outputs for inspection and reproducibility.

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    pan_cancer_dir = "results/unsupervised/pan_cancer",
    output_dir = "derived_results/model_prioritisation",
    feature_table = NA_character_,
    gene_annotation = "resources/ensembl_to_symbol.tsv"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be provided as --key value pairs.", call. = FALSE)
    }
    value <- args[[i + 1L]]
    key <- sub("^--", "", key)
    if (!key %in% names(out)) {
      stop("Unknown argument: --", key, call. = FALSE)
    }
    out[[key]] <- value
    i <- i + 2L
  }
  out
}

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
  path
}

lineage_full <- function(x) {
  map <- c(
    BRCA = "breast cancer",
    NBL = "neuroblastoma",
    RBL = "retinoblastoma"
  )
  y <- unname(map[as.character(x)])
  y[is.na(y)] <- as.character(x)[is.na(y)]
  y
}

infer_cohort <- function(sample_id) {
  fifelse(
    grepl("^TCGA-", sample_id), "TCGA-BRCA",
    fifelse(
      grepl("^TARGET-", sample_id), "TARGET-NBL",
      fifelse(grepl("_", sample_id), sub("_.*$", "", sample_id), "unknown")
    )
  )
}

make_display_labels <- function(sample_ids) {
  labels <- gsub("^NG-[^_]+_", "", sample_ids)
  labels <- gsub("_lib[^_]+.*$", "", labels)
  labels <- gsub("_", "-", labels)

  dup_mask <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  if (any(dup_mask)) {
    bio_suffix <- sub(".*_(\\d+)$", "\\1", sample_ids)
    labels[dup_mask] <- paste0(labels[dup_mask], "-", bio_suffix[dup_mask])
  }
  labels
}

clean_latex_cell <- function(x) {
  x <- gsub("\\\\allowbreak\\{\\}", "", x)
  x <- gsub("\\\\_", "_", x)
  x <- gsub("\\\\&", "&", x)
  trimws(x)
}

read_feature_table <- function(path) {
  if (grepl("\\.tsv$", path)) {
    dt <- fread(path)
    required <- c(
      "rank", "ensembl_gene_id", "source_category",
      "recurrent_disease", "isolate_source_contrasts"
    )
    missing <- setdiff(required, names(dt))
    if (length(missing)) {
      stop("Feature table missing columns: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    return(dt[, ..required])
  }

  lines <- readLines(path, warn = FALSE)
  start <- grep("\\\\label\\{tab:pancancer_257_genes\\}", lines)
  if (!length(start)) {
    stop("Could not find tab:pancancer_257_genes in ", path, call. = FALSE)
  }
  end <- grep("\\\\end\\{longtable\\}", lines)
  end <- end[end > start[1L]][1L]
  if (is.na(end)) {
    stop("Could not find end of pan-cancer feature longtable.", call. = FALSE)
  }
  block <- lines[start[1L]:end]
  rows <- grep("^\\s*[0-9]+\\s*&", block, value = TRUE)
  rows <- gsub("\\\\\\\\.*$", "", rows)
  parsed <- rbindlist(lapply(rows, function(row) {
    fields <- strsplit(row, "\\s*&\\s*", perl = TRUE)[[1L]]
    if (length(fields) != 5L) {
      stop("Unexpected feature table row: ", row, call. = FALSE)
    }
    data.table(
      rank = as.integer(clean_latex_cell(fields[1L])),
      ensembl_gene_id = clean_latex_cell(fields[2L]),
      source_category = clean_latex_cell(fields[3L]),
      recurrent_disease = clean_latex_cell(fields[4L]),
      isolate_source_contrasts = clean_latex_cell(fields[5L])
    )
  }))
  parsed[]
}

write_tsv <- function(dt, path) {
  fwrite(dt, path, sep = "\t", quote = FALSE, na = "")
}

validate_outputs <- function(
  patient_rankings,
  ecdf_input,
  ecdf_summary,
  retrieval,
  group_metadata,
  feature_set
) {
  model_counts <- patient_rankings[, .N, by = tumour_sample_id]
  if (!all(model_counts$N == 58L)) {
    stop("patient_to_cellline_rankings.tsv does not contain 58 models per tumour.",
         call. = FALSE)
  }
  if (!all(patient_rankings[, any(rank == 1), by = tumour_sample_id]$V1)) {
    stop("At least one tumour sample lacks rank 1.", call. = FALSE)
  }
  if (!all(patient_rankings$is_focal_lineage !=
           patient_rankings$is_non_focal_lineage)) {
    stop("Focal and non-focal indicators are not logical opposites.",
         call. = FALSE)
  }
  if (!all(ecdf_input$top10_indicator == (ecdf_input$rank <= 10))) {
    stop("ECDF top10_indicator does not match rank <= 10.", call. = FALSE)
  }
  if (!all(ecdf_summary$non_focal_median_ecdf >= 0 &
           ecdf_summary$non_focal_median_ecdf <= 1 &
           ecdf_summary$non_focal_q25_ecdf >= 0 &
           ecdf_summary$non_focal_q25_ecdf <= 1 &
           ecdf_summary$non_focal_q75_ecdf >= 0 &
           ecdf_summary$non_focal_q75_ecdf <= 1)) {
    stop("Non-focal ECDF summaries must be between 0 and 1.", call. = FALSE)
  }
  if (!all(ecdf_summary$non_focal_q25_ecdf <=
           ecdf_summary$non_focal_median_ecdf &
           ecdf_summary$non_focal_median_ecdf <=
           ecdf_summary$non_focal_q75_ecdf)) {
    stop("Non-focal ECDF quantiles are not ordered.", call. = FALSE)
  }
  if (uniqueN(retrieval$cell_line_group_id) != 56L) {
    stop("Cell line-centred retrieval does not contain 56 biological groups.",
         call. = FALSE)
  }
  top1 <- unique(retrieval[, .(
    cell_line_group_id,
    cell_line_lineage,
    top1_tumour_lineage,
    top1_correct_indicator
  )])
  if (sum(top1$top1_correct_indicator) != 54L) {
    stop("Cell line-centred top-1 accuracy does not reproduce 54/56.",
         call. = FALSE)
  }
  if (top1[cell_line_lineage == "retinoblastoma", uniqueN(cell_line_group_id)] !=
      9L) {
    stop("Retinoblastoma retrieval does not contain 9 biological groups.",
         call. = FALSE)
  }
  if (nrow(group_metadata) != 58L ||
      uniqueN(group_metadata$cell_line_group_id) != 56L) {
    stop("Cell line metadata must expose 58 profiles and 56 groups.",
         call. = FALSE)
  }
  if (uniqueN(feature_set$gene_id) != 257L || nrow(feature_set) != 257L) {
    stop("Feature set export must contain 257 unique genes.", call. = FALSE)
  }
}

args <- parse_args()

pan_dir <- normalizePath(args$pan_cancer_dir, mustWork = TRUE)
out_dir <- args$output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_rds <- require_file(file.path(
  pan_dir,
  "tumour_mapping/tumour_to_cellline_similarity/tumour_cellline_scores.rds"
))
meta_rds <- require_file(file.path(pan_dir, "pan_cancer_expr.rds"))
mapping_tsv <- require_file(file.path(
  pan_dir,
  "tumour_mapping/cellline_to_tumour_similarity/publication_figures/replicate_collapse_mapping.tsv"
))
components_tsv <- require_file(file.path(pan_dir, "tumour_components.tsv"))

feature_table <- args$feature_table
if (is.na(feature_table) || !nzchar(feature_table)) {
  candidate <- file.path(dirname(pan_dir), "Table_S5_pan_cancer_feature_set.tsv")
  feature_table <- require_file(candidate)
} else {
  feature_table <- require_file(feature_table)
}

scores <- readRDS(score_rds)
pan_obj <- readRDS(meta_rds)
meta <- as.data.table(pan_obj$meta)
mapping <- fread(mapping_tsv)

setnames(mapping, c("cell_line_raw", "cell_line_group", "cell_lineage"),
         c("cell_line_model_id", "cell_line_group_id", "cell_line_lineage_code"))
mapping[, cell_line_lineage := lineage_full(cell_line_lineage_code)]

meta_tum <- meta[type == "tumour"]
meta_cell <- meta[type == "cell_line"]

if (!all(rownames(scores) %in% meta_tum$sample_id)) {
  stop("Score matrix row names do not match tumour metadata.", call. = FALSE)
}
if (!all(colnames(scores) %in% meta_cell$sample_id)) {
  stop("Score matrix column names do not match cell line metadata.", call. = FALSE)
}
if (!all(colnames(scores) %in% mapping$cell_line_model_id)) {
  stop("Not all score-matrix cell line profiles are present in replicate mapping.",
       call. = FALSE)
}

score_dt <- as.data.table(as.table(scores))
setnames(score_dt, c("tumour_sample_id", "cell_line_model_id", "similarity_score"))
score_dt[, tumour_sample_id := as.character(tumour_sample_id)]
score_dt[, cell_line_model_id := as.character(cell_line_model_id)]
score_dt[, similarity_score := as.numeric(similarity_score)]

tumour_lineage_code <- setNames(meta_tum$lineage, meta_tum$sample_id)
cell_lineage_code <- setNames(meta_cell$lineage, meta_cell$sample_id)
score_dt[, tumour_lineage_code := tumour_lineage_code[tumour_sample_id]]
score_dt[, cell_line_lineage_code := cell_lineage_code[cell_line_model_id]]
score_dt <- merge(
  score_dt,
  mapping[, .(cell_line_model_id, cell_line_group_id)],
  by = "cell_line_model_id",
  all.x = TRUE,
  sort = FALSE
)

score_dt[, tumour_lineage := lineage_full(tumour_lineage_code)]
score_dt[, cell_line_lineage := lineage_full(cell_line_lineage_code)]
score_dt[, tumour_project_or_cohort := infer_cohort(tumour_sample_id)]
score_dt[, similarity_metric := "Spearman correlation"]
score_dt[, rank := frank(-similarity_score, ties.method = "average"),
         by = tumour_sample_id]
score_dt[, is_focal_lineage := cell_line_lineage == tumour_lineage]
score_dt[, is_non_focal_lineage := cell_line_lineage != tumour_lineage]

patient_rankings <- score_dt[, .(
  tumour_sample_id,
  tumour_lineage,
  tumour_project_or_cohort,
  cell_line_model_id,
  cell_line_group_id,
  cell_line_lineage,
  similarity_metric,
  similarity_score,
  rank,
  is_focal_lineage,
  is_non_focal_lineage
)]
setorder(patient_rankings, tumour_sample_id, rank, cell_line_model_id)

display_labels <- make_display_labels(meta_cell$sample_id)
names(display_labels) <- meta_cell$sample_id
displayed <- patient_rankings[is_focal_lineage == TRUE,
  .(
    median_rank = median(rank, na.rm = TRUE),
    median_score = median(similarity_score, na.rm = TRUE)
  ),
  by = .(tumour_lineage, cell_line_model_id)
][order(tumour_lineage, median_rank, -median_score, cell_line_model_id)]
displayed <- displayed[, head(.SD, 4L), by = tumour_lineage]

patient_rankings[, figure_label := unname(display_labels[cell_line_model_id])]
patient_rankings[, panel := fifelse(
  tumour_lineage == "breast cancer", "A",
  fifelse(tumour_lineage == "neuroblastoma", "B",
          fifelse(tumour_lineage == "retinoblastoma", "C", NA_character_))
)]
patient_rankings[, is_displayed_model :=
  paste(tumour_lineage, cell_line_model_id) %in%
  paste(displayed$tumour_lineage, displayed$cell_line_model_id)]

ecdf_input <- patient_rankings[, .(
  tumour_lineage,
  tumour_sample_id,
  cell_line_model_id,
  cell_line_group_id,
  cell_line_lineage,
  rank,
  top10_indicator = rank <= 10,
  is_displayed_model,
  is_focal_lineage,
  is_non_focal_lineage,
  similarity_score,
  panel,
  figure_label
)]
setorder(ecdf_input, tumour_lineage, tumour_sample_id, rank, cell_line_model_id)

thresholds <- data.table(rank_threshold = seq_len(ncol(scores)))
ecdf_by_model <- patient_rankings[is_non_focal_lineage == TRUE,
  {
    ranks <- rank
    data.table(
      rank_threshold = thresholds$rank_threshold,
      ecdf = vapply(thresholds$rank_threshold, function(k) {
        mean(ranks <= k)
      }, numeric(1L))
    )
  },
  by = .(tumour_lineage, cell_line_model_id)
]
ecdf_summary <- ecdf_by_model[, .(
  non_focal_median_ecdf = median(ecdf),
  non_focal_q25_ecdf = as.numeric(quantile(ecdf, 0.25, names = FALSE)),
  non_focal_q75_ecdf = as.numeric(quantile(ecdf, 0.75, names = FALSE)),
  n_non_focal_models = uniqueN(cell_line_model_id)
), by = .(tumour_lineage, rank_threshold)]
setorder(ecdf_summary, tumour_lineage, rank_threshold)

group_scores <- merge(
  score_dt[, .(
    tumour_sample_id,
    tumour_lineage,
    cell_line_model_id,
    similarity_score
  )],
  mapping[, .(cell_line_model_id, cell_line_group_id, cell_line_lineage)],
  by = "cell_line_model_id",
  all.x = TRUE,
  sort = FALSE
)[, .(
  similarity_score = median(similarity_score)
), by = .(
  cell_line_group_id,
  cell_line_lineage,
  tumour_sample_id,
  tumour_lineage
)]

setorder(group_scores, cell_line_group_id, -similarity_score, tumour_sample_id)
group_scores[, rank := seq_len(.N), by = cell_line_group_id]
group_scores[, same_lineage_indicator := cell_line_lineage == tumour_lineage]
group_scores[, similarity_metric := "Spearman correlation"]
group_scores[, rank_percentile := rank / .N, by = cell_line_group_id]

top1 <- group_scores[rank == 1L, .(
  cell_line_group_id,
  top1_tumour_lineage = tumour_lineage,
  top1_correct_indicator = same_lineage_indicator
)]
retrieval <- merge(group_scores, top1, by = "cell_line_group_id", sort = FALSE)
retrieval <- retrieval[, .(
  cell_line_group_id,
  cell_line_lineage,
  tumour_sample_id,
  tumour_lineage,
  similarity_metric,
  similarity_score,
  rank,
  same_lineage_indicator,
  rank_percentile,
  top1_tumour_lineage,
  top1_correct_indicator
)]
setorder(retrieval, cell_line_group_id, rank, tumour_sample_id)

components <- fread(components_tsv)
setnames(components, c("tumour", "component"),
         c("tumour_sample_id", "consensus_component_id"))
components <- merge(
  components,
  unique(score_dt[, .(tumour_sample_id, tumour_lineage)]),
  by = "tumour_sample_id",
  all.x = TRUE,
  sort = FALSE
)
component_labels <- components[, {
  counts <- table(tumour_lineage)
  dominant <- names(which.max(counts))
  purity <- as.numeric(max(counts)) / .N
  .(
    consensus_component_label = sprintf(
      "Component %s - %s dominant (n = %d, purity = %.2f)",
      unique(consensus_component_id), dominant, .N, purity
    )
  )
}, by = consensus_component_id]

top50 <- merge(
  retrieval[rank <= 50L, .(
    cell_line_group_id,
    cell_line_lineage,
    tumour_sample_id,
    tumour_lineage,
    similarity_score,
    rank,
    same_lineage_indicator
  )],
  components[, .(tumour_sample_id, consensus_component_id)],
  by = "tumour_sample_id",
  all.x = TRUE,
  sort = FALSE
)
top50 <- merge(top50, component_labels, by = "consensus_component_id",
               all.x = TRUE, sort = FALSE)
top50 <- top50[, .(
  cell_line_group_id,
  cell_line_lineage,
  tumour_sample_id,
  tumour_lineage,
  similarity_score,
  rank,
  consensus_component_id,
  consensus_component_label,
  same_lineage_indicator
)]
setorder(top50, cell_line_group_id, rank, tumour_sample_id)

group_metadata <- copy(mapping)
group_metadata[, n_profiles_in_group := .N, by = cell_line_group_id]
group_metadata[, canonical_cell_line_name := gsub("_", "-", cell_line_group_id)]
group_metadata[, original_profile_id := cell_line_model_id]
group_metadata[, cancer_lineage := cell_line_lineage]
group_metadata[, included_in_patient_to_cellline_ranking := TRUE]
group_metadata[, included_in_cell_line_centred_retrieval := TRUE]
group_metadata <- group_metadata[, .(
  cell_line_group_id,
  canonical_cell_line_name,
  original_profile_id,
  cancer_lineage,
  n_profiles_in_group,
  replicate_status,
  included_in_patient_to_cellline_ranking,
  included_in_cell_line_centred_retrieval,
  collapse_rule
)]
setorder(group_metadata, cell_line_group_id, original_profile_id)

features <- read_feature_table(feature_table)
features <- features[order(rank)]
features[, gene_id := sub("\\..*$", "", ensembl_gene_id)]
features <- unique(features, by = "gene_id")
features <- head(features, 257L)

annotation <- fread(require_file(args$gene_annotation))
setnames(annotation, "ensembl_id", "gene_id", skip_absent = TRUE)
features <- merge(
  features,
  annotation[, .(gene_id, gene_symbol = symbol)],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)
features[is.na(gene_symbol), gene_symbol := ""]
features[, source_lineage_or_contrast := fifelse(
  nzchar(recurrent_disease) & nzchar(isolate_source_contrasts),
  paste(recurrent_disease, isolate_source_contrasts, sep = "; "),
  fifelse(nzchar(recurrent_disease), recurrent_disease, isolate_source_contrasts)
)]
features[, selection_reason := fifelse(
  source_category == "recurrence-supported core" & nzchar(isolate_source_contrasts),
  "recurrence-supported core marker retained in pan-cancer feature set; also supported by isolate rescue contrast",
  fifelse(
    source_category == "recurrence-supported core",
    "recurrence-supported core marker retained in pan-cancer feature set",
    "isolate rescue marker retained in pan-cancer feature set"
  )
)]
feature_set <- features[order(rank), .(
  gene_id,
  gene_symbol,
  source_category,
  source_lineage_or_contrast,
  selection_reason
)]

validate_outputs(
  patient_rankings,
  ecdf_input,
  ecdf_summary,
  retrieval,
  group_metadata,
  feature_set
)

write_tsv(patient_rankings[, .(
  tumour_sample_id,
  tumour_lineage,
  tumour_project_or_cohort,
  cell_line_model_id,
  cell_line_group_id,
  cell_line_lineage,
  similarity_metric,
  similarity_score,
  rank,
  is_focal_lineage,
  is_non_focal_lineage
)], file.path(out_dir, "patient_to_cellline_rankings.tsv"))
write_tsv(ecdf_input, file.path(out_dir, "ecdf_input_table.tsv"))
write_tsv(ecdf_summary, file.path(out_dir, "ecdf_non_focal_summary.tsv"))
write_tsv(retrieval, file.path(out_dir, "cell_line_centred_retrieval_results.tsv"))
write_tsv(top50, file.path(out_dir, "top50_tumour_neighbourhoods_per_cell_line.tsv"))
write_tsv(group_metadata, file.path(out_dir, "cell_line_group_metadata.tsv"))
write_tsv(feature_set, file.path(out_dir, "feature_set_257_genes.tsv"))

summary_dt <- data.table(
  check = c(
    "patient tumour samples",
    "ranked cell line models per tumour",
    "biological cell line groups",
    "top-1 correct biological cell line groups",
    "retinoblastoma biological cell line groups",
    "feature genes exported"
  ),
  value = c(
    uniqueN(patient_rankings$tumour_sample_id),
    unique(patient_rankings[, .N, by = tumour_sample_id]$N),
    uniqueN(retrieval$cell_line_group_id),
    sum(unique(retrieval[, .(cell_line_group_id, top1_correct_indicator)])$top1_correct_indicator),
    unique(retrieval[cell_line_lineage == "retinoblastoma"]$cell_line_group_id) |> length(),
    nrow(feature_set)
  )
)
write_tsv(summary_dt, file.path(out_dir, "validation_summary.tsv"))

cat("Export complete: ", out_dir, "\n", sep = "")
print(summary_dt)
