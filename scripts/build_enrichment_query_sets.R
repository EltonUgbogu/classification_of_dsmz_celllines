#!/usr/bin/env Rscript

# Build g:Profiler-ready enrichment inputs from DESeq2 marker outputs.
#
# The script converts isolate-vs-rest and component-vs-rest DESeq2 results into
# query gene sets, matched custom backgrounds, and a manifest consumed by
# run_gprofiler_from_manifest.R. Query sets use DESeq2-derived markers; VST
# expression values are not used for enrichment testing.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--profile", type = "character"),
  make_option("--isolate-dir", type = "character"),
  make_option("--component-dir", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--min-query-genes", type = "integer", default = 5),
  make_option("--recurrence-k", type = "integer", default = 2),
  make_option("--component-marker-pooling", type = "character", default = "union"),
  make_option("--component-recurrence-min", type = "integer", default = 1),
  make_option("--component-rank-stat", type = "character", default = "median_abs_log2fc"),
  make_option("--operative-feature-set-gene-file", type = "character", default = ""),
  make_option("--operative-feature-set-rank-tsv", type = "character", default = ""),
  make_option("--recurrent-union-primary", type = "logical", default = TRUE),
  make_option("--operative-feature-set-primary", type = "logical", default = TRUE),
  make_option("--ordered-recurrent-union", type = "logical", default = TRUE),
  make_option("--ordered-operative-feature-set", type = "logical", default = FALSE),
  make_option("--disease-profiles", type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$profile) || is.null(opt$`isolate-dir`) ||
    is.null(opt$`component-dir`) || is.null(opt$outdir)) {
  stop("--profile, --isolate-dir, --component-dir, and --outdir are required", call. = FALSE)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
gene_root <- file.path(opt$outdir, "gene_sets")
dir.create(gene_root, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) message(sprintf(...))

# Remove Ensembl version suffixes so marker and background IDs match g:Profiler.
clean_gene <- function(x) sub("\\.[0-9]+$", "", trimws(as.character(x)))
present <- function(x) !is.null(x) && length(x) == 1 && !is.na(x) && nzchar(x)
sanitize_id <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
empty_dt <- function() data.table(gene_id = character(), log2FoldChange = numeric(),
                                  padj = numeric(), normalised_count_in_test_sample = numeric())

# Read a DESeq2 table while tolerating older column naming variants.
read_table <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(empty_dt())
  x <- fread(path)
  if (!"gene_id" %in% names(x)) names(x)[1] <- "gene_id"
  if (!"normalised_count_in_test_sample" %in% names(x)) {
    alt <- intersect(names(x), c("normalized_count_in_test_sample", "test_sample_normalised_count",
                                 "test_sample_normalized_count"))
    if (length(alt) > 0) {
      setnames(x, alt[1], "normalised_count_in_test_sample")
    } else {
      x[, normalised_count_in_test_sample := NA_real_]
    }
  }
  x[, gene_clean := clean_gene(gene_id)]
  x
}

# Derive a query-specific custom background from genes tested by DESeq2 and
# sufficiently expressed in the test sample. Per-contrast query sets use this
# contrast-specific DESeq2-tested background; higher-level query sets combine
# eligible marker-selection backgrounds across contributing contrasts.
background_from_table <- function(tbl) {
  if (nrow(tbl) == 0 || !"normalised_count_in_test_sample" %in% names(tbl)) return(character())
  unique(tbl[!is.na(padj) & !is.na(normalised_count_in_test_sample) &
               normalised_count_in_test_sample >= 10, gene_clean])
}

read_gene_list <- function(path) {
  if (!present(path) || !file.exists(path) || file.info(path)$size == 0) return(character())
  unique(clean_gene(readLines(path, warn = FALSE)))
}

resolve_path <- function(path, base_dir) {
  if (!present(path)) return("")
  if (grepl("^/", path)) return(path)
  file.path(base_dir, path)
}

find_marker_file <- function(markers_dir, contrast) {
  hits <- list.files(markers_dir, pattern = paste0("^", gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", contrast), ".*\\.txt$"),
                     full.names = TRUE)
  if (length(hits) == 0) return("")
  hits[1]
}

# Write one query directory containing the submitted genes, custom background,
# and optional rank statistics for ordered enrichment.
write_query_files <- function(query_id, genes, background, ranks = NULL) {
  qdir <- file.path(gene_root, query_id)
  dir.create(qdir, recursive = TRUE, showWarnings = FALSE)
  genes_dt <- data.table(gene_id = unique(genes))
  bg_dt <- data.table(gene_id = unique(background))
  genes_path <- file.path(qdir, "genes.tsv")
  bg_path <- file.path(qdir, "background.tsv")
  fwrite(genes_dt, genes_path, sep = "\t")
  fwrite(bg_dt, bg_path, sep = "\t")
  ranked_path <- file.path(qdir, "ranked_genes.tsv")
  if (!is.null(ranks) && nrow(ranks) > 0) {
    ranks <- unique(ranks[order(-rank_stat)], by = "gene_id")
    fwrite(ranks[, .(gene_id, rank_stat)], ranked_path, sep = "\t")
  } else {
    fwrite(data.table(gene_id = character(), rank_stat = numeric()), ranked_path, sep = "\t")
  }
  list(genes = genes_path, background = bg_path, ranked = ranked_path)
}

rows <- list()
configured_final_feature_genes <- read_gene_list(opt$`operative-feature-set-gene-file`)
final_feature_set_label <- "final pan-cancer feature set"

query_family_label <- function(family) {
  labels <- c(
    per_contrast = "per-contrast marker set",
    component_or_isolate = "component or isolate marker set",
    disease_recurrent = "disease recurrence-filtered marker set",
    recurrence_filtered_marker_union = "pan-cancer recurrence-filtered marker union",
    ordered_recurrence_filtered_marker_union = "ordered recurrence-filtered marker union",
    final_pan_cancer_feature_set = final_feature_set_label,
    ordered_final_pan_cancer_feature_set = paste("ordered", final_feature_set_label)
  )
  label <- labels[[family]]
  if (is.null(label) || is.na(label)) return(gsub("_", " ", family, fixed = TRUE))
  label
}

query_display_label <- function(family, disease, direction, ordered) {
  label <- query_family_label(family)
  if (family %in% c(
    "recurrence_filtered_marker_union",
    "ordered_recurrence_filtered_marker_union",
    "disease_recurrent"
  ) && direction != "all") {
    direction_label <- sub("_ordered$", "", direction)
    label <- paste0(label, " (", direction_label, ")")
  }
  if (family %in% c("per_contrast", "component_or_isolate") && direction != "all") {
    label <- paste0(label, " (", direction, ")")
  }
  label
}

# Add one row to the query manifest and record whether it is analyzable. Skipped
# rows are still written so downstream reports can explain absent results.
add_manifest_row <- function(query_id, family, profile, disease, group_id, contrast_id,
                             direction, ordered, rank_source, genes, background,
                             ranks = NULL, skip_reason = "") {
  query_id <- sanitize_id(query_id)
  genes <- unique(genes[nzchar(genes)])
  background <- unique(background[nzchar(background)])
  skip <- FALSE
  if (!present(skip_reason)) {
    if (length(genes) < opt$`min-query-genes`) {
      skip <- TRUE
      skip_reason <- "fewer_than_min_genes"
    } else if (length(background) == 0) {
      skip <- TRUE
      skip_reason <- "missing_background"
    } else if (ordered && (is.null(ranks) || nrow(ranks) == 0)) {
      skip <- TRUE
      skip_reason <- "missing_rank_statistics"
    }
  } else {
    skip <- TRUE
  }
  files <- write_query_files(query_id, genes, background, ranks)
  family_label <- query_family_label(family)
  label <- query_display_label(family, disease, direction, ordered)
  rows[[length(rows) + 1]] <<- data.table(
    query_id = query_id,
    query_family = family,
    query_family_label = family_label,
    query_label = label,
    profile = profile,
    disease = disease,
    group_id = group_id,
    contrast_id = contrast_id,
    direction = direction,
    ordered = ordered,
    rank_source = rank_source,
    genes_tsv = files$genes,
    background_tsv = files$background,
    ranked_genes_tsv = files$ranked,
    gene_count = length(genes),
    background_count = length(background),
    interpretation_hint = paste(
      label,
      ifelse(ordered, "ordered enrichment; rank_stat is used", "over-representation query"),
      sep = " | "
    ),
    skip = skip,
    skip_reason = ifelse(skip, skip_reason, "")
  )
}

# Collect all DESeq2-derived query ingredients for one disease profile.
collect_profile <- function(profile, profile_root) {
  isolate_dir <- file.path(profile_root, "deseq2_markers")
  component_dir <- file.path(profile_root, "deseq2", "component_vs_rest")
  out <- list(contrasts = list(), recurrent = character(), backgrounds = character(),
              up = character(), down = character(), ranks_up = data.table(), ranks_down = data.table())

  manifest_path <- file.path(isolate_dir, "markers", "marker_sets_manifest.tsv")
  if (file.exists(manifest_path)) {
    manifest <- fread(manifest_path, fill = TRUE)
    markers_dir <- file.path(isolate_dir, "markers")
    tables_dir <- file.path(isolate_dir, "tables")

    # Isolate queries preserve directionality: positive log2FC means enriched
    # in the isolate cell line relative to the remaining DSMZ samples.
    for (i in seq_len(nrow(manifest))) {
      contrast <- manifest$contrast[i]
      table_path <- if ("table_file" %in% names(manifest) && present(manifest$table_file[i])) {
        resolve_path(manifest$table_file[i], isolate_dir)
      } else {
        file.path(tables_dir, paste0(contrast, ".tsv"))
      }
      marker_path <- if ("marker_file" %in% names(manifest) && present(manifest$marker_file[i])) {
        resolve_path(manifest$marker_file[i], isolate_dir)
      } else {
        find_marker_file(markers_dir, contrast)
      }
      tbl <- read_table(table_path)
      bg <- background_from_table(tbl)
      genes <- read_gene_list(marker_path)
      signed <- tbl[gene_clean %in% genes]
      up <- unique(signed[log2FoldChange > 0, gene_clean])
      down <- unique(signed[log2FoldChange < 0, gene_clean])
      out$contrasts[[length(out$contrasts) + 1]] <- list(
        profile = profile, family = "per_contrast", group = "isolate",
        contrast = contrast, background = bg, all = genes, up = up, down = down,
        ranks_up = signed[log2FoldChange > 0, .(gene_id = gene_clean, rank_stat = abs(log2FoldChange))],
        ranks_down = signed[log2FoldChange < 0, .(gene_id = gene_clean, rank_stat = abs(log2FoldChange))]
      )
    }
  }

  comp_marker_dir <- file.path(component_dir, "markers")
  comp_files <- list.files(comp_marker_dir, pattern = "^component_[^/]+_vs_rest_(UP|DOWN)_top[0-9]+\\.tsv$",
                           full.names = TRUE)
  if (length(comp_files) > 0) {
    # Component queries pool markers across UP/DOWN component marker files. The
    # recurrence threshold controls how often a marker must appear to be retained.
    comp_groups <- split(comp_files, sub("^component_([^_]+)_vs_rest_.*$", "\\1", basename(comp_files)))
    for (comp in names(comp_groups)) {
      bg <- character()
      comp_up <- character()
      comp_down <- character()
      ranks_up <- data.table()
      ranks_down <- data.table()
      for (fp in comp_groups[[comp]]) {
        table_path <- file.path(component_dir, "tables", paste0("component_", comp, "_vs_rest.tsv"))
        tbl <- read_table(table_path)
        bg <- union(bg, background_from_table(tbl))
        markers <- fread(fp)
        if (!"gene_id" %in% names(markers)) names(markers)[1] <- "gene_id"
        markers[, gene_clean := clean_gene(gene_id)]
        markers <- merge(markers[, .(gene_clean)], tbl, by = "gene_clean", all.x = TRUE)
        if (grepl("_UP_", basename(fp))) {
          comp_up <- c(comp_up, markers$gene_clean)
          ranks_up <- rbind(ranks_up, markers[, .(gene_id = gene_clean, rank_stat = abs(log2FoldChange))],
                            fill = TRUE)
        } else {
          comp_down <- c(comp_down, markers$gene_clean)
          ranks_down <- rbind(ranks_down, markers[, .(gene_id = gene_clean, rank_stat = abs(log2FoldChange))],
                              fill = TRUE)
        }
      }
      summarise_ranks <- function(x) {
        if (nrow(x) == 0) return(x)
        x[!is.na(rank_stat), .(rank_stat = median(rank_stat, na.rm = TRUE)), by = gene_id]
      }
      up_counts <- table(comp_up)
      down_counts <- table(comp_down)
      comp_up <- names(up_counts)[up_counts >= opt$`component-recurrence-min`]
      comp_down <- names(down_counts)[down_counts >= opt$`component-recurrence-min`]
      out$contrasts[[length(out$contrasts) + 1]] <- list(
        profile = profile, family = "component_or_isolate", group = paste0("component_", comp),
        contrast = paste0("component_", comp, "_vs_rest"), background = bg,
        all = union(comp_up, comp_down), up = unique(comp_up), down = unique(comp_down),
        ranks_up = summarise_ranks(ranks_up[gene_id %in% comp_up]),
        ranks_down = summarise_ranks(ranks_down[gene_id %in% comp_down])
      )
    }
  }

  recurrent_file <- file.path(isolate_dir, "markers", paste0("unique_feature_set_recurrence_ge_", opt$`recurrence-k`, ".txt"))
  out$recurrent <- read_gene_list(recurrent_file)

  # Disease-level summaries combine contrast-level query sets and backgrounds.
  for (x in out$contrasts) {
    out$backgrounds <- union(out$backgrounds, x$background)
    out$up <- union(out$up, x$up)
    out$down <- union(out$down, x$down)
    out$ranks_up <- rbind(out$ranks_up, x$ranks_up, fill = TRUE)
    out$ranks_down <- rbind(out$ranks_down, x$ranks_down, fill = TRUE)
  }
  rank_collapse <- function(x) {
    if (nrow(x) == 0) return(x)
    x[!is.na(rank_stat), .(rank_stat = median(rank_stat, na.rm = TRUE)), by = gene_id]
  }
  out$ranks_up <- rank_collapse(out$ranks_up)
  out$ranks_down <- rank_collapse(out$ranks_down)
  out
}

# Disease profiles are resolved relative to the unsupervised-results root.
profile_root <- dirname(opt$`isolate-dir`)
unsup_root <- dirname(profile_root)
disease_profiles <- strsplit(opt$`disease-profiles`, ",", fixed = TRUE)[[1]]
disease_profiles <- disease_profiles[nzchar(disease_profiles)]
if (length(disease_profiles) == 0) disease_profiles <- opt$profile

profiles <- list()
for (p in disease_profiles) {
  p_root <- file.path(unsup_root, p)
  if (dir.exists(p_root)) profiles[[p]] <- collect_profile(p, p_root)
}
if (!opt$profile %in% names(profiles)) {
  profiles[[opt$profile]] <- collect_profile(opt$profile, profile_root)
}

# Per-contrast queries for the active profile support cell-line- or
# component-specific biological interpretation.
current <- profiles[[opt$profile]]
for (x in current$contrasts) {
  base_id <- paste(opt$profile, x$family, x$contrast, sep = "__")
  add_manifest_row(paste(base_id, "all", sep = "__"), x$family, opt$profile, opt$profile,
                   x$group, x$contrast, "all", FALSE, "", x$all, x$background)
  add_manifest_row(paste(base_id, "up", sep = "__"), x$family, opt$profile, opt$profile,
                   x$group, x$contrast, "up", FALSE, "", x$up, x$background)
  add_manifest_row(paste(base_id, "down", sep = "__"), x$family, opt$profile, opt$profile,
                   x$group, x$contrast, "down", FALSE, "", x$down, x$background)
}

# Disease-level recurrent queries summarise markers repeatedly observed across
# isolate contrasts within each disease profile.
for (p in names(profiles)) {
  pr <- profiles[[p]]
  add_manifest_row(paste(p, "disease_recurrent_ge", opt$`recurrence-k`, "all", sep = "__"),
                   "disease_recurrent", opt$profile, p, paste0("recurrence_ge_", opt$`recurrence-k`),
                   "", "all", FALSE, "", pr$recurrent, pr$backgrounds)
  add_manifest_row(paste(p, "disease_recurrent_ge", opt$`recurrence-k`, "up", sep = "__"),
                   "disease_recurrent", opt$profile, p, paste0("recurrence_ge_", opt$`recurrence-k`),
                   "", "up", FALSE, "", intersect(pr$recurrent, pr$up), pr$backgrounds)
  add_manifest_row(paste(p, "disease_recurrent_ge", opt$`recurrence-k`, "down", sep = "__"),
                   "disease_recurrent", opt$profile, p, paste0("recurrence_ge_", opt$`recurrence-k`),
                   "", "down", FALSE, "", intersect(pr$recurrent, pr$down), pr$backgrounds)
}

recurrent_union <- Reduce(union, lapply(profiles, `[[`, "recurrent"))
recurrent_bg <- Reduce(union, lapply(profiles, `[[`, "backgrounds"))
recurrent_up <- Reduce(union, lapply(profiles, `[[`, "up"))
recurrent_down <- Reduce(union, lapply(profiles, `[[`, "down"))
recurrent_ranks_up <- rbindlist(lapply(profiles, `[[`, "ranks_up"), fill = TRUE)
recurrent_ranks_down <- rbindlist(lapply(profiles, `[[`, "ranks_down"), fill = TRUE)
collapse_ranks <- function(x, genes) {
  if (nrow(x) == 0) return(data.table(gene_id = character(), rank_stat = numeric()))
  x[gene_id %in% genes & !is.na(rank_stat), .(rank_stat = median(rank_stat, na.rm = TRUE)), by = gene_id]
}
recurrent_ranks_up <- collapse_ranks(recurrent_ranks_up, intersect(recurrent_union, recurrent_up))
recurrent_ranks_down <- collapse_ranks(recurrent_ranks_down, intersect(recurrent_union, recurrent_down))

# Pan-cancer recurrence-filtered marker-union queries test recurrent disease
# marker sets against the combined DESeq2-tested background.
if (isTRUE(opt$`recurrent-union-primary`)) {
  add_manifest_row("recurrence_filtered_marker_union__all", "recurrence_filtered_marker_union",
                   opt$profile, "pan_cancer", "recurrence_filtered_marker_union", "",
                   "all", FALSE, "", recurrent_union, recurrent_bg)
  add_manifest_row("recurrence_filtered_marker_union__up", "recurrence_filtered_marker_union",
                   opt$profile, "pan_cancer", "recurrence_filtered_marker_union", "",
                   "up", FALSE, "", intersect(recurrent_union, recurrent_up), recurrent_bg)
  add_manifest_row("recurrence_filtered_marker_union__down", "recurrence_filtered_marker_union",
                   opt$profile, "pan_cancer", "recurrence_filtered_marker_union", "",
                   "down", FALSE, "", intersect(recurrent_union, recurrent_down), recurrent_bg)
}
if (isTRUE(opt$`ordered-recurrent-union`) && isTRUE(opt$`recurrent-union-primary`)) {
  add_manifest_row("ordered_recurrence_filtered_marker_union__up_ordered",
                   "ordered_recurrence_filtered_marker_union", opt$profile,
                   "pan_cancer", "recurrence_filtered_marker_union", "", "up_ordered", TRUE,
                   opt$`component-rank-stat`, recurrent_ranks_up$gene_id, recurrent_bg, recurrent_ranks_up)
  add_manifest_row("ordered_recurrence_filtered_marker_union__down_ordered",
                   "ordered_recurrence_filtered_marker_union", opt$profile,
                   "pan_cancer", "recurrence_filtered_marker_union", "", "down_ordered", TRUE,
                   opt$`component-rank-stat`, recurrent_ranks_down$gene_id, recurrent_bg, recurrent_ranks_down)
}

# The final pan-cancer feature set is the downstream selected gene panel, tested
# separately from the recurrence-filtered recurrent-marker union.
final_feature_genes <- configured_final_feature_genes
if (length(final_feature_genes) > 0) {
  if (isTRUE(opt$`operative-feature-set-primary`)) {
    add_manifest_row("final_pan_cancer_feature_set__all", "final_pan_cancer_feature_set",
                     opt$profile, "pan_cancer", "final_pan_cancer_feature_set", "",
                     "all", FALSE, "", final_feature_genes, recurrent_bg)
  }
} else if (isTRUE(opt$`operative-feature-set-primary`)) {
  add_manifest_row("final_pan_cancer_feature_set__all", "final_pan_cancer_feature_set",
                   opt$profile, "pan_cancer", "final_pan_cancer_feature_set", "",
                   "all", FALSE, "", character(), recurrent_bg,
                   skip_reason = "missing_final_pan_cancer_feature_set_gene_file")
}

if (isTRUE(opt$`ordered-operative-feature-set`) && isTRUE(opt$`operative-feature-set-primary`)) {
  if (present(opt$`operative-feature-set-rank-tsv`) && file.exists(opt$`operative-feature-set-rank-tsv`)) {
    ranks <- fread(opt$`operative-feature-set-rank-tsv`)
    if (!all(c("gene_id", "rank_stat", "direction") %in% names(ranks))) {
      add_manifest_row("ordered_final_pan_cancer_feature_set__up_ordered",
                       "ordered_final_pan_cancer_feature_set",
                       opt$profile, "pan_cancer", "final_pan_cancer_feature_set", "", "up_ordered",
                       TRUE, "initial_run_log2fc", character(), recurrent_bg,
                       skip_reason = "missing_rank_statistics")
    } else {
      ranks[, gene_id := clean_gene(gene_id)]
      add_manifest_row("ordered_final_pan_cancer_feature_set__up_ordered",
                       "ordered_final_pan_cancer_feature_set",
                       opt$profile, "pan_cancer", "final_pan_cancer_feature_set", "", "up_ordered",
                       TRUE, "initial_run_log2fc", ranks[direction == "up", gene_id], recurrent_bg,
                       ranks[direction == "up", .(gene_id, rank_stat)])
      add_manifest_row("ordered_final_pan_cancer_feature_set__down_ordered",
                       "ordered_final_pan_cancer_feature_set",
                       opt$profile, "pan_cancer", "final_pan_cancer_feature_set", "", "down_ordered",
                       TRUE, "initial_run_log2fc", ranks[direction == "down", gene_id], recurrent_bg,
                       ranks[direction == "down", .(gene_id, rank_stat)])
    }
  } else {
    add_manifest_row("ordered_final_pan_cancer_feature_set__up_ordered",
                     "ordered_final_pan_cancer_feature_set",
                     opt$profile, "pan_cancer", "final_pan_cancer_feature_set", "", "up_ordered",
                     TRUE, "initial_run_log2fc", character(), recurrent_bg,
                     skip_reason = "missing_rank_statistics")
    add_manifest_row("ordered_final_pan_cancer_feature_set__down_ordered",
                     "ordered_final_pan_cancer_feature_set",
                     opt$profile, "pan_cancer", "final_pan_cancer_feature_set", "", "down_ordered",
                     TRUE, "initial_run_log2fc", character(), recurrent_bg,
                     skip_reason = "missing_rank_statistics")
  }
}

manifest <- rbindlist(rows, fill = TRUE)
setorder(manifest, query_family, disease, group_id, contrast_id, direction)
manifest_path <- file.path(opt$outdir, "query_manifest.tsv")
skipped_path <- file.path(opt$outdir, "skipped_queries.tsv")
summary_path <- file.path(opt$outdir, "query_summary.tsv")
readme_path <- file.path(opt$outdir, "README_enrichment_queries.txt")

query_summary <- manifest[, .(
  n_queries = .N,
  n_run = sum(skip == FALSE),
  n_skipped = sum(skip == TRUE),
  min_gene_count = min(gene_count, na.rm = TRUE),
  median_gene_count = as.numeric(stats::median(gene_count, na.rm = TRUE)),
  max_gene_count = max(gene_count, na.rm = TRUE),
  gene_counts = paste(sort(unique(gene_count)), collapse = ","),
  median_background_count = as.numeric(stats::median(background_count, na.rm = TRUE))
), by = .(query_family, query_family_label, disease, direction, ordered)]
setorder(query_summary, query_family, disease, direction, ordered)

fwrite(manifest, manifest_path, sep = "\t")
fwrite(manifest[skip == TRUE], skipped_path, sep = "\t")
fwrite(query_summary, summary_path, sep = "\t")
writeLines(c(
  "Enrichment query-set build report",
  "",
  "Purpose:",
  "  Builds manifest-driven functional-enrichment query sets from DESeq2 marker outputs.",
  "",
  "Key files:",
  paste0("  query_manifest.tsv: functional-enrichment job manifest with runnable and skipped query records."),
  paste0("  skipped_queries.tsv: skipped query records and reasons."),
  paste0("  query_summary.tsv: runnable/skipped query counts by family, disease, and direction."),
  paste0("  gene_sets/<query_id>/: row-level genes.tsv, background.tsv, and ranked_genes.tsv inputs."),
  "",
  "Interpretation notes:",
  "  direction=up means markers with positive log2FC in the tested group.",
  "  direction=down means markers with negative log2FC in the tested group.",
  "  ordered=TRUE means g:Profiler receives genes ranked by rank_stat.",
  "  Per-contrast backgrounds are DESeq2-tested genes with sufficient normalised expression in the test sample.",
  "  Higher-level query sets use combined eligible marker-selection backgrounds."
), readme_path)
msg("Wrote %d query manifest rows to %s", nrow(manifest), manifest_path)
msg("Wrote %d skipped query rows to %s", nrow(manifest[skip == TRUE]), skipped_path)
msg("Wrote query summary to %s", summary_path)
