#!/usr/bin/env Rscript

# hema_summary_purity.R

# HEMA-wide *summary only*:

# - Reads per-project metadata + purity scores

# - Produces:

#   1) Before vs after barplot (CLL / BL / MPN)

#   2) Purity histogram faceted by cohort

#

# IMPORTANT:

# - Does NOT create any matrices or inputs for the main pipeline.

# - Only consumes per-project outputs and writes summary plots.



suppressPackageStartupMessages({

  library(dplyr)

  library(tibble)

  library(ggplot2)

  library(readr)

  library(RColorBrewer)

  library(purrr)

})



`%||%` <- function(x, y) if (!is.null(x)) x else y



# -------------------------------------------------------------------

# Config: per-project directories (modular, project-first)

# -------------------------------------------------------------------

HEMA_PROJECTS <- tibble::tribble(

  ~cohort, ~project_id, ~base_dir,

  "CLL",   "GSE237907", "/home/chu25/data/cll/GSE237907",

  "BL",    "GSE199413", "/home/chu25/data/pri_bl/GSE199413",

  "MPN",   "GSE228758", "/home/chu25/data/mpn/GSE228758"

)



PURITY_THRESHOLD_DEFAULT <- 0.7

OUTDIR_DEFAULT <- "/home/chu25/data/hema_summary"



# -------------------------------------------------------------------

# CLI parser: --key value

#   --outdir    where to write summary PDFs

#   --threshold purity threshold (default 0.7)

# -------------------------------------------------------------------

parse_cli <- function(args) {

  res <- list()

  i <- 1

  while (i <= length(args)) {

    key <- args[i]

    if (!startsWith(key, "--")) stop("Expected named arg starting with '--', got: ", key)

    if (i == length(args)) stop("No value provided for ", key)

    val <- args[i + 1]

    res[[sub("^--", "", key)]] <- val

    i <- i + 2

  }

  res

}



# -------------------------------------------------------------------

# Helpers: load metadata + purity per project

# -------------------------------------------------------------------

load_project_meta <- function(cohort, project_id, base_dir) {

  agg_dir  <- file.path(base_dir, "results", "aggregated")

  meta_csv <- file.path(agg_dir, "sample_metadata.csv")



  if (!file.exists(meta_csv)) {

    stop("[ERROR] sample_metadata.csv not found for ", cohort,

         " (", project_id, "): ", meta_csv)

  }



  meta <- readr::read_csv(meta_csv, show_col_types = FALSE)

  if (!("sample" %in% names(meta))) {

    stop("[ERROR] 'sample' column missing in ", meta_csv)

  }



  meta %>%

    mutate(

      cohort    = cohort,

      project   = project_id,

      base_dir  = base_dir

    )

}



load_project_purity <- function(cohort, project_id, base_dir) {

  agg_dir    <- file.path(base_dir, "results", "aggregated")

  purity_dir <- file.path(agg_dir, "tumour_purity")

  scores_csv <- file.path(purity_dir, "hema_estimate_scores.csv")



  if (!file.exists(scores_csv)) {

    stop("[ERROR] hema_estimate_scores.csv not found for ", cohort,

         " (", project_id, "): ", scores_csv)

  }



  # We wrote this with row.names = TRUE → first column is sample ID

  scores <- read.csv(scores_csv, row.names = 1, check.names = FALSE)

  scores <- scores %>%

    as.data.frame() %>%

    tibble::rownames_to_column("sample") %>%

    mutate(

      cohort  = cohort,

      project = project_id,

      base_dir = base_dir

    )



  if (!("tumourPurity" %in% colnames(scores))) {

    stop("[ERROR] 'tumourPurity' column missing in ", scores_csv,

         " (cohort ", cohort, ", project ", project_id, ")")

  }



  scores

}



# -------------------------------------------------------------------

# Main summary logic

# -------------------------------------------------------------------

run_hema_summary <- function(outdir = OUTDIR_DEFAULT,

                             purity_threshold = PURITY_THRESHOLD_DEFAULT) {

  message("[INFO] HEMA summary starting...")

  message("[INFO] Output directory: ", outdir)

  message("[INFO] Purity threshold: ", purity_threshold)



  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  outdir <- normalizePath(outdir, mustWork = FALSE)



  # 1) Load per-project metadata + purity

  meta_list <- purrr::pmap(

    HEMA_PROJECTS,

    ~ load_project_meta(..1, ..2, ..3)

  )

  purity_list <- purrr::pmap(

    HEMA_PROJECTS,

    ~ load_project_purity(..1, ..2, ..3)

  )



  meta_all   <- dplyr::bind_rows(meta_list)

  purity_all <- dplyr::bind_rows(purity_list)



  # 2) Join metadata + purity on (cohort, project, sample)

  joined <- meta_all %>%

    select(cohort, project, sample, everything()) %>%

    left_join(

      purity_all %>% select(cohort, project, sample, tumourPurity),

      by = c("cohort", "project", "sample")

    )



  if (any(is.na(joined$tumourPurity))) {

    n_na <- sum(is.na(joined$tumourPurity))

    message("[WARN] ", n_na, " samples have NA tumourPurity after join.")

  }



  # -----------------------------------------------------------------

  # A) BEFORE vs AFTER barplot (by cohort)

  # -----------------------------------------------------------------

  retention_by_cohort <- joined %>%

    group_by(cohort) %>%

    summarise(

      n_before = n(),

      n_after  = sum(!is.na(tumourPurity) & tumourPurity >= purity_threshold),

      .groups = "drop"

    ) %>%

    tidyr::pivot_longer(

      cols = c(n_before, n_after),

      names_to = "stage",

      values_to = "count"

    ) %>%

    mutate(

      stage = dplyr::recode(

        stage,

        n_before = "Before filtering",

        n_after  = sprintf("After purity ≥ %.2f", purity_threshold)

      ),

      stage = factor(stage, levels = c("Before filtering",

                                       sprintf("After purity ≥ %.2f", purity_threshold)))

    )



  message("[INFO] Sample counts by cohort (before / after):")

  print(retention_by_cohort)



  n_cohorts <- length(unique(retention_by_cohort$cohort))

  pal <- RColorBrewer::brewer.pal(max(3, n_cohorts), "Set2")[seq_len(n_cohorts)]

  names(pal) <- sort(unique(retention_by_cohort$cohort))



  p_ret <- ggplot(retention_by_cohort,

                  aes(x = cohort, y = count, fill = stage)) +

    geom_col(position = position_dodge(width = 0.8),

             width = 0.7, colour = "white", linewidth = 0.6) +

    geom_text(

      aes(label = count),

      position = position_dodge(width = 0.8),

      vjust = -0.6,

      size = 4.0,

      fontface = "bold"

    ) +

    scale_fill_manual(values = c("Before filtering" = "#D3D3D3",

                                 sprintf("After purity ≥ %.2f", purity_threshold) = "#4DAF4A"),

                      name = NULL) +

    labs(

      title    = "HEMA: Sample retention by cohort",

      subtitle = sprintf("Purity threshold: %.2f", purity_threshold),

      x        = "Cohort",

      y        = "Number of samples"

    ) +

    theme_minimal(base_size = 14) +

    theme(

      plot.title    = element_text(face = "bold", hjust = 0.5),

      plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),

      axis.title.x  = element_text(face = "bold"),

      axis.title.y  = element_text(face = "bold"),

      legend.position = "top",

      panel.grid.major.x = element_blank(),

      panel.grid.minor   = element_blank(),

      panel.border       = element_rect(colour = "grey80", fill = NA, linewidth = 0.5)

    )



  bar_pdf <- file.path(outdir, "hema_sample_retention_by_cohort.pdf")

  ggsave(bar_pdf, p_ret, width = 8, height = 6, dpi = 300)

  message("[SUCCESS] Retention barplot saved: ", bar_pdf)



  # -----------------------------------------------------------------

  # B) Purity histogram faceted by cohort

  # -----------------------------------------------------------------

  purity_df <- joined %>%

    filter(!is.na(tumourPurity)) %>%

    mutate(

      PurityGroup = cut(

        tumourPurity,

        breaks = c(-Inf, purity_threshold, Inf),

        labels = c(paste0("< ", purity_threshold),

                   paste0("≥ ", purity_threshold)),

        include.lowest = TRUE

      )

    )



  p_hist <- ggplot(purity_df, aes(x = tumourPurity, fill = cohort)) +

    geom_histogram(bins = 30, alpha = 0.85, colour = "white") +

    geom_vline(xintercept = purity_threshold,

               linetype = "dashed",

               colour   = "firebrick",

               linewidth = 1) +

    scale_fill_manual(values = pal, name = "Cohort") +

    labs(

      title = "HEMA: Distribution of tumour purity by cohort",

      x     = "Tumour purity",

      y     = "Number of samples"

    ) +

    facet_wrap(~ cohort, ncol = 1, scales = "free_y") +

    theme_minimal(base_size = 13) +

    theme(

      plot.title       = element_text(face = "bold", hjust = 0.5),

      legend.position  = "none",

      panel.grid.minor = element_blank(),

      panel.border     = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),

      strip.text       = element_text(face = "bold")

    )



  hist_pdf <- file.path(outdir, "hema_purity_histogram_by_cohort.pdf")

  ggsave(hist_pdf, p_hist, width = 9, height = 9, dpi = 300)

  message("[SUCCESS] Purity histogram saved: ", hist_pdf)



  message("[INFO] HEMA summary finished.")

  invisible(list(

    retention = retention_by_cohort,

    purity_df = purity_df

  ))

}



# -------------------------------------------------------------------

# MAIN

# -------------------------------------------------------------------

if (sys.nframe() == 0) {

  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {

    # Defaults

    run_hema_summary()

  } else {

    opts <- parse_cli(args)

    outdir <- opts[["outdir"]]    %||% OUTDIR_DEFAULT

    thr    <- opts[["threshold"]] %||% as.character(PURITY_THRESHOLD_DEFAULT)

    thr    <- as.numeric(thr)

    run_hema_summary(outdir = outdir, purity_threshold = thr)

  }

}
