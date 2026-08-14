## Batch-effect quantification and paired before/after PCA figures.
##
## One script for BRCA, NBL and RBL. It consumes the joint pre- and
## post-ComBat-seq VST matrices that the cohort preprocessing rules already
## produced and never rebuilds them, so the derived statistics and figures can
## be deleted and regenerated on their own without touching a count or VST
## matrix. Everything numerical lives in the shared module; this file is the
## Snakemake interface.
##
## Reproducibility contract:
##   * only the pinned in-repo GENCODE GTF is read (declared as a rule input);
##     no BioMart call, no download;
##   * one ordered gene-ID manifest per feature space is written and asserted to
##     index the pre- and post-correction matrices identically;
##   * every permutation test draws a seed derived from the master seed and the
##     (cohort, feature_space, predictor, stage, method) tuple, and that seed is
##     written into the result row;
##   * every output is promoted by rename, so an interrupted job leaves nothing
##     that looks complete;
##   * batch_effect_quantification_provenance.tsv is deterministic; volatile run
##     facts go to batch_effect_run_metadata.tsv.

options(stringsAsFactors = FALSE)

started_at <- Sys.time()

inputs <- snakemake@input
outputs <- snakemake@output
params <- snakemake@params

for (path in unname(outputs)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

if (length(snakemake@log) > 0) {
  log_file <- snakemake@log[[1]]
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)
}

source(inputs[["shared_figure_module"]])

COHORT <- as.character(params[["cohort"]])
repo_root <- as.character(params[["repo_root"]])
master_seed <- as.integer(params[["master_seed"]])
qc_top_genes <- as.integer(params[["qc_top_genes"]])
permutations <- as.integer(params[["permanova_permutations"]])
primary_predictor <- as.character(params[["primary_predictor"]])

cat(sprintf("[INFO] %s batch-effect quantification\n", COHORT))
cat(sprintf("[INFO] master seed %d; %s permutations per test\n",
            master_seed, figure_count(permutations)))
cat(sprintf("[INFO] seed derivation: %s\n", SEED_DERIVATION_RULE))

## ---------------------------------------------------------------------------
## Inputs. The VST matrices are read only; this rule never writes them.
## ---------------------------------------------------------------------------
joint_vst_pre <- readRDS(inputs[["joint_vst_pre_bc_rds"]])
joint_vst_post <- readRDS(inputs[["joint_vst_post_bc_rds"]])
if (!identical(colnames(joint_vst_pre), colnames(joint_vst_post))) {
  stop("[ERROR] Pre- and post-correction VST matrices differ in sample order")
}
if (!identical(rownames(joint_vst_pre), rownames(joint_vst_post))) {
  stop("[ERROR] Pre- and post-correction VST matrices differ in gene order")
}

coldata <- utils::read.delim(
  inputs[["coldata_tsv"]], stringsAsFactors = FALSE, check.names = FALSE
)
if (!"sample_id" %in% colnames(coldata)) {
  stop("[ERROR] coldata must carry a sample_id column: ", inputs[["coldata_tsv"]])
}
rownames(coldata) <- coldata$sample_id
sample_order <- colnames(joint_vst_pre)
if (!all(sample_order %in% rownames(coldata))) {
  stop("[ERROR] coldata does not cover every column of the joint VST matrices")
}
coldata <- coldata[sample_order, , drop = FALSE]

tumour_label <- as.character(params[["tumour_label"]])
dsmz_label <- as.character(params[["dsmz_label"]])
expected_tumours <- as.integer(params[["expected_tumours"]])
expected_collapsed_dsmz <- as.integer(params[["expected_collapsed_dsmz"]])
include_tumour_source_analysis <- isTRUE(params[["include_tumour_source_analysis"]])

manifest_output_key <- c(
  top3000 = "feature_manifest_top3000_tsv",
  protein_coding = "feature_manifest_protein_coding_tsv",
  full_expression = "feature_manifest_full_expression_tsv"
)
audit_manifest_output_key <- c(
  top3000 = "audit_legacy_feature_manifest_top3000_tsv",
  protein_coding = "audit_legacy_feature_manifest_protein_coding_tsv",
  full_expression = "audit_legacy_feature_manifest_full_expression_tsv"
)
write_manifest_block <- function(feature_spaces, output_key_map) {
  manifest_paths <- list()
  for (space_name in names(feature_spaces$spaces)) {
    manifest_path <- outputs[[output_key_map[[space_name]]]]
    write_feature_manifest(feature_spaces$spaces[[space_name]]$ids, manifest_path)
    manifest_paths[[space_name]] <- manifest_path
    cat(sprintf(
      "[INFO] manifest written | %-20s %s genes -> %s\n",
      space_name,
      figure_count(length(feature_spaces$spaces[[space_name]]$ids)),
      basename(manifest_path)
    ))
  }
  manifest_paths
}
verify_written_manifests <- function(feature_spaces, manifest_paths, context_label) {
  for (space_name in names(manifest_paths)) {
    written <- utils::read.delim(manifest_paths[[space_name]], stringsAsFactors = FALSE)
    if (!identical(as.character(written$ensembl_gene_id),
                   as.character(feature_spaces$spaces[[space_name]]$ids))) {
      stop("[ERROR] Written manifest does not match the analysed feature space: ", context_label, " / ", space_name)
    }
    if (!identical(written$order_index, seq_len(nrow(written)))) {
      stop("[ERROR] Manifest order_index is not a dense ascending sequence: ", context_label, " / ", space_name)
    }
  }
  cat(sprintf("[INFO] %s manifests re-read and confirmed identical to the analysed spaces\n", context_label))
}

plotting_pre <- joint_vst_pre
plotting_post <- joint_vst_post
plot_labels <- factor(coldata[sample_order, primary_predictor])
plotting_summary <- "DSMZ cell-line profiles"
active_predictors <- resolve_predictors(coldata, sample_order, primary_predictor, NULL)
active_tumour_ids <- sample_order[coldata[sample_order, primary_predictor] == tumour_label]
active_dsmz_ids <- sample_order[coldata[sample_order, primary_predictor] == dsmz_label]
active_analysis_note <- "joint pre-ComBat-seq VST and joint post-ComBat-seq VST -> three fixed feature spaces -> PCA -> PERMANOVA + betadisper"

if (all(c(
  "tumour_vst_pre_bc_rds",
  "tumour_vst_post_bc_rds",
  "dsmz_vst_pre_bc_collapsed_cellline_rds",
  "dsmz_vst_post_bc_collapsed_cellline_rds",
  "retained_tumour_ids_tsv"
) %in% names(inputs))) {
  tumour_plot_pre <- readRDS(inputs[["tumour_vst_pre_bc_rds"]])
  tumour_plot_post <- readRDS(inputs[["tumour_vst_post_bc_rds"]])
  dsmz_plot_pre <- readRDS(inputs[["dsmz_vst_pre_bc_collapsed_cellline_rds"]])
  dsmz_plot_post <- readRDS(inputs[["dsmz_vst_post_bc_collapsed_cellline_rds"]])
  retained_manifest <- utils::read.delim(
    inputs[["retained_tumour_ids_tsv"]], stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!"sample_id" %in% colnames(retained_manifest)) {
    stop("[ERROR] retained tumour manifest is missing sample_id: ", inputs[["retained_tumour_ids_tsv"]])
  }
  purity_passed_tumour_ids <- as.character(retained_manifest$sample_id)
  pre_tumour_ids <- colnames(tumour_plot_pre)
  post_tumour_ids <- colnames(tumour_plot_post)
  if (!identical(rownames(tumour_plot_pre), rownames(dsmz_plot_pre)) ||
      !identical(rownames(tumour_plot_post), rownames(dsmz_plot_post)) ||
      !identical(rownames(tumour_plot_pre), rownames(tumour_plot_post))) {
    stop("[ERROR] RBL plotting VST objects do not share an identical gene order")
  }
  if (!identical(pre_tumour_ids, post_tumour_ids)) {
    stop("[ERROR] identical(pre_tumour_ids, post_tumour_ids) = FALSE")
  }
  if (!identical(pre_tumour_ids, purity_passed_tumour_ids)) {
    stop("[ERROR] identical(plot_tumour_ids, purity_passed_tumour_ids) = FALSE")
  }
  if (length(pre_tumour_ids) != expected_tumours) {
    stop(sprintf("[ERROR] n_retained_tumours = %d; expected %d", length(pre_tumour_ids), expected_tumours))
  }
  if (!all(pre_tumour_ids %in% colnames(joint_vst_pre)) || !all(post_tumour_ids %in% colnames(joint_vst_post))) {
    stop("[ERROR] At least one retained tumour ID is absent from a joint VST matrix")
  }
  joint_pre_tumour_ids <- colnames(joint_vst_pre)[colnames(joint_vst_pre) %in% pre_tumour_ids]
  joint_post_tumour_ids <- colnames(joint_vst_post)[colnames(joint_vst_post) %in% post_tumour_ids]
  if (!identical(joint_pre_tumour_ids, pre_tumour_ids) || !identical(joint_post_tumour_ids, post_tumour_ids)) {
    stop("[ERROR] Retained tumour IDs are not present in the expected order within the joint VST matrices")
  }
  if (!identical(colnames(dsmz_plot_pre), colnames(dsmz_plot_post))) {
    stop("[ERROR] RBL collapsed DSMZ plotting VST objects differ in group order before and after correction")
  }
  if (ncol(dsmz_plot_pre) != expected_collapsed_dsmz || ncol(dsmz_plot_post) != expected_collapsed_dsmz) {
    stop(sprintf(
      "[ERROR] RBL plotting VST objects contain %d pre-BC and %d post-BC DSMZ groups; expected %d",
      ncol(dsmz_plot_pre), ncol(dsmz_plot_post), expected_collapsed_dsmz
    ))
  }
  plotting_pre <- cbind(tumour_plot_pre, dsmz_plot_pre)
  plotting_post <- cbind(tumour_plot_post, dsmz_plot_post)
  plot_labels <- factor(
    c(rep(tumour_label, ncol(tumour_plot_pre)), rep(dsmz_label, ncol(dsmz_plot_pre))),
    levels = c(dsmz_label, tumour_label)
  )
  plotting_summary <- "DSMZ cell-line groups"
  active_predictors <- list()
  active_predictors[[primary_predictor]] <- plot_labels
  active_tumour_ids <- pre_tumour_ids
  active_dsmz_ids <- colnames(dsmz_plot_pre)
  active_analysis_note <- paste(
    "purity-retained tumour pre/post VST matrices plus arithmetic-mean-collapsed DSMZ biological-group",
    "pre/post VST matrices -> three fixed feature spaces -> PCA -> PERMANOVA + betadisper"
  )
  if ("final_matrix_validation_path" %in% names(params) &&
      nzchar(as.character(params[["final_matrix_validation_path"]])) &&
      file.exists(as.character(params[["final_matrix_validation_path"]]))) {
    validation_report <- utils::read.delim(
      as.character(params[["final_matrix_validation_path"]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    retained_row <- validation_report[
      validation_report$check == "retained_raw_columns_equal_manifest",
      , drop = FALSE
    ]
    if (nrow(retained_row) == 1L && nzchar(retained_row$expected[[1L]])) {
      expected_from_validation <- strsplit(retained_row$expected[[1L]], ",", fixed = TRUE)[[1L]]
      if (!identical(pre_tumour_ids, expected_from_validation)) {
        stop("[ERROR] Plot tumour IDs do not match the existing final-matrix validation manifest order")
      }
    }
  }
  cat(sprintf(
    "[INFO] active analysis population for %s: %s tumours + %s collapsed DSMZ groups = %s plotted profiles\n",
    COHORT,
    figure_count(length(active_tumour_ids)),
    figure_count(length(active_dsmz_ids)),
    figure_count(ncol(plotting_pre))
  ))
}

label_sizes <- table(plot_labels)
embedding_subtitle <- sprintf(
  "%s patient tumours and %s %s",
  figure_count(sum(plot_labels == tumour_label)),
  figure_count(sum(plot_labels == dsmz_label)),
  plotting_summary
)
cat(sprintf(
  "[INFO] active analysis samples: %s; %s\n",
  figure_count(length(plot_labels)),
  paste(sprintf("%s=%d", names(label_sizes), as.integer(label_sizes)), collapse = " ")
))

## ---------------------------------------------------------------------------
## Feature spaces. Pinned annotation only; ranked and filtered on the active
## pre-correction analysis matrix, then reused unchanged for both stages.
## ---------------------------------------------------------------------------
gene_annotation <- read_gene_biotypes(inputs[["gene_annotation_gtf"]])
protein_coding <- select_protein_coding_features(
  plotting_pre, plotting_post, gene_annotation
)
write_tsv_atomic(
  data.frame(ensembl_gene_id = protein_coding$ids, stringsAsFactors = FALSE),
  outputs[["protein_coding_ids_tsv"]]
)
feature_spaces <- build_feature_spaces(
  plotting_pre, plotting_post, qc_top_genes, protein_coding$ids
)
manifest_paths <- write_manifest_block(feature_spaces, manifest_output_key)

## ---------------------------------------------------------------------------
## Active statistics and figures: always derived from the same objects the
## paired figures plot.
## ---------------------------------------------------------------------------
batch_effect <- quantify_batch_effect(
  COHORT,
  plotting_pre,
  plotting_post,
  feature_spaces$spaces,
  active_predictors,
  master_seed,
  permutations,
  tumour_count = length(active_tumour_ids),
  dsmz_count = length(active_dsmz_ids)
)
result_table <- batch_effect$table
if (include_tumour_source_analysis) {
  tumour_source_column <- as.character(params[["tumour_source_column"]])
  tumour_source_levels <- if (tumour_source_column %in% colnames(coldata)) {
    sort(unique(as.character(coldata[active_tumour_ids, tumour_source_column])))
  } else {
    character(0)
  }
  tumour_source_status <- NA_character_
  if (length(tumour_source_levels) < 2L) {
    tumour_source_status <- sprintf(
      "not applicable: %d tumour source level(s) present (%s); a one-way model needs at least two",
      length(tumour_source_levels), paste(tumour_source_levels, collapse = ", ")
    )
    cat(sprintf("[INFO] tumour_source analysis skipped -- %s\n", tumour_source_status))
  } else {
    tumour_predictors <- list(
      tumour_source = factor(coldata[active_tumour_ids, tumour_source_column])
    )
    tumour_source_effect <- quantify_batch_effect(
      COHORT,
      plotting_pre[, active_tumour_ids, drop = FALSE],
      plotting_post[, active_tumour_ids, drop = FALSE],
      feature_spaces$spaces,
      tumour_predictors,
      master_seed,
      permutations,
      tumour_count = length(active_tumour_ids),
      dsmz_count = 0L
    )
    result_table <- rbind(result_table, tumour_source_effect$table)
    tumour_source_status <- sprintf(
      "tumour-only subset: %d tumours, %d cell-line profiles excluded; levels %s",
      length(active_tumour_ids), length(active_dsmz_ids),
      paste(tumour_source_levels, collapse = ", ")
    )
  }
} else {
  tumour_source_status <- "not run in active RBL branch; all active rows are the 68-tumour + 9-DSMZ mixed-class analysis population"
  tumour_source_levels <- character(0)
}
result_table <- result_table[
  order(result_table$predictor, result_table$feature_space,
        result_table$method, result_table$stage), , drop = FALSE
]
rownames(result_table) <- NULL
write_tsv_atomic(result_table, outputs[["batch_effect_tsv"]])
verify_written_manifests(feature_spaces, manifest_paths, "active")

## ---------------------------------------------------------------------------
## Managed audit copy of the legacy joint branch. Only RBL declares these
## outputs; BRCA and NBL run the active branch alone.
## ---------------------------------------------------------------------------
legacy_outputs_requested <- all(c(
  "audit_legacy_batch_effect_tsv",
  "audit_legacy_batch_effect_provenance_tsv",
  "audit_legacy_batch_effect_run_metadata_tsv",
  "audit_legacy_protein_coding_ids_tsv",
  "audit_legacy_feature_manifest_top3000_tsv",
  "audit_legacy_feature_manifest_protein_coding_tsv",
  "audit_legacy_feature_manifest_full_expression_tsv"
) %in% names(outputs))

if (legacy_outputs_requested) {
  legacy_predictors <- resolve_predictors(coldata, sample_order, primary_predictor, NULL)
  legacy_protein_coding <- select_protein_coding_features(
    joint_vst_pre, joint_vst_post, gene_annotation
  )
  write_tsv_atomic(
    data.frame(ensembl_gene_id = legacy_protein_coding$ids, stringsAsFactors = FALSE),
    outputs[["audit_legacy_protein_coding_ids_tsv"]]
  )
  legacy_feature_spaces <- build_feature_spaces(
    joint_vst_pre, joint_vst_post, qc_top_genes, legacy_protein_coding$ids
  )
  legacy_manifest_paths <- write_manifest_block(legacy_feature_spaces, audit_manifest_output_key)
  legacy_batch_effect <- quantify_batch_effect(
    COHORT,
    joint_vst_pre,
    joint_vst_post,
    legacy_feature_spaces$spaces,
    legacy_predictors,
    master_seed,
    permutations,
    tumour_count = sum(coldata[sample_order, primary_predictor] == tumour_label),
    dsmz_count = sum(coldata[sample_order, primary_predictor] == dsmz_label)
  )
  legacy_result_table <- legacy_batch_effect$table[
    order(legacy_batch_effect$table$predictor, legacy_batch_effect$table$feature_space,
          legacy_batch_effect$table$method, legacy_batch_effect$table$stage), , drop = FALSE
  ]
  rownames(legacy_result_table) <- NULL
  write_tsv_atomic(legacy_result_table, outputs[["audit_legacy_batch_effect_tsv"]])
  verify_written_manifests(legacy_feature_spaces, legacy_manifest_paths, "legacy audit")
}

## ---------------------------------------------------------------------------
## Paired before/after PCA, one figure per feature space. PC1/PC2 only.
## ---------------------------------------------------------------------------
figure_output_key <- c(
  top3000 = "pca_before_after_top3000_pdf",
  protein_coding = "pca_before_after_protein_coding_pdf",
  full_expression = "pca_before_after_full_expression_pdf"
)
## The paired figures show the full joint set under the ComBat-seq batch factor;
## the tumour-only tumour_source model is a statistics-only companion.
for (space_name in names(batch_effect$embeddings)) {
  plot_stage_pca <- batch_effect$embeddings[[space_name]]
  primary <- batch_effect$table[
    batch_effect$table$feature_space == space_name &
      batch_effect$table$method == "permanova" &
      batch_effect$table$predictor == names(active_predictors)[1L], , drop = FALSE
  ]
  annotation <- sprintf(
    "Source/specimen-associated variation (%s)  |  PERMANOVA R2 %.3f -> %.3f",
    names(active_predictors)[1L],
    primary$r_squared[primary$stage == "before"],
    primary$r_squared[primary$stage == "after"]
  )
  plot_feature_space_pca(
    plot_stage_pca,
    plot_labels,
    outputs[[figure_output_key[[space_name]]]],
    sprintf("%s profiles before and after ComBat-seq", COHORT),
    embedding_subtitle,
    feature_spaces$spaces[[space_name]]$label,
    annotation
  )
}

## Legacy unsuffixed path, where the cohort still declares it. This is a byte
## copy of the top3000 figure, not a second rendering: copying rather than
## re-plotting is what makes "compatibility copy" verifiable by checksum and
## stops the two files drifting apart.
if ("pca_before_after_pdf" %in% names(outputs)) {
  legacy_target <- outputs[["pca_before_after_pdf"]]
  legacy_temp <- file.path(
    dirname(legacy_target),
    sprintf(".partial-%d-%s", Sys.getpid(), basename(legacy_target))
  )
  if (!file.copy(outputs[["pca_before_after_top3000_pdf"]], legacy_temp, overwrite = TRUE)) {
    stop("[ERROR] Could not stage the compatibility copy of the top3000 figure")
  }
  if (!file.rename(legacy_temp, legacy_target)) {
    stop("[ERROR] Could not promote the compatibility copy: ", legacy_target)
  }
  stopifnot(identical(
    tools::md5sum(outputs[["pca_before_after_top3000_pdf"]])[[1L]],
    tools::md5sum(legacy_target)[[1L]]
  ))
  cat("[INFO] pca_before_after_BC.pdf written as a byte-identical copy of the top3000 figure\n")
}

## ---------------------------------------------------------------------------
## Provenance.
## ---------------------------------------------------------------------------
code_paths <- list(
  shared_module = inputs[["shared_figure_module"]],
  quantification_script = inputs[["quantification_script"]],
  snakefile = inputs[["snakefile"]]
)
input_paths <- list(
  joint_vst_pre_bc = inputs[["joint_vst_pre_bc_rds"]],
  joint_vst_post_bc = inputs[["joint_vst_post_bc_rds"]],
  coldata = inputs[["coldata_tsv"]],
  retained_tumour_ids = if ("retained_tumour_ids_tsv" %in% names(inputs)) inputs[["retained_tumour_ids_tsv"]] else NA_character_
)
optional_plotting_inputs <- c(
  tumour_vst_pre_bc = "tumour_vst_pre_bc_rds",
  tumour_vst_post_bc = "tumour_vst_post_bc_rds",
  dsmz_vst_pre_bc_collapsed_cellline = "dsmz_vst_pre_bc_collapsed_cellline_rds",
  dsmz_vst_post_bc_collapsed_cellline = "dsmz_vst_post_bc_collapsed_cellline_rds"
)
for (field_name in names(optional_plotting_inputs)) {
  input_key <- optional_plotting_inputs[[field_name]]
  if (input_key %in% names(inputs)) {
    input_paths[[field_name]] <- inputs[[input_key]]
  }
}

provenance <- batch_effect_provenance_table(
  cohort = COHORT,
  repo_root = repo_root,
  snakemake_rule = snakemake@rule,
  snakemake_scriptdir = snakemake@scriptdir,
  snakemake_command = as.character(params[["snakemake_command"]]),
  snakemake_version = as.character(params[["snakemake_version"]]),
  code_paths = code_paths,
  input_paths = input_paths,
  gene_annotation_path = inputs[["gene_annotation_gtf"]],
  gene_annotation_release = as.character(params[["gene_annotation_release"]]),
  gene_annotation_source = as.character(params[["gene_annotation_source"]]),
  protein_coding_ids_path = outputs[["protein_coding_ids_tsv"]],
  feature_manifest_paths = manifest_paths,
  protein_coding_counts = protein_coding$counts,
  feature_space_counts = feature_spaces$counts,
  feature_space_names = names(feature_spaces$spaces),
  predictor_names = unique(result_table$predictor),
  group_sizes_label = paste(
    names(label_sizes), as.integer(label_sizes), sep = ":", collapse = ","
  ),
  n_samples = length(plot_labels),
  tumour_count = length(active_tumour_ids),
  dsmz_count = length(active_dsmz_ids),
  master_seed = master_seed,
  permutations = permutations,
  batch_effect_table = result_table,
  tumour_source_status = tumour_source_status,
  tumour_source_levels = tumour_source_levels,
  analysis_population_note = active_analysis_note
)
write_tsv_atomic(provenance, outputs[["batch_effect_provenance_tsv"]])

if (legacy_outputs_requested) {
  legacy_provenance <- batch_effect_provenance_table(
    cohort = COHORT,
    repo_root = repo_root,
    snakemake_rule = snakemake@rule,
    snakemake_scriptdir = snakemake@scriptdir,
    snakemake_command = as.character(params[["snakemake_command"]]),
    snakemake_version = as.character(params[["snakemake_version"]]),
    code_paths = code_paths,
    input_paths = list(
      joint_vst_pre_bc = inputs[["joint_vst_pre_bc_rds"]],
      joint_vst_post_bc = inputs[["joint_vst_post_bc_rds"]],
      coldata = inputs[["coldata_tsv"]]
    ),
    gene_annotation_path = inputs[["gene_annotation_gtf"]],
    gene_annotation_release = as.character(params[["gene_annotation_release"]]),
    gene_annotation_source = as.character(params[["gene_annotation_source"]]),
    protein_coding_ids_path = outputs[["audit_legacy_protein_coding_ids_tsv"]],
    feature_manifest_paths = legacy_manifest_paths,
    protein_coding_counts = legacy_protein_coding$counts,
    feature_space_counts = legacy_feature_spaces$counts,
    feature_space_names = names(legacy_feature_spaces$spaces),
    predictor_names = unique(legacy_result_table$predictor),
    group_sizes_label = paste(
      names(table(coldata[sample_order, primary_predictor])),
      as.integer(table(coldata[sample_order, primary_predictor])),
      sep = ":",
      collapse = ","
    ),
    n_samples = length(sample_order),
    tumour_count = sum(coldata[sample_order, primary_predictor] == tumour_label),
    dsmz_count = sum(coldata[sample_order, primary_predictor] == dsmz_label),
    master_seed = master_seed,
    permutations = permutations,
    batch_effect_table = legacy_result_table,
    tumour_source_status = "legacy audit branch only",
    tumour_source_levels = character(0),
    analysis_population_note = "legacy audit branch on uncollapsed joint pre/post VST matrices (68 tumours + 11 DSMZ libraries = 79 profiles)"
  )
  write_tsv_atomic(legacy_provenance, outputs[["audit_legacy_batch_effect_provenance_tsv"]])
}

write_tsv_atomic(
  run_metadata_table(COHORT, started_at, Sys.time()),
  outputs[["batch_effect_run_metadata_tsv"]]
)
if (legacy_outputs_requested) {
  write_tsv_atomic(
    run_metadata_table(paste0(COHORT, "_legacy_joint_79_profile"), started_at, Sys.time()),
    outputs[["audit_legacy_batch_effect_run_metadata_tsv"]]
  )
}

cat(sprintf(
  "[SUCCESS] %s batch-effect quantification complete: %d result rows across %d feature spaces and %d predictor(s) [%s]\n",
  COHORT, nrow(result_table), length(feature_spaces$spaces),
  length(unique(result_table$predictor)),
  paste(unique(result_table$predictor), collapse = ", ")
))
