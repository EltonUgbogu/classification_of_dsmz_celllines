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

plot_labels <- coldata[sample_order, primary_predictor]
label_sizes <- table(plot_labels)
tumour_label <- as.character(params[["tumour_label"]])
dsmz_label <- as.character(params[["dsmz_label"]])
embedding_subtitle <- sprintf(
  "%s patient tumours and %s DSMZ cell-line profiles",
  figure_count(sum(plot_labels == tumour_label)),
  figure_count(sum(plot_labels == dsmz_label))
)

cat(sprintf(
  "[INFO] %s samples; %s\n",
  figure_count(length(sample_order)),
  paste(sprintf("%s=%d", names(label_sizes), as.integer(label_sizes)), collapse = " ")
))

## ---------------------------------------------------------------------------
## Feature spaces. Pinned annotation only; ranked and filtered on the
## pre-correction matrix, then reused unchanged for both stages.
## ---------------------------------------------------------------------------
gene_annotation <- read_gene_biotypes(inputs[["gene_annotation_gtf"]])
protein_coding <- select_protein_coding_features(
  joint_vst_pre, joint_vst_post, gene_annotation
)
write_tsv_atomic(
  data.frame(ensembl_gene_id = protein_coding$ids, stringsAsFactors = FALSE),
  outputs[["protein_coding_ids_tsv"]]
)

feature_spaces <- build_feature_spaces(
  joint_vst_pre, joint_vst_post, qc_top_genes, protein_coding$ids
)

## Ordered manifests, one per feature space, written before the statistics so a
## failed run still shows exactly which genes were about to be used.
manifest_output_key <- c(
  top3000 = "feature_manifest_top3000_tsv",
  protein_coding = "feature_manifest_protein_coding_tsv",
  full_expression = "feature_manifest_full_expression_tsv"
)
manifest_paths <- list()
for (space_name in names(feature_spaces$spaces)) {
  manifest_path <- outputs[[manifest_output_key[[space_name]]]]
  write_feature_manifest(feature_spaces$spaces[[space_name]]$ids, manifest_path)
  manifest_paths[[space_name]] <- manifest_path
  cat(sprintf(
    "[INFO] manifest written | %-20s %s genes -> %s\n",
    space_name,
    figure_count(length(feature_spaces$spaces[[space_name]]$ids)),
    basename(manifest_path)
  ))
}

## ---------------------------------------------------------------------------
## Statistics. quantify_batch_effect() re-asserts, per feature space, that the
## manifest indexes both matrices in the identical order before testing.
## ---------------------------------------------------------------------------
## Analysis 1: the ComBat-seq batch factor, over the full joint set.
predictors <- resolve_predictors(coldata, sample_order, primary_predictor, NULL)

batch_effect <- quantify_batch_effect(
  COHORT,
  joint_vst_pre,
  joint_vst_post,
  feature_spaces$spaces,
  predictors,
  master_seed,
  permutations
)
result_table <- batch_effect$table

## ---------------------------------------------------------------------------
## Analysis 2: tumour-source structure, tumours only.
##
## A separate one-way model on the tumour subset -- never a second term beside
## `dataset` in one PERMANOVA. The cell lines are dropped before the PCA is
## recomputed, so this asks only "do the tumour studies differ from each other",
## with no cell-line level in the design. Cohorts with a single tumour source
## (BRCA) are skipped and the reason is recorded in the provenance.
## ---------------------------------------------------------------------------
tumour_source_column <- as.character(params[["tumour_source_column"]])
tumour_ids <- sample_order[coldata[sample_order, primary_predictor] == tumour_label]
tumour_source_levels <- if (tumour_source_column %in% colnames(coldata)) {
  sort(unique(as.character(coldata[tumour_ids, tumour_source_column])))
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
  dropped <- length(sample_order) - length(tumour_ids)
  cat(sprintf(
    "[INFO] tumour_source analysis | %s tumours retained, %s cell-line profiles excluded; levels: %s\n",
    figure_count(length(tumour_ids)), figure_count(dropped),
    paste(tumour_source_levels, collapse = ", ")
  ))
  ## Fail closed: no cell line may survive into the tumour-only design.
  if (any(coldata[tumour_ids, primary_predictor] == dsmz_label)) {
    stop("[ERROR] A DSMZ profile survived into the tumour-only tumour_source subset")
  }
  if (dsmz_label %in% tumour_source_levels) {
    stop("[ERROR] The cell-line label appears as a tumour_source level")
  }

  tumour_predictors <- list(
    tumour_source = factor(coldata[tumour_ids, tumour_source_column])
  )
  cat(sprintf(
    "[INFO] Predictor 'tumour_source': %s\n", describe_levels(tumour_predictors$tumour_source)
  ))
  tumour_source_effect <- quantify_batch_effect(
    COHORT,
    joint_vst_pre[, tumour_ids, drop = FALSE],
    joint_vst_post[, tumour_ids, drop = FALSE],
    feature_spaces$spaces,
    tumour_predictors,
    master_seed,
    permutations
  )
  result_table <- rbind(result_table, tumour_source_effect$table)
  tumour_source_status <- sprintf(
    "tumour-only subset: %d tumours, %d cell-line profiles excluded; levels %s",
    length(tumour_ids), dropped, paste(tumour_source_levels, collapse = ", ")
  )
}

result_table <- result_table[
  order(result_table$predictor, result_table$feature_space,
        result_table$method, result_table$stage), , drop = FALSE
]
rownames(result_table) <- NULL
write_tsv_atomic(result_table, outputs[["batch_effect_tsv"]])

## Independent re-read check: the manifests on disk must still be exactly the
## row sets the statistics used.
for (space_name in names(manifest_paths)) {
  written <- utils::read.delim(manifest_paths[[space_name]], stringsAsFactors = FALSE)
  if (!identical(as.character(written$ensembl_gene_id),
                 as.character(feature_spaces$spaces[[space_name]]$ids))) {
    stop("[ERROR] Written manifest does not match the analysed feature space: ", space_name)
  }
  if (!identical(written$order_index, seq_len(nrow(written)))) {
    stop("[ERROR] Manifest order_index is not a dense ascending sequence: ", space_name)
  }
}
cat("[INFO] all feature manifests re-read and confirmed identical to the analysed spaces\n")

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
  primary <- batch_effect$table[
    batch_effect$table$feature_space == space_name &
      batch_effect$table$method == "permanova" &
      batch_effect$table$predictor == names(predictors)[1L], , drop = FALSE
  ]
  annotation <- sprintf(
    "Source/specimen-associated variation (%s)  |  PERMANOVA R2 %.3f -> %.3f",
    names(predictors)[1L],
    primary$r_squared[primary$stage == "before"],
    primary$r_squared[primary$stage == "after"]
  )
  plot_feature_space_pca(
    batch_effect$embeddings[[space_name]],
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
  coldata = inputs[["coldata_tsv"]]
)

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
  n_samples = length(sample_order),
  master_seed = master_seed,
  permutations = permutations,
  batch_effect_table = result_table,
  tumour_source_status = tumour_source_status,
  tumour_source_levels = tumour_source_levels
)
write_tsv_atomic(provenance, outputs[["batch_effect_provenance_tsv"]])

write_tsv_atomic(
  run_metadata_table(COHORT, started_at, Sys.time()),
  outputs[["batch_effect_run_metadata_tsv"]]
)

cat(sprintf(
  "[SUCCESS] %s batch-effect quantification complete: %d result rows across %d feature spaces and %d predictor(s) [%s]\n",
  COHORT, nrow(result_table), length(feature_spaces$spaces),
  length(unique(result_table$predictor)),
  paste(unique(result_table$predictor), collapse = ", ")
))
