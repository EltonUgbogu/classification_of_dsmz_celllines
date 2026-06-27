#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(optparse)})
opts <- parse_args(OptionParser(option_list=list(
  make_option("--clean-features", type="character"),
  make_option("--cell-line-rds", type="character"),
  make_option("--multicohort-rds", type="character"),
  make_option("--dsmz-outdir", type="character"),
  make_option("--multicohort-outdir", type="character"),
  make_option("--k", type="integer", default=20)
)))
required <- c("clean-features","cell-line-rds","multicohort-rds","dsmz-outdir","multicohort-outdir")
for (key in required) if (is.null(opts[[key]])) stop("Missing --", key)
dir.create(opts[["dsmz-outdir"]], recursive=TRUE, showWarnings=FALSE)
dir.create(opts[["multicohort-outdir"]], recursive=TRUE, showWarnings=FALSE)
genes <- scan(opts[["clean-features"]], what="character", quiet=TRUE)
if (length(genes) != 379L || anyDuplicated(genes)) stop("Active clean list must contain 379 unique genes")

clean_ids <- function(x) sub("\\..*$", "", x)
coverage_table <- function(ids, label) {
  clean <- clean_ids(ids)
  data.table(
    method="ranked_marker_source_pan_cancer_panel",
    matrix=label,
    clean_gene_id=genes,
    found=genes %in% clean,
    matched_matrix_row=match(genes, clean),
    matrix_duplicate_clean_id_count=sum(duplicated(clean)),
    requested_gene_count=length(genes),
    found_gene_count=sum(genes %in% clean),
    missing_gene_count=sum(!genes %in% clean)
  )
}
subset_matrix <- function(x) {
  ids <- clean_ids(rownames(x))
  idx <- match(genes, ids)
  if (anyNA(idx)) stop("Missing active features: ", paste(genes[is.na(idx)], collapse=","))
  out <- x[idx,,drop=FALSE]
  rownames(out) <- genes
  out
}
write_wide <- function(x, path) {
  dt <- as.data.table(x, keep.rownames="clean_gene_id")
  fwrite(dt, path, sep="\t")
}

cl_obj <- readRDS(opts[["cell-line-rds"]])
if (!all(c("expr","meta") %in% names(cl_obj))) stop("cell-line RDS lacks expr/meta")
cl_expr <- subset_matrix(cl_obj$expr)
if (nrow(cl_expr) != 379L) stop("Cell-line matrix does not contain the expected current feature-panel genes")
write_wide(cl_expr, file.path(opts[["dsmz-outdir"]], "ranked_marker_source_panel_dsmz_expression_matrix.tsv"))
fwrite(coverage_table(rownames(cl_obj$expr), "pan_cancer_feature_expr_cell_lines_only"),
       file.path(opts[["dsmz-outdir"]], "ranked_marker_source_panel_dsmz_expression_feature_coverage.tsv"), sep="\t")
sim <- cor(cl_expr, method="spearman", use="pairwise.complete.obs")
write_wide(sim, file.path(opts[["dsmz-outdir"]], "ranked_marker_source_panel_dsmz_similarity_matrix.tsv"))
knn <- rbindlist(lapply(colnames(sim), function(sample_id) {
  vals <- sim[,sample_id]
  vals <- vals[names(vals) != sample_id]
  ord <- order(-vals, names(vals))
  keep <- ord[seq_len(min(opts$k, length(ord)))]
  data.table(from=sample_id, to=names(vals)[keep], neighbour_rank=seq_along(keep), similarity=as.numeric(vals[keep]))
}))
fwrite(knn, file.path(opts[["dsmz-outdir"]], "ranked_marker_source_panel_dsmz_knn_edges.tsv"), sep="\t")

mc <- readRDS(opts[["multicohort-rds"]])
if (!is.matrix(mc) && !is.data.frame(mc)) stop("multicohort RDS is not a matrix")
mc_expr <- subset_matrix(as.matrix(mc))
write_wide(mc_expr, file.path(opts[["multicohort-outdir"]], "ranked_marker_source_panel_multicohort_expression_matrix.tsv"))
fwrite(coverage_table(rownames(mc), "multicohort_joint_expr_matrix"),
       file.path(opts[["multicohort-outdir"]], "ranked_marker_source_panel_multicohort_feature_coverage.tsv"), sep="\t")
cat("DSMZ matrix:", nrow(cl_expr), "x", ncol(cl_expr), "\n")
cat("Multicohort matrix:", nrow(mc_expr), "x", ncol(mc_expr), "\n")
