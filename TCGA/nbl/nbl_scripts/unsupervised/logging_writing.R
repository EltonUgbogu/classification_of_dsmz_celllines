
write_dimension_report <- function(outdir, dims_list) {
  dim_report <- file.path(outdir, "dimension_report.txt")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  con <- file(dim_report, open = "wt")
  on.exit(close(con), add = TRUE)
  lapply(names(dims_list), function(nm) write_dims_line(con, nm, dims_list[[nm]]))
}

safe_rds_save <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
}

log_dims <- function(tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  cat(sprintf("[DIM] %-20s: %8s genes x %6s samples\n", tag, format(ng, big.mark=","), format(ns, big.mark=",")))
}

write_dims_line <- function(con, tag, mat) {
  ng <- if (!is.null(mat)) nrow(mat) else NA_integer_
  ns <- if (!is.null(mat)) ncol(mat) else NA_integer_
  writeLines(sprintf("%-24s\tgenes=%d\tsamples=%d", tag, ng, ns), con)
}

## form clustering.R
ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

safe_rds_save <- function(object, file) {
  ensure_dir(dirname(file))
  saveRDS(object, file = file)
}

run_pca_and_save <- function(V_adj, outdir, n_hvg = 3000, max_pc = 30) {
  ensure_dir(outdir)
  pcs <- make_pcs_matrix(V_adj, n_hvg = n_hvg, max_pc = max_pc)
  safe_rds_save(pcs, file.path(outdir, "pca_objects.rds"))
  pcs
}