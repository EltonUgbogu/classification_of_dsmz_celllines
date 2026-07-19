#!/usr/bin/env Rscript

set.seed(20260719)

make_profile <- function(profile, outdir, n_genes = 40, n_cell_lines = 3, n_tumours = 4) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  genes <- sprintf("ENSG%011d", seq_len(n_genes))

  cell_lines <- sprintf("%s_CELL_LINE_SYN_%03d", toupper(profile), seq_len(n_cell_lines))
  tumours <- sprintf("%s_TUMOUR_SYN_%03d", toupper(profile), seq_len(n_tumours))
  samples <- c(cell_lines, tumours)

  expr <- matrix(
    rnorm(n_genes * length(samples), mean = 8, sd = 1.5),
    nrow = n_genes,
    ncol = length(samples),
    dimnames = list(genes, samples)
  )

  metadata <- data.frame(
    sample_id = samples,
    profile = profile,
    cancer_type = toupper(profile),
    sample_type = c(rep("cell_line", n_cell_lines), rep("tumour", n_tumours)),
    source = "synthetic_minimal",
    stringsAsFactors = FALSE
  )

  saveRDS(expr, file.path(outdir, paste0(profile, "_vst_joint.rds")))

  write.table(
    metadata,
    file.path(outdir, "metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

base <- "examples/minimal"

make_profile("brca", file.path(base, "brca"))
make_profile("nbl",  file.path(base, "nbl"))
make_profile("rbl",  file.path(base, "rbl"))

cat("Synthetic minimal inputs written to examples/minimal\n")
