

#/usr/bin/env Rscript

path <- Sys.getenv(
  "JOINT_EXPR_MATRIX_RDS",
  unset = file.path(
    "results",
    "unsupervised",
    "multicohort_cancer",
    "inputs",
    "joint_expr_matrix.rds"
  )
)

readRDS(path)
