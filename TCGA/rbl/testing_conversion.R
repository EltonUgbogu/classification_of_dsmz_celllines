#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  # No heavy deps needed; keep it minimal
})

RDS_PATH <- "/home/chu25/data/rbl/rbl_raw_count_matrix.rds"

cat("[INFO] Reading RDS from:", RDS_PATH, "\n")
mat_test <- readRDS(RDS_PATH)

if (!is.matrix(mat_test)) {
  cat("[WARN] Object is not a matrix. Class:\n")
  print(class(mat_test))
} else {
  cat("[INFO] Loaded a matrix.\n")
}

cat("[CHECK] Matrix dims: ", paste(dim(mat_test), collapse = " × "), "\n")

rn <- rownames(mat_test)
cn <- colnames(mat_test)

cat("[CHECK] Rowname class: "); print(class(rn))
cat("[CHECK] Colname class: "); print(class(cn))

cat("[CHECK] First 5 rownames:\n")
print(head(rn, 5))

cat("[CHECK] First 5 colnames:\n")
print(head(cn, 5))

if (!is.null(rn) && !is.null(cn) &&
    is.character(rn) && is.character(cn)) {
  message("[CHECK SUCCESS] Row/Col names are non-NULL plain character vectors")
} else {
  warning("[CHECK FAIL] Row/Col names are NULL or not plain characters")
}
