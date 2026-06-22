#!/usr/bin/env bash
set -euo pipefail

# tidyestimate is archived on CRAN and is not packaged by conda-forge or
# bioconda. Snakemake runs this hook after creating tcga-r-env so the package
# is installed reproducibly from the archived source tarball.
Rscript --vanilla - <<'RSCRIPT'
pkg <- "tidyestimate"
expected_version <- "1.1.0"
src <- sprintf(
  "https://cran.r-project.org/src/contrib/Archive/tidyestimate/tidyestimate_%s.tar.gz",
  expected_version
)

if (requireNamespace(pkg, quietly = TRUE)) {
  installed_version <- as.character(utils::packageVersion(pkg))
  if (identical(installed_version, expected_version)) {
    message(sprintf("[OK] %s %s already installed.", pkg, installed_version))
    quit(save = "no", status = 0)
  }
}

message(sprintf("[INFO] Installing %s %s from CRAN archive.", pkg, expected_version))
utils::install.packages(src, repos = NULL, type = "source")

if (!requireNamespace(pkg, quietly = TRUE)) {
  stop(sprintf("%s installation did not produce a loadable namespace.", pkg))
}

installed_version <- as.character(utils::packageVersion(pkg))
if (!identical(installed_version, expected_version)) {
  stop(sprintf(
    "%s version mismatch: expected %s, found %s.",
    pkg, expected_version, installed_version
  ))
}

if (!requireNamespace(pkg, quietly = TRUE)) {
  stop(sprintf("%s cannot be loaded after installation.", pkg))
}

message(sprintf("[OK] %s %s installed and loads.", pkg, installed_version))
RSCRIPT
