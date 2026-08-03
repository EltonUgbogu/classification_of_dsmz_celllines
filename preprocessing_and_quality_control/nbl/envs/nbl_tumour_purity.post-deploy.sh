#!/usr/bin/env bash
set -euo pipefail

# tidyestimate is not packaged in the configured Conda channels. Install the
# exact audited upstream revision into the newly deployed Snakemake environment.
"${CONDA_PREFIX}/bin/Rscript" -e \
  'remotes::install_github("KaiAragaki/tidyestimate@10455270b97eabeadbe707009a0f5dc9f440c0e0", upgrade = "never", dependencies = FALSE)'

"${CONDA_PREFIX}/bin/Rscript" -e \
  'stopifnot(as.character(utils::packageVersion("tidyestimate")) == "1.1.1")'
