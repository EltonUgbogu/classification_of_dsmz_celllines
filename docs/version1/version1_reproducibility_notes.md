# Version 1 reproducibility notes

The protected Version 1 thesis workflow remains controlled by `Snakefile` and `config/config.yaml`.

The Version 2 development workflow uses `Snakefile.v2` and `config/config_v2.yaml`. Version 2 dry-runs or runs should not be used to regenerate Version 1 outputs.

If a Version 1 input required by Version 2 is missing, run the corresponding protected Version 1 target explicitly through `Snakefile` after inspecting the dry-run.

