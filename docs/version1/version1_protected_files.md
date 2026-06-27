# Version 1 protected files

The following files and folders are treated as protected Version 1 thesis workflow material during Version 2 development:

- `Snakefile`
- `config/config.yaml`
- `results/unsupervised/`
- `results/unsupervised/pan_cancer/`
- `results/unsupervised/multicohort_cancer/`

Version 2 reads selected Version 1 result files as inputs. Version 2 output rules write only to `results/version2/` and `docs/version2/`.

