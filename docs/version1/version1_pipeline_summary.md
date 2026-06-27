# Version 1 pipeline summary

Version 1 is the protected thesis workflow controlled by the repository-root `Snakefile`.

Version 1 produces the submitted thesis pipeline outputs under existing result namespaces such as:

- `results/unsupervised/`
- `results/unsupervised/pan_cancer/`
- `results/unsupervised/multicohort_cancer/`

Version 2 may read selected Version 1 outputs as protected inputs, but Version 2 must not regenerate or overwrite Version 1 results.

