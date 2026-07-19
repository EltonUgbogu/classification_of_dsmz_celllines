# Synthetic minimal example inputs

This directory contains tiny synthetic BRCA, NBL, and RBL example inputs for testing the repository structure.

The files use fake sample identifiers and fake numeric expression values. They are not biologically meaningful and must not be used for scientific interpretation.

The RDS files mimic the object type expected by the workflow:

- R object class: `matrix` / `array`
- Rows: Ensembl-like gene IDs
- Columns: synthetic sample IDs
- Values: numeric VST-like expression values

## Files

    examples/minimal/
    ├── brca/
    │   ├── brca_vst_joint.rds
    │   └── metadata.tsv
    ├── nbl/
    │   ├── nbl_vst_joint.rds
    │   └── metadata.tsv
    ├── rbl/
    │   ├── rbl_vst_joint.rds
    │   └── metadata.tsv
    └── make_synthetic_minimal_inputs.R

## Regenerate the synthetic inputs

    Rscript examples/minimal/make_synthetic_minimal_inputs.R

## Validate the minimal example

    snakemake -n --use-conda test_minimal_inputs --config pipeline_profile=brca
    snakemake --use-conda --cores 1 test_minimal_inputs --config pipeline_profile=brca

The validation rule writes:

    results/tests/minimal_inputs.ok

This output file is generated during testing and should not be committed.
