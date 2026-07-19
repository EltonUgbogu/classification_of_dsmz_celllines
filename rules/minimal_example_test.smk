rule test_minimal_inputs:
    """
    Validate synthetic minimal BRCA/NBL/RBL example inputs.

    This rule checks that the repository ships tiny schema-compatible
    example inputs without redistributing real patient or cell-line data.
    """
    input:
        script="scripts/validate_minimal_example_inputs.R",
        brca_rds="examples/minimal/brca/brca_vst_joint.rds",
        brca_meta="examples/minimal/brca/metadata.tsv",
        nbl_rds="examples/minimal/nbl/nbl_vst_joint.rds",
        nbl_meta="examples/minimal/nbl/metadata.tsv",
        rbl_rds="examples/minimal/rbl/rbl_vst_joint.rds",
        rbl_meta="examples/minimal/rbl/metadata.tsv"
    output:
        touch("results/tests/minimal_inputs.ok")
    shell:
        """
        Rscript {input.script} examples/minimal
        """
