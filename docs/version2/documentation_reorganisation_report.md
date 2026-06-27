# Documentation reorganisation report

## Scope

This report records the documentation movement plan for separating protected Version 1 thesis-pipeline documentation from Version 2 development documentation.

## Movement plan recorded before moving files

The existing root-level `docs/` files were inspected before reorganisation. They describe the submitted thesis pipeline, previous debugging, or Version 1 reproducibility/support-threshold work rather than the new Version 2 development workflow.

## Files moved to `docs/version1/`

- `docs/AGENT_TASK_PIPELINE_ANALYSIS.md`
- `docs/DSMZ_GRAPH_RULE_DEBUG.md`
- `docs/PIPELINE_DOCUMENTATION.md`
- `docs/PIPELINE_PROFILES_AND_OUTPUT_NAMESPACES.md`
- `docs/POST_RESOLUTION_EDGE_SUPPORT_STRATIFICATION.md`
- `docs/SUPPORT_THRESHOLD_PASS_REPORT.md`
- `docs/SUPPORT_THRESHOLD_PASS_diff_Snakefile.patch`
- `docs/SUPPORT_THRESHOLD_PASS_diff_config.patch`
- `docs/concepts_and_interpretation.md`
- `docs/method_notes.md`

These documents are treated as Version 1 preservation or thesis-pipeline documentation.

## Files moved to `docs/version2/`

No pre-existing root-level documentation was classified as Version 2 development documentation. New Version 2 documentation is written directly to `docs/version2/`.

## Files left in place

After the movement, the root `docs/` directory is intended to contain only documentation namespace folders:

- `docs/version1/`
- `docs/version2/`

## Files requiring user inspection

- `docs/version1/AGENT_TASK_PIPELINE_ANALYSIS.md`: task-oriented wording may be useful historically but should not be treated as final thesis prose.
- `docs/version1/DSMZ_GRAPH_RULE_DEBUG.md`: debugging context should remain Version 1 preservation material unless the user wants it archived elsewhere.
- `docs/version1/SUPPORT_THRESHOLD_PASS_diff_Snakefile.patch`: patch artefact retained for reproducibility history.
- `docs/version1/SUPPORT_THRESHOLD_PASS_diff_config.patch`: patch artefact retained for reproducibility history.

## Git move note

`git mv` should be used where the files are tracked by Git. If a file is untracked in the current working tree, a plain filesystem move preserves the file content but cannot preserve Git history that does not yet exist.

