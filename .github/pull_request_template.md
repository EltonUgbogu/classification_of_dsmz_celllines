## Summary

Describe what changed and why.

## Type of change

* [ ] Bug fix
* [ ] Refactoring with no intended scientific change
* [ ] New analysis or workflow stage
* [ ] Snakemake workflow change
* [ ] Configuration change
* [ ] Data or manifest schema change
* [ ] Documentation-only change
* [ ] Test or validation improvement

## Scientific and workflow impact

* [ ] The scientific interpretation is unchanged.
* [ ] The scientific interpretation has changed and is explained below.
* [ ] Cohort-specific effects for BRCA, NBL, RBL or HEME were considered.
* [ ] Pan-cancer downstream effects were considered.
* [ ] Random seeds and deterministic behaviour remain controlled where applicable.

### Impact description

Describe the affected analysis stages, cohorts, outputs and assumptions.

## Workflow and data contracts

* [ ] Relevant Snakemake rules were inspected.
* [ ] Rule-to-script relationships remain correct.
* [ ] Input and output paths remain compatible.
* [ ] Manifest schemas remain compatible.
* [ ] Required columns, identifiers and file formats remain compatible.
* [ ] Downstream consumers of changed outputs were inspected.
* [ ] Configuration keys and profile behaviour remain compatible.
* [ ] Not applicable — explain below.

### Contract changes

List any changed inputs, outputs, columns, identifiers, filenames or configuration keys.

## Architecture documentation

* [ ] `docs/pipeline_architecture.md` was updated.
* [ ] `reports/agent_report/pipeline_architecture_map.tsv` was updated.
* [ ] Documentation remains accurate and no update was required.
* [ ] Not applicable — this change has no architectural or data-contract effect.

### Documentation decision

Explain briefly why documentation was updated or why no update was required.

## Code Review Graph

* [ ] The Code Review Graph index was updated or rebuilt.
* [ ] CRG change-impact or minimal-context analysis was reviewed.
* [ ] CRG findings were verified against the source code.
* [ ] Snakemake, YAML, manifests and filename-mediated dependencies were checked manually.
* [ ] Not applicable — documentation-only or trivial isolated change.

CRG is a navigation aid. Absence of a graph edge is not evidence that no workflow or data dependency exists.

## Validation

List the exact checks performed.

```text
# Commands and checks
```

* [ ] Relevant syntax checks passed.
* [ ] Relevant unit or contract tests passed.
* [ ] Relevant Snakemake dry run passed.
* [ ] Expected output schemas were checked.
* [ ] No full pipeline run was necessary.
* [ ] Full or targeted pipeline execution was performed where appropriate.
* [ ] Tests could not be run — explain below.

### Validation limitations

Describe tests not run, unavailable dependencies or unresolved uncertainty.

## Security and repository hygiene

* [ ] No credentials, secrets or controlled data were added.
* [ ] No local `.code-review-graph/` database or machine-specific export was added.
* [ ] Generated outputs were not committed unintentionally.
* [ ] Shell and subprocess inputs are safely handled where relevant.
* [ ] `git diff --cached` was reviewed before committing.

## Files requiring particular reviewer attention

List the most important files and explain why.

## Remaining risks or follow-up work

List anything deliberately left unresolved.
