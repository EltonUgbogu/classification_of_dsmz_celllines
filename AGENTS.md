# Repository agent guidance

## Code Review Graph workflow

For repository-wide reviews, debugging and change-impact analysis:

1. Confirm that `code-review-graph status` is current for `HEAD`.
2. Query minimal context before broad file searches, then query the smallest relevant file or function neighbourhood.
3. Use graph results to identify likely R, Python or shell source, callers, dependants, hubs, bridges and affected flows.
4. Treat graph output as navigation evidence and preliminary blast-radius guidance, not proof of correctness.
5. Open and inspect every returned source file before making a claim or change.
6. Inspect the root `Snakefile`, included `rules/*.smk`, active profiles, YAML configuration, manifests, rule inputs/outputs, wildcards and filename-based dependencies manually.
7. Never infer that no dependency exists merely because the graph does not show one.
8. For multi-file reviews, run change-impact analysis against a verified comparison base and corroborate it with the Git diff and Snakemake/data-contract inspection.
9. Use an incremental update after ordinary tracked R, Python or shell changes; rebuild after substantial changes, rebases or parser/configuration changes.
10. Keep cloud embeddings disabled unless the user explicitly authorises repository-text transmission to a named provider.
11. Do not commit `.code-review-graph/` databases or HTML, JSON, GraphML or SVG exports; inspect exports for absolute paths and source-derived identifiers before sharing.
12. Record graph omissions, false positives and unsupported workflow/data relationships in review reports.
13. Maintain a resumable task handoff under `reports/handsoff/` for long-running work.
14. Write completed evaluations and supporting evidence under `reports/agent_report/`, not `reports/handsoff/`.
15. Do not enable watch mode, persistent daemons, Git hooks or pull-request automation without explicit approval.
