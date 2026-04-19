# Method Notes

## Final selection

Clustering creates candidate tumour-cell-line neighbourhoods independently across the configured feature-selection and distance directions. Consensus then aggregates support across those analytical choices via `p_consensus`, and `summarize_p_consensus_all.R` ranks cell lines using the emitted composite score (`CPI` / `composite_score`) together with `frac_ge_thr`, `median_p`, and `max_p`.

`validation/04_model_selection_summary.R` does not invent a new rule. It reformats the existing `p_consensus_best_cell_lines_ranked.tsv` output into a shortlist that makes the final candidate models explicit, and it appends graph/community context only when that context is already present in the profile outputs.

Graph and community outputs are therefore supportive interpretation layers. They are not a separate model-selection system unless a downstream analysis explicitly uses them.
