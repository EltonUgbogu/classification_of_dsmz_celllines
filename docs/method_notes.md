# Method Notes

## Final selection

Clustering creates candidate tumour-cell-line neighbourhoods independently across the configured feature-selection and distance directions. Consensus then aggregates support across those analytical choices via `p_consensus`, and `summarize_p_consensus_all.R` ranks cell lines using the emitted composite score (`CPI` / `composite_score`) together with `frac_ge_thr`, `median_p`, and `max_p`.

`validation/04_model_selection_summary.R` does not invent a new rule. It reformats the existing `p_consensus_best_cell_lines_ranked.tsv` output into a shortlist that makes the final candidate models explicit, and it appends graph/community context only when that context is already present in the profile outputs.

Graph and community outputs are therefore supportive interpretation layers. They are not a separate model-selection system unless a downstream analysis explicitly uses them.

## Graph interpretation and post-resolution edge support

The union-supported graph is a permissive graph layer that maximises sensitivity
to any edge supported by at least one configured feature--distance
representation. The support-threshold consensus network applies the same
majority-style support rule across cohorts. The resolved-neighbour graph is a
stability-based edge filter. Unlike the union graph, it retains only edges that
meet the predefined support/intersection criteria across feature–distance
representations. This reduces the risk of overfitting to
representation-specific edges, but does not prove that overfitting is absent. It
should not be changed to force isolated cell lines into a component, and it
should not be called optimal unless a formal optimisation criterion has been
defined.

Fragmentation in the RBL-only resolved graph is interpreted under the same
majority-threshold, intersection, and resolved-neighbour rules applied across the
focal cohorts. It indicates that some RBL-only union-supported edges are not
sufficiently recurrent under those rules. An isolate in the resolved graph means
that the cell line lacks a stable resolved cell-line neighbour under the current
rule. It does not mean that the cell line lacks tumour similarity, cancer
relevance, or all possible similarity evidence.

Post-resolution edge-support stratification is a downstream audit layer. It does
not replace the resolved neighbour graph. It distinguishes resolved isolates
that are completely unsupported from those that retain weaker union-supported
or recurrent support outside the resolved graph.
