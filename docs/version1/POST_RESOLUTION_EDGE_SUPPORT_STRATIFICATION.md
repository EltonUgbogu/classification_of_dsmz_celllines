# Post-Resolution Edge-Support Stratification

This note records the downstream audit layer added for the focal NBL, BRCA, and
RBL graph outputs. The rule reads existing graph sidecars and does not alter the
union-supported graph, support-threshold consensus network, or resolved
neighbour graph.

## Outputs

- `results/tables/NBL_post_resolution_edge_support_stratification.tsv`
- `results/tables/BRCA_post_resolution_edge_support_stratification.tsv`
- `results/tables/RBL_post_resolution_edge_support_stratification.tsv`
- `results/tables/combined_post_resolution_edge_support_stratification.tsv`

## Results/Discussion-Ready Text

The union-supported graph is a permissive layer that maximises sensitivity to
any edge supported by the configured feature--distance representations. The
resolved-neighbour graph is a stability-based edge filter. Unlike the union
graph, it retains only edges that meet the predefined support/intersection
criteria across feature–distance representations. This reduces the risk of
overfitting to representation-specific edges, but does not prove that overfitting
is absent. Post-resolution edge-support stratification does not replace the
resolved graph. It audits edges present in the permissive union graph but absent
from the resolved graph, clarifying whether resolved isolates are completely
unsupported or retain weaker or recurrent support outside the resolved-neighbour
rule. This distinction is important for interpreting NBL, BRCA, and RBL
cell-line neighbourhoods without overstating biological conclusions from
union-only edges.

## Interpretation Boundary

The `MULTICOHORT_CANCER` graph is unsupervised with respect to cancer type:
cancer-type labels are used only after graph construction to interpret connected
components and neighbourhoods. RBL-labelled cell lines should therefore be
described as emerging as a connected group in the unsupervised multicohort
similarity space, not as a component produced by cancer-type information given
to the graph construction procedure.

RBL-only fragmentation is interpreted under the same majority-threshold,
intersection, and resolved-neighbour rules applied across cohorts. It indicates
that some RBL-only union-supported edges are not sufficiently recurrent under
the stability-based edge filter. An isolate in the resolved graph lacks a
stable resolved cell-line neighbour under that rule; it is not evidence that the
cell line lacks tumour similarity, cancer relevance, or all possible similarity
evidence.
