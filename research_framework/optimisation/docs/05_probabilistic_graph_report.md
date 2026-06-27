# Probabilistic graph report

## Scope

The probabilistic graph step estimated cell-line graph edges from expected shared probabilistic tumour-neighbourhood evidence.

## Model

- The initial implementation focuses on latent neighbourhood events `N_ct` and graph-edge events `E_cd`.
- Expected shared neighbourhood evidence was computed as the sum over tumours of products of posterior neighbourhood probabilities for two cell lines.
- Null model: within_cancer_tumour_permutation.
- Null iterations: 1000.
- Posterior edge probability was estimated as one minus the empirical upper-tail null probability.

## Edge probability classes

- high_probability: 923 edge(s).
- low_probability: 1154 edge(s).
- moderate_probability: 17 edge(s).
- uncertain_probability: 41 edge(s).

## Outputs

- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_cellline_edges.tsv`
- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_edge_null_summary.tsv`

## Manual inspection required

- Inspect high-probability edges that are absent from deterministic thesis graphs.
- Inspect deterministic edges with low posterior edge probabilities.
- Inspect cancer types where the null mean is close to observed shared neighbourhood evidence because threshold selection may be unstable.

## Probabilistic graph resolution

Selected threshold: 0.5.
Expected false edge rate at selected threshold: 0.0176021735654763.
Threshold note: configured expected false edge-rate criterion satisfied.

## Resolution outputs

- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_resolved_edges.tsv`
- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_components.tsv`
- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_isolates.tsv`
- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_bridge_anchors.tsv`
- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_threshold_sweep.tsv`

## Manual inspection required

- Inspect thresholds where the expected false edge rate is close to the configured maximum.
- Inspect probabilistic isolates that differ from deterministic isolate status.
- Inspect bridge-like anchors with high posterior bridge scores before interpreting them biologically.
