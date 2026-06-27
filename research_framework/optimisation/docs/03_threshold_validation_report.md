# Threshold validation report

## Scope

The threshold validation step evaluated posterior edge-probability thresholds for the probabilistic graph.

## Threshold selection rule

For each candidate threshold, the expected false edge rate was computed as the sum of one minus posterior edge probability across retained edges, divided by the number of retained edges.

Configured maximum expected false edge rate: 0.05.
Selected threshold: 0.5.
Expected false edge rate at selected threshold: 0.0176021735654763.
Selection note: configured expected false edge-rate criterion satisfied.

## Output

- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_threshold_sweep.tsv`

## Manual inspection required

- Inspect candidate thresholds with very few retained edges.
- Inspect thresholds near the configured expected false edge-rate limit.
- Inspect whether the selected threshold removes deterministic thesis edges that remain scientifically plausible.
