# Probabilistic neighbourhood report

## Scope

The probabilistic neighbourhood step treated thesis-framework p_consensus values as observations and estimated posterior neighbourhood probabilities with a Beta-binomial approximation.

## Probability interpretation

`p_consensus_observed` is an observed support proportion from the thesis framework, not a posterior probability. Posterior quantities in this report are produced only after modelling representation evidence with priors and reliability weights.

## Model

- Prior: Beta(1, 1).
- Representation evidence threshold for support calls: 0.7.
- Weighted support equals the sum of representation reliability weights for representations with observed support calls.

## Outputs

- `research_framework/optimisation/results/probabilistic_neighbourhoods/posterior_neighbourhood_probabilities.tsv.gz`
- `research_framework/optimisation/results/probabilistic_neighbourhoods/neighbourhood_probability_summary.tsv`

## Manual inspection required

- Inspect pairs with wide posterior intervals before treating them as stable tumour-neighbourhood relations.
- Inspect pairs where high observed p_consensus values do not translate into high posterior probabilities because the contributing representations were downweighted.
