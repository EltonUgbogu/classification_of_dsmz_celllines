# Optimisation framework overview

## Scope

The optimisation framework is a research layer above the thesis framework. It consumes thesis-framework outputs and estimates representation diagnostics, reliability weights, posterior neighbourhood probabilities, probabilistic graph edges, probabilistic graph-resolution outputs, and ranking uncertainty.

## Separation from the thesis framework

The thesis framework remains the empirical baseline. It is not rewritten, replaced, or reinterpreted as a probabilistic model. The optimisation framework treats thesis outputs as observations.

## Probability interpretation

Existing `p_consensus` values are observed values from the thesis framework. They are not posterior probabilities and are not statistical p-values. Posterior quantities are created only after the optimisation framework combines representation evidence with priors and reliability weights.

## PDG-inspired layer

The graph model in this directory is PDG-inspired. It represents conditional probabilistic claims over stable tumour-neighbourhood events and cell-line graph edges, and it records disagreement between evidence sources. It is not a full exact probabilistic dependency graph solver.

## Manual inspection

The user must inspect input completeness, downweighted representations, threshold sweeps, probabilistic edges that disagree with deterministic graph edges, probabilistic isolates, bridge-like anchors, and ranking decisions classified as uncertain.

