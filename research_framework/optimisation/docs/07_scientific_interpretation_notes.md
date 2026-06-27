# Scientific interpretation notes

## Scope statement

The optimisation framework is a research layer above the thesis framework. It does not replace the thesis framework.

## Probability interpretation

`p_consensus` is an observed support proportion from the thesis framework and must not be interpreted as a posterior probability. Posterior probabilities are generated only by the optimisation model.

## Representation diagnostics

- Representations scored for redundancy: 76.
- Representations scored for stability: 76.
- High-weight representations: 19.

## Threshold validation

Probabilistic graph thresholds were selected by expected false edge-rate control. The threshold sweep table reports the expected false edge rate for each candidate threshold.

## Probabilistic neighbourhoods

- Neighbourhood probability classes reported: high_posterior_probability, low_posterior_probability, moderate_posterior_probability, weak_posterior_probability.

## Probabilistic graph model

- High-probability graph edges: 923.
- Uncertain graph edges: 41.
The current implementation is PDG-inspired and focuses on probabilistic neighbourhood variables and graph-edge variables.

## Isolates and bridge-like anchors

- Probabilistic isolates: 0.
- Probabilistic bridge-like anchors: 16.

## Ranking uncertainty

- high_uncertainty: 1 cell line(s).
- moderate_uncertainty: 55 cell line(s).

## Questions the user must be prepared to answer

### Is p_consensus a probability?

No. It is an observed thesis-framework support proportion. Posterior probabilities are estimated only in the optimisation layer.

### Why was a probabilistic graph model not used as the main thesis framework?

The thesis framework was designed as a deterministic and consensus-based empirical baseline. The probabilistic graph model is a research extension layered above those outputs.

### How would the deterministic thesis framework be extended probabilistically?

The deterministic outputs can be treated as observations, representation reliability can be estimated, and posterior neighbourhood and graph-edge probabilities can be inferred from weighted evidence.

### How are thresholds selected without being arbitrary?

The probabilistic graph threshold is selected by controlling the expected false edge rate across retained edges.

### How does the probabilistic model handle disagreement between feature-distance representations?

Disagreement is represented through representation-specific evidence, redundancy diagnostics, reliability weights, and posterior uncertainty intervals.

### What is the difference between representation support and posterior support?

Representation support is an observed call or value from a thesis representation. Posterior support is inferred after modelling representation evidence with priors and reliability weights.

### What does a probabilistic isolate mean?

A probabilistic isolate is a cell line with no retained probabilistic graph edge at the selected threshold and an estimated probability of being disconnected under the model.

### What does a probabilistic bridge-like anchor mean?

A probabilistic bridge-like anchor is a node with high bridge score in the retained probabilistic graph topology. It is a graph uncertainty result and requires manual scientific inspection.

## Scientific limitations

- The current graph model is PDG-inspired, not a full formal probabilistic dependency graph solver.
- Representation reliability weights depend on available thesis outputs and metadata completeness.
- The null model currently uses within-cancer tumour-label permutation; other null models should be compared before publication-level claims.
- Ranking uncertainty is limited by the structure of available thesis ranking tables.
- A full PDG implementation would require explicit variable domains, conditional dependency claims, consistency objectives, and formal uncertainty propagation across all graph variables.

## Public-release portability

- The public README contains no detected local or HPC-specific absolute paths.
- The optimisation README contains no detected local or HPC-specific absolute paths.
- All config paths are repository-relative.
- The optimisation Snakefile uses repository-relative paths.
- Scripts accept paths from config or command-line arguments.
- The workflow is intended to run from a fresh clone using repository-relative paths.
- The thesis framework was not modified by this research layer.

## Machine-readable portability audit

- `research_framework/optimisation/results/input_audit/public_portability_audit.tsv`
