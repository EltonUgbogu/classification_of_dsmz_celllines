# Research optimisation framework

This directory contains a research optimisation layer above the thesis framework. The thesis framework remains the empirical baseline and is not reinterpreted as a probabilistic model.

The optimisation layer treats thesis outputs as observations. Existing `p_consensus` values are not posterior probabilities. The probabilistic layer estimates posterior neighbourhood probabilities from representation evidence, reliability weights, and explicit priors.

The PDG-inspired graph model integrates uncertain and potentially conflicting evidence sources. It is a tractable research implementation inspired by probabilistic dependency graphs, not a full formal PDG solver.

Run all commands from the repository root.

## Public clone workflow

```bash
git clone <repository-url>
cd classification_of_dsmz_celllines
git checkout v2.0-optimisation
```

## Dry run

```bash
snakemake -n \
  -s research_framework/optimisation/workflow/Snakefile.optimisation \
  --configfile research_framework/optimisation/config/optimisation.yaml
```

## Execution

```bash
snakemake -j 8 \
  -s research_framework/optimisation/workflow/Snakefile.optimisation \
  --configfile research_framework/optimisation/config/optimisation.yaml
```

## Output layout

The workflow writes tables and reports under:

- `research_framework/optimisation/results/input_audit/`
- `research_framework/optimisation/results/representation_diagnostics/`
- `research_framework/optimisation/results/representation_weights/`
- `research_framework/optimisation/results/probabilistic_neighbourhoods/`
- `research_framework/optimisation/results/probabilistic_graphs/`
- `research_framework/optimisation/results/ranking_uncertainty/`
- `research_framework/optimisation/results/logs/`
- `research_framework/optimisation/docs/`

## Reports to read first

Read these reports before interpreting probabilistic outputs:

1. `research_framework/optimisation/docs/01_input_audit_report.md`
2. `research_framework/optimisation/docs/02_representation_diagnostics_report.md`
3. `research_framework/optimisation/docs/03_threshold_validation_report.md`
4. `research_framework/optimisation/docs/04_probabilistic_neighbourhood_report.md`
5. `research_framework/optimisation/docs/05_probabilistic_graph_report.md`
6. `research_framework/optimisation/docs/07_scientific_interpretation_notes.md`

## Interpretation guardrails

- The thesis framework produces deterministic and consensus-based empirical outputs.
- The optimisation framework evaluates stability, redundancy, threshold behaviour, probabilistic neighbourhoods, probabilistic graph edges, and ranking uncertainty.
- `p_consensus` is an observed thesis-framework value and must not be called a posterior probability.
- Posterior quantities are generated only by the optimisation model.
- Disagreement between representations is an output to inspect, not a nuisance to hide.
