# Ranking uncertainty report

## Scope

The ranking uncertainty step consumed available thesis ranking tables and combined rank margins with posterior neighbourhood and probabilistic graph context when available.

## Assumptions

- Rank margin was computed as the difference between the mean score for the top-ranked tumour lineage and the second-ranked tumour lineage for each cell line.
- Graph-ranking agreement was recorded as lineage agreement where probabilistic component information was available.

## Outputs

- `research_framework/optimisation/results/ranking_uncertainty/ranking_uncertainty_summary.tsv`
- `research_framework/optimisation/results/ranking_uncertainty/cellline_rank_confidence.tsv`

## Manual inspection required

- Inspect high-uncertainty cell lines before making cancer-type matching claims.
- Inspect cases where rank margin is low despite high posterior neighbourhood evidence.
