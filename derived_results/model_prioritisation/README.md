# Model Prioritisation Derived Results

This directory contains derived analysis outputs for inspecting the patient-referenced model-prioritisation and cell line-centred retrieval analyses.

## Key Terms

Patient-to-cell-line ranking: Each patient tumour sample was compared against ranked cell line models using Spearman similarity in the curated 257-gene pan-cancer feature space.

Cell line-centred retrieval: Each biological cell line group was compared against clinical tumour samples to assess whether its top-ranked tumour belonged to the same lineage.

Non-focal lineage models: For each patient tumour cohort, non-focal lineage models are all ranked cell line models whose annotated cancer lineage differs from the queried tumour lineage. The median and interquartile range across these models define the internal negative-reference distribution in the ECDF plots.

Cell line models and biological cell line groups: Figure 15 uses 58 ranked cell line model profiles. Cell line-centred retrieval collapses replicate libraries into 56 biological cell line groups; RBL-15 and RBL-20 each have two RNA-seq profiles collapsed by median per-tumour Spearman similarity before retrieval ranking.

## Files

- `patient_to_cellline_rankings.tsv`: all patient-to-cell-line similarity scores and ranks for 1,128 patient tumour samples against 58 cell line model profiles.
- `ecdf_input_table.tsv`: the per-sample, per-model rank table used to draw the Figure 15 ECDFs.
- `ecdf_non_focal_summary.tsv`: median, lower-quartile, and upper-quartile ECDF values for non-focal lineage models at each rank threshold.
- `cell_line_centred_retrieval_results.tsv`: full cell line-centred tumour retrieval ranks for 56 biological cell line groups.
- `top50_tumour_neighbourhoods_per_cell_line.tsv`: top-50 tumour-neighbourhood assignments per biological cell line group with consensus component annotations.
- `cell_line_group_metadata.tsv`: mapping between original RNA-seq profile IDs, canonical biological cell line groups, replicate status, and inclusion flags.
- `feature_set_257_genes.tsv`: exported pan-cancer feature set annotations.
- `validation_summary.tsv`: compact validation summary produced by the export script.

## Data Availability

Raw controlled-access sequencing data are not redistributed. These TSV files are derived analysis outputs intended to support inspection and reproducibility of ranking and retrieval analyses.
