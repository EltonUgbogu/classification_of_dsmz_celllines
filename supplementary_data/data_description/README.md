# Data description

This folder documents processed analysis inputs, retained sample identifiers, exclusions, and processing notes for the active pan-cancer analysis.

At the time this documentation was prepared, the raw/pre-filter patient tumour data used upstream of tumour-purity estimation were not available in the working project directory. Therefore, the patient tumour component of the reported analysis is documented from the available post-threshold analysis set.

The tumour samples listed here are samples that had already passed the tumour-purity threshold used for the reported analysis. They should therefore be interpreted as the tumour-purity-passed analysis cohort, not as the complete raw/pre-filter patient tumour cohort originally considered before tumour-purity filtering.

The tumour sample identifiers used in the reported analysis are provided in:

- `tumour_purity_passed_sample_ids.tsv`

Cell-line sample/profile identifiers used by each disease analysis are provided in:

- `cell_line_sample_ids_by_cohort.tsv`

DSMZ cell-line identifiers and metadata are provided separately. The full DSMZ metadata table contains 167 cell-line entries and is documented in:

- `dsmz_167_cell_line_ids.tsv`
- `dsmz_167_cell_line_metadata.tsv`

Where retinoblastoma-specific DSMZ cell-line profiles are used, the relevant subset is documented in:

- `dsmz_rbl_cell_line_subset.tsv`

Ongoing raw-data reconstruction workflows for RBL and NBL are tracked separately and should not be interpreted as the input dataset for the current reported analysis unless the reported figures and tables were regenerated from those outputs.

## File status vocabulary

- `used_in_reported_analysis`: processed object or metadata used to document the reported analysis.
- `not_available_at_time_of_documentation`: raw/pre-filter or excluded material not present in the working project directory.
- `under_reconstruction_for_future_rerun`: workflow output being rebuilt for possible future reruns, not a current reported-analysis input.
- `not_used_in_reported_analysis`: available object not used for the reported analysis.
- `legacy_reference_only`: older reference object retained for traceability, not a current analysis input.

## Files

- `analysis_input_files.tsv`: processed analysis inputs confirmed or flagged by the audit.
- `tumour_purity_passed_sample_ids.tsv`: tumour-purity-passed tumour sample/profile IDs by cohort.
- `cell_line_sample_ids_by_cohort.tsv`: cell-line sample/profile IDs by cohort and source object.
- `dsmz_167_cell_line_ids.tsv`: selected DSMZ identifier columns for all 167 metadata rows.
- `dsmz_167_cell_line_metadata.tsv`: full DSMZ metadata table.
- `dsmz_rbl_cell_line_subset.tsv`: retinoblastoma-relevant DSMZ rows.
- `excluded_or_not_available_sources.tsv`: unavailable or not-used sources and reasons.
- `processing_summary.md`: concise processing and availability summary.
- `reproducibility_checksums.tsv`: SHA256 checksums for curated data-description files in this directory.
