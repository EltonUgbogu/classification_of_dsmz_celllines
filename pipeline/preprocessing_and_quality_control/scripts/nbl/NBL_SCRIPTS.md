# NBL Data Acquisition and QC Scripts

Scripts for downloading, processing, and validating neuroblastoma (NBL) RNA-seq data from public repositories (GEO, SRA, GDC/TARGET).

## Scripts

### SRR ID Extraction

- **`run_download_nbl_srr_ids.sh`** -- SLURM wrapper that activates the `sra3` conda environment and runs the SRR extraction script for all NBL datasets.
- **`download_nbl_tumour_sample_srr_ids.py`** -- Retrieves SRR accession lists for NBL RNA-seq datasets (GSE100148, GSE189367, SRP409177) from NCBI SRA/GEO. Supports filtering by GEO sample titles or SRA RunInfo metadata. Outputs paired-end SRR ID lists and cached RunInfo CSVs per dataset.

### FASTQ Download

- **`download_nbl_primary_pe.sh`** -- SLURM job that downloads paired-end FASTQ files from SRA using `prefetch`, `fasterq-dump`, and `pigz`. Skips already-downloaded accessions and logs failures with structured reasons.

### TARGET-NBL STAR Counts

- **`generate_target_nbl_star_manifest.py`** -- Queries the GDC API to generate a download manifest for TARGET-NBL primary-tumour STAR count quantification files.
- **`download_target_nbl_star_counts.sh`** -- SLURM job that downloads TARGET-NBL STAR count files via `gdc-client` and builds file-to-aliquot ID mappings.

### Sample ID Organisation

- **`list_nbl_sample_ids_by_cohort.R`** -- Reads a preprocessed VST-normalised expression matrix (RDS) and exports sample IDs stratified by cohort (GSE100148, GSE189367, GSE49711, SRP409177, TARGET).

### Validation and Comparison

- **`generate_GSE189367_comparison_nbl.py`** -- Compares GSE189367 sample-export SRR IDs against FASTQ files on disk, producing an intersection/difference report.
- **`generate_nbl_target_id_comparison_report.py`** -- Compares TARGET-NBL sample IDs from downstream analysis against raw aliquot IDs from metadata to identify samples filtered during preprocessing.
- **`compare_target_nbl_id_sets.py`** -- Alternative implementation of the TARGET-NBL ID set comparison above.

## Archive

The `archieve/` subdirectory contains deprecated or earlier versions of these scripts:

| File | Description |
|------|-------------|
| `download_target_nbl.sh` | Older TARGET-NBL download script using a direct `gdc-client` path |
| `extract_target_submitter_ids.py` | Extracts TARGET submitter IDs from a GDC JSON metadata cart |
| `verify_target_nbl_paired_rnaseq.py` | Verifies paired-end sequencing layout for TARGET-NBL files via the GDC API |
| `build_target_expression_matrix.py` | Combines individual STAR count files into a unified expression matrix |
| `target_nbl_file_to_sample.tsv` | Python script (misnamed) that maps GDC file UUIDs to sample IDs |

## Data Flow

```
SRR ID extraction  -->  FASTQ download  -->  ID validation / reconciliation
                                                      |
TARGET manifest generation  -->  STAR count download --+-->  Sample ID export by cohort
```
