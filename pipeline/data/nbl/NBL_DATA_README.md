# Neuroblastoma (NBL) RNA-seq Data — Documentation

**Location:** `/work/ugbogu/pipeline/data/nbl/`  
**Last updated:** 2026-04-08  
**Maintainer:** Elton Ugbogu

---

## Overview

This directory holds all RNA-seq input data and intermediate outputs for the
neuroblastoma cohort. It contains four independent data sources:

| Source | Cohort ID | Repository | Samples | Read length | Platform |
|--------|-----------|------------|---------|-------------|----------|
| GEO | GSE100148 | SRA / NCBI | 15 | 100 bp PE | Illumina HiSeq 2500 |
| GEO | GSE189367 | SRA / NCBI | 57 | ~290 bp PE | Illumina NovaSeq 6000 |
| GEO | SRP409177 | SRA / NCBI | 47 | 202 bp PE | Illumina HiSeq 4000 |
| GDC | TARGET-NBL | NCI GDC | 162 | — | Illumina (various) |

**Total samples across all cohorts: 281**

All samples are primary tumour, paired-end RNA-seq from *Homo sapiens* (taxid 9606).

---

## Directory Structure

```
data/nbl/
├── GSE100148/               # GEO cohort — ATRX wild-type NBL (HiSeq 2500)
│   ├── fastq/               # Paired-end FASTQ files (*_1.fastq.gz, *_2.fastq.gz)
│   ├── sra/                 # Intermediate SRA archives (post-download, pre-dump)
│   ├── metadata/            # SRA run tables and GSM mappings
│   ├── logs/                # Download and prefetch logs
│   ├── tmp/                 # Temporary download staging
│   └── removed_nbl_srr_id_primary_tumour_pe.txt  # Excluded SRR IDs
│
├── GSE189367/               # GEO cohort — NBL primary tumours (NovaSeq 6000)
│   ├── fastq/
│   ├── sra/
│   ├── metadata/
│   ├── logs/
│   ├── tmp/
│   ├── removed_nbl_srr_id_primary_tumour_pe.txt
│   └── retry_gse189367_missing_srrs.sh  # Script to re-download failed SRRs
│
├── SRP409177/               # SRA study — NBL (HiSeq 4000, CHLA)
│   ├── fastq/
│   ├── sra/
│   ├── metadata/
│   ├── logs/
│   ├── tmp/
│   └── removed_nbl_srr_id_primary_tumour_pe.txt
│
├── target_nbl/              # GDC TARGET-NBL — STAR augmented counts
│   ├── count_data/          # 162 UUID directories, one per sample
│   │   └── {file_uuid}/
│   │       └── *.rna_seq.augmented_star_gene_counts.tsv
│   ├── metadata/
│   │   ├── target_nbl_file_to_aliquot.tsv   # UUID → aliquot_submitter_id mapping
│   │   ├── target_nbl_star_counts_mapping.tsv  # Full GDC metadata per file
│   │   ├── target_nbl_aliquot_ids.txt
│   │   └── target_nbl_aliquot_ids_unique.txt
│   ├── manifests/
│   │   ├── gdc_manifest_target_nbl_star_counts.txt
│   │   └── payload_counts_manifest.json
│   ├── manifest.tsv         # GDC download manifest (id, filename, md5, size)
│   ├── logs/
│   └── merged/              # Output of merge_target_nbl_counts.R (generated)
│       ├── target_nbl_count.rds
│       ├── target_nbl_sample_metadata.csv
│       └── target_nbl_gene_list.txt
│
├── count_data/              # Combined count matrix output (generated)
│   ├── nbl_tumour_count.rds              # genes × 274 samples (all cohorts)
│   ├── nbl_tumour_sample_metadata.csv    # cohort, sample ID, strandedness
│   └── nbl_tumour_gene_list_hgnc.txt     # gene universe (Ensembl IDs, no version)
│
└── preprocessing_results/   # Downstream processed data
    ├── vst_NBL_joint_batch_corrected.rds  # VST-normalised, batch-corrected matrix
    ├── tumour_sample_ids.txt
    └── sample_id_exports/
        ├── GSE100148_sample_ids.tsv
        ├── GSE189367_sample_ids.tsv
        ├── SRP409177_sample_ids.tsv
        ├── TARGET_sample_ids.tsv
        ├── all_sample_ids_by_cohort.tsv
        └── sample_id_counts_by_cohort.tsv
```

---

## GEO Cohort Details

### GSE100148 — ATRX wild-type neuroblastoma

| Field | Value |
|-------|-------|
| GEO accession | GSE100148 |
| SRA study | SRP109627 |
| BioProject | PRJNA390893 |
| Samples (final) | 15 paired-end |
| Read length | 100 bp |
| Platform | Illumina HiSeq 2500 |
| Library layout | PAIRED |
| Tissue | Primary tumour |
| Sample IDs | SRR78285xx series |

**Metadata files:**
- `SraRunInfo.csv` — full SRA run table for all runs in the study
- `SraRunInfo_primary.csv` — filtered to primary tumour samples
- `SraRunInfo_primary_PE.csv` — further filtered to paired-end only (used for download)
- `gsm_titles.tsv` — GSM accession → sample title mapping
- `primary_gsms.tsv` — GSM accessions for primary tumour samples
- `removed_nbl_srr_id_primary_tumour_pe.txt` — SRR IDs explicitly excluded from analysis

---

### GSE189367 — NBL primary tumours (NovaSeq)

| Field | Value |
|-------|-------|
| GEO accession | GSE189367 |
| SRA study | SRP347311 |
| BioProject | PRJNA782665 |
| Samples (final) | 57 paired-end |
| Read length | ~290 bp |
| Platform | Illumina NovaSeq 6000 |
| Library layout | PAIRED |
| Tissue | Primary tumour |
| Sample IDs | SRR17010xxx / SRR17011xxx series |

**Metadata files:**
- `SraRunInfo.csv`, `SraRunInfo_primary.csv`, `SraRunInfo_primary_PE.csv`
- `removed_nbl_srr_id_primary_tumour_pe.txt` — excluded SRR IDs
- `retry_gse189367_missing_srrs.sh` — re-download script for failed SRRs

**Note:** Some SRR downloads required retry. See `logs/retry_fail_*.txt` and
`logs/retry_success_*.txt` for the retry history. `SRR17010986` previously had
a missing R2; check `removed_nbl_srr_id_primary_tumour_pe.txt` for current status.

---

### SRP409177 — NBL (HiSeq 4000, CHLA)

| Field | Value |
|-------|-------|
| SRA study | SRP409177 |
| BioProject | PRJNA904244 |
| Submitting centre | CHLA (Children's Hospital Los Angeles) |
| Samples (final) | 47 paired-end |
| Read length | 202 bp |
| Platform | Illumina HiSeq 4000 |
| Library layout | PAIRED |
| Tissue | Primary tumour |
| Sample IDs | SRR22373xxx series |

**Metadata files:**
- `SraRunInfo.csv`, `SraRunInfo_primary.csv`, `SraRunInfo_primary_PE.csv`
- `removed_nbl_srr_id_primary_tumour_pe.txt` — excluded SRR IDs

---

## TARGET-NBL (GDC)

| Field | Value |
|-------|-------|
| GDC project | TARGET-NBL |
| Samples | 162 |
| Data type | Gene Expression Quantification |
| Workflow | STAR - Counts |
| Annotation | GENCODE v36 |
| File format | Augmented STAR gene count TSV |
| Sample identifiers | `aliquot_submitter_id` (e.g. `TARGET-30-PARACS-01A-01R`) |
| Directory naming | GDC `file_id` (UUID) |
| Sample type | Primary Tumor |
| Tissue type | Tumor |

**Count file format** (`*.rna_seq.augmented_star_gene_counts.tsv`):

```
# gene-model: GENCODE v36
gene_id    gene_name    gene_type    unstranded    stranded_first    stranded_second    tpm_unstranded    fpkm_unstranded    fpkm_uq_unstranded
N_unmapped ...
N_multimapping ...
ENSG00000000003.15    TSPAN6    protein_coding    1373    679    694    ...
```

**Key mapping files:**
- `metadata/target_nbl_file_to_aliquot.tsv` — maps `file_id` (UUID directory name)
  to `aliquot_submitter_id`. Required by `merge_target_nbl_counts.R`.
- `metadata/target_nbl_star_counts_mapping.tsv` — full GDC metadata per file
  including case ID, sample type, submitter IDs.
- `manifest.tsv` — GDC download manifest with md5 checksums for integrity verification.

---

## Sample Filtering

Each GEO cohort has a `removed_nbl_srr_id_primary_tumour_pe.txt` file listing SRR IDs
that were explicitly excluded. Reasons include:
- Non-primary tumour samples
- Single-end libraries (pipeline requires paired-end)
- Failed or incomplete downloads

The Snakefile discovers valid samples by scanning `fastq/` directories for
`*_1.fastq.gz` / `*_2.fastq.gz` pairs at runtime, so excluded samples never
enter the pipeline regardless of these files.

---

## Generated Outputs

### Per-sample (Snakefile outputs)

| File | Location |
|------|----------|
| Aligned BAM | `{cohort}/results/star/{sample}.Aligned.sortedByCoord.out.bam` |
| BAM index | `{cohort}/results/star/{sample}.Aligned.sortedByCoord.out.bam.bai` |
| STAR log | `{cohort}/results/star/{sample}.Log.final.out` |
| Gene count table | `{cohort}/results/star/valid_count/{sample}.tab` |
| Picard metrics | `{cohort}/results/picard/{sample}_*.txt` |
| FastQC report | `{cohort}/results/fastqc/{sample_read}_fastqc.html` |
| MultiQC report | `{cohort}/results/multiqc/multiqc_report.html` |

The `.tab` files are 4-column STAR `ReadsPerGene.out.tab` format:
`gene_id | unstranded | stranded_first | stranded_second`.
The appropriate strandedness column is detected automatically per sample.

### Combined count matrix

Built by `merge_star_counts.R` after all cohorts are aligned:

| File | Description |
|------|-------------|
| `count_data/nbl_tumour_count.rds` | genes × 274 samples integer count matrix (R matrix) |
| `count_data/nbl_tumour_sample_metadata.csv` | cohort, sample ID, strandedness column used |
| `count_data/nbl_tumour_gene_list_hgnc.txt` | Ensembl gene IDs (version suffix stripped) |

**Gene universe:** intersection of GENCODE v44 (GEO cohorts) and GENCODE v36 (TARGET-NBL).
Ensembl version suffixes are stripped (`ENSG00000000003.15` → `ENSG00000000003`).

**Sample name conventions:**

| Cohort | Column name format | Example |
|--------|--------------------|---------|
| GSE100148 | SRR accession | `SRR7828575` |
| GSE189367 | SRR accession | `SRR17011027` |
| SRP409177 | SRR accession | `SRR22373272` |
| TARGET-NBL | aliquot_submitter_id | `TARGET-30-PARACS-01A-01R` |

### Preprocessed data

| File | Description |
|------|-------------|
| `preprocessing_results/vst_NBL_joint_batch_corrected.rds` | VST-normalised, batch-corrected expression matrix across all cohorts |
| `preprocessing_results/tumour_sample_ids.txt` | Final sample ID list after QC |
| `preprocessing_results/sample_id_exports/*.tsv` | Per-cohort sample ID exports |

---

## Reference Genome

GEO cohorts (GSE100148, GSE189367, SRP409177) are aligned against:

| Resource | Path | Version |
|----------|------|---------|
| STAR index | `/work/ugbogu/pipeline/data/reference/star_gencode_v44` | GENCODE v44 |
| GTF | `/work/ugbogu/pipeline/data/reference/gencode_v44/gencode.v44.basic.annotation.gtf` | GENCODE v44 |
| FASTA | `/work/ugbogu/pipeline/data/reference/gencode_v44/GRCh38.primary_assembly.genome.fa` | GRCh38 |

TARGET-NBL data was pre-aligned by GDC using GENCODE v36.

---

## Alignment Pipeline

GEO cohorts are processed by:
`pipeline/preprocessing_and_quality_control/Snakefile`

Steps per sample:
1. **FastQC** — per-read quality control (R1 and R2 separately)
2. **STAR 2-pass** — splice-aware alignment, coordinate-sorted BAM + gene counts
3. **samtools index** — BAM indexing
4. **Picard ValidateSamFile** — BAM integrity check
5. **Picard CollectAlignmentSummaryMetrics** — mapping statistics
6. **Picard CollectQualityYieldMetrics** — base quality summary
7. **Picard CollectInsertSizeMetrics** — insert size distribution
8. **MultiQC** — per-project aggregated QC report

Count matrix assembly:
- `merge_target_nbl_counts.R` — merges TARGET-NBL TSV files → `target_nbl/merged/`
- `merge_star_counts.R` — merges GEO `.tab` files + TARGET-NBL RDS → `count_data/`

---

## Data Download Scripts

| Script | Location | Purpose |
|--------|----------|---------|
| `download_nbl_primary_pe.sh` | `scripts/nbl/` | Download GEO cohort FASTQs via prefetch + fasterq-dump |
| `download_target_nbl_star_counts.sh` | `scripts/nbl/` | Download TARGET-NBL counts from GDC |
| `retry_gse189367_missing_srrs.sh` | `GSE189367/` | Retry failed GSE189367 downloads |

---

## Tumour Purity Analysis

**Script:** `preprocessing_and_quality_control/tumour_purity_analysis/nbl/nbl_tumour_purity_analysis.R`  
**Runs:** standalone, after the Snakefile has produced the count matrix

### Purpose

Estimates tumour microenvironment signals (stromal score, immune score, ESTIMATE
score) from bulk RNA-seq expression data using the ESTIMATE method, then filters
samples below a purity threshold to remove those dominated by non-tumour signal.
Low-purity samples can distort downstream analyses such as clustering,
differential expression, and cohort comparison.

### Inputs

| File | Source |
|------|--------|
| `count_data/nbl_tumour_count.rds` | Produced by `merge_star_counts.R` (Snakefile) |
| `count_data/nbl_tumour_sample_metadata.csv` | Produced by `merge_star_counts.R` (Snakefile) |
| `count_data/nbl_ensembl_to_hgnc.tsv` | Ensembl → HGNC symbol mapping table |

### Outputs (written to `results/tumour_purity_analysis/nbl/`)

| File | Description |
|------|-------------|
| `estimate_scores.csv` | Per-sample stromal, immune, and ESTIMATE scores + derived purity |
| `filtered_hgnc_counts_purity0.7.rds` | Filtered count matrix in HGNC gene symbol space |
| `nbl_filtered_count_matrix_ensembl_purity0.7.rds` | Same filtered samples in Ensembl ID space |
| `sample_retention_by_project_same_axis.pdf` | Bar plot of sample counts before/after filtering by cohort |
| `purity_diagnostics.pdf` | Scatter plots (purity vs stromal/immune score) and purity distribution histogram |
| `nbl_purity_pipeline_<timestamp>.log` | Full run log |

### Purity threshold

The threshold used for this project is **0.7**. Samples with estimated tumour
purity below 0.7 are excluded from the filtered matrix.

> **Note for future research:** The threshold is a configurable parameter.
> To explore the effect of stricter or more lenient filtering, pass a different
> value via the `PURITY_THRESHOLD` environment variable before running the script.
> The `purity_diagnostics.pdf` shows the full purity distribution and can guide
> threshold selection for other cohorts or research questions.

### Sample-to-cohort traceability

Project assignment for every sample is resolved by direct lookup in
`nbl_tumour_sample_metadata.csv`, not by parsing sample name conventions.
This ensures SRR accessions (GEO cohorts) and TARGET aliquot IDs are correctly
attributed to their source cohort in all plots and reports.

### Invocation

```bash
SE_PATH=/work/ugbogu/pipeline/data/nbl/count_data/nbl_tumour_count.rds \
META_CSV=/work/ugbogu/pipeline/data/nbl/count_data/nbl_tumour_sample_metadata.csv \
MAP_TSV=/work/ugbogu/pipeline/data/nbl/count_data/nbl_ensembl_to_hgnc.tsv \
OUTPUT_DIR=/work/ugbogu/pipeline/results/tumour_purity_analysis/nbl \
Rscript preprocessing_and_quality_control/tumour_purity_analysis/nbl/nbl_tumour_purity_analysis.R
```
