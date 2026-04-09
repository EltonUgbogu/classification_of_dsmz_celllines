# Neuroblastoma (NBL) RNA-seq Data — Documentation

**Location:** `/work/ugbogu/pipeline/data/nbl/`  
**Last updated:** 2026-04-09  
**Maintainer:** Elton Ugbogu

---

## Overview

This directory holds all RNA-seq input data for the neuroblastoma cohort pipeline.
Four independent data sources have been downloaded and are ready for alignment:

| Source | Cohort ID | Repository | FASTQs on disk | Read length | Platform |
|--------|-----------|------------|----------------|-------------|----------|
| GEO | GSE100148 | SRA / NCBI | 15 | 100 bp PE | Illumina HiSeq 2500 |
| GEO | GSE189367 | SRA / NCBI | 64 | ~290 bp PE | Illumina NovaSeq 6000 |
| GEO | SRP409177 | SRA / NCBI | 47 | 202 bp PE | Illumina HiSeq 4000 |
| GDC | TARGET-NBL | NCI GDC | 162 | — | Illumina (various) |

**Total samples available for alignment: 288**

All samples are primary tumour, paired-end RNA-seq from *Homo sapiens* (taxid 9606).

> **Pipeline status:** The alignment Snakefile has not yet been run on the current
> FASTQ data. `count_data/` does not exist. `preprocessing_results/` contains
> outputs from a **prior analysis run** using different sample sets and should not
> be treated as reflecting the current data. See the
> [Preprocessing Results](#preprocessing-results-prior-run) section for details.

---

## Directory Structure

```
data/nbl/
├── GSE100148/               # GEO cohort — ATRX wild-type NBL (HiSeq 2500)
│   ├── all_primary_tumour_pe.txt    # All primary-tumour PE SRR IDs (pre-exclusion)
│   ├── all_primary_tumour_srr.txt   # All primary-tumour SRR IDs (pre-exclusion)
│   ├── fastq/               # 15 paired-end samples (*_1.fastq.gz, *_2.fastq.gz)
│   ├── logs/
│   │   ├── GSE100148_88684.err      # SLURM job stderr (job 88684)
│   │   ├── GSE100148_88684.out      # SLURM job stdout (job 88684)
│   │   ├── fasterq_dump/
│   │   ├── pigz/
│   │   └── prefetch/
│   ├── metadata/            # SRA run tables and GSM mappings
│   ├── removed_nbl_srr_id_primary_tumour_pe.txt  # Excluded SRR IDs (empty — none excluded)
│   ├── sra/                 # Intermediate SRA archives (post-download, pre-dump)
│   └── tmp/
│       └── 88684/           # SLURM job staging directory
│
├── GSE189367/               # GEO cohort — NBL primary tumours (NovaSeq 6000)
│   ├── all_primary_tumour_pe.txt
│   ├── all_primary_tumour_srr.txt
│   ├── fastq/               # 64 paired-end samples
│   ├── gse189367_ids_in_sample_and_fastq.txt   # Diagnostic: SRR IDs present in both sets
│   ├── gse189367_ids_only_in_fastq.txt         # Diagnostic: SRRs on disk absent from prior-run sample export
│   ├── gse189367_ids_only_in_sample_export.txt # Diagnostic: (empty)
│   ├── gse189367_sample_vs_fastq_report.txt    # Diagnostic: full comparison report
│   ├── logs/
│   │   ├── GSE189367_88685.err
│   │   ├── GSE189367_88685.out
│   │   ├── fasterq_dump/
│   │   ├── pigz/
│   │   ├── prefetch/
│   │   ├── retry_fail_20260408_100150.txt
│   │   ├── retry_missing_88692.err / .out
│   │   ├── retry_missing_88705.err / .out
│   │   └── retry_success_*.txt
│   ├── metadata/
│   ├── removed_nbl_srr_id_primary_tumour_pe.txt  # 7 SRRs: prefetch_failed_unresolvable_after_retry
│   ├── retry_gse189367_missing_srrs.sh
│   ├── sra/
│   └── tmp/
│       ├── 88685/
│       └── manual_retry_SRR170109xx/  # 24 per-SRR staging dirs (retries completed)
│
├── SRP409177/               # SRA study — NBL (HiSeq 4000, CHLA)
│   ├── all_primary_tumour_pe.txt
│   ├── all_primary_tumour_srr.txt
│   ├── fastq/               # 47 paired-end samples
│   ├── logs/
│   │   ├── SRP409177_88686.err
│   │   ├── SRP409177_88686.out
│   │   ├── fasterq_dump/
│   │   ├── pigz/
│   │   └── prefetch/
│   ├── metadata/
│   ├── removed_nbl_srr_id_primary_tumour_pe.txt  # (empty — none excluded)
│   ├── sra/
│   └── tmp/
│       └── 88686/
│
├── target_nbl/              # GDC TARGET-NBL — STAR augmented counts
│   ├── count_data/          # 162 UUID directories, one per sample
│   │   └── {file_uuid}/
│   │       └── *.rna_seq.augmented_star_gene_counts.tsv
│   ├── download.log
│   ├── logs/
│   │   ├── target_nbl_counts_88683.err  # (empty — clean run)
│   │   └── target_nbl_counts_88683.out
│   ├── manifest.tsv         # GDC download manifest (id, filename, md5, size)
│   ├── manifests/
│   │   ├── gdc_manifest_target_nbl_star_counts.txt
│   │   └── payload_counts_manifest.json
│   ├── metadata/
│   │   ├── payload_counts.json
│   │   ├── target_nbl_aliquot_ids.txt
│   │   ├── target_nbl_aliquot_ids_unique.txt
│   │   ├── target_nbl_file_to_aliquot.tsv
│   │   └── target_nbl_star_counts_mapping.tsv
│   ├── target_nbl_id_set_comparison_report.txt
│   ├── target_nbl_ids_in_both.txt
│   ├── target_nbl_only_in_aliquot_list.txt
│   └── target_nbl_only_in_sample_list.txt  # (empty)
│
├── preprocessing_results/   # ⚠ PRIOR RUN — does not reflect current FASTQ data
│   ├── vst_NBL_joint_batch_corrected.rds  # Placeholder/stub (142 bytes — not a real matrix)
│   ├── tumour_sample_ids.txt              # 244 sample IDs from prior run
│   └── sample_id_exports/
│       ├── GSE49711_sample_ids.tsv        # Prior-run cohort not in current pipeline (5 samples)
│       ├── GSE100148_sample_ids.tsv       # Prior-run IDs (SRR5710xxx — differ from current FASTQs)
│       ├── GSE189367_sample_ids.tsv       # Prior-run IDs
│       ├── SRP409177_sample_ids.tsv       # Prior-run IDs
│       ├── TARGET_sample_ids.tsv          # Prior-run IDs (148 samples)
│       ├── all_sample_ids_by_cohort.tsv
│       └── sample_id_counts_by_cohort.tsv
│
└── count_data/              # ⚠ NOT YET GENERATED — pipeline output (Snakefile + merge scripts)
    ├── nbl_tumour_count.rds
    ├── nbl_tumour_sample_metadata.csv
    ├── nbl_tumour_gene_list_hgnc.txt
    └── nbl_ensembl_to_hgnc.tsv
```

---

## GEO Cohort Details

### GSE100148 — ATRX wild-type neuroblastoma

| Field | Value |
|-------|-------|
| GEO accession | GSE100148 |
| SRA study | SRP109627 |
| BioProject | PRJNA390893 |
| FASTQs on disk | 15 paired-end (SRR7828562–SRR7828576) |
| Read length | 100 bp |
| Platform | Illumina HiSeq 2500 |
| Library layout | PAIRED |
| Tissue | Primary tumour |

**Metadata files (`metadata/`):**
- `SraRunInfo.csv` — full SRA run table (includes both SRR5710xxx and SRR7828xxx series from SRP109627)
- `SraRunInfo_primary.csv` — filtered to primary tumour samples
- `SraRunInfo_primary_PE.csv` — further filtered to paired-end only (used for download)
- `gsm_titles.tsv` — GSM accession → sample title mapping
- `primary_gsms.tsv` — GSM accessions for primary tumour samples

**Root files:**
- `all_primary_tumour_srr.txt`, `all_primary_tumour_pe.txt` — pre-exclusion SRR ID lists (15 entries each)
- `removed_nbl_srr_id_primary_tumour_pe.txt` — (empty; no samples excluded)

---

### GSE189367 — NBL primary tumours (NovaSeq)

| Field | Value |
|-------|-------|
| GEO accession | GSE189367 |
| SRA study | SRP347311 |
| BioProject | PRJNA782665 |
| FASTQs on disk | 64 paired-end (SRR17010968–SRR17011031) |
| Read length | ~290 bp |
| Platform | Illumina NovaSeq 6000 |
| Library layout | PAIRED |
| Tissue | Primary tumour |

**Metadata files (`metadata/`):**
- `SraRunInfo.csv`, `SraRunInfo_primary.csv`, `SraRunInfo_primary_PE.csv` (64 samples each)

**Root files:**
- `all_primary_tumour_srr.txt`, `all_primary_tumour_pe.txt`
- `removed_nbl_srr_id_primary_tumour_pe.txt` — 7 SRRs marked `prefetch_failed_unresolvable_after_retry`
  (SRR17010968–70, 73–76); all were subsequently recovered via manual retry and are present in `fastq/`
- `retry_gse189367_missing_srrs.sh` — re-download script used for manual retry

**Diagnostic files:**
- `gse189367_sample_vs_fastq_report.txt` — comparison of prior-run sample export vs current FASTQ directory
- `gse189367_ids_only_in_fastq.txt` — 21 SRRs present on disk but absent from the prior-run sample export

**Download retry history (`logs/`):**
- SLURM retry jobs 88692 and 88705 logs; `retry_fail/success_*.txt` batch retry outcome records

---

### SRP409177 — NBL (HiSeq 4000, CHLA)

| Field | Value |
|-------|-------|
| SRA study | SRP409177 |
| BioProject | PRJNA904244 |
| Submitting centre | CHLA (Children's Hospital Los Angeles) |
| FASTQs on disk | 47 paired-end (SRR22373268–SRR22373314) |
| Read length | 202 bp |
| Platform | Illumina HiSeq 4000 |
| Library layout | PAIRED |
| Tissue | Primary tumour |

**Metadata files (`metadata/`):**
- `SraRunInfo.csv`, `SraRunInfo_primary.csv`, `SraRunInfo_primary_PE.csv`

**Root files:**
- `all_primary_tumour_srr.txt`, `all_primary_tumour_pe.txt`
- `removed_nbl_srr_id_primary_tumour_pe.txt` — (empty; no samples excluded)

---

## TARGET-NBL (GDC)

| Field | Value |
|-------|-------|
| GDC project | TARGET-NBL |
| Count files on disk | 162 |
| Data type | Gene Expression Quantification |
| Workflow | STAR - Counts |
| Annotation | GENCODE v36 |
| File format | Augmented STAR gene count TSV |
| Sample identifiers | `aliquot_submitter_id` (e.g. `TARGET-30-PARACS-01A-01R`) |
| Directory naming | GDC `file_id` (UUID) |
| Sample type | Primary Tumor |

**Count file format** (`*.rna_seq.augmented_star_gene_counts.tsv`):

```
# gene-model: GENCODE v36
gene_id    gene_name    gene_type    unstranded    stranded_first    stranded_second    tpm_unstranded    fpkm_unstranded    fpkm_uq_unstranded
N_unmapped ...
ENSG00000000003.15    TSPAN6    protein_coding    1373    679    694    ...
```

**Key mapping files (`metadata/`):**
- `target_nbl_file_to_aliquot.tsv` — maps `file_id` (UUID dir name) to `aliquot_submitter_id`; required by `merge_target_nbl_counts.R`
- `target_nbl_star_counts_mapping.tsv` — full GDC metadata per file (case ID, sample type, submitter IDs)
- `target_nbl_aliquot_ids.txt`, `target_nbl_aliquot_ids_unique.txt` — aliquot ID lists (162 total)
- `payload_counts.json` — GDC API payload used to query count files (reproducibility record)

**Diagnostic files (root-level):**
- `target_nbl_id_set_comparison_report.txt` — comparison of aliquot metadata list (162) vs a prior-run
  sample export (148); the 14-sample difference reflects prior analysis filtering, not current data
- `target_nbl_ids_in_both.txt` — 148 IDs shared between the two sets
- `target_nbl_only_in_aliquot_list.txt` — 14 aliquot IDs not in the prior-run sample export
- `target_nbl_only_in_sample_list.txt` — (empty)

---

## Sample Exclusion Files

Each GEO cohort has a `removed_nbl_srr_id_primary_tumour_pe.txt` at its root.
For GSE100148 and SRP409177 this file is empty (no exclusions). For GSE189367,
7 SRRs are listed as `prefetch_failed_unresolvable_after_retry` — these were later
recovered by manual retry and their FASTQs are present in `fastq/`.

The Snakefile discovers valid samples by scanning `fastq/` for `*_1.fastq.gz` /
`*_2.fastq.gz` pairs at runtime; this file does not gate pipeline execution.

---

## Expected Pipeline Outputs (Not Yet Generated)

### Per-sample (Snakefile outputs)

| File | Location |
|------|----------|
| Aligned BAM | `{cohort}/results/star/{sample}.Aligned.sortedByCoord.out.bam` |
| BAM index | `.bam.bai` |
| STAR log | `{cohort}/results/star/{sample}.Log.final.out` |
| Gene count table | `{cohort}/results/star/valid_count/{sample}.tab` |
| Picard metrics | `{cohort}/results/picard/{sample}_*.txt` |
| FastQC report | `{cohort}/results/fastqc/{sample_read}_fastqc.html` |
| MultiQC report | `{cohort}/results/multiqc/multiqc_report.html` |

The `.tab` files are 4-column STAR `ReadsPerGene.out.tab` format:
`gene_id | unstranded | stranded_first | stranded_second`.
Strandedness column is detected automatically per sample.

### Combined count matrix (merge scripts)

| File | Description |
|------|-------------|
| `count_data/nbl_tumour_count.rds` | genes × N samples integer count matrix |
| `count_data/nbl_tumour_sample_metadata.csv` | cohort, sample ID, strandedness |
| `count_data/nbl_tumour_gene_list_hgnc.txt` | Ensembl gene IDs (version suffix stripped) |
| `count_data/nbl_ensembl_to_hgnc.tsv` | Ensembl → HGNC symbol mapping |

**Gene universe:** intersection of GENCODE v44 (GEO cohorts) and GENCODE v36 (TARGET-NBL).

**Sample name conventions (once generated):**

| Cohort | Column name format | Example |
|--------|--------------------|---------|
| GSE100148 | SRR accession | `SRR7828562` |
| GSE189367 | SRR accession | `SRR17010968` |
| SRP409177 | SRR accession | `SRR22373268` |
| TARGET-NBL | aliquot_submitter_id | `TARGET-30-PARACS-01A-01R` |

---

## Preprocessing Results — Prior Run

`preprocessing_results/` contains outputs from a **previous pipeline run** on a
different sample set. These files do **not** correspond to the FASTQ data currently
in `fastq/` directories.

| File | Prior-run value | Note |
|------|-----------------|------|
| `sample_id_exports/GSE100148_sample_ids.tsv` | 4 samples (SRR5710xxx) | SRR5710xxx series; current FASTQs are SRR7828xxx |
| `sample_id_exports/GSE189367_sample_ids.tsv` | 43 samples | Subset of the 64 currently on disk |
| `sample_id_exports/SRP409177_sample_ids.tsv` | 44 samples | Vs 47 currently on disk |
| `sample_id_exports/TARGET_sample_ids.tsv` | 148 samples | Vs 162 currently on disk |
| `sample_id_exports/GSE49711_sample_ids.tsv` | 5 samples | Cohort not in current pipeline |
| `sample_id_exports/all_sample_ids_by_cohort.tsv` | 244 samples total | Includes GSE49711 |
| `tumour_sample_ids.txt` | 244 sample IDs | Mixed cohorts from prior run |
| `vst_NBL_joint_batch_corrected.rds` | 142 bytes | Placeholder/stub — not a real matrix |

These files should be regenerated once the current Snakefile run completes and
`merge_star_counts.R` has been executed.

---

## Reference Genome

GEO cohorts are aligned against:

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
1. **FastQC** — per-read QC (R1 and R2 separately)
2. **STAR 2-pass** — splice-aware alignment, coordinate-sorted BAM + gene counts
3. **samtools index** — BAM indexing
4. **Picard ValidateSamFile** — BAM integrity check
5. **Picard CollectAlignmentSummaryMetrics**
6. **Picard CollectQualityYieldMetrics**
7. **Picard CollectInsertSizeMetrics**
8. **MultiQC** — per-project aggregated QC report

Count matrix assembly:
- `merge_target_nbl_counts.R` — merges TARGET-NBL TSV files from `target_nbl/count_data/`
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
**Status:** Not yet run on current data. Depends on `count_data/` being generated first.

### Purpose

Estimates tumour microenvironment signals (stromal, immune, ESTIMATE scores) using
the ESTIMATE method, then filters samples below a purity threshold to remove those
dominated by non-tumour signal.

### Inputs (all pending Snakefile + merge scripts)

| File | Source |
|------|--------|
| `count_data/nbl_tumour_count.rds` | `merge_star_counts.R` |
| `count_data/nbl_tumour_sample_metadata.csv` | `merge_star_counts.R` |
| `count_data/nbl_ensembl_to_hgnc.tsv` | Ensembl → HGNC mapping table |

### Outputs (written to `results/tumour_purity_analysis/nbl/`)

| File | Description |
|------|-------------|
| `estimate_scores.csv` | Per-sample stromal, immune, ESTIMATE scores + derived purity |
| `filtered_hgnc_counts_purity0.7.rds` | Filtered count matrix in HGNC gene symbol space |
| `nbl_filtered_count_matrix_ensembl_purity0.7.rds` | Filtered samples in Ensembl ID space |
| `sample_retention_by_project_same_axis.pdf` | Sample counts before/after filtering by cohort |
| `purity_diagnostics.pdf` | Purity distribution and diagnostic scatter plots |
| `nbl_purity_pipeline_<timestamp>.log` | Full run log |

### Purity threshold

Default threshold: **0.7**. Configurable via `PURITY_THRESHOLD` environment variable.

### Invocation

```bash
SE_PATH=/work/ugbogu/pipeline/data/nbl/count_data/nbl_tumour_count.rds \
META_CSV=/work/ugbogu/pipeline/data/nbl/count_data/nbl_tumour_sample_metadata.csv \
MAP_TSV=/work/ugbogu/pipeline/data/nbl/count_data/nbl_ensembl_to_hgnc.tsv \
OUTPUT_DIR=/work/ugbogu/pipeline/results/tumour_purity_analysis/nbl \
Rscript preprocessing_and_quality_control/tumour_purity_analysis/nbl/nbl_tumour_purity_analysis.R
```
