# DSMZ–TCGA RNA-seq Correlation Workflow

This workflow processes DSMZ RNA-seq data, prepares expression matrices, and computes correlations with TCGA tumor cohorts.
All steps are implemented as **standalone** scripts (R, Python, and shell). Run in sequence or adapt to your scheduler.

---

## 1) Sample Table Generation
Script: `generate_samples_tsv.py` – builds `config/samples.tsv` from FASTQs.

## 2) Gene Length Extraction
Script: `gtf_exon_lengths.py` – computes total non-overlapping exon length per gene from GTF.

## 3) Alignment and Read Counting
Tool: STAR – aligns FASTQs and produces `ReadsPerGene.out.tab` per sample.

## 4) TPM Calculation
Script: `compute_tpm.py` – combines STAR counts and gene lengths to compute TPM.

## 5) BAM Coverage QC
Script: `bam_coverage.py` – summarizes mapped reads per BAM using samtools idxstats.

## 6) Strandedness Inference
Script: `run_infer_experiment.py` – wraps RSeQC infer_experiment.py.

## 7) STAR → Counts Summary
Rmd: `star_2_countdata.Rmd` – merges STAR count files, maps Ensembl IDs to gene symbols.

## 8) Preliminary Correlation
Rmd: `preliminary_analysis.Rmd` – harmonizes genes, optional purity adjustment, Spearman correlations.

## 9) Full DSMZ–TCGA Correlation
R: `dsmz_tcga_correlation.R` – organ-constrained mapping, Spearman correlations, UMAP visualization.

---
