# !/usr/bin/env snakemake
# Snakefile for breast cancer subtype analysis using TCGA and DSMZ data
# Author: Elton Ugbogu
# Description:
#   - Load TCGA BRCA and DSMZ breast cancer cell line data
#   - Perform joint VST normalization
#   - Unsupervised clustering using HVGs
#   - UMAP visualization and marker expression analysis
#   - Map DSMZ samples to TCGA subtypes

# -----------------------------
# Paths
# -----------------------------
TCGA_SE_RDS   = "/home/chu25/data/tcga/ALL_TCGA_STAR_Counts_SummarizedExperiment_filtered.rds"
DSMZ_RDS      = "/home/chu25/data/dsmz/DSMZ_count_gene.rds"
DSMZ_META_CSV = "/home/chu25/data/dsmz/DSMZ_metadata.csv"
OUTDIR        = "/home/chu25/dsmz/results/brca_vst_dual"
R_SCRIPT      = "breast_cancer_subtypes.R"
LOGDIR        = "logs"
INIT_MARKER   = f"{OUTDIR}/.init"   # <- marker file for directory creation

# Core outputs
VST_TCGA_RDS       = f"{OUTDIR}/VST_tcga_brca.rds"
VST_DSMZ_RDS       = f"{OUTDIR}/VST_dsmz_brca.rds"
DIM_REPORT         = f"{OUTDIR}/dimension_report.txt"
KMEANS_METRICS     = f"{OUTDIR}/kmeans_metrics_HVG.csv"
KMEANS_ELBOW_PDF   = f"{OUTDIR}/kmeans_elbow_HVG.pdf"
KMEANS_PANEL_PDF   = f"{OUTDIR}/kmeans_metrics_panel_HVG.pdf"
CLUSTERS_CSV       = f"{OUTDIR}/clusters_HVG_kmeans.csv"
UMAP_DS_PDF        = f"{OUTDIR}/umap_dataset.pdf"
UMAP_TCGA_PDF      = f"{OUTDIR}/umap_tcga_subtype.pdf"
UMAP_DSMZ_PDF      = f"{OUTDIR}/umap_dsmz_cluster.pdf"
UMAP_COMBINED_PDF  = f"{OUTDIR}/umap_combined.pdf"
UMAP_EMBED_CSV     = f"{OUTDIR}/umap_discovery_HVG_coords.csv"
HEATMAP_MARKER_PDF = f"{OUTDIR}/heatmap_marker_vst.pdf"
HVG_LIST_TXT       = f"{OUTDIR}/HVGs_used_for_discovery_3k.txt"
MAP_DSMZ_TCGA_CSV  = f"{OUTDIR}/mapping_dsmz_to_tcga.csv"
SENTINEL           = f"{OUTDIR}/_SUCCESS.txt"

rule all:
    input:
        VST_TCGA_RDS,
        VST_DSMZ_RDS,
        DIM_REPORT,
        KMEANS_METRICS,
        KMEANS_ELBOW_PDF,
        KMEANS_PANEL_PDF,
        CLUSTERS_CSV,
        UMAP_DS_PDF,
        UMAP_TCGA_PDF,
        UMAP_DSMZ_PDF,
        UMAP_COMBINED_PDF,
        UMAP_EMBED_CSV,
        HEATMAP_MARKER_PDF,
        HVG_LIST_TXT,
        MAP_DSMZ_TCGA_CSV,
        SENTINEL

# Create OUTDIR via a marker file instead of directory() output
rule make_outdir:
    output:
        INIT_MARKER
    shell:
        "mkdir -p {OUTDIR} && touch {output}"

rule breast_cancer_subtypes:
    input:
        tcga = TCGA_SE_RDS,
        dsmz = DSMZ_RDS,
        meta = DSMZ_META_CSV,
        init = INIT_MARKER         # <- require the dir marker
    output:
        VST_TCGA_RDS,
        VST_DSMZ_RDS,
        DIM_REPORT,
        KMEANS_METRICS,
        KMEANS_ELBOW_PDF,
        KMEANS_PANEL_PDF,
        CLUSTERS_CSV,
        UMAP_DS_PDF,
        UMAP_TCGA_PDF,
        UMAP_DSMZ_PDF,
        UMAP_COMBINED_PDF,
        UMAP_EMBED_CSV,
        HEATMAP_MARKER_PDF,
        HVG_LIST_TXT,
        MAP_DSMZ_TCGA_CSV,
        SENTINEL
    threads: 8
    resources:
        mem_mb = 64000
    log:
        out = f"{LOGDIR}/brca_subtypes.out",
        err = f"{LOGDIR}/brca_subtypes.err"
    conda:
        "deseq2-r-env"     # keep DESeq2-safe env for this rule
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{LOGDIR}"
        Rscript "{R_SCRIPT}" > "{log.out}" 2> "{log.err}"
        touch "{SENTINEL}"
        """