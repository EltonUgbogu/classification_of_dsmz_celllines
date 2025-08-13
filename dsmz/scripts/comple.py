#!/usr/bin/env python3
"""
Generate a PDF README that summarizes the full DSMZ–TCGA workflow.
Writes: README_full_workflow.pdf
"""

from reportlab.lib.pagesizes import LETTER
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, ListFlowable, ListItem
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import os
import sys

OUTFILE = "README_full_workflow.pdf"

def register_fonts():
    """
    Try to register DejaVu Sans for full Unicode support.
    Fallback to built-in Helvetica if not found.
    """
    try:
        # Common Linux path; adjust if you keep TTFs elsewhere.
        dejavu_path = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
        if os.path.exists(dejavu_path):
            pdfmetrics.registerFont(TTFont("DejaVuSans", dejavu_path))
            return "DejaVuSans"
    except Exception:
        pass
    return "Helvetica"

def styles_for(font_name: str):
    ss = getSampleStyleSheet()

    # Base body
    body = ParagraphStyle(
        "Body",
        parent=ss["BodyText"],
        fontName=font_name,
        fontSize=10.5,
        leading=14,
        spaceAfter=6,
    )

    title = ParagraphStyle(
        "Title",
        parent=ss["Title"],
        fontName=font_name,
        fontSize=20,
        leading=24,
        spaceAfter=14,
    )

    h1 = ParagraphStyle(
        "Heading1",
        parent=ss["Heading1"],
        fontName=font_name,
        fontSize=14.5,
        leading=18,
        spaceBefore=12,
        spaceAfter=6,
    )

    h2 = ParagraphStyle(
        "Heading2",
        parent=ss["Heading2"],
        fontName=font_name,
        fontSize=12.5,
        leading=16,
        spaceBefore=10,
        spaceAfter=4,
    )

    bullet = ParagraphStyle(
        "Bullet",
        parent=body,
        leftIndent=14,
        bulletIndent=6,
        spaceBefore=2,
        spaceAfter=2,
    )

    code = ParagraphStyle(
        "Code",
        parent=body,
        fontName=font_name,
        fontSize=9.5,
        leading=13,
        backColor="#F5F5F5",
        leftIndent=8,
        rightIndent=8,
        spaceBefore=4,
        spaceAfter=6,
    )

    return {"title": title, "h1": h1, "h2": h2, "body": body, "bullet": bullet, "code": code}

def bullet_list(items, style, level=0):
    return ListFlowable(
        [ListItem(Paragraph(txt, style), leftIndent=level * 12) for txt in items],
        bulletType="bullet",
        start="•",
        leftIndent=0,
        bulletFontName=style.fontName,
        bulletFontSize=style.fontSize,
    )

def main():
    font = register_fonts()
    st = styles_for(font)

    doc = SimpleDocTemplate(
        OUTFILE,
        pagesize=LETTER,
        leftMargin=0.9 * inch,
        rightMargin=0.9 * inch,
        topMargin=0.8 * inch,
        bottomMargin=0.8 * inch,
        title="DSMZ–TCGA Workflow README",
        author="",
    )

    story = []
    story.append(Paragraph("DSMZ–TCGA Correlation Workflow", st["title"]))
    story.append(Paragraph(
        "This README summarizes the end-to-end workflow for building RNA-seq count data "
        "from STAR outputs, performing quality checks, inferring strandedness, and computing "
        "DSMZ–TCGA expression correlations with optional tumor purity adjustment and UMAP visualization.",
        st["body"]
    ))

    # Inputs / Key Data
    story.append(Paragraph("Key Inputs", st["h1"]))
    story.append(bullet_list([
        "<b>FASTQ lists</b>: text file of FASTQ paths for sample discovery.",
        "<b>BAM files</b>: STAR-aligned, coordinate-sorted BAMs plus .bai.",
        "<b>STAR ReadsPerGene</b> outputs: one file per sample.",
        "<b>GENCODE GTF</b>: e.g., Homo_sapiens.GRCh38.107.gtf.",
        "<b>TCGA SummarizedExperiment</b>: gene counts for TCGA tumors (RDS).",
        "<b>DSMZ metadata</b>: CSV with sample_name and Origin (organ source).",
        "<b>(Optional) Purity table</b>: per-sample tumor purity for TCGA.",
    ], st["bullet"]))

    # Outputs
    story.append(Paragraph("Key Outputs", st["h1"]))
    story.append(bullet_list([
        "<b>samples.tsv</b>: paired-end sample manifest for downstream use.",
        "<b>bam_coverage.tsv</b>: mapped-read totals and coverage status.",
        "<b>*.infer.txt</b>: strandedness inference per BAM (RSeQC).",
        "<b>DSMZ_star2count.rds</b>: merged counts from STAR outputs.",
        "<b>DSMZ_count_gene.rds</b>: counts with Ensembl and gene_name columns.",
        "<b>spearman_all_wide.csv</b> and <b>spearman_all_long.csv</b>: all DSMZ vs. TCGA cohort Spearman correlations.",
        "<b>dsmz_tcga_mapping.csv</b> and <b>spearman_best_by_sample.csv</b>: organ-constrained and best overall matches.",
        "<b>diagnostic PDFs</b>: violin plots, heatmaps, UMAP embeddings.",
    ], st["bullet"]))

    # Workflow Steps
    story.append(Paragraph("Workflow Steps", st["h1"]))

    # 1
    story.append(Paragraph("1) Generate Sample Manifest", st["h2"]))
    story.append(Paragraph(
        "Script: <b>generate_samples_tsv.py</b>. Reads a text file listing FASTQ files, groups R1/R2 by prefix, "
        "and writes <b>config/samples.tsv</b> for downstream processes.", st["body"]))
    story.append(Paragraph("Example:", st["body"]))
    story.append(Paragraph(
        "python generate_samples_tsv.py  (paths configured in the script)",
        st["code"]))

    # 2
    story.append(Paragraph("2) BAM Coverage QC", st["h2"]))
    story.append(Paragraph(
        "Script: <b>bam_coverage.py</b>. Runs <i>samtools idxstats</i> and totals mapped reads per BAM to flag low-coverage samples.", st["body"]))
    story.append(Paragraph("Example:", st["body"]))
    story.append(Paragraph(
        "python bam_coverage.py  (adjust bam_dir, output_path, min_reads in the script)",
        st["code"]))

    # 3
    story.append(Paragraph("3) Strandedness Inference", st["h2"]))
    story.append(Paragraph(
        "Script: <b>run_infer_experiment.py</b>. Validates BAMs, runs RSeQC infer_experiment.py with BED12 annotation, "
        "and summarizes strandedness calls.", st["body"]))
    story.append(Paragraph("Example:", st["body"]))
    story.append(Paragraph(
        "python run_infer_experiment.py  (paths configured at top of script)",
        st["code"]))

    # 4
    story.append(Paragraph("4) Build DSMZ Count Tables From STAR", st["h2"]))
    story.append(Paragraph(
        "R Markdown: <b>star_2_countdata.Rmd</b>. Collects all STAR ReadsPerGene.out.tab files, extracts unstranded counts, "
        "merges across samples, maps Ensembl IDs to gene_name via GENCODE GTF, and writes:", st["body"]))
    story.append(bullet_list([
        "<b>DSMZ_star2count.rds</b> – merged raw counts with gene_id_ver.",
        "<b>DSMZ_count_gene.rds</b> – tidy table with Ensembl_ID_with_version, Ensembl_ID (no version), gene_name, then samples.",
    ], st["bullet"]))
    story.append(Paragraph("Example:", st["body"]))
    story.append(Paragraph(
        "Rscript -e 'rmarkdown::render(\"star_2_countdata.Rmd\", output_format=\"pdf_document\")'",
        st["code"]))

    # 5
    story.append(Paragraph("5) Compute Gene Lengths (Optional)", st["h2"]))
    story.append(Paragraph(
        "Script: <b>gtf_exon_lengths.py</b>. Parses GTF exon features, merges overlaps, and sums non-overlapping exon length per gene.", st["body"]))
    story.append(Paragraph("Example:", st["body"]))
    story.append(Paragraph(
        "python gtf_exon_lengths.py  (driven by Snakemake or direct args, depending on your variant)",
        st["code"]))

    # 6
    story.append(Paragraph("6) TPM Calculation (Optional)", st["h2"]))
    story.append(Paragraph(
        "Script: <b>compute_tpm.py</b>. Merges raw counts with gene lengths and computes TPM per sample.", st["body"]))
    story.append(Paragraph("Example:", st["body"]))
    story.append(Paragraph(
        "python compute_tpm.py  (driven by Snakemake or direct args, depending on your variant)",
        st["code"]))

    # 7
    story.append(Paragraph("7) DSMZ–TCGA Correlation Analysis", st["h2"]))
    story.append(Paragraph(
        "R script: <b>dsmz_tcga_correlation.R</b>. Harmonizes Ensembl IDs (version handling), aligns DSMZ counts to metadata, "
        "optionally filters and regresses TCGA by tumor purity, computes TCGA cohort means, selects variable genes if requested, "
        "computes Spearman correlations (DSMZ samples by TCGA cohorts), performs organ-constrained mapping, "
        "and generates plots (violin, heatmaps) and UMAP embeddings.", st["body"]))
    story.append(Paragraph("Key outputs:", st["body"]))
    story.append(bullet_list([
        "<b>spearman_all_wide.csv</b>, <b>spearman_all_long.csv</b> – all correlations.",
        "<b>dsmz_tcga_mapping.csv</b>, <b>spearman_best_by_sample.csv</b> – organ-constrained and best overall matches.",
        "<b>violin_rho_by_organ.pdf</b>, <b>heatmap_top25_dsmz_vs_tcga_cohorts.pdf</b>, "
        "<b>organ_tcga_assignment_counts.pdf</b>, <b>umap_*.pdf</b> – diagnostics and embeddings.",
    ], st["bullet"]))
    story.append(Paragraph("Example:", st["body"]))
    story.append(Paragraph(
        "Rscript dsmz_tcga_correlation.R",
        st["code"]))

    # Notes
    story.append(Paragraph("Notes and Tips", st["h1"]))
    story.append(bullet_list([
        "Keep BAMs coordinate-sorted and indexed. Coverage and strandedness steps assume valid .bam + .bai.",
        "Ensure genome build consistency between BAMs and BED/GTF (e.g., GRCh38 hg38).",
        "If you see very few shared genes between DSMZ and TCGA, inspect Ensembl ID versioning and mapping.",
        "On HPC, submit each step as a standalone job, or wire them with a workflow manager later.",
    ], st["bullet"]))

    doc.build(story)
    print(f"Wrote {OUTFILE}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
