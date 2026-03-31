#!/usr/bin/env Rscript

# =============================================================================
# enrichment_analysis.R
# =============================================================================
#
# DESCRIPTION
# -----------
# This script performs functional enrichment analysis on gene lists derived
# from differential expression or clustering analyses. Enrichment analysis
# identifies biological pathways, processes, and molecular functions that are
# statistically overrepresented in a gene list relative to a background
# genome, thereby providing biological interpretation of high-throughput
# genomics results.
#
# The script integrates three complementary knowledge bases:
#
#   (1) Gene Ontology (GO): Structured vocabulary describing gene function
#       across Biological Process (BP), Molecular Function (MF), and
#       Cellular Component (CC) domains.
#
#   (2) MSigDB Hallmark Gene Sets: Curated signatures representing well-defined
#       biological states and processes, designed to reduce redundancy whilst
#       retaining biological coherence.
#
#   (3) Reactome Pathways: Manually curated, peer-reviewed pathway database
#       providing detailed mechanistic representations of biological processes.
#
#
# BIOLOGICAL CONTEXT: FUNCTIONAL ENRICHMENT
# -----------------------------------------
# Differential expression analysis identifies genes with altered expression
# between conditions, but individual gene lists are difficult to interpret
# biologically. Functional enrichment analysis addresses this by asking:
#
#   "Are genes in my list associated with particular biological functions
#    more often than expected by chance?"
#
# This question is formalised using the hypergeometric test (Fisher's exact
# test), which evaluates whether the overlap between a gene list and a
# functional category exceeds random expectation.
#
#
# STATISTICAL FRAMEWORK
# ---------------------
# \subsection{Hypergeometric Test}
# For a gene set $S$ of size $k$ from a genome of $N$ genes, with $K$ genes
# annotated to a functional term and $x$ genes in $S$ also annotated to that
# term, the probability of observing $x$ or more overlapping genes by chance
# is:
#
# \begin{equation}
#     P(X \geq x) = 1 - \sum_{i=0}^{x-1}
#         \frac{\binom{K}{i} \binom{N-K}{k-i}}{\binom{N}{k}}
# \end{equation}
#
# This is equivalent to drawing $k$ balls without replacement from an urn
# containing $K$ white balls (term-annotated genes) and $N-K$ black balls
# (non-annotated genes), and asking for the probability of drawing $\geq x$
# white balls.
#
# \subsection{Multiple Testing Correction}
# Thousands of functional terms are tested simultaneously, necessitating
# correction for multiple comparisons. The script employs the
# Benjamini-Hochberg (BH) procedure to control the false discovery rate
# (FDR), providing the expected proportion of false positives among
# significant results:
#
# \begin{equation}
#     q_i = \min\left(1, \frac{p_i \cdot m}{\text{rank}(p_i)}\right)
# \end{equation}
#
# where $m$ is the total number of tests. An FDR threshold of 0.05 indicates
# that approximately 5\% of reported enrichments are expected to be false
# positives.
#
# \subsection{Gene Set Size Filtering}
# Very small gene sets ($< 10$ genes) provide insufficient statistical power,
# whilst very large gene sets ($> 500$ genes) are often too general to be
# biologically informative. The script applies size filters to focus on
# interpretable, well-powered results.
#
#
# KNOWLEDGE BASES
# ---------------
# \subsection{Gene Ontology (GO)}
# The Gene Ontology Consortium maintains a controlled vocabulary organised
# into three orthogonal domains:
#
# \begin{description}
#   \item[Biological Process (BP)] Multi-step molecular events with defined
#         beginning and end (e.g., "apoptotic process", "cell cycle")
#   \item[Molecular Function (MF)] Elemental activities at the molecular level
#         (e.g., "kinase activity", "DNA binding")
#   \item[Cellular Component (CC)] Locations within or outside the cell
#         (e.g., "nucleus", "plasma membrane")
# \end{description}
#
# GO terms are organised hierarchically; enrichment at a child term often
# propagates to parent terms, creating apparent redundancy. The
# \texttt{simplify()} function (not used here) can reduce this redundancy
# by removing highly similar terms.
#
# \subsection{MSigDB Hallmark Gene Sets}
# The Molecular Signatures Database (MSigDB) Hallmark collection comprises
# 50 gene sets representing well-defined biological states and processes.
# Key design principles:
#
# \begin{itemize}
#   \item Derived from multiple founder sets to reduce noise
#   \item Refined to minimise overlap between sets
#   \item Cover fundamental processes: proliferation, metabolism, signalling,
#         immune response, and development
# \end{itemize}
#
# Example Hallmark sets relevant to cancer:
#   - HALLMARK\_E2F\_TARGETS (proliferation)
#   - HALLMARK\_EPITHELIAL\_MESENCHYMAL\_TRANSITION
#   - HALLMARK\_INFLAMMATORY\_RESPONSE
#   - HALLMARK\_MYC\_TARGETS\_V1/V2
#
# \subsection{Reactome Pathways}
# Reactome provides manually curated, peer-reviewed pathway annotations
# organised hierarchically. Unlike GO, Reactome pathways represent
# mechanistic relationships between molecules (reactions, complexes,
# regulatory interactions).
#
# Reactome is particularly valuable for:
#   - Signalling pathway analysis (e.g., MAPK, PI3K-AKT, Wnt)
#   - Metabolic pathway analysis
#   - Identifying druggable pathway nodes
#
#
# IDENTIFIER CONVERSION
# ---------------------
# Enrichment databases typically use Entrez Gene IDs as the primary key.
# This script accepts gene symbols (HGNC nomenclature) and converts them
# to Entrez IDs using the \texttt{bitr()} function from clusterProfiler.
#
# Conversion failures may occur due to:
#   - Outdated or non-standard gene symbols
#   - Genes lacking Entrez annotations (e.g., novel transcripts)
#   - Ambiguous symbols mapping to multiple Entrez IDs
#
# The script reports conversion statistics to alert users to potential
# data loss.
#
#
# INPUT REQUIREMENTS
# ------------------
# The script expects a tab-separated file containing:
#
#   - A column with gene symbols (default: "symbol")
#   - Additional columns are ignored but preserved in memory
#
# Typical inputs include:
#   - Differentially expressed gene lists from DESeq2/edgeR/limma
#   - Cluster marker genes from Seurat/Scanpy
#   - Genes from custom filtering pipelines
#
#
# OUTPUT FILES
# ------------
# The script produces enrichment results in tab-separated format:
#
#   \begin{verbatim}
#   {outdir}/
#   +-- {prefix}_GO.tsv        # Gene Ontology enrichment results
#   +-- {prefix}_Hallmark.tsv  # MSigDB Hallmark enrichment (if available)
#   +-- {prefix}_Reactome.tsv  # Reactome pathway enrichment (if available)
#   \end{verbatim}
#
# Each output file contains columns:
#   - ID: Term/pathway identifier
#   - Description: Human-readable term name
#   - GeneRatio: Proportion of query genes in the term
#   - BgRatio: Proportion of background genes in the term
#   - pvalue: Hypergeometric test $p$-value
#   - p.adjust: BH-adjusted $p$-value (FDR)
#   - qvalue: Local FDR estimate
#   - geneID: List of query genes in the term
#   - Count: Number of query genes in the term
#
#
# USAGE EXAMPLE
# -------------
#   \begin{verbatim}
#   Rscript enrichment_analysis.R \
#     --markers component_1_vs_rest_UP_top500.annotated.tsv \
#     --gene_col symbol \
#     --outdir results/enrichment \
#     --prefix component_1_UP \
#     --pvalueCutoff 0.05 \
#     --qvalueCutoff 0.2
#   \end{verbatim}
#
#
# DEPENDENCIES
# ------------
# Required:
#   - clusterProfiler: Core enrichment analysis framework
#   - org.Hs.eg.db: Human gene annotation database
#   - enrichplot: Visualisation functions
#   - ggplot2, dplyr, readr: Data manipulation and I/O
#
# Optional (enhances functionality):
#   - msigdbr: Access to MSigDB gene sets including Hallmark
#   - ReactomePA: Reactome pathway enrichment
#
# Bioconductor installation:
#   \begin{verbatim}
#   BiocManager::install(c("clusterProfiler", "org.Hs.eg.db",
#                          "enrichplot", "ReactomePA"))
#   install.packages("msigdbr")
#   \end{verbatim}
#
#
# REFERENCES
# ----------
# Wu, T., et al. (2021). clusterProfiler 4.0: A universal enrichment tool
# for interpreting omics data. \textit{The Innovation}, 2(3), 100141.
#
# Subramanian, A., et al. (2005). Gene set enrichment analysis: A knowledge-
# based approach for interpreting genome-wide expression profiles.
# \textit{PNAS}, 102(43), 15545--15550.
#
# Liberzon, A., et al. (2015). The Molecular Signatures Database Hallmark
# Gene Set Collection. \textit{Cell Systems}, 1(6), 417--425.
#
# Jassal, B., et al. (2020). The Reactome Pathway Knowledgebase.
# \textit{Nucleic Acids Research}, 48(D1), D498--D503.
#
#
# AUTHOR
# ------
# Developed as part of the differential expression analysis pipeline for
# biological interpretation of marker gene lists in cancer genomics analyses.
#
# =============================================================================

# -----------------------------------------------------------------------------
# Package Loading
# -----------------------------------------------------------------------------
# The script loads required packages with startup messages suppressed.
# Core packages are mandatory; optional packages (ReactomePA, msigdbr)
# are loaded conditionally to allow the script to run in environments
# with partial installations.
#
# \paragraph{Package Purposes:}
# \begin{description}
#   \item[\texttt{clusterProfiler}] Core enrichment analysis framework
#         providing \texttt{enrichGO()}, \texttt{enricher()}, and
#         identifier conversion utilities
#   \item[\texttt{org.Hs.eg.db}] Human genome annotation database for
#         identifier mapping and GO annotations
#   \item[\texttt{enrichplot}] Publication-quality visualisation functions
#         for enrichment results (dot plots, gene-concept networks)
#   \item[\texttt{ggplot2}] Grammar of graphics plotting framework
#   \item[\texttt{dplyr, readr}] Data manipulation and I/O utilities
# \end{description}

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# -----------------------------------------------------------------------------
# Optional Package Loading: ReactomePA
# -----------------------------------------------------------------------------
# ReactomePA provides access to the Reactome pathway database. The package
# requires network access for initial setup and may not be available in
# all environments. The script checks availability and proceeds without
# Reactome enrichment if the package is missing.

has_reactome <- requireNamespace("ReactomePA", quietly = TRUE)
if (has_reactome) {
  suppressPackageStartupMessages(library(ReactomePA))
}

# -----------------------------------------------------------------------------
# Command-Line Argument Definitions
# -----------------------------------------------------------------------------
# The script accepts arguments organised into three categories:
#
# \paragraph{Input/Output:}
# \texttt{--markers}, \texttt{--outdir}, \texttt{--prefix}
#
# \paragraph{Gene Identifier:}
# \texttt{--gene\_col} specifies which column contains gene symbols
#
# \paragraph{Statistical Thresholds:}
# \begin{description}
#   \item[\texttt{--minGSSize}] Minimum gene set size (default: 10).
#         Sets smaller than this lack statistical power.
#   \item[\texttt{--maxGSSize}] Maximum gene set size (default: 500).
#         Very large sets are often too general to be informative.
#   \item[\texttt{--pvalueCutoff}] Nominal $p$-value threshold (default: 0.05).
#         Terms exceeding this are excluded from results.
#   \item[\texttt{--qvalueCutoff}] FDR threshold (default: 0.2).
#         More permissive than $p$-value to retain borderline findings
#         for exploratory analysis.
# \end{description}

opt_list <- list(
  make_option("--markers", type = "character",
              help = "Path to marker TSV file with gene symbols"),
  make_option("--gene_col", type = "character", default = "symbol",
              help = "Column name containing gene symbols"),
  make_option("--outdir", type = "character",
              help = "Output directory"),
  make_option("--prefix", type = "character",
              help = "Output file prefix"),
  make_option("--minGSSize", type = "integer", default = 10,
              help = "Minimum gene set size"),
  make_option("--maxGSSize", type = "integer", default = 500,
              help = "Maximum gene set size"),
  make_option("--pvalueCutoff", type = "numeric", default = 0.05,
              help = "P-value cutoff"),
  make_option("--qvalueCutoff", type = "numeric", default = 0.2,
              help = "Q-value cutoff")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# -----------------------------------------------------------------------------
# Input Validation
# -----------------------------------------------------------------------------
# Early validation ensures informative error messages for missing arguments.

if (is.null(opt$markers) || is.null(opt$outdir) || is.null(opt$prefix)) {
  stop("--markers, --outdir, and --prefix are required")
}

# -----------------------------------------------------------------------------
# Output Directory Initialisation
# -----------------------------------------------------------------------------

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Marker Gene Loading
# -----------------------------------------------------------------------------
# The script loads the marker gene file and extracts the specified gene
# symbol column. Missing values and empty strings are removed to prevent
# downstream errors in identifier conversion.

cat("[INFO] Loading markers...\n")
markers_df <- read_tsv(opt$markers, show_col_types = FALSE)

if (!opt$gene_col %in% colnames(markers_df)) {
  stop(paste0("Gene column '", opt$gene_col, "' not found"))
}

# Extract and clean gene symbols
gene_list <- markers_df[[opt$gene_col]]
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
gene_list <- unique(gene_list)

cat(sprintf("[INFO] Input genes: %d\n", length(gene_list)))

# -----------------------------------------------------------------------------
# Gene Symbol to Entrez ID Conversion
# -----------------------------------------------------------------------------
# Most enrichment databases use Entrez Gene IDs as the primary identifier.
# The \texttt{bitr()} (Biological ID Translator) function from clusterProfiler
# performs batch conversion between identifier types.
#
# \paragraph{Conversion Parameters:}
# \begin{description}
#   \item[\texttt{fromType = "SYMBOL"}] Input identifier type (HGNC symbols)
#   \item[\texttt{toType = "ENTREZID"}] Target identifier type
#   \item[\texttt{OrgDb = org.Hs.eg.db}] Annotation database for Homo sapiens
#   \item[\texttt{drop = TRUE}] Remove unmapped identifiers from output
# \end{description}
#
# Conversion rates below 80\% may indicate issues with gene symbol formatting
# or the use of non-standard nomenclature.

gene_map <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)

entrez_ids <- unique(gene_map$ENTREZID)
cat(sprintf("[INFO] Mapped to Entrez IDs: %d\n", length(entrez_ids)))

# -----------------------------------------------------------------------------
# Minimum Gene Count Validation
# -----------------------------------------------------------------------------
# Enrichment analysis requires a minimum number of genes for meaningful
# statistical inference. With fewer than 10 genes, the hypergeometric test
# has severely limited power, and results are unlikely to be reliable.

if (length(entrez_ids) < 10) {
  stop("Too few genes mapped to Entrez IDs for enrichment analysis")
}

# -----------------------------------------------------------------------------
# Gene Ontology Enrichment Analysis
# -----------------------------------------------------------------------------
# The \texttt{enrichGO()} function performs over-representation analysis
# against Gene Ontology terms. Setting \texttt{ont = "ALL"} tests all three
# GO domains (BP, MF, CC) simultaneously.
#
# \paragraph{Key Parameters:}
# \begin{description}
#   \item[\texttt{pAdjustMethod = "BH"}] Benjamini-Hochberg FDR correction
#   \item[\texttt{readable = TRUE}] Convert Entrez IDs back to gene symbols
#         in output for interpretability
# \end{description}
#
# The function returns an \texttt{enrichResult} S4 object containing:
#   - @result: Data frame of enrichment statistics
#   - @gene: Input gene list
#   - @keytype: Identifier type used
#   - @ontology: GO domain(s) tested

cat("[INFO] Running GO enrichment (BP, MF, CC)...\n")
go_enrich <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = opt$pvalueCutoff,
  qvalueCutoff = opt$qvalueCutoff,
  minGSSize = opt$minGSSize,
  maxGSSize = opt$maxGSSize,
  readable = TRUE
)

# -----------------------------------------------------------------------------
# MSigDB Hallmark Gene Set Enrichment (Optional)
# -----------------------------------------------------------------------------
# The msigdbr package provides programmatic access to the Molecular Signatures
# Database (MSigDB). The Hallmark collection (collection = "h") contains 50
# curated gene sets representing fundamental biological processes.
#
# Unlike GO enrichment, Hallmark analysis uses \texttt{enricher()}, a generic
# function that accepts custom gene set definitions via the TERM2GENE argument.
#
# The analysis is wrapped in tryCatch to handle network failures gracefully,
# as msigdbr requires internet access to download gene set definitions.

hallmark_enrich <- NULL
if (requireNamespace("msigdbr", quietly = TRUE)) {
  tryCatch({
    cat("[INFO] Attempting MSigDB Hallmark enrichment...\n")
    suppressPackageStartupMessages(library(msigdbr))

    # Retrieve Hallmark gene sets for Homo sapiens
    # The returned data frame contains columns: gs_name, entrez_gene, etc.
    hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "h")

    # Run enrichment using the generic enricher function
    # TERM2GENE maps gene set names to member gene IDs
    hallmark_enrich <- enricher(
      gene = as.character(entrez_ids),
      TERM2GENE = hallmark_sets[, c("gs_name", "entrez_gene")],
      pAdjustMethod = "BH",
      pvalueCutoff = opt$pvalueCutoff,
      qvalueCutoff = opt$qvalueCutoff,
      minGSSize = opt$minGSSize,
      maxGSSize = opt$maxGSSize
    )
    cat("[OK] MSigDB Hallmark enrichment completed\n")
  }, error = function(e) {
    cat(sprintf("[WARN] MSigDB Hallmark enrichment failed: %s\n", e$message))
    cat("      Continuing with GO enrichment only\n")
  })
} else {
  cat("[INFO] msigdbr not available, using GO enrichment only\n")
}

# -----------------------------------------------------------------------------
# Reactome Pathway Enrichment (Optional)
# -----------------------------------------------------------------------------
# ReactomePA provides enrichment analysis against the Reactome pathway
# database. Reactome offers mechanistic pathway representations that
# complement the functional annotations in GO.
#
# The \texttt{enrichPathway()} function is analogous to \texttt{enrichGO()}
# but queries the Reactome database instead of GO.
#
# \paragraph{Reactome-Specific Considerations:}
# \begin{itemize}
#   \item Pathways are hierarchically organised; related pathways may
#         show correlated enrichment
#   \item Signalling pathways are particularly well-annotated
#   \item The database is updated quarterly with peer review
# \end{itemize}

reactome_enrich <- NULL
if (has_reactome) {
  cat("[INFO] Running Reactome pathway enrichment...\n")
  reactome_enrich <- enrichPathway(
    gene = entrez_ids,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = opt$pvalueCutoff,
    qvalueCutoff = opt$qvalueCutoff,
    minGSSize = opt$minGSSize,
    maxGSSize = opt$maxGSSize,
    readable = TRUE
  )
} else {
  cat("[INFO] ReactomePA not available, skipping Reactome enrichment\n")
}

# -----------------------------------------------------------------------------
# Gene Ontology Results Export and Summary
# -----------------------------------------------------------------------------
# Enrichment results are exported as tab-separated files for downstream
# analysis and visualisation. The script also prints a summary of top
# enriched terms to the console for immediate interpretation.
#
# \paragraph{Interpreting Results:}
# \begin{itemize}
#   \item Terms with FDR $< 0.05$ are conventionally considered significant
#   \item Count indicates the number of query genes in the term; higher
#         counts provide more confidence
#   \item GeneRatio (query genes / term size) indicates specificity
#   \item Redundant terms (parent-child relationships) may appear; consider
#         using \texttt{simplify()} for downstream visualisation
# \end{itemize}

if (!is.null(go_enrich) && nrow(go_enrich@result) > 0) {
  go_file <- file.path(opt$outdir, paste0(opt$prefix, "_GO.tsv"))
  write_tsv(go_enrich@result, go_file)
  cat(sprintf("[OK] GO enrichment: %d terms (FDR<%.2f)\n",
              sum(go_enrich@result$p.adjust < opt$pvalueCutoff),
              opt$pvalueCutoff))
  cat(sprintf("     Written to: %s\n", go_file))

  # Display top 10 GO terms for immediate interpretation
  go_top <- head(go_enrich@result[go_enrich@result$p.adjust < opt$pvalueCutoff, ], 20)
  if (nrow(go_top) > 0) {
    cat("\nTop GO terms:\n")
    for (i in 1:min(10, nrow(go_top))) {
      cat(sprintf("  %d. %s [%s] (FDR=%.2e, Count=%d)\n",
                  i, go_top$Description[i], go_top$ONTOLOGY[i],
                  go_top$p.adjust[i], go_top$Count[i]))
    }
  }
} else {
  cat("[WARN] No significant GO terms found\n")
}

# -----------------------------------------------------------------------------
# Hallmark Results Export and Summary
# -----------------------------------------------------------------------------
# Hallmark results provide a high-level overview of biological themes
# represented in the gene list. The 50 Hallmark sets are designed to be
# non-redundant, making interpretation more straightforward than GO.
#
# \paragraph{Cancer-Relevant Hallmark Sets:}
# \begin{itemize}
#   \item Proliferation: E2F\_TARGETS, MYC\_TARGETS, G2M\_CHECKPOINT
#   \item Metabolism: OXIDATIVE\_PHOSPHORYLATION, GLYCOLYSIS, FATTY\_ACID\_METABOLISM
#   \item Immune: INFLAMMATORY\_RESPONSE, INTERFERON\_ALPHA/GAMMA\_RESPONSE
#   \item Signalling: KRAS\_SIGNALING, NOTCH\_SIGNALING, WNT\_BETA\_CATENIN
#   \item Stress: HYPOXIA, UNFOLDED\_PROTEIN\_RESPONSE, DNA\_REPAIR
# \end{itemize}

if (!is.null(hallmark_enrich) && nrow(hallmark_enrich@result) > 0) {
  hallmark_file <- file.path(opt$outdir, paste0(opt$prefix, "_Hallmark.tsv"))
  write_tsv(hallmark_enrich@result, hallmark_file)
  cat(sprintf("\n[OK] Hallmark enrichment: %d pathways (FDR<%.2f)\n",
              sum(hallmark_enrich@result$p.adjust < opt$pvalueCutoff),
              opt$pvalueCutoff))
  cat(sprintf("     Written to: %s\n", hallmark_file))

  # Display top Hallmark pathways
  hallmark_top <- head(hallmark_enrich@result[hallmark_enrich@result$p.adjust < opt$pvalueCutoff, ], 20)
  if (nrow(hallmark_top) > 0) {
    cat("\nTop Hallmark pathways:\n")
    for (i in 1:min(10, nrow(hallmark_top))) {
      cat(sprintf("  %d. %s (FDR=%.2e, Count=%d)\n",
                  i, hallmark_top$Description[i],
                  hallmark_top$p.adjust[i], hallmark_top$Count[i]))
    }
  }
}

# -----------------------------------------------------------------------------
# Reactome Results Export and Summary
# -----------------------------------------------------------------------------
# Reactome results complement GO and Hallmark by providing mechanistic
# pathway context. Enriched Reactome pathways may suggest specific
# molecular mechanisms or therapeutic targets.

if (!is.null(reactome_enrich) && nrow(reactome_enrich@result) > 0) {
  reactome_file <- file.path(opt$outdir, paste0(opt$prefix, "_Reactome.tsv"))
  write_tsv(reactome_enrich@result, reactome_file)
  cat(sprintf("\n[OK] Reactome enrichment: %d pathways (FDR<%.2f)\n",
              sum(reactome_enrich@result$p.adjust < opt$pvalueCutoff),
              opt$pvalueCutoff))
  cat(sprintf("     Written to: %s\n", reactome_file))

  # Display top Reactome pathways
  reactome_top <- head(reactome_enrich@result[reactome_enrich@result$p.adjust < opt$pvalueCutoff, ], 20)
  if (nrow(reactome_top) > 0) {
    cat("\nTop Reactome pathways:\n")
    for (i in 1:min(10, nrow(reactome_top))) {
      cat(sprintf("  %d. %s (FDR=%.2e, Count=%d)\n",
                  i, reactome_top$Description[i],
                  reactome_top$p.adjust[i], reactome_top$Count[i]))
    }
  }
} else {
  cat("[WARN] No significant Reactome pathways found\n")
}

# -----------------------------------------------------------------------------
# Execution Summary
# -----------------------------------------------------------------------------

cat("\n[OK] Enrichment analysis complete!\n")