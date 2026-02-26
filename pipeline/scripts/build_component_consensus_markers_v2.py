#!/usr/bin/env python3
"""
=============================================================================
build_component_consensus_markers.py
=============================================================================

DESCRIPTION
-----------
This script aggregates differential expression results from multiple anchor
samples to derive consensus marker genes for each graph component. By
requiring that marker genes show consistent differential expression across
multiple independent anchors, the analysis identifies robust transcriptional
signatures that characterise each component's biological identity.


BIOLOGICAL CONTEXT: CONSENSUS MARKER IDENTIFICATION
----------------------------------------------------
In graph-based clustering of cancer transcriptomes, each component represents
a group of samples with similar expression profiles. Characterising these
components requires identification of genes that distinguish component
members from non-members. However, single-sample differential expression
analyses are noisy and may identify spurious markers.

The consensus approach addresses this limitation by:

  (1) Using Multiple Anchors:
      Each component contains multiple samples (anchors) that can serve as
      representatives. By performing differential expression analysis for
      each anchor against samples outside the component, the method generates
      multiple independent assessments of component-defining genes.

  (2) Requiring Reproducibility:
      Genes must show consistent differential expression across multiple
      anchors to be considered robust markers. This requirement filters out
      sample-specific effects and technical noise, retaining only genes that
      reliably distinguish the component.

  (3) Separating Core and Extended Markers:
      Genes supported by many anchors (core markers) represent the most
      confident component signatures, while genes supported by fewer anchors
      (extended markers) may capture subtype-specific or context-dependent
      expression patterns.


DIFFERENTIAL EXPRESSION CONTEXT
-------------------------------
The input files are DESeq2 results comparing each anchor sample against
samples outside its component. DESeq2 performs negative binomial modelling
of RNA-seq count data, accounting for library size differences and
biological variability to estimate log2 fold changes and statistical
significance.

Key DESeq2 Output Columns:
  - gene_id:        Gene identifier (typically Ensembl ID with version)
  - baseMean:       Mean normalised count across all samples
  - log2FoldChange: Estimated effect size (positive = upregulated in anchor)
  - pvalue:         Raw p-value from the Wald test
  - padj:           Benjamini-Hochberg adjusted p-value (FDR)


CONSENSUS AGGREGATION STRATEGY
------------------------------
For each component, the script processes all anchor-specific DESeq2 results:

  Step 1: Filter Significant Genes
  Each anchor's results are filtered to retain genes meeting significance
  thresholds (adjusted p-value, minimum expression, fold change magnitude).
  Genes are separated into upregulated (log2FC > 0) and downregulated
  (log2FC < 0) sets.

  Step 2: Fallback for Low-Signal Anchors
  If an anchor yields too few significant genes (below --min_anchor_genes),
  a fallback strategy retains the top N genes ranked by significance and
  effect size. This ensures all anchors contribute to consensus building
  even when stringent thresholds would exclude them.

  Step 3: Remove Flip Genes
  Genes that appear as upregulated in some anchors but downregulated in
  others within the same component are removed. These "flip genes" likely
  reflect anchor-specific rather than component-defining expression patterns.

  Step 4: Aggregate Across Anchors
  For each gene appearing in any anchor's filtered results, the script
  computes:
    - support_n:       Number of anchors showing this gene as significant
    - median_log2FC:   Median fold change across supporting anchors
    - median_baseMean: Median expression level
    - best_padj:       Most significant adjusted p-value
    - rank_score:      Composite score (support_n * median_abs_log2FC)

  Step 5: Stratify by Support Level
  Genes are classified as core markers (supported by >= threshold anchors)
  or extended markers (supported by fewer anchors). Core markers represent
  the most robust component signatures.


QUALITY CONTROL: JACCARD SIMILARITY
-----------------------------------
The script computes pairwise Jaccard similarity between anchor gene sets
to assess reproducibility. Jaccard similarity measures the overlap between
two sets:

  J(A, B) = |A intersection B| / |A union B|

Values range from 0 (no overlap) to 1 (identical sets). High Jaccard
similarity between anchors from the same component indicates consistent
identification of component markers, validating the biological coherence
of the clustering.


INPUT REQUIREMENTS
------------------
The script expects a directory structure where each disease has DESeq2
results organised as:

  {root}/{disease}/deseq2_markers/tables/anchor_{name}_vs_outside_component_{N}.tsv

Each TSV file must contain columns: gene_id, baseMean, log2FoldChange,
pvalue, padj.


OUTPUT FILES
------------
Per-Component Outputs (in {outdir}/{disease}/component_{N}/):

  {anchor}.UP.sig.tsv / {anchor}.DOWN.sig.tsv
      Filtered significant genes for each anchor, separated by direction.

  component_consensus.UP.all.tsv / component_consensus.DOWN.all.tsv
      All consensus genes with aggregation statistics.

  component_consensus.UP.core.tsv / component_consensus.DOWN.core.tsv
      Core markers meeting minimum support threshold.

  component_consensus.UP.extended.tsv / component_consensus.DOWN.extended.tsv
      Extended markers below core support threshold.

Global Outputs (in {outdir}/):

  anchor_agreement_qc.tsv
      Pairwise Jaccard similarity between anchors for quality control.

  all_anchor_sig_long.tsv
      Long-format table of all significant genes across all anchors.


USAGE EXAMPLE
-------------
  python build_component_consensus_markers.py \
    --root results/unsupervised \
    --outdir results/consensus_markers \
    --padj 0.05 \
    --lfc 0.5 \
    --core_support_min 2


DEPENDENCIES
------------
Python packages: pandas, numpy, argparse, pathlib, re

=============================================================================
"""

import argparse
import re
from pathlib import Path
import pandas as pd
import numpy as np


# =============================================================================
# CONSTANTS AND PATTERNS
# =============================================================================
# Regular expression pattern for parsing anchor DESeq2 result filenames.
# The pattern extracts the anchor sample name and component number from
# filenames following the convention:
#   anchor_{sample_name}_vs_outside_component_{component_number}.tsv
#
# Capturing groups:
#   Group 1: Anchor sample name (any characters, non-greedy)
#   Group 2: Component number (digits)
# -----------------------------------------------------------------------------

ANCHOR_RE = re.compile(r"^anchor_(.+?)_vs_outside_component_(\d+)\.tsv$")


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def strip_ensg_version(gene_id: str) -> str:
    """
    Remove version suffix from Ensembl gene identifiers.
    
    Ensembl gene IDs include version numbers (e.g., ENSG00000141510.16) that
    change when gene annotations are updated. For cross-dataset comparisons
    and database lookups, the version-stripped ID (ENSG00000141510) provides
    more stable matching.
    
    Parameters
    ----------
    gene_id : str
        Ensembl gene identifier, potentially with version suffix.
    
    Returns
    -------
    str
        Gene identifier with version suffix removed.
    
    Examples
    --------
    >>> strip_ensg_version("ENSG00000141510.16")
    'ENSG00000141510'
    >>> strip_ensg_version("ENSG00000141510")
    'ENSG00000141510'
    """
    if pd.isna(gene_id):
        return gene_id
    return str(gene_id).split(".", 1)[0]


# =============================================================================
# DATA LOADING AND PREPROCESSING
# =============================================================================

def read_deseq2_table(path: Path) -> pd.DataFrame:
    """
    Load and preprocess a DESeq2 differential expression results table.
    
    This function performs several preprocessing steps to prepare DESeq2
    output for downstream consensus analysis:
    
      1. Validates that required columns are present
      2. Strips version suffixes from Ensembl gene IDs
      3. Converts numeric columns to appropriate types
      4. Removes genes with zero expression or missing fold changes
    
    The function preserves the original gene ID in a separate column
    (gene_id_raw) to enable tracing back to source data when needed.
    
    Parameters
    ----------
    path : Path
        Path to the DESeq2 results TSV file.
    
    Returns
    -------
    pd.DataFrame
        Preprocessed DESeq2 results with columns:
        - gene_id:        Version-stripped Ensembl ID
        - gene_id_raw:    Original gene ID with version
        - baseMean:       Mean normalised expression
        - log2FoldChange: Effect size estimate
        - pvalue:         Raw p-value
        - padj:           Adjusted p-value (FDR)
    
    Raises
    ------
    ValueError
        If required columns are missing from the input file.
    
    Notes
    -----
    Genes with baseMean = 0 or missing log2FoldChange are excluded as they
    represent untested or problematic cases in the DESeq2 analysis.
    """
    df = pd.read_csv(path, sep="\t")
    
    # Validate required columns are present
    required = ["gene_id", "baseMean", "log2FoldChange", "pvalue", "padj"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")

    df = df.copy()
    
    # Preserve original gene ID before stripping version
    df["gene_id_raw"] = df["gene_id"]
    df["gene_id"] = df["gene_id"].map(strip_ensg_version)

    # Ensure numeric columns have correct types
    # errors="coerce" converts non-numeric values to NaN
    for c in ["baseMean", "log2FoldChange", "pvalue", "padj"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Remove genes that were not meaningfully tested:
    # - baseMean = 0 indicates no detectable expression
    # - Missing log2FoldChange indicates the gene was filtered by DESeq2
    df = df[df["baseMean"].fillna(0) > 0].copy()
    df = df[~df["log2FoldChange"].isna()].copy()
    
    return df


# =============================================================================
# DIFFERENTIAL EXPRESSION FILTERING
# =============================================================================

def split_up_down(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Separate genes into upregulated and downregulated sets.
    
    In differential expression analysis, the sign of log2FoldChange indicates
    the direction of expression change:
      - Positive log2FC: Gene is upregulated in the test condition (anchor)
        relative to the reference (outside component)
      - Negative log2FC: Gene is downregulated in anchor relative to outside
    
    Separating by direction enables independent consensus analysis for
    upregulated and downregulated markers, which may have distinct biological
    interpretations.
    
    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results table with log2FoldChange column.
    
    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Tuple of (upregulated_genes, downregulated_genes) DataFrames.
    """
    up = df[df["log2FoldChange"] > 0].copy()
    down = df[df["log2FoldChange"] < 0].copy()
    return up, down


def filter_sig(df: pd.DataFrame, padj_max: float, baseMean_min: float, 
               lfc_min: float) -> pd.DataFrame:
    """
    Filter genes to retain only statistically significant results.
    
    This function applies three filtering criteria commonly used in
    differential expression analysis:
    
      1. Statistical Significance (padj):
         Adjusted p-value threshold controls the false discovery rate.
         A threshold of 0.05 means accepting up to 5% false positives
         among called significant genes.
    
      2. Expression Level (baseMean):
         Minimum expression threshold excludes lowly-expressed genes
         where fold change estimates are unreliable. A threshold of 1.0
         is conservative; some analyses use 10 or higher.
    
      3. Effect Size (log2FoldChange):
         Minimum fold change threshold ensures biological relevance.
         A threshold of 0.5 corresponds to approximately 1.4-fold change,
         while 1.0 corresponds to 2-fold change.
    
    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results table (typically already split by direction).
    padj_max : float
        Maximum adjusted p-value (e.g., 0.05 for 5% FDR).
    baseMean_min : float
        Minimum mean expression level.
    lfc_min : float
        Minimum absolute log2 fold change.
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing only significant genes.
    """
    return df[
        (df["padj"].notna()) &
        (df["padj"] <= padj_max) &
        (df["baseMean"] >= baseMean_min) &
        (df["log2FoldChange"].abs() >= lfc_min)
    ].copy()


def topn_fallback(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """
    Fallback strategy when stringent filtering yields too few genes.
    
    Some anchors may yield very few significant genes under standard
    thresholds due to low statistical power, high biological variability,
    or genuine lack of differential expression. To ensure all anchors
    contribute to consensus building, this fallback retains the top N
    genes ranked by statistical evidence.
    
    Ranking Strategy:
    Genes are sorted by:
      1. Adjusted p-value (ascending - most significant first)
      2. Absolute log2 fold change (descending - largest effects first)
    
    This prioritises genes with the strongest statistical support while
    using effect size as a tiebreaker.
    
    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results table (typically one direction: up or down).
    n : int
        Number of top genes to retain.
    
    Returns
    -------
    pd.DataFrame
        Top N genes by significance and effect size.
    
    Notes
    -----
    The fallback is applied separately to upregulated and downregulated
    genes to maintain directional balance in the consensus analysis.
    """
    # Restrict to genes with valid adjusted p-values (actually tested)
    df2 = df[df["padj"].notna()].copy()
    if df2.empty:
        return df2
    
    # Create absolute fold change for ranking
    df2["abs_lfc"] = df2["log2FoldChange"].abs()
    
    # Sort by significance (primary) and effect size (secondary)
    df2 = df2.sort_values(
        ["padj", "abs_lfc"], 
        ascending=[True, False]
    ).head(n).drop(columns=["abs_lfc"])
    
    return df2


# =============================================================================
# CONSENSUS AGGREGATION
# =============================================================================

def aggregate_consensus(anchor_dfs: list[pd.DataFrame], direction: str) -> pd.DataFrame:
    """
    Aggregate differential expression results across multiple anchors.
    
    This function combines filtered gene lists from multiple anchor samples
    to identify consensus markers. For each gene appearing in any anchor's
    results, the function computes summary statistics that characterise
    the gene's behaviour across anchors.
    
    Aggregation Statistics:
    
      support_n:
          Number of anchors where this gene was significant. Higher support
          indicates more reproducible differential expression.
    
      median_log2FC:
          Median fold change across supporting anchors. The median is robust
          to outliers from individual anchors with unusual effect sizes.
    
      median_abs_log2FC:
          Median absolute fold change, used for ranking regardless of
          direction (since input is already split by direction).
    
      median_baseMean:
          Median expression level, indicating whether the gene is lowly
          or highly expressed.
    
      best_padj:
          Minimum adjusted p-value across anchors, representing the
          strongest statistical evidence for this gene.
    
      rank_score:
          Composite score combining support and effect size:
          rank_score = support_n * median_abs_log2FC
          This prioritises genes that are both reproducible and have
          large expression differences.
    
    Parameters
    ----------
    anchor_dfs : list[pd.DataFrame]
        List of filtered DESeq2 results DataFrames, one per anchor.
    direction : str
        Direction label ("UP" or "DOWN") for output annotation.
    
    Returns
    -------
    pd.DataFrame
        Consensus table with one row per gene and aggregation statistics,
        sorted by support_n (primary), rank_score (secondary), and
        best_padj (tertiary).
    """
    if not anchor_dfs:
        return pd.DataFrame()

    # Combine all anchor results into long format with anchor index
    long = []
    for i, df in enumerate(anchor_dfs, start=1):
        tmp = df[["gene_id", "baseMean", "log2FoldChange", "padj"]].copy()
        tmp["anchor_idx"] = i
        long.append(tmp)
    long = pd.concat(long, ignore_index=True)

    # Aggregate statistics per gene across all anchors
    g = long.groupby("gene_id", as_index=False).agg(
        # Number of unique anchors supporting this gene
        support_n=("anchor_idx", "nunique"),
        # Median fold change (preserves direction)
        median_log2FC=("log2FoldChange", "median"),
        # Median absolute fold change (for ranking)
        median_abs_log2FC=("log2FoldChange", lambda x: float(np.median(np.abs(x)))),
        # Median expression level
        median_baseMean=("baseMean", "median"),
        # Best (minimum) adjusted p-value
        best_padj=("padj", "min"),
        # Total observations (for debugging)
        n_obs=("log2FoldChange", "size"),
    )
    
    # Compute composite ranking score
    g["rank_score"] = g["support_n"] * g["median_abs_log2FC"]
    
    # Add direction annotation
    g["direction"] = direction
    
    # Sort by support (most reproducible), then rank score, then significance
    g = g.sort_values(
        ["support_n", "rank_score", "best_padj"], 
        ascending=[False, False, True]
    ).reset_index(drop=True)
    
    return g


# =============================================================================
# QUALITY CONTROL METRICS
# =============================================================================

def jaccard(a: set, b: set) -> float:
    """
    Compute Jaccard similarity coefficient between two sets.
    
    The Jaccard index measures the similarity between two sets as the
    size of their intersection divided by the size of their union:
    
      J(A, B) = |A ∩ B| / |A ∪ B|
    
    Values range from 0 (completely disjoint) to 1 (identical sets).
    
    In the context of differential expression analysis, Jaccard similarity
    between anchor gene sets indicates reproducibility:
      - High Jaccard (> 0.3): Good agreement between anchors
      - Moderate Jaccard (0.1-0.3): Partial overlap, some shared markers
      - Low Jaccard (< 0.1): Poor agreement, few shared markers
    
    Parameters
    ----------
    a : set
        First set of gene identifiers.
    b : set
        Second set of gene identifiers.
    
    Returns
    -------
    float
        Jaccard similarity coefficient in range [0, 1].
    
    Notes
    -----
    Edge cases:
      - Both sets empty: Returns 1.0 (trivial identity)
      - One set empty: Returns 0.0 (no possible overlap)
    """
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def comp_sort_key(c: str):
    """
    Generate sort key for component identifiers.
    
    Component identifiers may be numeric strings ("1", "2", "10") or
    alphanumeric. This function enables natural numeric sorting where
    "2" comes before "10", while handling non-numeric identifiers by
    placing them at the end.
    
    Parameters
    ----------
    c : str
        Component identifier string.
    
    Returns
    -------
    int
        Integer sort key (parsed value for numeric IDs, large value for others).
    """
    try:
        return int(c)
    except ValueError:
        return 10**9


# =============================================================================
# MAIN ANALYSIS WORKFLOW
# =============================================================================

def main():
    """
    Main entry point for consensus marker identification.
    
    This function orchestrates the complete analysis workflow:
    
      1. Parse command-line arguments for threshold configuration
      2. Discover disease directories containing DESeq2 results
      3. For each component in each disease:
         a. Load all anchor DESeq2 result files
         b. Filter significant genes (with fallback for low-signal anchors)
         c. Remove flip genes showing inconsistent direction
         d. Aggregate consensus statistics across anchors
         e. Stratify into core and extended marker sets
         f. Save per-anchor and consensus output files
      4. Compute pairwise Jaccard similarity for quality control
      5. Compile global summary tables
    """
    
    # -------------------------------------------------------------------------
    # Argument Parsing
    # Define command-line interface with configurable thresholds for
    # significance filtering and consensus building.
    # -------------------------------------------------------------------------
    
    ap = argparse.ArgumentParser(
        description="Build component consensus markers from DESeq2 anchor tables."
    )
    ap.add_argument("--root", required=True, 
                    help="Pipeline results root, e.g. results/unsupervised")
    ap.add_argument("--outdir", required=True, 
                    help="Output directory for consensus results")
    ap.add_argument("--padj", type=float, default=0.05, 
                    help="Adjusted p-value threshold (default: 0.05)")
    ap.add_argument("--basemean", type=float, default=1.0, 
                    help="Minimum baseMean expression threshold (default: 1.0)")
    ap.add_argument("--lfc", type=float, default=0.5, 
                    help="Minimum absolute log2 fold change (default: 0.5)")
    ap.add_argument("--core_support_min", type=int, default=2, 
                    help="Minimum anchor support for core markers (default: 2)")
    ap.add_argument("--min_anchor_genes", type=int, default=25, 
                    help="Minimum significant genes before fallback (default: 25)")
    ap.add_argument("--fallback_topn", type=int, default=200, 
                    help="Top N genes to retain in fallback mode (default: 200)")
    args = ap.parse_args()

    root = Path(args.root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # Disease Directory Discovery
    # Scan the root directory for subdirectories containing DESeq2 marker
    # results, indicating diseases with completed differential expression
    # analyses.
    # -------------------------------------------------------------------------
    
    diseases = []
    for d in root.iterdir():
        if d.is_dir() and (d / "deseq2_markers").exists():
            diseases.append(d.name)
    
    if not diseases:
        raise SystemExit(f"No disease directories found under {root} containing deseq2_markers/")

    # Accumulators for global outputs
    qc_rows = []          # Pairwise Jaccard similarity records
    all_sig_long = []     # All significant genes in long format

    # -------------------------------------------------------------------------
    # Main Processing Loop: Iterate Over Diseases
    # -------------------------------------------------------------------------
    
    for disease in sorted(diseases):
        tables_dir = root / disease / "deseq2_markers" / "tables"
        if not tables_dir.exists():
            print(f"[WARN] Missing tables dir for {disease}: {tables_dir}")
            continue

        # ---------------------------------------------------------------------
        # Discover Anchor Files and Organise by Component
        # Parse filenames to extract anchor names and component assignments,
        # grouping files by component for joint processing.
        # ---------------------------------------------------------------------
        
        anchors_by_comp: dict[str, list[tuple[str, Path]]] = {}
        for f in tables_dir.glob("anchor_*_vs_outside_component_*.tsv"):
            m = ANCHOR_RE.match(f.name)
            if not m:
                continue
            anchor_name, comp = m.group(1), m.group(2)
            anchors_by_comp.setdefault(comp, []).append((anchor_name, f))

        # ---------------------------------------------------------------------
        # Process Each Component
        # For each component, load all anchor results, filter significant
        # genes, remove inconsistent genes, and aggregate consensus markers.
        # ---------------------------------------------------------------------
        
        for comp, anchor_files in sorted(anchors_by_comp.items(), 
                                          key=lambda kv: comp_sort_key(kv[0])):
            # Create output directory for this component
            comp_out = outdir / disease / f"component_{comp}"
            comp_out.mkdir(parents=True, exist_ok=True)

            # Track gene sets for Jaccard computation and flip gene removal
            up_sets = {}      # anchor_name -> set of upregulated gene IDs
            down_sets = {}    # anchor_name -> set of downregulated gene IDs
            up_filtered = []  # List of filtered UP DataFrames
            down_filtered = [] # List of filtered DOWN DataFrames

            # -----------------------------------------------------------------
            # Process Each Anchor
            # Load DESeq2 results, apply significance filters (with fallback),
            # and track gene sets for downstream aggregation.
            # -----------------------------------------------------------------
            
            for anchor_name, path in sorted(anchor_files, key=lambda x: x[0]):
                # Load and preprocess DESeq2 results
                df = read_deseq2_table(path)
                
                # Separate by expression direction
                up, down = split_up_down(df)

                # Apply significance filters
                up_sig = filter_sig(up, args.padj, args.basemean, args.lfc)
                down_sig = filter_sig(down, args.padj, args.basemean, args.lfc)

                # Apply fallback if too few genes pass filters
                # This ensures low-signal anchors still contribute to consensus
                fallback_used = False
                if len(up_sig) < args.min_anchor_genes:
                    up_sig = topn_fallback(up, args.fallback_topn)
                    fallback_used = True
                if len(down_sig) < args.min_anchor_genes:
                    down_sig = topn_fallback(down, args.fallback_topn)
                    fallback_used = True

                # Store gene sets for Jaccard computation
                up_sets[anchor_name] = set(up_sig["gene_id"].tolist())
                down_sets[anchor_name] = set(down_sig["gene_id"].tolist())

                # Add metadata columns for tracking
                up_sig = up_sig.assign(
                    anchor=anchor_name, 
                    component=comp, 
                    disease=disease, 
                    direction="UP", 
                    fallback_used=fallback_used
                )
                down_sig = down_sig.assign(
                    anchor=anchor_name, 
                    component=comp, 
                    disease=disease, 
                    direction="DOWN", 
                    fallback_used=fallback_used
                )

                # Accumulate for consensus aggregation
                up_filtered.append(up_sig)
                down_filtered.append(down_sig)

                # Accumulate for global long-format output
                all_sig_long.append(up_sig[[
                    "disease", "component", "anchor", "direction", "gene_id",
                    "gene_id_raw", "baseMean", "log2FoldChange", "padj", "fallback_used"
                ]])
                all_sig_long.append(down_sig[[
                    "disease", "component", "anchor", "direction", "gene_id",
                    "gene_id_raw", "baseMean", "log2FoldChange", "padj", "fallback_used"
                ]])

                # Save per-anchor filtered results
                up_sig.to_csv(comp_out / f"{anchor_name}.UP.sig.tsv", sep="\t", index=False)
                down_sig.to_csv(comp_out / f"{anchor_name}.DOWN.sig.tsv", sep="\t", index=False)

            # -----------------------------------------------------------------
            # Remove Flip Genes
            # Genes appearing in both UP and DOWN sets across different anchors
            # within the same component show inconsistent direction and are
            # excluded from consensus analysis. These may represent:
            #   - Genes with anchor-specific regulation
            #   - Genes near the fold change threshold fluctuating by chance
            #   - Technical artefacts in specific samples
            # -----------------------------------------------------------------
            
            # Identify genes appearing in both directions across any anchors
            all_up_genes = set().union(*up_sets.values()) if up_sets else set()
            all_down_genes = set().union(*down_sets.values()) if down_sets else set()
            flip_genes = all_up_genes & all_down_genes
            
            if flip_genes:
                # Remove flip genes from all gene sets
                for k in up_sets:
                    up_sets[k] -= flip_genes
                for k in down_sets:
                    down_sets[k] -= flip_genes
                
                # Remove flip genes from filtered DataFrames
                up_filtered = [df[~df["gene_id"].isin(flip_genes)].copy() 
                              for df in up_filtered]
                down_filtered = [df[~df["gene_id"].isin(flip_genes)].copy() 
                                for df in down_filtered]

            # -----------------------------------------------------------------
            # Aggregate Consensus Markers
            # Combine results across anchors to identify reproducibly
            # differential genes with support from multiple anchors.
            # -----------------------------------------------------------------
            
            cons_up = aggregate_consensus(up_filtered, direction="UP")
            cons_down = aggregate_consensus(down_filtered, direction="DOWN")

            # -----------------------------------------------------------------
            # Stratify by Support Level
            # Core markers: Supported by >= core_support_min anchors
            # Extended markers: Supported by fewer anchors
            #
            # This stratification enables downstream analyses to focus on
            # the most robust markers (core) while retaining potentially
            # interesting genes with lower support (extended).
            # -----------------------------------------------------------------
            
            core_up = cons_up[cons_up["support_n"] >= args.core_support_min].copy()
            core_down = cons_down[cons_down["support_n"] >= args.core_support_min].copy()
            ext_up = cons_up[cons_up["support_n"] < args.core_support_min].copy()
            ext_down = cons_down[cons_down["support_n"] < args.core_support_min].copy()

            # Save consensus outputs
            cons_up.to_csv(comp_out / "component_consensus.UP.all.tsv", sep="\t", index=False)
            cons_down.to_csv(comp_out / "component_consensus.DOWN.all.tsv", sep="\t", index=False)
            core_up.to_csv(comp_out / "component_consensus.UP.core.tsv", sep="\t", index=False)
            core_down.to_csv(comp_out / "component_consensus.DOWN.core.tsv", sep="\t", index=False)
            ext_up.to_csv(comp_out / "component_consensus.UP.extended.tsv", sep="\t", index=False)
            ext_down.to_csv(comp_out / "component_consensus.DOWN.extended.tsv", sep="\t", index=False)

            # -----------------------------------------------------------------
            # Compute Pairwise Jaccard Similarity
            # For each pair of anchors within the component, compute the
            # Jaccard similarity of their significant gene sets. This
            # quality control metric indicates reproducibility.
            # -----------------------------------------------------------------
            
            anchors = sorted(up_sets.keys())
            for i in range(len(anchors)):
                for j in range(i + 1, len(anchors)):
                    a, b = anchors[i], anchors[j]
                    qc_rows.append({
                        "disease": disease,
                        "component": comp,
                        "anchor_a": a,
                        "anchor_b": b,
                        "jaccard_up": jaccard(up_sets[a], up_sets[b]),
                        "jaccard_down": jaccard(down_sets[a], down_sets[b]),
                        "n_up_a": len(up_sets[a]),
                        "n_up_b": len(up_sets[b]),
                        "n_down_a": len(down_sets[a]),
                        "n_down_b": len(down_sets[b]),
                        "n_flip_genes": len(flip_genes),
                    })

    # -------------------------------------------------------------------------
    # Write Global Output Files
    # -------------------------------------------------------------------------
    
    # Quality control table with pairwise anchor agreement
    qc = pd.DataFrame(qc_rows)
    qc_out = outdir / "anchor_agreement_qc.tsv"
    qc.to_csv(qc_out, sep="\t", index=False)

    # Long-format table of all significant genes across all anchors
    if all_sig_long:
        all_long = pd.concat(all_sig_long, ignore_index=True)
        all_long.to_csv(outdir / "all_anchor_sig_long.tsv", sep="\t", index=False)

    print(f"[OK] Wrote QC: {qc_out}")
    print(f"[OK] Consensus outputs in: {outdir}")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    main()