#!/usr/bin/env python3
"""
confirm_nbl_gene_list.py

Confirm whether an input NBL (neuroblastoma) gene list is informative for:
- neuroblastoma biology (enrichment themes)
- neuroblastoma lineage markers (ADR vs MES, sympathetic neuron/crest programs)

Works with gprofiler-official==1.0.0 (uses GProfiler.profile()).

Install deps:
  pip install gprofiler-official pandas scipy
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd
from scipy.stats import fisher_exact

from gprofiler import GProfiler


# -----------------------------
# Curated NBL marker panels
# -----------------------------
# NOTE: These are pragmatic, literature-inspired panels used as "sanity checks".
# Replace/extend with your preferred signatures (e.g., from a specific paper or dataset).

# Core neuroblastoma markers including oncogenes, lineage factors, and differentiation markers
NBL_CORE: Set[str] = {
    # commonly referenced NBL lineage / drivers / differentiation markers
    "MYCN", "ALK", "PHOX2B", "PHOX2A", "HAND2", "GATA3", "ISL1", "DBH",
    "TH", "DDC", "SLC6A2", "SLC18A1", "SLC18A2", "CHGA", "CHGB",
    "NCAM1", "NRCAM", "GD2A",  # GD2A is not a standard gene symbol; keep if your mapping uses it
    "BCL2", "LIN28B", "LMO1", "NTRK1", "NTRK2", "NTRK3",
    "TERT", "ATRX", "DUSP26"
}

# Adrenergic (ADR) state markers - sympathetic neuron-like phenotype with high MYCN/catecholamine pathway
NBL_ADR: Set[str] = {
    "PHOX2B", "PHOX2A", "HAND2", "GATA3", "ISL1", "ASCL1",
    "TH", "DDC", "DBH", "SLC6A2", "SLC18A1", "SLC18A2",
    "CHGA", "CHGB", "DLK1"
}

# Mesenchymal (MES) state markers - invasive, therapy-resistant phenotype with EMT features
NBL_MES: Set[str] = {
    "PRRX1", "FOSL1", "FOSL2", "VIM", "FN1", "ITGA5", "COL1A1", "COL1A2",
    "COL3A1", "CXCL8", "IL6", "TAGLN", "ACTA2", "SNAI2", "ZEB1", "ZEB2",
    "TGFBI"
}

# Neural crest / Schwann-like / developmental program markers - reflecting developmental origins
NBL_DEV: Set[str] = {
    "SOX10", "S100B", "MPZ", "PLP1", "PMP22", "ERBB3", "NRG1",
    "GFAP", "FABP7", "NES", "TUBB3", "MAP2"
}

# Proliferation markers to identify if gene list only captures cell cycle rather than NBL-specific biology
PROLIFERATION: Set[str] = {"MKI67", "TOP2A", "CCNB1", "CDC20", "BIRC5", "AURKA", "AURKB", "UBE2C", "MELK"}

# Dictionary organizing all marker panels for easy iteration
PANELS: Dict[str, Set[str]] = {
    "NBL_core_markers": NBL_CORE,
    "NBL_ADR_program": NBL_ADR,
    "NBL_MES_program": NBL_MES,
    "NBL_development_program": NBL_DEV,
    "Proliferation_program": PROLIFERATION,
}

# g:Profiler pathway databases to query for functional enrichment analysis
GPROF_SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"]


def read_gene_list(path: str) -> List[str]:
    """
    Read gene list from a text file, expecting one gene ID per line.
    
    Handles:
    - Comment lines starting with #
    - Whitespace stripping
    - Ensembl version suffix removal (e.g., ENSG00000000003.15 -> ENSG00000000003)
    - Deduplication while preserving order
    
    Args:
        path: File path to gene list
        
    Returns:
        List of unique gene IDs
    """
    genes: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            # Skip empty lines and comments
            if not g or g.startswith("#"):
                continue
            g = g.split()[0]  # Take only first token if multiple columns
            g = g.split(".")[0]  # Strip Ensembl version suffix
            genes.append(g)
    
    # De-duplicate while preserving order
    seen = set()
    out = []
    for g in genes:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out


def ensure_outdir(outdir: str) -> Path:
    """
    Create output directory if it doesn't exist.
    
    Args:
        outdir: Path to output directory
        
    Returns:
        Path object for the directory
    """
    p = Path(outdir)
    p.mkdir(parents=True, exist_ok=True)
    return p


def convert_ensembl_to_symbols(gp: GProfiler, genes: List[str], organism: str) -> pd.DataFrame:
    """
    Convert Ensembl gene IDs to gene symbols using g:Profiler's convert endpoint.
    
    Args:
        gp: GProfiler instance
        genes: List of Ensembl gene IDs
        organism: Organism code (e.g., 'hsapiens')
        
    Returns:
        DataFrame with columns 'input_id' and 'symbol'
    """
    # Query g:Profiler conversion service
    conv = gp.convert(organism=organism, query=genes, target_namespace="ENSG")
    df = pd.DataFrame(conv)
    
    # Handle empty result
    if df.empty:
        return pd.DataFrame({"input_id": genes, "symbol": [None] * len(genes)})

    # Identify the incoming gene ID column (API may vary)
    incoming = "incoming" if "incoming" in df.columns else df.columns[0]
    
    # Find the symbol/name column using heuristics
    symbol_col = None
    for c in df.columns:
        if c.lower() in {"name", "converted_name", "symbol"}:
            symbol_col = c
            break
    if symbol_col is None:
        for c in df.columns:
            if c.lower().endswith("name"):
                symbol_col = c
                break

    # Standardize column names
    df = df.rename(columns={incoming: "input_id"})
    # Strip version suffix from input IDs
    df["input_id"] = df["input_id"].astype(str).str.split(".").str[0]

    # Handle case where no symbol column found
    if symbol_col is None:
        df["symbol"] = None
        return df[["input_id", "symbol"]].drop_duplicates()

    # Clean up symbol column
    df = df.rename(columns={symbol_col: "symbol"})
    df["symbol"] = df["symbol"].where(df["symbol"].notna(), None)
    # Convert string 'nan' to None
    df.loc[df["symbol"].astype(str).str.lower() == "nan", "symbol"] = None
    
    return df[["input_id", "symbol"]].drop_duplicates()


def fisher_enrichment(gene_symbols: Set[str], signature: Set[str], universe_size: int) -> Tuple[float, float, int]:
    """
    Perform one-sided Fisher exact test for gene set overlap enrichment.
    
    Tests whether the overlap between gene_symbols and signature is greater than expected by chance.
    
    Contingency table:
                     In signature    Not in signature
    In gene list         a                  c
    Not in gene list     b                  d
    
    Args:
        gene_symbols: Set of gene symbols from input list
        signature: Reference gene signature to test against
        universe_size: Total number of genes in the universe (background)
        
    Returns:
        Tuple of (odds_ratio, p_value, overlap_count)
    """
    overlap = len(gene_symbols & signature)
    a = overlap  # In both
    b = len(signature) - overlap  # In signature only
    c = len(gene_symbols) - overlap  # In gene list only
    d = universe_size - (a + b + c)  # In neither
    
    # Sanity check for valid contingency table
    if d < 0:
        return float("nan"), 1.0, overlap
    
    # One-sided test (greater) - is overlap enriched?
    oddsratio, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return oddsratio, p, overlap


def run_enrichment(gp: GProfiler, genes: List[str], organism: str, fdr: float, sources: List[str]) -> pd.DataFrame:
    """
    Run functional enrichment analysis using g:Profiler.
    
    Args:
        gp: GProfiler instance
        genes: List of gene IDs to analyze
        organism: Organism code
        fdr: FDR threshold for significance
        sources: List of pathway databases to query (e.g., GO:BP, KEGG)
        
    Returns:
        DataFrame of significant enriched terms, sorted by source and p-value
    """
    # Call g:Profiler API (gprofiler-official==1.0.0 uses profile() method)
    df = gp.profile(
        organism=organism,
        query=genes,
        sources=sources,
        user_threshold=fdr,
        significance_threshold_method="fdr",
        no_evidences=False,  # Include evidence genes
        ordered=False,  # Treat query as unordered gene set
    )
    df = pd.DataFrame(df)
    
    if df.empty:
        return df

    # Select relevant columns for output
    keep = [c for c in [
        "source", "native", "name", "p_value", "significant",
        "term_size", "query_size", "intersection_size",
        "effective_domain_size", "intersections"
    ] if c in df.columns]

    return df[keep].sort_values(["source", "p_value"], ascending=[True, True])


def main() -> None:
    """
    Main function to orchestrate NBL gene list validation workflow:
    
    1. Read input gene list (Ensembl IDs)
    2. Convert Ensembl IDs to gene symbols
    3. Test overlap with curated NBL marker panels using Fisher exact test
    4. Run functional enrichment analysis at multiple FDR thresholds
    5. Generate summary report
    """
    # Parse command-line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True, help="Text file: one Ensembl gene ID per line (ENSG...)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--organism", default="hsapiens", help="g:Profiler organism (default: hsapiens)")
    ap.add_argument("--fdrs", default="0.01,0.05", help="Comma-separated FDR thresholds to test")
    ap.add_argument("--universe_size", type=int, default=20000, help="Universe size for Fisher tests")
    args = ap.parse_args()

    # Set up output directory
    outdir = ensure_outdir(args.outdir)
    
    # Read and parse gene list
    genes = read_gene_list(args.genes)

    # Initialize g:Profiler client
    gp = GProfiler()

    # Convert Ensembl IDs to gene symbols
    mapping = convert_ensembl_to_symbols(gp, genes, args.organism)
    mapping.to_csv(outdir / "ensembl_to_symbol.tsv", sep="\t", index=False)

    # Extract valid gene symbols (non-null, non-empty)
    symbols = set([s for s in mapping["symbol"].dropna().astype(str) if s and s.lower() != "none"])

    # Test overlap with each curated NBL marker panel
    overlap_rows = []
    for panel_name, panel_genes in PANELS.items():
        # Run Fisher exact test
        odds, p, k = fisher_enrichment(symbols, panel_genes, args.universe_size)
        # Get list of overlapping genes
        overlap_list = sorted(symbols & panel_genes)
        overlap_rows.append({
            "panel": panel_name,
            "panel_size": len(panel_genes),
            "list_symbol_count": len(symbols),
            "overlap_count": k,
            "overlap_genes": ",".join(overlap_list),
            "fisher_oddsratio": odds,
            "fisher_p_greater": p,
        })

    # Save overlap statistics
    overlap_df = pd.DataFrame(overlap_rows).sort_values("fisher_p_greater", ascending=True)
    overlap_df.to_csv(outdir / "overlap_report.tsv", sep="\t", index=False)

    # Run functional enrichment at multiple FDR thresholds
    fdrs = [float(x.strip()) for x in args.fdrs.split(",") if x.strip()]
    enrich_summaries = []
    for fdr in fdrs:
        enr = run_enrichment(gp, genes, args.organism, fdr, GPROF_SOURCES)
        outpath = outdir / f"enrichment_FDR{fdr:.2f}.tsv"
        enr.to_csv(outpath, sep="\t", index=False)

        # Summarize enrichment results for this FDR
        if enr.empty:
            enrich_summaries.append(f"FDR {fdr:.2f}: no significant terms (sources={','.join(GPROF_SOURCES)})")
        else:
            counts = enr["source"].value_counts().to_dict() if "source" in enr.columns else {}
            enrich_summaries.append(f"FDR {fdr:.2f}: {len(enr)} terms; per-source: {counts}")

    # Generate human-readable summary report
    with open(outdir / "summary.txt", "w", encoding="utf-8") as f:
        f.write("NBL gene list confirmation summary\n")
        f.write("=================================\n\n")
        f.write(f"Input genes: {len(genes)} (unique)\n")
        f.write(f"Converted gene symbols: {len(symbols)}\n\n")

        f.write("Panel overlap (one-sided Fisher exact, greater):\n")
        # Show top 10 most significant panel overlaps
        top = overlap_df.head(10)
        for _, r in top.iterrows():
            f.write(
                f"  - {r['panel']}: overlap={r['overlap_count']}/{r['panel_size']} "
                f"(p={r['fisher_p_greater']:.3e}, OR={r['fisher_oddsratio']})\n"
            )

        f.write("\nEnrichment summaries:\n")
        for line in enrich_summaries:
            f.write(f"  - {line}\n")

        # Provide interpretation guidance
        f.write("\nNotes:\n")
        f.write("- Strong signal for NBL_ADR or NBL_MES supports state-specific informativeness.\n")
        f.write("- Strong NBL_core + developmental programs supports lineage specificity.\n")
        f.write("- If only Proliferation_program hits and enrichment is generic, the list may reflect general variability.\n")

    # Print completion message with output file listing
    print(f"[OK] Wrote outputs to: {outdir}")
    print(" - ensembl_to_symbol.tsv")
    print(" - overlap_report.tsv")
    for fdr in fdrs:
        print(f" - enrichment_FDR{fdr:.2f}.tsv")
    print(" - summary.txt")


if __name__ == "__main__":
    main()