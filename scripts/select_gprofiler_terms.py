#!/usr/bin/env python3
"""
Filter and select top enrichment terms from g:Profiler results.

This script processes g:Profiler output to select the most significant Gene Ontology (GO)
terms and pathway enrichments based on statistical significance and intersection size.
"""

import argparse
import pandas as pd
from typing import Optional, Tuple

# Recognized pathway databases in g:Profiler output
PATHWAY_SOURCES = ["KEGG", "REAC", "WP"]


def resolve_go_source_and_top(args) -> Tuple[str, int]:
    """
    Determine which GO namespace to use and how many terms to select.
    
    Handles backwards compatibility with legacy flags while supporting new unified interface.
    
    Priority order:
      1) --go {MF,BP,CC} with optional --top_go (new interface)
      2) --top_bp or --top_cc if provided (legacy specific namespace flags)
      3) Default: GO:MF with --top_mf value (legacy default)
    
    Args:
        args: Parsed command-line arguments
        
    Returns:
        Tuple of (GO source string like "GO:MF", number of top terms to select)
    """
    # Priority 1: New unified --go flag
    if args.go is not None:
        go = args.go.upper()
        # Use --top_go if provided, otherwise default to 10
        top = args.top_go if args.top_go is not None else 10
        return f"GO:{go}", int(top)
    
    # Priority 2: Legacy --top_bp flag
    if args.top_bp is not None:
        return "GO:BP", int(args.top_bp)
    
    # Priority 3: Legacy --top_cc flag
    if args.top_cc is not None:
        return "GO:CC", int(args.top_cc)
    
    # Default: GO:MF with --top_mf value
    return "GO:MF", int(args.top_mf)


def pick_terms(
    df: pd.DataFrame,
    fdr: float,
    min_is: int,
    go_source: str,
    top_go: int,
    top_pathways: int,
    max_per_db: Optional[int],
) -> pd.DataFrame:
    """
    Filter and select top enrichment terms from g:Profiler results.
    
    Applies FDR and intersection size filters, then selects top GO terms
    from specified namespace and top pathway terms across databases.
    
    Args:
        df: g:Profiler results DataFrame
        fdr: FDR (adjusted p-value) threshold for significance
        min_is: Minimum intersection size (genes in query overlapping with term)
        go_source: GO namespace to select from (e.g., "GO:MF", "GO:BP", "GO:CC")
        top_go: Number of top GO terms to select
        top_pathways: Total number of pathway terms to select
        max_per_db: Maximum terms per pathway database (None to disable cap)
        
    Returns:
        DataFrame with selected terms, sorted by source and significance
    """
    df = df.copy()
    
    # Filter by statistical significance and minimum gene overlap
    df = df[(df["adjusted_p_value"] < fdr) & (df["intersection_size"] >= min_is)]
    
    # Sort by significance (p-value) first, then by intersection size
    df = df.sort_values(["adjusted_p_value", "intersection_size"], ascending=[True, False])
    
    # Select top N terms from the specified GO namespace
    go_df = df[df["source"] == go_source].head(top_go)
    
    # Process pathway databases (KEGG, Reactome, WikiPathways)
    pw = df[df["source"].isin(PATHWAY_SOURCES)]
    
    if max_per_db is None:
        # Simple case: take top N pathways regardless of database
        pw = pw.head(top_pathways)
    else:
        # Cap per-database: ensure diversity across pathway sources
        picked = []
        per_db = {s: 0 for s in PATHWAY_SOURCES}  # Track count per database
        
        for _, row in pw.iterrows():
            s = row["source"]
            # Add term if we haven't hit the per-database cap
            if per_db[s] < max_per_db:
                picked.append(row)
                per_db[s] += 1
            # Stop when we've selected enough total pathway terms
            if len(picked) >= top_pathways:
                break
        
        pw = pd.DataFrame(picked)
    
    # Combine GO terms and pathway terms
    out = pd.concat([go_df, pw], ignore_index=True)
    
    # Select relevant columns for output
    keep_cols = [
        "source",              # Database (GO:MF, KEGG, etc.)
        "term_name",           # Human-readable term name
        "term_id",             # Database-specific ID
        "adjusted_p_value",    # FDR-corrected p-value
        "intersection_size",   # Number of query genes in this term
        "term_size",           # Total genes annotated to this term
        "query_size",          # Total genes in query
        "intersections"        # Specific genes in overlap
    ]
    # Only keep columns that exist in the dataframe
    keep_cols = [c for c in keep_cols if c in out.columns]
    
    # Sort final output: by source, then significance, then intersection size
    return out[keep_cols].sort_values(
        ["source", "adjusted_p_value", "intersection_size"],
        ascending=[True, True, False],
    )


def main():
    """Parse arguments and execute term selection workflow."""
    ap = argparse.ArgumentParser()
    
    # Input/output files
    ap.add_argument("--infile", required=True, help="g:Profiler intersections CSV")
    ap.add_argument("--outfile", required=True, help="output CSV")
    
    # Filtering thresholds
    ap.add_argument("--fdr", type=float, default=0.01)
    ap.add_argument("--min_is", type=int, default=3)
    
    # GO namespace selection (backwards compatible interface)
    ap.add_argument("--top_mf", type=int, default=10, help="(legacy) top GO:MF terms")
    ap.add_argument("--top_bp", type=int, default=None, help="top GO:BP terms")
    ap.add_argument("--top_cc", type=int, default=None, help="top GO:CC terms")
    
    # GO namespace selection (new unified interface)
    ap.add_argument("--go", choices=["MF", "BP", "CC"], default=None,
                    help="alternative selector for GO namespace")
    ap.add_argument("--top_go", type=int, default=None, help="used with --go (default 10)")
    
    # Pathway selection
    ap.add_argument("--top_pathways", type=int, default=5)
    ap.add_argument("--max_per_db", type=int, default=2,
                    help="cap per pathway DB (KEGG/REAC/WP). Set 0 to disable.")
    
    args = ap.parse_args()
    
    # Load input data
    df = pd.read_csv(args.infile)
    
    # Convert max_per_db=0 to None (disables capping)
    max_per_db = None if args.max_per_db == 0 else args.max_per_db
    
    # Resolve which GO namespace and how many terms to select
    go_source, top_go = resolve_go_source_and_top(args)
    
    # Filter and select top terms
    out = pick_terms(df, args.fdr, args.min_is, go_source, top_go, args.top_pathways, max_per_db)
    
    # Write output
    out.to_csv(args.outfile, index=False)
    print(f"Wrote {len(out)} rows -> {args.outfile}")


if __name__ == "__main__":
    main()