#!/usr/bin/env python3
"""
confirm_brca_gene_list.py

Confirm whether an input BRCA gene list is informative for:
- breast cancer biology (enrichment themes)
- breast cancer subtypes (Luminal / HER2-enriched / Basal-like signals)

Requires internet access (uses g:Profiler public API via gprofiler-official).
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd
from scipy.stats import fisher_exact

try:
    from gprofiler import GProfiler
except ImportError as e:
    raise SystemExit(
        "Missing dependency: gprofiler-official\n"
        "Install with:\n"
        "  pip install gprofiler-official pandas scipy\n"
    ) from e


# -----------------------------
# Curated panels / marker sets
# -----------------------------

PAM50: Set[str] = {
    "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CDC20","CDC6","CENPF",
    "CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXC1","FOXA1","GPR160",
    "GRB7","KIF2C","MELK","MIA","MKI67","MLPH","MMP11","MYBL2","NAT1","ORC6",
    "PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T",
    "BASAL1","KRT5","KRT14","KRT17","MDM2","NDC80","NUF2","RACGAP1","STON2","TRIM29"
}
# Note: Some published PAM50 lists vary slightly by source; this is a commonly used set.
# If you have an "official" list file, you can replace this set.

SUBTYPE_MARKERS: Dict[str, Set[str]] = {
    # “Proxies” at RNA level for clinical markers & subtype programs
    "Clinical_markers_ER_PR_HER2": {"ESR1", "PGR", "ERBB2", "GRB7"},
    "Luminal_program": {"ESR1","PGR","GATA3","FOXA1","XBP1","NAT1","BCL2","SLC39A6","TFF1","TFF3"},
    "HER2_enriched_program": {"ERBB2","GRB7","FGFR4","PTK6","STARD3","TCAP"},
    "Basal_like_program": {"KRT5","KRT14","KRT17","EGFR","VIM","ITGA6","ITGB4","LAMC2","FOXC1"},
    "Proliferation_program": {"MKI67","TOP2A","CCNB1","CDC20","BIRC5","AURKA","AURKB","UBE2C","MELK"},
}

# g:Profiler sources to query
GPROF_SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"]


def read_gene_list(path: str) -> List[str]:
    genes = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if not g or g.startswith("#"):
                continue
            # strip Ensembl version suffix ENSG... .8
            g = g.split()[0]
            g = g.split(".")[0]
            genes.append(g)
    # de-duplicate while preserving order
    seen = set()
    out = []
    for g in genes:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out


def ensure_outdir(outdir: str) -> Path:
    p = Path(outdir)
    p.mkdir(parents=True, exist_ok=True)
    return p


def convert_ensembl_to_symbols(gp: GProfiler, genes: List[str], organism: str) -> pd.DataFrame:
    """
    Uses g:Profiler convert endpoint.
    Returns a dataframe mapping input -> converted_name (symbol) where possible.
    """
    conv = gp.convert(
        organism=organism,
        query=genes,
        target_namespace="ENSG",
    )
    # conv includes columns like "incoming", "converted", "name", etc depending on API response
    # We’ll keep it flexible and standardize.
    df = pd.DataFrame(conv)

    # Heuristic column picks
    incoming = "incoming" if "incoming" in df.columns else df.columns[0]
    name_col_candidates = [c for c in df.columns if c.lower() in {"name", "converted_name", "symbol"}]
    symbol_col = name_col_candidates[0] if name_col_candidates else None

    if symbol_col is None:
        # sometimes g:Profiler returns 'name' for gene symbol
        for c in df.columns:
            if c.lower().endswith("name"):
                symbol_col = c
                break

    if symbol_col is None:
        # fallback: keep raw table, user can inspect
        df = df.rename(columns={incoming: "input_id"})
        df["symbol"] = None
        return df[["input_id", "symbol"]].drop_duplicates()

    df = df.rename(columns={incoming: "input_id", symbol_col: "symbol"})
    # Clean
    df["input_id"] = df["input_id"].astype(str).str.split(".").str[0]
    df["symbol"] = df["symbol"].astype(str).where(df["symbol"].notna(), None)
    # Some rows may have "nan" as string
    df.loc[df["symbol"].str.lower() == "nan", "symbol"] = None
    return df[["input_id", "symbol"]].drop_duplicates()


def fisher_enrichment(
    gene_symbols: Set[str],
    signature: Set[str],
    universe_size: int,
) -> Tuple[float, float, int]:
    """
    One-sided Fisher exact test (greater) for overlap enrichment.
    Returns (oddsratio, pvalue, overlap_count).
    """
    overlap = len(gene_symbols & signature)
    a = overlap
    b = len(signature) - overlap
    c = len(gene_symbols) - overlap
    d = universe_size - (a + b + c)
    if d < 0:
        # Universe too small; clip to avoid crash, but warn via p=1
        return float("nan"), 1.0, overlap
    oddsratio, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return oddsratio, p, overlap


def run_enrichment(
    gp: GProfiler,
    genes: List[str],
    organism: str,
    fdr: float,
    sources: List[str],
) -> pd.DataFrame:
    # gprofiler-official 1.x uses .profile(); older/other clients may use .gost()
    if hasattr(gp, "profile"):
        df = gp.profile(
            organism=organism,
            query=genes,
            sources=sources,
            user_threshold=fdr,
            significance_threshold_method="fdr",
            no_evidences=False,
            ordered=False,
        )
        df = pd.DataFrame(df)
    elif hasattr(gp, "gost"):
        res = gp.gost(
            organism=organism,
            query=genes,
            sources=sources,
            user_threshold=fdr,
            correction_method="fdr",
            ordered=False,
        )
        df = pd.DataFrame(res["result"]) if isinstance(res, dict) and "result" in res else pd.DataFrame(res)
    else:
        raise RuntimeError("No g:Profiler method found. Install/upgrade gprofiler-official.")

    if df.empty:
        return df

    keep = [
        c for c in [
            "source", "native", "name", "p_value", "significant",
            "term_size", "query_size", "intersection_size",
            "effective_domain_size", "intersections",
        ]
        if c in df.columns
    ]
    return df[keep].sort_values(["source", "p_value"], ascending=[True, True])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True, help="Text file: one Ensembl gene ID per line (ENSG...)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--organism", default="hsapiens", help="g:Profiler organism (default: hsapiens)")
    ap.add_argument("--fdrs", default="0.01,0.05", help="Comma-separated FDR thresholds to test (default: 0.01,0.05)")
    ap.add_argument("--universe_size", type=int, default=20000, help="Universe size for Fisher tests (default 20000)")
    args = ap.parse_args()

    outdir = ensure_outdir(args.outdir)
    genes = read_gene_list(args.genes)

    gp = GProfiler()

    # Convert Ensembl -> symbols
    mapping = convert_ensembl_to_symbols(gp, genes, args.organism)
    mapping.to_csv(outdir / "ensembl_to_symbol.tsv", sep="\t", index=False)

    symbols = set([s for s in mapping["symbol"].dropna().astype(str) if s and s.lower() != "none"])
    # Some conversions may return multiple mappings; we treat as set for overlap stats

    # Overlap reports
    overlap_rows = []

    panels = {"PAM50": PAM50, **SUBTYPE_MARKERS}
    for panel_name, panel_genes in panels.items():
        odds, p, k = fisher_enrichment(symbols, panel_genes, args.universe_size)
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

    overlap_df = pd.DataFrame(overlap_rows).sort_values("fisher_p_greater", ascending=True)
    overlap_df.to_csv(outdir / "overlap_report.tsv", sep="\t", index=False)

    # Enrichment at requested FDRs
    fdrs = []
    for x in args.fdrs.split(","):
        x = x.strip()
        if not x:
            continue
        fdrs.append(float(x))

    enrich_summaries = []
    for fdr in fdrs:
        enr = run_enrichment(gp, genes, args.organism, fdr, GPROF_SOURCES)
        outpath = outdir / f"enrichment_FDR{fdr:.2f}.tsv"
        enr.to_csv(outpath, sep="\t", index=False)

        # Quick source summary
        if enr.empty:
            enrich_summaries.append(f"FDR {fdr:.2f}: no significant terms (sources={','.join(GPROF_SOURCES)})")
        else:
            counts = enr["source"].value_counts().to_dict() if "source" in enr.columns else {}
            enrich_summaries.append(f"FDR {fdr:.2f}: {len(enr)} terms; per-source: {counts}")

    # Write a short human summary
    with open(outdir / "summary.txt", "w", encoding="utf-8") as f:
        f.write("BRCA gene list confirmation summary\n")
        f.write("=================================\n\n")
        f.write(f"Input genes: {len(genes)} (unique)\n")
        f.write(f"Converted gene symbols: {len(symbols)}\n\n")

        f.write("Panel overlap (one-sided Fisher exact, greater):\n")
        top = overlap_df.head(10)
        for _, r in top.iterrows():
            f.write(
                f"  - {r['panel']}: overlap={r['overlap_count']}/{r['panel_size']} "
                f"(p={r['fisher_p_greater']:.3e}, OR={r['fisher_oddsratio']})\n"
            )
        f.write("\nEnrichment summaries:\n")
        for line in enrich_summaries:
            f.write(f"  - {line}\n")

        f.write("\nNotes:\n")
        f.write("- If PAM50 / Luminal / HER2 / Basal overlaps are weak and enrichment is mostly generic\n")
        f.write("  (translation/ribosome/cell cycle only), the list may reflect general variability rather than subtype programs.\n")
        f.write("- Try comparing FDR=0.01 vs 0.05 outputs to see if breast-specific pathways appear at the relaxed threshold.\n")

    print(f"[OK] Wrote outputs to: {outdir}")
    print(f" - ensembl_to_symbol.tsv")
    print(f" - overlap_report.tsv")
    for fdr in fdrs:
        print(f" - enrichment_FDR{fdr:.2f}.tsv")
    print(" - summary.txt")


if __name__ == "__main__":
    main()
