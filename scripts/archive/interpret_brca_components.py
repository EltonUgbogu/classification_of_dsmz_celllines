#!/usr/bin/env python3
"""
Interpret all BRCA components at once using gene signature scoring.

This script:
1. Loads component consensus markers (UP/DOWN core)
2. Annotates them with gene symbols
3. Scores each component against known BRCA gene signatures
4. Generates a summary interpretation table
"""

from pathlib import Path
import pandas as pd
import sys

# ---- Configuration ----
CONS_DIR = Path("results/unsupervised/component_consensus_markers_v2/brca")
QC_PATH = Path("results/unsupervised/component_consensus_markers_v2/anchor_agreement_qc.tsv")
MAP_PATH = Path("resources/ensembl_to_symbol.tsv")

# Simple signature gene sets (symbols)
SIG = {
    "Luminal_HR": {"ESR1", "PGR", "GATA3", "FOXA1", "KRT8", "KRT18", "KRT19", "EPCAM", "SCUBE2", "MLPH", "XBP1"},
    "HER2_Amp": {"ERBB2", "GRB7", "PGAP3", "STARD3", "MIEN1", "PERLD1", "TCAP"},
    "Basal": {"KRT5", "KRT14", "KRT15", "KRT17", "TP63", "EGFR", "ITGA6", "LAMC2"},
    "EMT": {"VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1", "COL1A1", "COL1A2"},
    "Prolif": {"MKI67", "TOP2A", "CDK1", "CCNB1", "CCNB2", "BIRC5", "UBE2C", "AURKA"},
}


def load_map():
    """Load Ensembl ID to symbol mapping."""
    if not MAP_PATH.exists():
        raise FileNotFoundError(f"Mapping file not found: {MAP_PATH}")
    m = pd.read_csv(MAP_PATH, sep="\t")
    need = {"ensembl_id", "symbol"}
    if not need.issubset(set(m.columns)):
        raise ValueError(f"{MAP_PATH} must contain columns: {need}")
    # Remove duplicates, keep first
    m = m.drop_duplicates("ensembl_id", keep="first")
    return m


def annotate(df, m):
    """Annotate gene_id column with symbols."""
    df = df.copy()
    # Strip version from Ensembl IDs
    df["ensembl_id"] = df["gene_id"].astype(str).str.split(".", n=1).str[0]
    # Merge with mapping
    df_annotated = df.merge(m[["ensembl_id", "symbol"]], on="ensembl_id", how="left")
    return df_annotated


def score_component(up_core, down_core):
    """Score component against gene signatures."""
    up_syms = set(up_core["symbol"].dropna().str.upper())
    down_syms = set(down_core["symbol"].dropna().str.upper())
    
    # Convert signature sets to uppercase for matching
    scores = {}
    for name, genes in SIG.items():
        genes_upper = {g.upper() for g in genes}
        scores[f"{name}_UP_hits"] = len(up_syms & genes_upper)
        scores[f"{name}_DOWN_hits"] = len(down_syms & genes_upper)
        scores[f"{name}_UP_genes"] = ";".join(sorted(up_syms & genes_upper))
        scores[f"{name}_DOWN_genes"] = ";".join(sorted(down_syms & genes_upper))
    return scores


def pick_label(scores):
    """Assign a simple label based on signature scores."""
    lum_up = scores["Luminal_HR_UP_hits"]
    lum_down = scores["Luminal_HR_DOWN_hits"]
    her2_up = scores["HER2_Amp_UP_hits"]
    her2_down = scores["HER2_Amp_DOWN_hits"]
    bas_up = scores["Basal_UP_hits"]
    bas_down = scores["Basal_DOWN_hits"]
    emt_up = scores["EMT_UP_hits"]
    emt_down = scores["EMT_DOWN_hits"]
    prol_up = scores["Prolif_UP_hits"]
    prol_down = scores["Prolif_DOWN_hits"]
    
    # Strong luminal signature (UP)
    if lum_up >= 2:
        if her2_up >= 1:
            return "Luminal + HER2 axis"
        return "Luminal / HR+"
    
    # HER2 signature (UP)
    if her2_up >= 2:
        return "HER2-enriched-like"
    
    # Basal signature (UP)
    if bas_up >= 2:
        return "Basal-like"
    
    # EMT signature (UP)
    if emt_up >= 2:
        return "EMT/Mesenchymal-like"
    
    # Luminal suppression (DOWN) suggests non-luminal
    if lum_down >= 3:
        if bas_up >= 1 or emt_up >= 1:
            return "Basal/EMT-like (Luminal-suppressed)"
        if her2_up >= 1:
            return "HER2-like (Luminal-suppressed)"
        return "Non-luminal (Luminal-suppressed)"
    
    # EMT suppression suggests epithelial
    if emt_down >= 2:
        return "Epithelial (EMT-suppressed)"
    
    # Proliferation
    if prol_up >= 2:
        return "Proliferation-high"
    
    return "Unclear / mixed (needs enrichment)"


def main():
    """Main execution."""
    print("[INFO] Loading Ensembl to symbol mapping...")
    m = load_map()
    print(f"[INFO] Loaded {len(m)} gene mappings")
    
    print("[INFO] Loading anchor agreement QC...")
    if not QC_PATH.exists():
        print(f"[WARN] QC file not found: {QC_PATH}, skipping QC metrics")
        qc = pd.DataFrame()
    else:
        qc = pd.read_csv(QC_PATH, sep="\t")
        if "disease" in qc.columns:
            qc = qc[qc["disease"].str.upper() == "BRCA"]
        print(f"[INFO] Loaded QC data for {len(qc)} anchor pairs")
    
    print("[INFO] Processing components...")
    rows = []
    
    # Find all component directories
    comp_dirs = sorted([d for d in CONS_DIR.glob("component_*") if d.is_dir()], 
                       key=lambda p: int(p.name.split("_")[1]) if p.name.split("_")[1].isdigit() else 999)
    
    for comp_dir in comp_dirs:
        comp = comp_dir.name.split("_")[1]
        print(f"  Processing component {comp}...")
        
        # Try core first, fall back to extended if core is empty
        up_core_path = comp_dir / "component_consensus.UP.core.tsv"
        down_core_path = comp_dir / "component_consensus.DOWN.core.tsv"
        up_extended_path = comp_dir / "component_consensus.UP.extended.tsv"
        down_extended_path = comp_dir / "component_consensus.DOWN.extended.tsv"
        
        # Load core or extended
        if up_core_path.exists():
            up_core = pd.read_csv(up_core_path, sep="\t")
            if len(up_core) == 0 and up_extended_path.exists():
                print(f"    [INFO] Component {comp}: core UP empty, using extended")
                up_core = pd.read_csv(up_extended_path, sep="\t")
        elif up_extended_path.exists():
            print(f"    [INFO] Component {comp}: no core UP, using extended")
            up_core = pd.read_csv(up_extended_path, sep="\t")
        else:
            print(f"    [WARN] Component {comp}: Missing UP marker files, skipping")
            continue
        
        if down_core_path.exists():
            down_core = pd.read_csv(down_core_path, sep="\t")
            if len(down_core) == 0 and down_extended_path.exists():
                print(f"    [INFO] Component {comp}: core DOWN empty, using extended")
                down_core = pd.read_csv(down_extended_path, sep="\t")
        elif down_extended_path.exists():
            print(f"    [INFO] Component {comp}: no core DOWN, using extended")
            down_core = pd.read_csv(down_extended_path, sep="\t")
        else:
            print(f"    [WARN] Component {comp}: Missing DOWN marker files, skipping")
            continue
        
        up_core = annotate(up_core, m)
        down_core = annotate(down_core, m)
        
        # Score component
        scores = score_component(up_core, down_core)
        label = pick_label(scores)
        
        # QC summary for this component
        if len(qc) > 0 and "component" in qc.columns:
            qc_comp = qc[qc["component"].astype(str) == str(comp)]
            j_up = qc_comp["jaccard_up"].mean() if len(qc_comp) > 0 else float("nan")
            j_dn = qc_comp["jaccard_down"].mean() if len(qc_comp) > 0 else float("nan")
            # Count unique anchors (handle both anchor_a/anchor_b and anchor1/anchor2)
            anchor_cols = [c for c in ["anchor_a", "anchor_b", "anchor1", "anchor2"] if c in qc_comp.columns]
            if len(anchor_cols) > 0:
                n_anchors = len(set(qc_comp[anchor_cols].values.flatten()))
            else:
                n_anchors = 0
        else:
            j_up = float("nan")
            j_dn = float("nan")
            n_anchors = 0
        
        # Top marker symbols (ranked by rank_score)
        top_up = ";".join(
            up_core.sort_values("rank_score", ascending=False)["symbol"]
            .dropna()
            .head(10)
            .tolist()
        )
        top_dn = ";".join(
            down_core.sort_values("rank_score", ascending=False)["symbol"]
            .dropna()
            .head(10)
            .tolist()
        )
        
        rows.append({
            "component": int(comp),
            "label": label,
            "n_up_core": len(up_core),
            "n_down_core": len(down_core),
            "n_anchors": n_anchors,
            "mean_jaccard_up": j_up,
            "mean_jaccard_down": j_dn,
            "top10_up_symbols": top_up,
            "top10_down_symbols": top_dn,
            **scores
        })
    
    # Create output DataFrame
    out = pd.DataFrame(rows).sort_values("component")
    
    # Write output
    out_path = CONS_DIR / "BRCA_all_components_interpretation.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\n[OK] Wrote interpretation to: {out_path}")
    
    # Print summary
    print("\n" + "="*80)
    print("BRCA Component Interpretation Summary")
    print("="*80)
    print(out.to_string(index=False))
    print("="*80)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
