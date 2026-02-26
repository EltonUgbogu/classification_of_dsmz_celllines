#!/usr/bin/env python3
"""Fail if consensus edges TSV contains ambiguous bare base IDs (e.g. RBL_15, RBL_20)."""
import pandas as pd
import sys
from pathlib import Path

# Default path for RBL
DEFAULT_EDGES = "results/unsupervised/rbl/tumour_neighbourhoods/final_consensus_all/dsmz_cellline_consensus_edges.tsv"

AMBIGUOUS_BASES = {"RBL_15", "RBL_20"}

def main():
    p = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(DEFAULT_EDGES)
    if not p.exists():
        print(f"[SKIP] File not found: {p}")
        sys.exit(0)
    df = pd.read_csv(p, sep="\t")
    nodes = set(df["node1"]).union(df["node2"])
    bad = [n for n in nodes if n in AMBIGUOUS_BASES]
    if bad:
        print("[FAIL] Ambiguous bare IDs present:", bad)
        sys.exit(1)
    print("[OK] No ambiguous bare IDs in consensus edges.")

if __name__ == "__main__":
    main()
