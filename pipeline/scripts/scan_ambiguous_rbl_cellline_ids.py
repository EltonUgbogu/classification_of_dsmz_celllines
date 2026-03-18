#!/usr/bin/env python3
"""
Scan direction edge files to find which ones contain ambiguous base IDs
(e.g., RBL_15 instead of RBL_15_2 or RBL_15_4).
"""
import pandas as pd
from pathlib import Path
from collections import defaultdict

tumour_nh_dir = Path("results/unsupervised/rbl/tumour_neighbourhoods")
name_map_path = Path("results/unsupervised/rbl/tumour_neighbourhoods/final_consensus_all/dsmz_cellline_shortnames.tsv")

# Load mapping: long_id -> short_id
m = pd.read_csv(name_map_path, sep="\t", dtype=str)
long2short = dict(zip(m["long_id"], m["short_id"]))

# Build: base -> variants (only for bases that have >1 variant)
base2vars = defaultdict(set)
for short in long2short.values():
    if short.startswith("RBL_") and short.count("_") == 2:
        base = "_".join(short.split("_")[:2])   # RBL_15_2 -> RBL_15
        base2vars[base].add(short)

ambiguous_bases = {b: sorted(v) for b, v in base2vars.items() if len(v) > 1}
if not ambiguous_bases:
    print("[INFO] No ambiguous bases found in mapping (unexpected for your RBL).")
else:
    print("[INFO] Ambiguous bases (must NOT appear as bare IDs):")
    for b, v in sorted(ambiguous_bases.items()):
        print(f"  - {b} -> {v}")

# Find all direction edge files
edge_files = sorted(tumour_nh_dir.glob("*/final_consensus/DSMZ_DSMZ_graph_edges_*.tsv"))
if not edge_files:
    raise SystemExit(f"[ERROR] No direction edge files found under: {tumour_nh_dir}")

print(f"\n[INFO] Scanning {len(edge_files)} direction edge files...")

hits = []  # (file, base, count, examples)

for f in edge_files:
    try:
        df = pd.read_csv(f, sep="\t", dtype=str)
    except Exception as e:
        print(f"[WARN] failed to read {f}: {e}")
        continue

    # guess columns (support both naming conventions)
    if {"cell_line1", "cell_line2"}.issubset(df.columns):
        a, b = "cell_line1", "cell_line2"
    elif {"node1", "node2"}.issubset(df.columns):
        a, b = "node1", "node2"
    else:
        print(f"[WARN] Skipping {f}: missing edge columns. Found: {list(df.columns)}")
        continue

    # apply long->short mapping *only* for long IDs present
    df[a] = df[a].map(lambda x: long2short.get(x, x))
    df[b] = df[b].map(lambda x: long2short.get(x, x))

    nodes = pd.concat([df[a], df[b]]).dropna().astype(str)

    for base, variants in ambiguous_bases.items():
        mask = (df[a] == base) | (df[b] == base)
        c = int(mask.sum())
        if c > 0:
            ex = df.loc[mask, [a, b]].head(5).to_dict("records")
            hits.append((str(f), base, c, variants, ex))

if not hits:
    print("\n[OK] No offending bare base IDs found in any direction file.")
else:
    print("\n[FAIL] Found offending bare base IDs in direction files:\n")
    for f, base, c, variants, ex in sorted(hits, key=lambda x: (-x[2], x[1], x[0])):
        print(f"- {f}")
        print(f"    base: {base}  (expected one of: {variants})")
        print(f"    rows: {c}")
        print(f"    examples:")
        for r in ex:
            print(f"      {r}")
        print()
