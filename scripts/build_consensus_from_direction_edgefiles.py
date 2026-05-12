#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import re
from collections import defaultdict
from pathlib import Path

def looks_long(x: str) -> bool:
    return bool(re.match(r"^NG-\d+_", str(x)))

def load_name_map(path):
    m = pd.read_csv(path, sep="\t")
    if not {"long_id", "short_id"}.issubset(m.columns):
        raise SystemExit("[ERROR] name_map must have columns: long_id, short_id")
    return dict(zip(m["long_id"].astype(str), m["short_id"].astype(str)))

def ambiguous_bases_from_name_map(path):
    # Bases with >1 short_id must never appear as bare IDs
    m = pd.read_csv(path, sep="\t", dtype=str)
    if "short_id" not in m.columns:
        return set()
    bases = m["short_id"].str.replace(r"_\d+$", "", regex=True)
    counts = bases.value_counts()
    return set(counts[counts > 1].index.tolist())

def canonicalise_edge_df(df, long2short):
    if {"node1", "node2"}.issubset(df.columns):
        a, b = "node1", "node2"
    elif {"cell_line1", "cell_line2"}.issubset(df.columns):
        a, b = "cell_line1", "cell_line2"
    else:
        raise SystemExit(f"[ERROR] edge df missing node columns. Found: {list(df.columns)}")
    df = df.copy()
    df[a] = df[a].astype(str).map(lambda x: long2short.get(x, x))
    df[b] = df[b].astype(str).map(lambda x: long2short.get(x, x))
    return df, a, b

def main():
    ap = argparse.ArgumentParser(
        description="Aggregate cell_line_similarity_graph_edges_<direction>.tsv across all directions into a consensus graph."
    )
    ap.add_argument("--tumour_nh_dir", required=True)
    ap.add_argument("--out_edges", required=True)
    ap.add_argument("--min_support", type=int, default=1)
    ap.add_argument("--min_mean_sim", type=float, default=None)
    ap.add_argument("--name_map", type=str, default=None)
    ap.add_argument("--require_short", action="store_true")
    args = ap.parse_args()

    root = Path(args.tumour_nh_dir)
    if not root.exists():
        raise SystemExit(f"Not found: {root}")

    long2short = None
    ambiguous_bases = set()
    if args.name_map:
        name_map_path = Path(args.name_map)
        if not name_map_path.exists():
            raise SystemExit(f"[ERROR] name_map file not found: {name_map_path}")
        long2short = load_name_map(name_map_path)
        ambiguous_bases = ambiguous_bases_from_name_map(name_map_path)
        if ambiguous_bases:
            print(f"[INFO] Ambiguous bases (must NOT appear as bare IDs): {sorted(ambiguous_bases)}")
        print(f"[INFO] Loaded {len(long2short)} ID mappings from {name_map_path}")

    dir_calls = defaultdict(set)
    dir_sims  = defaultdict(list)

    edge_files = sorted(root.glob("*/final_consensus/cell_line_similarity_graph_edges_*.tsv"))
    if not edge_files:
        raise SystemExit(f"No edge files found under {root}/*/final_consensus/")

    for ef in edge_files:
        direction = ef.name.replace("cell_line_similarity_graph_edges_", "").replace(".tsv", "")
        df = pd.read_csv(ef, sep="\t")

        if long2short:
            df, c1, c2 = canonicalise_edge_df(df, long2short)
        else:
            if {"cell_line1", "cell_line2"}.issubset(df.columns):
                c1, c2 = "cell_line1", "cell_line2"
            elif {"node1", "node2"}.issubset(df.columns):
                c1, c2 = "node1", "node2"
            else:
                raise SystemExit(f"{ef}: cannot find edge columns. Columns={list(df.columns)}")

        if args.require_short:
            nodes = set(df[c1]).union(df[c2])
            bad = sorted([n for n in nodes if looks_long(str(n))])
            if bad:
                raise SystemExit(f"[ERROR] Long IDs still present after mapping in {ef.name}: {bad[:10]}")

        if ambiguous_bases:
            nodes = set(df[c1]).union(df[c2])
            bare = sorted([n for n in nodes if str(n) in ambiguous_bases])
            if bare:
                offenders = df[(df[c1].isin(bare)) | (df[c2].isin(bare))]
                ex = offenders.head(5)[[c1, c2]].to_dict("records")
                raise SystemExit(
                    f"[ERROR] Ambiguous bare base IDs in {ef.name} (cannot be resolved safely). "
                    f"Bases: {bare}. Expected replicate-specific IDs. Examples: {ex}"
                )

        has_sim = "similarity" in df.columns

        for _, r in df.iterrows():
            u = str(r[c1]).strip()
            v = str(r[c2]).strip()
            if not u or not v or u == v:
                continue
            dir_calls[(u, v)].add(direction)
            if has_sim:
                try:
                    dir_sims[(u, v)].append(float(r["similarity"]))
                except Exception:
                    pass

    und_rows = []
    seen = set()

    for (u, v), dirs_uv in dir_calls.items():
        key = tuple(sorted((u, v)))
        if key in seen:
            continue
        seen.add(key)

        dirs_vu = dir_calls.get((v, u), set())
        reciprocal = 1 if len(dirs_uv | dirs_vu) >= 2 else 0

        dirs_union = sorted(dirs_uv.union(dirs_vu))
        support = len(dirs_union)

        sims = dir_sims.get((u, v), []) + dir_sims.get((v, u), [])
        mean_sim = float(np.mean(sims)) if sims else np.nan
        sum_sim  = float(np.sum(sims)) if sims else np.nan
        max_sim  = float(np.max(sims)) if sims else np.nan

        if support < args.min_support:
            continue
        if args.min_mean_sim is not None and (np.isnan(mean_sim) or mean_sim < args.min_mean_sim):
            continue

        # support_type labels whether an edge has multi direction consensus or not
        support_type = "multi supported" if support >= 2 else "single supported"

        und_rows.append({
            "node1": key[0],
            "node2": key[1],
            "reciprocal": reciprocal,
            "support_directions": support,
            "support_weight_mean": mean_sim,
            "support_weight_sum": sum_sim,
            "support_weight_max": max_sim,
            "methods_union": ";".join(dirs_union),
            "direction_list": ";".join(dirs_union),
            "support_type": support_type,
        })

    out = pd.DataFrame(und_rows).sort_values(
        ["support_directions", "reciprocal", "support_weight_mean"],
        ascending=[False, False, False]
    )

    out_path = Path(args.out_edges)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {out_path} with {len(out)} edges from {len(edge_files)} directions")

if __name__ == "__main__":
    main()
