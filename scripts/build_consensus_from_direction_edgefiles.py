#!/usr/bin/env python3
import argparse
import itertools
import pandas as pd
import numpy as np
import re
from collections import defaultdict
from pathlib import Path

def looks_long(x: str) -> bool:
    """Check if an ID looks like a profile-level NG identifier."""
    return bool(re.match(r"^NG[-_][A-Za-z0-9]+_", str(x)))


def canonical_cell_line_id(x: str) -> str:
    """Map multicohort profile-level IDs to biological cell-line/group IDs."""
    node = str(x).strip()
    m = re.match(r"^NG[-_][^_]+_(.+?)_lib", node)
    if m:
        node = m.group(1)
    node = re.sub(r"^(RBL_\d+)_\d+$", r"\1", node)
    return node


def map_and_canonicalise(x: str, long2short=None) -> str:
    node = str(x).strip()
    if long2short is not None:
        node = long2short.get(node, node)
    node = canonical_cell_line_id(node)
    if long2short is not None:
        node = long2short.get(node, node)
    return node

def load_name_map(path):
    """Load long_id -> short_id mapping from TSV file."""
    m = pd.read_csv(path, sep="\t")
    if not {"long_id", "short_id"}.issubset(m.columns):
        raise SystemExit("[ERROR] name_map must have columns: long_id, short_id")
    return dict(zip(m["long_id"].astype(str), m["short_id"].astype(str)))


def ambiguous_bases_from_name_map(path):
    """Bases that have >1 short_id (e.g. RBL_15 -> RBL_15_2, RBL_15_4). These must never appear as bare IDs."""
    m = pd.read_csv(path, sep="\t", dtype=str)
    if "short_id" not in m.columns:
        return set()
    bases = m["short_id"].str.replace(r"_\d+$", "", regex=True)
    counts = bases.value_counts()
    return set(counts[counts > 1].index.tolist())

def canonicalise_edge_df(df, long2short):
    """
    Canonicalise edge dataframe by mapping long IDs to short IDs.
    Supports either (node1,node2) or (cell_line1,cell_line2) columns.
    Returns: (canonicalised_df, col1_name, col2_name)
    """
    if {"node1", "node2"}.issubset(df.columns):
        a, b = "node1", "node2"
    elif {"cell_line1", "cell_line2"}.issubset(df.columns):
        a, b = "cell_line1", "cell_line2"
    else:
        raise SystemExit(f"[ERROR] edge df missing node columns. Found: {list(df.columns)}")

    df = df.copy()
    df[a] = df[a].astype(str).map(lambda x: map_and_canonicalise(x, long2short))
    df[b] = df[b].astype(str).map(lambda x: map_and_canonicalise(x, long2short))
    return df, a, b


def load_node_universe(paths, long2short=None):
    """Load a complete cell-line node universe from one or more TSV files."""
    if not paths:
        return []

    preferred_cols = (
        "short_id",
        "cell_line",
        "cell_line_display",
        "cell_line1",
        "node1",
    )
    nodes = []
    seen = set()
    for raw_path in paths:
        path = Path(raw_path)
        if not path.exists():
            raise SystemExit(f"[ERROR] node universe file not found: {path}")
        df = pd.read_csv(path, sep="\t", dtype=str)
        col = next((c for c in preferred_cols if c in df.columns), None)
        if col is None:
            raise SystemExit(
                f"[ERROR] cannot infer node column for {path}. "
                f"Expected one of {preferred_cols}; found {list(df.columns)}"
            )
        for value in df[col].dropna().astype(str):
            node = value.strip()
            node = map_and_canonicalise(node, long2short=long2short)
            if node and node not in seen:
                seen.add(node)
                nodes.append(node)
    return sorted(nodes)


def validate_node_names(nodes, ambiguous_bases, require_short=False):
    if require_short:
        bad = sorted(n for n in nodes if looks_long(str(n)))
        if bad:
            raise SystemExit(
                f"[ERROR] Long IDs still present in node universe after mapping: {bad[:10]}"
            )
    if ambiguous_bases:
        bare = sorted(n for n in nodes if str(n) in ambiguous_bases)
        if bare:
            raise SystemExit(
                "[ERROR] Ambiguous bare base IDs in node universe "
                f"(cannot be resolved safely): {bare}"
            )

def main():
    ap = argparse.ArgumentParser(
        description="Aggregate cell_line_similarity_graph_edges_<direction>.tsv across all directions into a consensus graph."
    )
    ap.add_argument("--tumour_nh_dir", required=True,
                    help="tumour_neighbourhoods directory (contains direction folders + final_consensus_all)")
    ap.add_argument("--out_edges", required=True,
                    help="Output TSV for consensus edges")
    ap.add_argument("--min_support", type=int, default=1,
                    help="Minimum number of directions supporting an edge (default 1)")
    ap.add_argument("--min_mean_sim", type=float, default=None,
                    help="Optional: minimum mean similarity across supporting directions")
    ap.add_argument("--name_map", type=str, default=None,
                    help="Path to TSV mapping file with columns: long_id, short_id")
    ap.add_argument("--require_short", action="store_true",
                    help="Fail if long IDs remain after mapping")
    ap.add_argument("--directions", type=str, default=None,
                    help="Optional comma-separated direction list. When set, aggregate only these directions.")
    ap.add_argument("--mode",
                    choices=("majority_threshold", "union_supported_edges"),
                    default="majority_threshold",
                    help="Edge-retention mode for --out_edges.")
    ap.add_argument("--cohort", type=str, default="",
                    help="Cohort/profile label written to support tables.")
    ap.add_argument("--out_support", type=str, default=None,
                    help="Optional TSV containing all unordered node pairs and support annotations.")
    ap.add_argument("--node_universe", type=str, default=None,
                    help="Optional comma-separated TSV files containing the full node universe.")
    args = ap.parse_args()

    root = Path(args.tumour_nh_dir)
    if not root.exists():
        raise SystemExit(f"Not found: {root}")

    # Load name mapping and ambiguous bases if provided
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

    # (u,v) directed calls -> set(directions), list(similarities)
    dir_calls = defaultdict(set)
    dir_sims  = defaultdict(list)
    observed_nodes = set()
    canonicalisation_rows = []

    # Discover all per-direction edge files, or restrict to the configured
    # representation list passed by Snakemake.
    if args.directions:
        directions = [d.strip() for d in args.directions.split(",") if d.strip()]
        edge_files = [
            root / d / "final_consensus" / f"cell_line_similarity_graph_edges_{d}.tsv"
            for d in directions
        ]
        missing = [str(p) for p in edge_files if not p.exists()]
        if missing:
            raise SystemExit("[ERROR] Missing configured direction edge files:\n" + "\n".join(missing))
    else:
        edge_files = sorted(root.glob("*/final_consensus/cell_line_similarity_graph_edges_*.tsv"))
    if not edge_files:
        raise SystemExit(f"No edge files found under {root}/*/final_consensus/")

    for ef in edge_files:
        direction = ef.name.replace("cell_line_similarity_graph_edges_", "").replace(".tsv", "")
        df = pd.read_csv(ef, sep="\t")

        # Canonicalise IDs if mapping provided
        if long2short:
            df, c1, c2 = canonicalise_edge_df(df, long2short)
        else:
            # expect either (cell_line1, cell_line2) or (node1, node2)
            if {"cell_line1", "cell_line2"}.issubset(df.columns):
                c1, c2 = "cell_line1", "cell_line2"
            elif {"node1", "node2"}.issubset(df.columns):
                c1, c2 = "node1", "node2"
            else:
                raise SystemExit(f"{ef}: cannot find edge columns. Columns={list(df.columns)}")

        # Check for remaining long IDs if required
        if args.require_short:
            nodes = set(df[c1]).union(df[c2])
            bad = sorted([n for n in nodes if looks_long(str(n))])
            if bad:
                raise SystemExit(f"[ERROR] Long IDs still present after mapping in {ef.name}: {bad[:10]}")

        # Fail if any ambiguous bare base IDs present (e.g. RBL_15, RBL_20)
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
            u = map_and_canonicalise(str(r[c1]).strip(), long2short=long2short)
            v = map_and_canonicalise(str(r[c2]).strip(), long2short=long2short)
            if not u or not v or u == v:
                continue

            a, b = sorted((u, v))
            observed_nodes.update((a, b))
            dir_calls[(a, b)].add(direction)
            if has_sim:
                try:
                    dir_sims[(a, b)].append(float(r["similarity"]))
                except Exception:
                    pass
        canonical_pairs = {
            tuple(sorted((
                map_and_canonicalise(str(r[c1]).strip(), long2short=long2short),
                map_and_canonicalise(str(r[c2]).strip(), long2short=long2short),
            )))
            for _, r in df.iterrows()
            if map_and_canonicalise(str(r[c1]).strip(), long2short=long2short)
            and map_and_canonicalise(str(r[c2]).strip(), long2short=long2short)
            and map_and_canonicalise(str(r[c1]).strip(), long2short=long2short)
            != map_and_canonicalise(str(r[c2]).strip(), long2short=long2short)
        }
        canonicalisation_rows.append((direction, len(df), len(canonical_pairs)))

    node_paths = [p.strip() for p in (args.node_universe or "").split(",") if p.strip()]
    node_universe = load_node_universe(node_paths, long2short=long2short)
    if node_universe:
        validate_node_names(
            node_universe,
            ambiguous_bases=ambiguous_bases,
            require_short=args.require_short,
        )
        missing_observed = sorted(observed_nodes.difference(node_universe))
        if missing_observed:
            raise SystemExit(
                "[ERROR] Observed edge endpoints are absent from the node universe: "
                + ", ".join(missing_observed[:20])
            )
    else:
        node_universe = sorted(observed_nodes)
        print("[WARN] --node_universe not supplied; node set inferred from edge endpoints only.")

    raw_profile_nodes = sorted(n for n in node_universe if looks_long(str(n)))
    if raw_profile_nodes:
        raise SystemExit(
            "[ERROR] Profile-level nodes remain after canonicalisation: "
            + ", ".join(raw_profile_nodes[:20])
        )

    for direction, input_rows, canonical_pairs in canonicalisation_rows:
        if input_rows != canonical_pairs:
            print(
                f"[INFO] {direction}: canonical edge rows collapsed "
                f"{input_rows} -> {canonical_pairs} unordered pairs"
            )

    n_configured = len(edge_files)
    support_threshold = int(args.min_support)
    support_rows = []
    for u, v in itertools.combinations(node_universe, 2):
        dirs_uv = dir_calls.get((u, v), set())
        dirs_vu = dir_calls.get((v, u), set())
        dirs_union = sorted(dirs_uv.union(dirs_vu))
        support = len(dirs_union)
        reciprocal = 1 if dirs_uv and dirs_vu else 0
        sims = dir_sims.get((u, v), []) + dir_sims.get((v, u), [])
        mean_sim = float(np.mean(sims)) if sims else np.nan
        sum_sim = float(np.sum(sims)) if sims else np.nan
        max_sim = float(np.max(sims)) if sims else np.nan
        majority_retained = support >= support_threshold
        union_retained = support >= 1
        union_edge_style = "dashed" if support == 1 else ("solid" if support >= 2 else "none")
        support_rows.append({
            "cohort": args.cohort,
            "node_a": u,
            "node_b": v,
            "support": support,
            "supporting_directions": ";".join(dirs_union),
            "n_configured_directions": n_configured,
            "support_threshold": support_threshold,
            "majority_threshold_retained": majority_retained,
            "union_supported_edges_retained": union_retained,
            "union_edge_style": union_edge_style,
            "reciprocal": reciprocal,
            "support_weight_mean": mean_sim,
            "support_weight_sum": sum_sim,
            "support_weight_max": max_sim,
        })

    support_df = pd.DataFrame(support_rows)
    if args.out_support:
        support_path = Path(args.out_support)
        support_path.parent.mkdir(parents=True, exist_ok=True)
        support_df.to_csv(support_path, sep="\t", index=False)

    if args.mode == "majority_threshold":
        retain_col = "majority_threshold_retained"
    else:
        retain_col = "union_supported_edges_retained"

    out = support_df[support_df[retain_col]].copy()
    if args.min_mean_sim is not None:
        out = out[
            out["support_weight_mean"].notna()
            & (out["support_weight_mean"] >= args.min_mean_sim)
        ].copy()
    out["node1"] = out["node_a"]
    out["node2"] = out["node_b"]
    out["support_directions"] = out["support"]
    out["methods_union"] = out["supporting_directions"]
    out["edge_style"] = np.where(out["support"] == 1, "dashed", "solid")
    out["graph_product"] = args.mode
    out = out.sort_values(
        ["support_directions", "reciprocal", "support_weight_mean", "node1", "node2"],
        ascending=[False, False, False, True, True],
    )
    out = out[
        [
            "cohort",
            "graph_product",
            "node1",
            "node2",
            "node_a",
            "node_b",
            "support",
            "support_directions",
            "supporting_directions",
            "methods_union",
            "n_configured_directions",
            "support_threshold",
            "majority_threshold_retained",
            "union_supported_edges_retained",
            "union_edge_style",
            "edge_style",
            "reciprocal",
            "support_weight_mean",
            "support_weight_sum",
            "support_weight_max",
        ]
    ]

    out_path = Path(args.out_edges)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(
        f"[OK] wrote {out_path} with {len(out)} {args.mode} edges "
        f"from {len(edge_files)} feature--distance representations"
    )
    if args.out_support:
        print(f"[OK] wrote support table {args.out_support} with {len(support_df)} unordered node pairs")

if __name__ == "__main__":
    main()
