#!/usr/bin/env python3
import argparse
import math
from pathlib import Path
import pandas as pd


def read_edges(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["cell_line1", "cell_line2", "similarity"])
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["cell_line1", "cell_line2", "similarity"])
    keep = [c for c in ["cell_line1", "cell_line2", "similarity"] if c in df.columns]
    return df[keep].copy()


def edge_set(df: pd.DataFrame):
    if df.empty:
      return set()
    return {tuple(sorted((str(a), str(b)))) for a, b in zip(df["cell_line1"], df["cell_line2"])}


def parse_resolved(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    records = []
    for _, row in df.iterrows():
        node = str(row["cell_line"])
        raw = "" if pd.isna(row.get("final_neighbors")) else str(row.get("final_neighbors"))
        neighbours = sorted({x.strip() for x in raw.split(";") if x.strip()})
        records.append({"cell_line": node, "neighbours": neighbours, "n_final": len(neighbours)})
    return pd.DataFrame(records)


def resolved_edge_set(df: pd.DataFrame):
    edges = set()
    for _, row in df.iterrows():
        node = row["cell_line"]
        for nb in row["neighbours"]:
            if nb != node:
                edges.add(tuple(sorted((node, nb))))
    return edges


def degree_map(edges):
    deg = {}
    for a, b in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    return deg


def pearson_from_dicts(a, b):
    keys = sorted(set(a) | set(b))
    if len(keys) < 2:
        return math.nan
    xa = [a.get(k, 0) for k in keys]
    xb = [b.get(k, 0) for k in keys]
    ma = sum(xa) / len(xa)
    mb = sum(xb) / len(xb)
    va = sum((x - ma) ** 2 for x in xa)
    vb = sum((x - mb) ** 2 for x in xb)
    if va == 0 or vb == 0:
        return math.nan
    cov = sum((x - ma) * (y - mb) for x, y in zip(xa, xb))
    return cov / math.sqrt(va * vb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metric_root", required=True)
    ap.add_argument("--metrics", default="pearson,jaccard")
    ap.add_argument("--directions", required=True)
    ap.add_argument("--resolved_pearson", required=True)
    ap.add_argument("--resolved_jaccard", required=True)
    ap.add_argument("--majority_pearson", required=True)
    ap.add_argument("--majority_jaccard", required=True)
    ap.add_argument("--union_pearson", required=True)
    ap.add_argument("--union_jaccard", required=True)
    ap.add_argument("--out_representation", required=True)
    ap.add_argument("--out_pairwise", required=True)
    ap.add_argument("--out_edge_agreement", required=True)
    ap.add_argument("--out_resolved_graph", required=True)
    ap.add_argument("--out_resolved_neighbours", required=True)
    ap.add_argument("--out_provenance", required=True)
    args = ap.parse_args()

    root = Path(args.metric_root)
    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    directions = [d.strip() for d in args.directions.split(",") if d.strip()]

    provenance_rows = []
    rep_rows = []
    pair_frames = {}
    threshold_by_metric_direction = {}
    edge_agreement_rows = []

    for metric in metrics:
        for direction in directions:
            base = root / metric / direction / "final_consensus"
            pairs_path = base / f"cell_line_similarity_pairs_{direction}.tsv"
            prov_path = base / f"cell_line_similarity_graph_provenance_{direction}.tsv"
            edges_path = base / f"cell_line_similarity_graph_edges_{direction}.tsv"
            node_path = base / f"cell_line_similarity_graph_node_summary_{direction}.tsv"

            pair_df = pd.read_csv(pairs_path, sep="\t")
            pair_frames[(metric, direction)] = pair_df
            prov_df = pd.read_csv(prov_path, sep="\t")
            provenance_rows.append(prov_df)
            threshold_by_metric_direction[(metric, direction)] = {
                "similarity_quantile": prov_df["similarity_quantile"].iloc[0],
                "computed_similarity_threshold": prov_df["computed_similarity_threshold"].iloc[0],
                "similarity_definition": prov_df["similarity_definition"].iloc[0],
                "similarity_question": prov_df["similarity_question"].iloc[0],
            }

            node_df = pd.read_csv(node_path, sep="\t")
            edge_df = read_edges(edges_path)
            rep_rows.append({
                "direction": direction,
                "similarity_metric": metric,
                "similarity_definition": prov_df["similarity_definition"].iloc[0],
                "n_nodes": int(prov_df["n_nodes"].iloc[0]),
                "n_candidate_pairs": int(prov_df["n_candidate_pairs"].iloc[0]),
                "n_defined_pairs": int(prov_df["n_defined_pairs"].iloc[0]),
                "n_na_pairs": int(prov_df["n_na_pairs"].iloc[0]),
                "similarity_quantile": prov_df["similarity_quantile"].iloc[0],
                "computed_similarity_threshold": prov_df["computed_similarity_threshold"].iloc[0],
                "n_pairs_gt_similarity_threshold": int(prov_df["n_pairs_gt_similarity_threshold"].iloc[0]),
                "n_pairs_eq_similarity_threshold": int(prov_df["n_pairs_eq_similarity_threshold"].iloc[0]),
                "n_edges": int(prov_df["selected_edge_count"].iloc[0]),
                "edge_fraction": prov_df["edge_fraction"].iloc[0],
                "density": prov_df["graph_density"].iloc[0],
                "n_isolates": int((node_df["degree"] == 0).sum()),
                "min_pairwise_active_tumours": prov_df["min_pairwise_active_tumours"].iloc[0],
                "pairwise_active_tumour_count_q1": prov_df["pairwise_active_tumour_count_q1"].iloc[0],
                "median_pairwise_active_tumours": prov_df["median_pairwise_active_tumours"].iloc[0],
                "pairwise_active_tumour_count_q3": prov_df["pairwise_active_tumour_count_q3"].iloc[0],
                "max_pairwise_active_tumours": prov_df["max_pairwise_active_tumours"].iloc[0],
                "n_pairs_lt_3_active_tumours": int(prov_df["n_pairs_lt_3_active_tumours"].iloc[0]),
                "n_pairs_lt_5_active_tumours": int(prov_df["n_pairs_lt_5_active_tumours"].iloc[0]),
                "n_pairs_lt_10_active_tumours": int(prov_df["n_pairs_lt_10_active_tumours"].iloc[0]),
            })

    for direction in directions:
        pearson_pairs = pair_frames[("pearson", direction)].rename(columns={"similarity": "pearson_similarity"})
        jaccard_pairs = pair_frames[("jaccard", direction)].rename(columns={"similarity": "jaccard_similarity"})
        merged = pearson_pairs.merge(
            jaccard_pairs[
                [
                    "cell_line1",
                    "cell_line2",
                    "jaccard_similarity",
                    "n_pairwise_active_tumours",
                    "n_shared_selected_tumours",
                    "undefined_similarity_reason",
                ]
            ].rename(
                columns={
                    "n_pairwise_active_tumours": "jaccard_n_pairwise_active_tumours",
                    "n_shared_selected_tumours": "jaccard_n_shared_selected_tumours",
                    "undefined_similarity_reason": "jaccard_undefined_similarity_reason",
                }
            ),
            on=["cell_line1", "cell_line2"],
            how="outer",
        )
        merged = merged.rename(
            columns={
                "n_pairwise_active_tumours": "pearson_n_pairwise_active_tumours",
                "n_shared_selected_tumours": "pearson_n_shared_selected_tumours",
                "undefined_similarity_reason": "pearson_undefined_similarity_reason",
            }
        )
        merged["direction"] = direction
        if "jaccard_n_pairwise_active_tumours" in merged.columns:
            merged["n_pairwise_active_tumours"] = merged["pearson_n_pairwise_active_tumours"].combine_first(
                merged["jaccard_n_pairwise_active_tumours"]
            )
        else:
            merged["n_pairwise_active_tumours"] = merged["pearson_n_pairwise_active_tumours"]
        if "jaccard_n_shared_selected_tumours" in merged.columns:
            merged["n_shared_selected_tumours"] = merged["pearson_n_shared_selected_tumours"].combine_first(
                merged["jaccard_n_shared_selected_tumours"]
            )
        else:
            merged["n_shared_selected_tumours"] = merged["pearson_n_shared_selected_tumours"]
        pair_frames[("merged", direction)] = merged

        pe = edge_set(read_edges(root / "pearson" / direction / "final_consensus" / f"cell_line_similarity_graph_edges_{direction}.tsv"))
        je = edge_set(read_edges(root / "jaccard" / direction / "final_consensus" / f"cell_line_similarity_graph_edges_{direction}.tsv"))
        union = pe | je
        edge_agreement_rows.append({
            "graph_level": "representation",
            "direction": direction,
            "pearson_edges": len(pe),
            "jaccard_edges": len(je),
            "shared_edges": len(pe & je),
            "pearson_only_edges": len(pe - je),
            "jaccard_only_edges": len(je - pe),
            "graph_edge_set_jaccard": (len(pe & je) / len(union)) if union else math.nan,
            "exact_edge_set_equality": pe == je,
        })

    pd.concat(provenance_rows, ignore_index=True).to_csv(args.out_provenance, sep="\t", index=False)
    pd.DataFrame(rep_rows).to_csv(args.out_representation, sep="\t", index=False)
    pairwise_frames = []
    for direction in directions:
        merged = pair_frames[("merged", direction)].copy()
        merged["pearson_similarity_quantile"] = threshold_by_metric_direction[("pearson", direction)]["similarity_quantile"]
        merged["jaccard_similarity_quantile"] = threshold_by_metric_direction[("jaccard", direction)]["similarity_quantile"]
        merged["pearson_computed_similarity_threshold"] = threshold_by_metric_direction[("pearson", direction)]["computed_similarity_threshold"]
        merged["jaccard_computed_similarity_threshold"] = threshold_by_metric_direction[("jaccard", direction)]["computed_similarity_threshold"]
        merged["pearson_selected"] = merged["pearson_similarity"] >= merged["pearson_computed_similarity_threshold"]
        merged.loc[merged["pearson_similarity"].isna(), "pearson_selected"] = False
        merged["jaccard_selected"] = merged["jaccard_similarity"] >= merged["jaccard_computed_similarity_threshold"]
        merged.loc[merged["jaccard_similarity"].isna(), "jaccard_selected"] = False
        merged["pearson_similarity_definition"] = threshold_by_metric_direction[("pearson", direction)]["similarity_definition"]
        merged["jaccard_similarity_definition"] = threshold_by_metric_direction[("jaccard", direction)]["similarity_definition"]
        merged["pearson_similarity_question"] = threshold_by_metric_direction[("pearson", direction)]["similarity_question"]
        merged["jaccard_similarity_question"] = threshold_by_metric_direction[("jaccard", direction)]["similarity_question"]
        pairwise_frames.append(merged)
    pd.concat(pairwise_frames, ignore_index=True).to_csv(args.out_pairwise, sep="\t", index=False)

    resolved_p = parse_resolved(Path(args.resolved_pearson))
    resolved_j = parse_resolved(Path(args.resolved_jaccard))
    resolved_merged = resolved_p.merge(resolved_j, on="cell_line", suffixes=("_pearson", "_jaccard"))
    resolved_merged["shared_neighbours"] = resolved_merged.apply(
        lambda r: ";".join(sorted(set(r["neighbours_pearson"]) & set(r["neighbours_jaccard"]))), axis=1
    )
    resolved_merged["pearson_only_neighbours"] = resolved_merged.apply(
        lambda r: ";".join(sorted(set(r["neighbours_pearson"]) - set(r["neighbours_jaccard"]))), axis=1
    )
    resolved_merged["jaccard_only_neighbours"] = resolved_merged.apply(
        lambda r: ";".join(sorted(set(r["neighbours_jaccard"]) - set(r["neighbours_pearson"]))), axis=1
    )
    resolved_merged["neighbour_set_jaccard"] = resolved_merged.apply(
        lambda r: (
            len(set(r["neighbours_pearson"]) & set(r["neighbours_jaccard"])) /
            len(set(r["neighbours_pearson"]) | set(r["neighbours_jaccard"]))
        ) if (set(r["neighbours_pearson"]) | set(r["neighbours_jaccard"])) else math.nan,
        axis=1,
    )
    resolved_merged.to_csv(args.out_resolved_neighbours, sep="\t", index=False)

    resolved_edges_p = resolved_edge_set(resolved_p)
    resolved_edges_j = resolved_edge_set(resolved_j)
    deg_p = degree_map(resolved_edges_p)
    deg_j = degree_map(resolved_edges_j)

    graph_levels = {
        "majority_threshold": (
            edge_set(read_edges(Path(args.majority_pearson))),
            edge_set(read_edges(Path(args.majority_jaccard))),
        ),
        "union_supported_edges": (
            edge_set(read_edges(Path(args.union_pearson))),
            edge_set(read_edges(Path(args.union_jaccard))),
        ),
        "resolved": (resolved_edges_p, resolved_edges_j),
    }

    resolved_rows = []
    for level, (pe, je) in graph_levels.items():
        union = pe | je
        resolved_rows.append({
            "graph_level": level,
            "pearson_edges": len(pe),
            "jaccard_edges": len(je),
            "shared_edges": len(pe & je),
            "pearson_only_edges": len(pe - je),
            "jaccard_only_edges": len(je - pe),
            "graph_edge_set_jaccard": (len(pe & je) / len(union)) if union else math.nan,
            "exact_edge_set_equality": pe == je,
            "degree_correlation": pearson_from_dicts(deg_p if level == "resolved" else degree_map(pe),
                                                    deg_j if level == "resolved" else degree_map(je)),
        })
    pd.DataFrame(resolved_rows).to_csv(args.out_resolved_graph, sep="\t", index=False)

    edge_agreement_rows.extend(resolved_rows)
    pd.DataFrame(edge_agreement_rows).to_csv(args.out_edge_agreement, sep="\t", index=False)


if __name__ == "__main__":
    main()
