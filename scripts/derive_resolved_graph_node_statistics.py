#!/usr/bin/env python3
"""
Derive analytical node statistics for the resolved cell-line graph.

The unit of analysis is a biological cell-line node in the resolved graph.
The script converts resolved-neighbour assignments into deterministic graph
edge, node-statistics, and anchor-centrality audit tables. It does not calculate
node layout, label positions, figure dimensions, or legend settings.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from visualize_resolved_dsmz_graph import (  # noqa: E402
    build_shortname_map,
    canonical_cell_line_id,
    canonical_neighbour_string,
    compute_component_centrality,
    edges_from_resolved_neighbours,
    parse_neighbours,
)


def require_columns(table: pd.DataFrame, path: Path, columns: list[str]) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise SystemExit(
            f"[Graph node statistics] {path} is missing required column(s): "
            + ", ".join(missing)
        )


def read_resolved_neighbour_table(path: Path) -> tuple[pd.DataFrame, str]:
    table = pd.read_csv(path, sep="\t", dtype=str)
    require_columns(table, path, ["cell_line"])
    if "final_neighbors" in table.columns:
        neighbour_column = "final_neighbors"
    elif "final_neighbours" in table.columns:
        neighbour_column = "final_neighbours"
    else:
        raise SystemExit(
            "[Graph node statistics] resolved-neighbour table must contain "
            "final_neighbors or final_neighbours"
        )
    table["cell_line"] = table["cell_line"].astype(str).str.strip().map(canonical_cell_line_id)
    table[neighbour_column] = (
        table[neighbour_column]
        .fillna("")
        .astype(str)
        .map(canonical_neighbour_string)
    )
    table = table[table["cell_line"].astype(bool)].copy()
    if table.empty:
        raise SystemExit("[Graph node statistics] resolved-neighbour table has no cell-line nodes")
    if table["cell_line"].duplicated().any():
        aggregation = {
            neighbour_column: lambda values: canonical_neighbour_string(
                ";".join(str(value) for value in values if str(value).strip())
            )
        }
        for column in ("reps_used", "best_overall_dir", "winner_dir"):
            if column in table.columns:
                aggregation[column] = lambda values: ";".join(
                    sorted({str(value).strip() for value in values if str(value).strip()})
                )
        table = table.groupby("cell_line", as_index=False).agg(aggregation)
    return table, neighbour_column


def build_resolved_graph(table: pd.DataFrame, neighbour_column: str) -> tuple[nx.Graph, int, set[str]]:
    node_set = set(table["cell_line"].unique())
    graph = nx.Graph()
    graph.add_nodes_from(sorted(node_set))
    self_loops = 0
    ignored_neighbours: set[str] = set()
    for row in table.itertuples(index=False):
        source = str(getattr(row, "cell_line")).strip()
        for target in parse_neighbours(str(getattr(row, neighbour_column) or "")):
            if target == source:
                self_loops += 1
                continue
            if target not in node_set:
                ignored_neighbours.add(target)
                continue
            graph.add_edge(source, target)
    return graph, self_loops, ignored_neighbours


def aggregate_resolved_edges(
    table: pd.DataFrame,
    neighbour_column: str,
    graph: nx.Graph,
    short_map: dict[str, str],
    resolved_table_path: Path,
    output_edges_path: Path,
) -> None:
    node_set = set(graph.nodes())
    tumour_neighbourhood_root = resolved_table_path.parent.parent
    edge_records: list[dict[str, object]] = []
    if "best_overall_dir" in table.columns and "winner_dir" in table.columns:
        directions_used: set[str] = set()
        for row in table.itertuples(index=False):
            best_overall = getattr(row, "best_overall_dir", "")
            winner = getattr(row, "winner_dir", "")
            if pd.notna(best_overall) and str(best_overall).strip():
                directions_used.add(str(best_overall).strip())
            if pd.notna(winner) and str(winner).strip():
                directions_used.add(str(winner).strip())
        if tumour_neighbourhood_root.exists():
            for direction_dir in tumour_neighbourhood_root.iterdir():
                if direction_dir.is_dir() and direction_dir.name != "final_consensus_all":
                    direction = direction_dir.name
                    edge_file = (
                        direction_dir
                        / "final_consensus"
                        / f"cell_line_similarity_graph_edges_{direction}.tsv"
                    )
                    if edge_file.exists():
                        directions_used.add(direction)
        for direction in sorted(directions_used):
            edge_file = (
                tumour_neighbourhood_root
                / direction
                / "final_consensus"
                / f"cell_line_similarity_graph_edges_{direction}.tsv"
            )
            if not edge_file.exists():
                continue
            direction_edges = pd.read_csv(edge_file, sep="\t")
            if not {"cell_line1", "cell_line2"}.issubset(direction_edges.columns):
                continue
            for edge_row in direction_edges.itertuples(index=False):
                node1 = canonical_cell_line_id(str(getattr(edge_row, "cell_line1")).strip())
                node2 = canonical_cell_line_id(str(getattr(edge_row, "cell_line2")).strip())
                if node1 == node2 or node1 not in node_set or node2 not in node_set:
                    continue
                if not graph.has_edge(node1, node2):
                    continue
                if node1 > node2:
                    node1, node2 = node2, node1
                similarity = getattr(edge_row, "similarity", 1.0)
                try:
                    similarity = float(similarity)
                except (TypeError, ValueError):
                    similarity = 1.0
                edge_records.append(
                    {
                        "node1": node1,
                        "node2": node2,
                        "direction": direction,
                        "similarity": similarity,
                    }
                )

    if edge_records:
        edges_table = pd.DataFrame(edge_records)
        methods_map: dict[str, str] = {}
        for row in table.itertuples(index=False):
            cell_line = str(getattr(row, "cell_line")).strip()
            methods = str(getattr(row, "reps_used", "") or "").strip()
            if not methods:
                best_overall = str(getattr(row, "best_overall_dir", "") or "").strip()
                winner = str(getattr(row, "winner_dir", "") or "").strip()
                if best_overall or winner:
                    methods = f"best_overall={best_overall};winner={winner}"
            methods_map[cell_line] = methods
        aggregated_rows: list[dict[str, object]] = []
        for (node1, node2), group in edges_table.groupby(["node1", "node2"], sort=True):
            directions = sorted(group["direction"].astype(str).unique())
            similarities = pd.to_numeric(group["similarity"], errors="coerce").dropna()
            aggregated_rows.append(
                {
                    "node1": node1,
                    "node2": node2,
                    "support_directions": len(directions),
                    "support_weight_mean": float(similarities.mean()) if len(similarities) else 1.0,
                    "support_weight_sum": float(similarities.sum()) if len(similarities) else 1.0,
                    "support_weight_max": float(similarities.max()) if len(similarities) else 1.0,
                    "methods_union": methods_map.get(node1, methods_map.get(node2, "")),
                    "node1_short": short_map.get(node1, node1),
                    "node2_short": short_map.get(node2, node2),
                }
            )
        output_edges = pd.DataFrame(aggregated_rows)
    else:
        neighbour_column_edges = edges_from_resolved_neighbours(table, neighbour_column, node_set)
        output_edges = pd.DataFrame(
            neighbour_column_edges,
            columns=[
                "node1",
                "node2",
                "support_directions",
                "support_weight_mean",
                "support_weight_sum",
                "support_weight_max",
                "methods_union",
            ],
        )
        if not output_edges.empty:
            output_edges["node1_short"] = output_edges["node1"].map(lambda value: short_map.get(value, value))
            output_edges["node2_short"] = output_edges["node2"].map(lambda value: short_map.get(value, value))

    output_edges = output_edges.sort_values(["node1", "node2"], kind="mergesort")
    output_edges_path.parent.mkdir(parents=True, exist_ok=True)
    output_edges.to_csv(output_edges_path, sep="\t", index=False)
    print(
        f"[Graph node statistics] Wrote {len(output_edges)} resolved graph edges "
        f"from {resolved_table_path}"
    )


def write_node_statistics_and_anchor_audit(
    graph: nx.Graph,
    short_map: dict[str, str],
    output_node_stats_path: Path,
    output_anchor_audit_path: Path,
    cohort: str,
    center_top_n: int,
) -> None:
    isolates = sorted(node for node in graph.nodes() if graph.degree(node) == 0)
    graph_non_isolates = graph.subgraph([node for node in graph.nodes() if node not in set(isolates)]).copy()
    (
        components,
        anchor_nodes,
        betweenness_unnormalised,
        betweenness_normalised,
        selected_by_unnormalised,
        selected_by_normalised,
        selected_by_degree,
    ) = compute_component_centrality(graph_non_isolates, center_top_n=center_top_n)
    for node in graph.nodes():
        betweenness_unnormalised.setdefault(node, 0.0)
        betweenness_normalised.setdefault(node, 0.0)

    component_by_node = {node: component_index for component_index, component in enumerate(components) for node in component}
    for isolate_index, node in enumerate(isolates, start=1):
        component_by_node[node] = -isolate_index

    anchor_set = set(anchor_nodes)
    rows: list[dict[str, object]] = []
    for node in sorted(graph.nodes()):
        component_id = component_by_node.get(node, -1)
        selected_unnormalised = selected_by_unnormalised.get(component_id)
        selected_normalised = selected_by_normalised.get(component_id)
        selected_degree = selected_by_degree.get(component_id)
        rows.append(
            {
                "cell_line": node,
                "cell_line_short": short_map.get(node, node),
                "degree": int(graph.degree(node)),
                "betweenness": float(betweenness_unnormalised.get(node, 0.0)),
                "betweenness_unnormalised": float(betweenness_unnormalised.get(node, 0.0)),
                "betweenness_normalised": float(betweenness_normalised.get(node, 0.0)),
                "component": int(component_id),
                "is_isolate": bool(node in set(isolates)),
                "is_central": bool(node in anchor_set),
                "selected_by_unnormalised": bool(node == selected_unnormalised) if selected_unnormalised else False,
                "selected_by_normalised": bool(node == selected_normalised) if selected_normalised else False,
                "selected_by_degree": bool(node == selected_degree) if selected_degree else False,
                "canonical_selected": bool(node == selected_unnormalised) if selected_unnormalised else False,
                "canonical_bridge_selected": bool(node == selected_unnormalised) if selected_unnormalised else False,
                "most_connected_selected": bool(node == selected_degree) if selected_degree else False,
                "centrality_metric": "degree_and_betweenness",
                "centrality_normalised": False,
                "centrality_weighted": False,
                "centrality_scope": "within_component",
                "centrality_tie_break": "highest_value_then_node_id",
            }
        )

    node_stats = pd.DataFrame(rows)
    output_node_stats_path.parent.mkdir(parents=True, exist_ok=True)
    node_stats.to_csv(output_node_stats_path, sep="\t", index=False)

    anchor_audit = node_stats.rename(
        columns={
            "cell_line": "node_id",
            "cell_line_short": "display_label",
            "component": "component_id",
        }
    )[
        [
            "component_id",
            "node_id",
            "display_label",
            "degree",
            "betweenness_normalised",
            "betweenness_unnormalised",
            "selected_by_normalised",
            "selected_by_unnormalised",
            "selected_by_degree",
            "canonical_selected",
            "canonical_bridge_selected",
            "most_connected_selected",
        ]
    ]
    anchor_audit.insert(0, "cohort", cohort)
    anchor_audit["changed_relative_to_legacy"] = (
        anchor_audit["selected_by_normalised"] != anchor_audit["selected_by_unnormalised"]
    )
    anchor_audit["degree_vs_betweenness_node_differs"] = (
        anchor_audit["selected_by_degree"] != anchor_audit["selected_by_unnormalised"]
    )
    anchor_audit["canonical_metric"] = "degree_and_unnormalised_betweenness"
    anchor_audit["alternative_metric"] = "normalised_betweenness"
    anchor_audit["centrality_scope"] = "within_component"
    anchor_audit["centrality_weighted"] = False
    anchor_audit["centrality_tie_break"] = "highest_value_then_node_id"
    anchor_audit["degree_anchor_selected"] = anchor_audit["most_connected_selected"]
    anchor_audit["bridge_betweenness_selected"] = anchor_audit["canonical_bridge_selected"]
    anchor_audit["anchor_selected"] = (
        anchor_audit["degree_anchor_selected"] | anchor_audit["bridge_betweenness_selected"]
    )

    def anchor_reason(row: pd.Series) -> str:
        by_degree = bool(row["degree_anchor_selected"])
        by_betweenness = bool(row["bridge_betweenness_selected"])
        if by_degree and by_betweenness:
            return "degree;betweenness_unnormalised"
        if by_degree:
            return "degree"
        if by_betweenness:
            return "betweenness_unnormalised"
        return ""

    anchor_audit["anchor_selection_reason"] = anchor_audit.apply(anchor_reason, axis=1)
    output_anchor_audit_path.parent.mkdir(parents=True, exist_ok=True)
    anchor_audit.to_csv(output_anchor_audit_path, sep="\t", index=False)
    print(
        f"[Graph node statistics] Wrote {len(node_stats)} node rows and "
        f"{int(anchor_audit['anchor_selected'].sum())} selected anchors"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolved-neighbours", required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--out-shortnames", required=True)
    parser.add_argument("--out-edges", required=True)
    parser.add_argument("--out-node-stats", required=True)
    parser.add_argument("--out-anchor-audit", required=True)
    parser.add_argument("--center-top-n", type=int, default=1)
    args = parser.parse_args()

    resolved_path = Path(args.resolved_neighbours)
    resolved_table, neighbour_column = read_resolved_neighbour_table(resolved_path)
    graph, self_loops, ignored_neighbours = build_resolved_graph(resolved_table, neighbour_column)
    short_map = build_shortname_map(graph.nodes())
    if len(set(short_map.values())) != len(short_map.values()):
        raise SystemExit("[Graph node statistics] non-unique display-name identifiers")

    shortname_rows = [
        {"long_id": long_id, "short_id": short_id}
        for long_id, short_id in sorted(short_map.items())
    ]
    out_shortnames = Path(args.out_shortnames)
    out_shortnames.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(shortname_rows).to_csv(out_shortnames, sep="\t", index=False)
    aggregate_resolved_edges(
        resolved_table,
        neighbour_column,
        graph,
        short_map,
        resolved_path,
        Path(args.out_edges),
    )
    write_node_statistics_and_anchor_audit(
        graph,
        short_map,
        Path(args.out_node_stats),
        Path(args.out_anchor_audit),
        args.cohort,
        args.center_top_n,
    )
    print(
        "[Graph node statistics] "
        f"nodes={graph.number_of_nodes()} edges={graph.number_of_edges()} "
        f"components={nx.number_connected_components(graph)} "
        f"isolates={sum(1 for node in graph.nodes() if graph.degree(node) == 0)} "
        f"self_loops_ignored={self_loops} "
        f"external_neighbours_ignored={len(ignored_neighbours)}"
    )


if __name__ == "__main__":
    main()
