#!/usr/bin/env python3
"""Static component panels and lightweight HTML inspection for pan-cancer graphs."""

from __future__ import annotations

import argparse
import html
import json
import math
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd


LINEAGE_COLOURS = {
    "BRCA": "#0072B2",
    "NBL": "#E69F00",
    "RBL": "#CC79A7",
    "HEME": "#009E73",
    "Normal": "#7F7F7F",
    "Unknown": "#999999",
}


def norm_id(value: object) -> str:
    return "".join(ch for ch in str(value).upper() if ch.isalnum())


def read_tsv(path: str | None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def build_display_maps(display_df: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    long_to_short: dict[str, str] = {}
    short_to_long: dict[str, str] = {}
    if {"long_id", "short_id"}.issubset(display_df.columns):
        for _, row in display_df.iterrows():
            long_id = str(row["long_id"])
            short_id = str(row["short_id"])
            long_to_short[long_id] = short_id
            short_to_long[short_id] = long_id
    return long_to_short, short_to_long


def infer_lineage(node: str, label: str, metadata: pd.DataFrame, short_to_long: dict[str, str]) -> str:
    candidates = [node, label, short_to_long.get(node, ""), short_to_long.get(label, "")]
    candidates = [norm_id(x) for x in candidates if str(x).strip()]
    if not metadata.empty and {"sample_id", "cancer_type"}.issubset(metadata.columns):
        rows = []
        norm_samples = metadata["_norm_sample_id"].tolist()
        for candidate in candidates:
            if len(candidate) < 3:
                continue
            rows.extend(i for i, sample in enumerate(norm_samples) if candidate in sample)
        if rows:
            values = [str(metadata.iloc[i]["cancer_type"]) for i in sorted(set(rows))]
            counts = Counter(values)
            if len(counts) == 1:
                return next(iter(counts))
            return counts.most_common(1)[0][0]
    joined = " ".join([node, label]).upper()
    if any(token in joined for token in ("RBL", "WERI", "Y_79", "Y79")):
        return "RBL"
    if any(token in joined for token in ("CHP", "LAN", "KELLY", "SIMA", "NBL", "GIME", "GI_ME", "SKNBE", "SHSY", "MHHNB", "CCLF", "NGP", "IMR32", "NMB")):
        return "NBL"
    if any(token in joined for token in ("BT_474", "CAL_", "EFM", "MDA", "HCC", "MCF", "T_47D", "SK_BR", "JIMT", "DU_4475", "COLO_824", "HS_578T", "EVSA", "KPL", "MFM", "IPH")):
        return "BRCA"
    if "NC_NC" in joined or "NC-NC" in joined:
        return "Normal"
    return "Unknown"


def read_graph(args: argparse.Namespace) -> tuple[nx.Graph, pd.DataFrame]:
    edges = read_tsv(args.edges)
    node_stats = read_tsv(args.node_stats)
    display = read_tsv(args.display_names)
    metadata = read_tsv(args.metadata)
    if not metadata.empty and "sample_id" in metadata.columns:
        metadata["_norm_sample_id"] = metadata["sample_id"].map(norm_id)

    long_to_short, short_to_long = build_display_maps(display)
    G = nx.Graph()
    for _, row in edges.iterrows():
        n1 = str(row["node1"])
        n2 = str(row["node2"])
        attrs = {k: row[k] for k in edges.columns if k not in {"node1", "node2"}}
        G.add_edge(n1, n2, **attrs)

    known_nodes: set[str] = set(G.nodes())
    if "cell_line" in node_stats.columns:
        known_nodes.update(node_stats["cell_line"].astype(str))
    if "short_id" in display.columns:
        known_nodes.update(display["short_id"].astype(str))
    for node in sorted(known_nodes):
        G.add_node(node)

    stats_by_node = {}
    if "cell_line" in node_stats.columns:
        stats_by_node = {str(row["cell_line"]): row for _, row in node_stats.iterrows()}

    rows = []
    components = sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)[0]))
    component_by_node = {}
    for idx, comp in enumerate(components, start=1):
        for node in comp:
            component_by_node[node] = idx

    for node in sorted(G.nodes()):
        row = stats_by_node.get(node)
        display_label = str(row["cell_line_short"]) if row is not None and "cell_line_short" in node_stats.columns else long_to_short.get(node, node)
        lineage = infer_lineage(node, display_label, metadata, short_to_long)
        degree = int(G.degree(node))
        betweenness = float(row["betweenness_unnormalised"]) if row is not None and "betweenness_unnormalised" in node_stats.columns and pd.notna(row["betweenness_unnormalised"]) else 0.0
        anchor = False
        if row is not None:
            for col in ("anchor_selected", "canonical_selected", "canonical_bridge_selected", "most_connected_selected", "is_central"):
                if col in node_stats.columns and parse_bool(row[col]):
                    anchor = True
        rows.append({
            "node_id": node,
            "display_label": display_label,
            "component_id": component_by_node[node],
            "component_size": len(components[component_by_node[node] - 1]),
            "lineage": lineage,
            "degree": degree,
            "betweenness_unnormalised": betweenness,
            "is_isolate": degree == 0,
            "anchor_or_focal": anchor,
        })
    node_df = pd.DataFrame(rows)
    nx.set_node_attributes(G, node_df.set_index("node_id").to_dict("index"))
    return G, node_df


def lineage_composition(labels: list[str]) -> str:
    counts = Counter(labels)
    order = sorted(counts, key=lambda x: (-counts[x], x))
    return ";".join(f"{k}={counts[k]}" for k in order)


def write_tables(G: nx.Graph, node_df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    component_rows = []
    for comp_id, group in node_df.groupby("component_id", sort=True):
        nodes = group["node_id"].tolist()
        sub = G.subgraph(nodes)
        lineage_counts = Counter(group["lineage"])
        majority_lineage, majority_n = lineage_counts.most_common(1)[0]
        component_rows.append({
            "graph_type": args.graph_type,
            "component_id": comp_id,
            "component_size": len(nodes),
            "n_edges": sub.number_of_edges(),
            "density": nx.density(sub) if len(nodes) > 1 else 0.0,
            "lineage_composition": lineage_composition(group["lineage"].tolist()),
            "majority_lineage": majority_lineage,
            "purity": majority_n / len(nodes),
            "n_isolates": int(group["is_isolate"].sum()),
            "n_anchor_or_focal": int(group["anchor_or_focal"].sum()),
            "nodes": ";".join(group.sort_values("display_label")["display_label"]),
        })
    comp_df = pd.DataFrame(component_rows)
    node_df.to_csv(args.node_labels_out, sep="\t", index=False)
    comp_df.to_csv(args.component_summary_out, sep="\t", index=False)
    return comp_df


def component_layout(sub: nx.Graph, seed: int) -> dict[str, np.ndarray]:
    n = sub.number_of_nodes()
    if n == 1:
        node = next(iter(sub.nodes()))
        return {node: np.array([0.0, 0.0])}
    if n == 2:
        nodes = sorted(sub.nodes())
        return {nodes[0]: np.array([-0.8, 0.0]), nodes[1]: np.array([0.8, 0.0])}
    k = 1.4 / math.sqrt(max(n, 1))
    pos = nx.spring_layout(sub, seed=seed, k=k, iterations=400, weight=None)
    coords = np.vstack([pos[n] for n in sub.nodes()])
    span = np.maximum(coords.max(axis=0) - coords.min(axis=0), 1e-6)
    scale = np.array([2.8 / span[0], 2.2 / span[1]])
    centre = coords.mean(axis=0)
    return {node: (np.array(p) - centre) * scale for node, p in pos.items()}


def draw_component_panel(ax, G: nx.Graph, node_df: pd.DataFrame, comp_id: int, title_prefix: str) -> None:
    group = node_df[node_df["component_id"] == comp_id].copy()
    nodes = group["node_id"].tolist()
    sub = G.subgraph(nodes).copy()
    pos = component_layout(sub, seed=1000 + int(comp_id))
    edge_support = [float(sub.edges[e].get("support_directions", 1) or 1) for e in sub.edges()]
    max_support = max(edge_support) if edge_support else 1.0
    for (u, v), support in zip(sub.edges(), edge_support):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        lw = 0.8 + 2.8 * support / max_support
        ax.plot([x0, x1], [y0, y1], color="#444444", lw=lw, alpha=0.65, zorder=1)
    for lineage, rows in group.groupby("lineage", sort=True):
        xs = [pos[n][0] for n in rows["node_id"]]
        ys = [pos[n][1] for n in rows["node_id"]]
        ax.scatter(xs, ys, s=380, c=LINEAGE_COLOURS.get(lineage, LINEAGE_COLOURS["Unknown"]),
                   edgecolors="#1A1A2E", linewidths=1.0,
                   marker="D" if bool(rows["is_isolate"].all()) else "o",
                   label=lineage, zorder=3)
    for _, row in group.iterrows():
        x, y = pos[row["node_id"]]
        if row["anchor_or_focal"]:
            ax.scatter([x], [y], s=620, facecolors="none", edgecolors="#111111", linewidths=2.1, zorder=4)
        ax.text(x, y, row["display_label"], ha="center", va="center",
                fontsize=8.5, fontweight="bold", color="#111827", zorder=5,
                bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.80))
    comp_lineage = lineage_composition(group["lineage"].tolist())
    ax.set_title(f"{title_prefix} C{comp_id}: n={len(nodes)}; {comp_lineage}", fontsize=10, loc="left")
    ax.axis("off")
    ax.set_aspect("equal")


def save_component_panels(G: nx.Graph, node_df: pd.DataFrame, args: argparse.Namespace) -> None:
    component_ids = sorted(node_df["component_id"].unique(), key=lambda c: (-int((node_df["component_id"] == c).sum()), int(c)))
    pages = [component_ids[i:i + 6] for i in range(0, len(component_ids), 6)]
    first_fig = None
    with PdfPages(args.panels_pdf) as pdf:
        for page_i, comps in enumerate(pages, start=1):
            fig, axes = plt.subplots(3, 2, figsize=(14, 18))
            axes = axes.ravel()
            for ax, comp_id in zip(axes, comps):
                draw_component_panel(ax, G, node_df, int(comp_id), args.graph_type)
            for ax in axes[len(comps):]:
                ax.axis("off")
            handles = [
                Line2D([], [], marker="o", color="w", markerfacecolor=colour,
                       markeredgecolor="#1A1A2E", markersize=8, label=lineage)
                for lineage, colour in LINEAGE_COLOURS.items()
                if lineage in set(node_df["lineage"])
            ]
            fig.legend(handles=handles, loc="lower center", ncol=min(6, len(handles)), frameon=False)
            fig.suptitle(f"{args.label}: component panels, page {page_i}/{len(pages)}",
                         fontsize=14, fontweight="bold")
            fig.tight_layout(rect=[0.02, 0.04, 0.98, 0.96])
            pdf.savefig(fig)
            if page_i == 1:
                first_fig = fig
                fig.savefig(args.panels_png, dpi=220)
                fig.savefig(args.panels_svg)
            else:
                plt.close(fig)
        if first_fig is not None:
            plt.close(first_fig)


def html_colour(lineage: str) -> str:
    return LINEAGE_COLOURS.get(lineage, LINEAGE_COLOURS["Unknown"])


def save_interactive_html(G: nx.Graph, node_df: pd.DataFrame, comp_df: pd.DataFrame, args: argparse.Namespace) -> None:
    components = sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)[0]))
    positions = {}
    x_cursor = 80.0
    y_cursor = 80.0
    max_row_h = 0.0
    for idx, comp in enumerate(components, start=1):
        sub = G.subgraph(comp).copy()
        local = component_layout(sub, seed=2000 + idx)
        n = len(comp)
        box_w = max(220.0, 95.0 * math.sqrt(n))
        box_h = max(160.0, 75.0 * math.sqrt(n))
        if x_cursor + box_w > 1800:
            x_cursor = 80.0
            y_cursor += max_row_h + 120.0
            max_row_h = 0.0
        for node, p in local.items():
            positions[node] = (x_cursor + box_w / 2.0 + p[0] * box_w / 4.0,
                               y_cursor + box_h / 2.0 + p[1] * box_h / 4.0)
        x_cursor += box_w + 90.0
        max_row_h = max(max_row_h, box_h)
    width = 1900
    height = int(max(900, y_cursor + max_row_h + 140))
    node_records = node_df.set_index("node_id").to_dict("index")
    edge_support = [float(G.edges[e].get("support_directions", 1) or 1) for e in G.edges()]
    max_support = max(edge_support) if edge_support else 1.0
    edge_elems = []
    for u, v in G.edges():
        x0, y0 = positions[u]
        x1, y1 = positions[v]
        support = float(G.edges[u, v].get("support_directions", 1) or 1)
        edge_elems.append(
            f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" '
            f'stroke="#555" stroke-opacity="0.48" stroke-width="{0.8 + 4.0 * support / max_support:.2f}">'
            f'<title>{html.escape(str(u))} -- {html.escape(str(v))}; support={support:g}</title></line>'
        )
    node_elems = []
    label_elems = []
    for node, rec in node_records.items():
        x, y = positions[node]
        radius = 9 if rec["is_isolate"] else 11
        tooltip = "; ".join([
            f"node={rec['display_label']}",
            f"lineage={rec['lineage']}",
            f"component=C{rec['component_id']}",
            f"degree={rec['degree']}",
            f"betweenness={float(rec['betweenness_unnormalised']):.3g}",
            f"anchor={rec['anchor_or_focal']}",
            f"isolate={rec['is_isolate']}",
        ])
        node_elems.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius}" fill="{html_colour(rec["lineage"])}" '
            f'stroke="#1A1A2E" stroke-width="{2.4 if rec["anchor_or_focal"] else 1.2}">'
            f'<title>{html.escape(tooltip)}</title></circle>'
        )
        label_elems.append(
            f'<text x="{x:.2f}" y="{y - radius - 3:.2f}" text-anchor="middle" '
            f'font-size="10" font-weight="700" paint-order="stroke" stroke="white" stroke-width="3" '
            f'fill="#111827">{html.escape(str(rec["display_label"]))}</text>'
        )
    comp_labels = []
    for _, row in comp_df.iterrows():
        nodes = node_df[node_df["component_id"] == row["component_id"]]["node_id"]
        xs = [positions[n][0] for n in nodes]
        ys = [positions[n][1] for n in nodes]
        comp_labels.append(
            f'<text x="{np.mean(xs):.2f}" y="{min(ys) - 36:.2f}" text-anchor="middle" '
            f'font-size="15" font-weight="800" fill="#111827" '
            f'paint-order="stroke" stroke="white" stroke-width="4">C{int(row["component_id"])} '
            f'n={int(row["component_size"])} {html.escape(str(row["lineage_composition"]))}</text>'
        )
    html_text = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>{html.escape(args.label)}</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 0; }}
#bar {{ padding: 10px 14px; border-bottom: 1px solid #ddd; background: #f7f7f7; position: sticky; top: 0; z-index: 2; }}
svg {{ width: 100vw; height: calc(100vh - 48px); display: block; background: white; cursor: grab; }}
svg:active {{ cursor: grabbing; }}
</style></head><body>
<div id="bar"><strong>{html.escape(args.label)}</strong>
 | wheel to zoom, drag to pan, hover nodes/edges for metadata. Communities/components describe graph organisation, not ground-truth lineage labels.</div>
<svg id="canvas" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">
<g id="viewport">
{''.join(edge_elems)}
{''.join(node_elems)}
{''.join(label_elems)}
{''.join(comp_labels)}
</g></svg>
<script>
const svg = document.getElementById('canvas');
const vp = document.getElementById('viewport');
let scale = 1, tx = 0, ty = 0, dragging = false, last = null;
function update() {{ vp.setAttribute('transform', `translate(${{tx}},${{ty}}) scale(${{scale}})`); }}
svg.addEventListener('wheel', (e) => {{
  e.preventDefault();
  const delta = e.deltaY < 0 ? 1.12 : 0.89;
  scale = Math.min(8, Math.max(0.25, scale * delta));
  update();
}});
svg.addEventListener('mousedown', (e) => {{ dragging = true; last = [e.clientX, e.clientY]; }});
window.addEventListener('mouseup', () => {{ dragging = false; }});
window.addEventListener('mousemove', (e) => {{
  if (!dragging || !last) return;
  tx += (e.clientX - last[0]);
  ty += (e.clientY - last[1]);
  last = [e.clientX, e.clientY];
  update();
}});
</script></body></html>
"""
    Path(args.interactive_html).write_text(html_text, encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph-type", required=True, choices=["resolved", "support_threshold"])
    ap.add_argument("--label", default="MULTICOHORT_CANCER")
    ap.add_argument("--edges", required=True)
    ap.add_argument("--node-stats", default=None)
    ap.add_argument("--display-names", default=None)
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--panels-pdf", required=True)
    ap.add_argument("--panels-png", required=True)
    ap.add_argument("--panels-svg", required=True)
    ap.add_argument("--interactive-html", required=True)
    ap.add_argument("--component-summary-out", required=True)
    ap.add_argument("--node-labels-out", required=True)
    args = ap.parse_args()

    for path in [
        args.panels_pdf, args.panels_png, args.panels_svg,
        args.interactive_html, args.component_summary_out, args.node_labels_out,
    ]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    G, node_df = read_graph(args)
    comp_df = write_tables(G, node_df, args)
    save_component_panels(G, node_df, args)
    save_interactive_html(G, node_df, comp_df, args)
    print(
        f"[OK] {args.graph_type}: nodes={G.number_of_nodes()} edges={G.number_of_edges()} "
        f"components={nx.number_connected_components(G)} outputs={args.panels_pdf}, {args.interactive_html}"
    )


if __name__ == "__main__":
    main()
