#!/usr/bin/env python3
"""
visualize_resolved_dsmz_graph.py
================================
Thesis figure for the resolved cell-line similarity / neighbourhood graph
for one cohort.

Visual language (single source of truth: scripts/graph_plot_style.py):
  * Node fill encodes connected component using a colour-blind-friendly
    palette with no green hue; isolates share one neutral grey.
    Component colours are structural identifiers only.
  * Isolates (degree 0) are drawn as diamonds inside a compact labelled
    "Isolates (degree 0)" strip below the rest of the layout.
  * Bridge-like anchor nodes (maximum unnormalised betweenness within each
    component) are drawn with non-fill rings.
  * Most-connected nodes (maximum degree within each component) are drawn
    with non-fill square outlines. Nodes satisfying both criteria carry both
    overlays.
  * Edge width is proportional to support count (support_directions).
    Dashed style flags single-supported (support == 1) or fallback edges
    that have no direction-based support; solid for multi-supported edges.
  * Labels are short names with a thin white halo (no boxed background).
  * A legend documents every encoding actually present, placed outside the
    axes so it never occludes nodes.

The script preserves all declared Snakemake rule outputs:
  * Figure trio: {out_prefix}.{png,pdf,svg}
  * patient_referenced_cell_line_display_names.tsv, or the profile-specific
    multicohort_cancer_* equivalent for the multicohort_cancer profile
  * patient_referenced_resolved_cell_line_neighbourhood_graph_edges.tsv
  * patient_referenced_resolved_cell_line_neighbourhood_graph_node_stats.tsv
  * patient_referenced_resolved_cell_line_neighbourhood_anchor_centrality_audit.tsv

Centrality selection is unchanged (highest unnormalised betweenness and
highest degree within each component, ties broken by lexicographically
largest node id). This script changes only plotting layout and overlay
rendering.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

mpl.rcParams["text.usetex"] = False
mpl.rcParams["mathtext.default"] = "regular"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import graph_plot_style as gps  # noqa: E402


# =============================================================================
# Helpers (unchanged in spirit; preserved for side-output schemas)
# =============================================================================

def parse_neighbours(s):
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return []
    s = str(s).strip()
    if not s:
        return []
    return [x.strip() for x in s.split(";") if x.strip()]


def canonical_cell_line_id(value: str) -> str:
    """Map multicohort profile-level IDs to biological cell-line/group IDs."""
    node = str(value).strip()
    m = re.match(r"^NG[-_][^_]+_(.+?)_lib", node)
    if m:
        node = m.group(1)
    node = re.sub(r"^(RBL_\d+)_\d+$", r"\1", node)
    return node


def canonical_neighbour_string(value: str) -> str:
    nodes = [canonical_cell_line_id(x) for x in parse_neighbours(value)]
    nodes = sorted({x for x in nodes if x})
    return ";".join(nodes)


def build_shortname_map(node_ids):
    """
    RBL naming convention:
      base = RBL_<num> from long_id (or WERI_RB1 / Y_79)
      unique base -> short_id = base; otherwise base_<trailing_digits>.
    For non-RBL cohorts the long_id is returned unchanged (typical for
    BRCA/NBL where IDs are already short).
    """
    node_ids = sorted(set(node_ids))
    base_for: dict[str, str] = {}
    for s in node_ids:
        m = re.search(r"_RBL_(\d+)_", s)
        if m:
            base_for[s] = f"RBL_{m.group(1)}"
            continue
        if "_WERI_RB1_" in s:
            base_for[s] = "WERI_RB1"
            continue
        if "_Y_79_" in s:
            base_for[s] = "Y_79"
            continue
        base_for[s] = s

    by_base = defaultdict(list)
    for long_id, base in base_for.items():
        by_base[base].append(long_id)

    short_map: dict[str, str] = {}
    for base, longs in by_base.items():
        longs = sorted(longs)
        if len(longs) == 1:
            short_map[longs[0]] = base
            continue
        for long_id in longs:
            last = str(long_id).split("_")[-1]
            if not last.isdigit():
                raise ValueError(
                    f"[shortname] Duplicate base {base} but long_id has "
                    f"no numeric trailing token: {long_id}"
                )
            short_map[long_id] = f"{base}_{last}"
    return short_map


def edges_from_resolved_neighbours(df, nb_col, node_set):
    edges = []
    seen = set()
    for _, row in df.iterrows():
        src = str(row.get("cell_line", "")).strip()
        if not src:
            continue
        for tgt in parse_neighbours(str(row.get(nb_col, "") or "")):
            if not tgt or tgt.upper() in {"NA", "NAN", "."}:
                continue
            if tgt == src or tgt not in node_set:
                continue
            n1, n2 = sorted([src, tgt])
            if (n1, n2) in seen:
                continue
            seen.add((n1, n2))
            edges.append({
                "node1": n1, "node2": n2,
                "support_directions": "",
                "support_weight_mean": "",
                "support_weight_sum": "",
                "support_weight_max": "",
                "methods_union": "fallback=resolved_neighbors",
            })
    return edges


def write_tsv_if_changed(df: pd.DataFrame, output_path: str | Path, suffix: str = "") -> bool:
    """
    Write a TSV only when its content changes.

    The resolved-graph Snakemake rule declares provenance TSVs alongside the
    figure outputs. Preserving identical files keeps plot-only reruns from
    changing graph/centrality audit artifacts.
    """
    output_path = Path(output_path)
    content = df.to_csv(sep="\t", index=False).encode("utf-8")
    if output_path.exists():
        if output_path.read_bytes() == content:
            print(f"[OK] Unchanged: {output_path}{suffix}")
            return False
    output_path.write_bytes(content)
    print(f"[OK] Saved: {output_path}{suffix}")
    return True


def sidecar_output_path(output_prefix: str, suffix: str) -> Path:
    """Return provenance TSV path with a prefix matching the figure name."""
    output_prefix_path = Path(output_prefix)
    figure_stem = output_prefix_path.name.lower()
    if figure_stem.startswith("fig_multicohort_cancer_"):
        provenance_prefix = "multicohort_cancer_"
    elif "_pan_cancer_" in figure_stem:
        provenance_prefix = "pan_cancer_"
    else:
        provenance_prefix = "patient_referenced_"
    out_dir = output_prefix_path.parent
    if str(out_dir) == ".":
        return Path(f"{provenance_prefix}{suffix}")
    return out_dir / f"{provenance_prefix}{suffix}"


def legacy_sidecar_output_path(output_prefix: str, suffix: str) -> Path | None:
    """Return legacy patient-referenced TSV path for multicohort compatibility."""
    output_prefix_path = Path(output_prefix)
    if "_pan_cancer_" not in output_prefix_path.name:
        return None
    out_dir = output_prefix_path.parent
    if str(out_dir) == ".":
        return Path(f"patient_referenced_{suffix}")
    return out_dir / f"patient_referenced_{suffix}"


def write_legacy_sidecar_if_needed(
    df: pd.DataFrame,
    output_prefix: str,
    suffix: str,
    message_suffix: str = "",
) -> None:
    legacy_output = legacy_sidecar_output_path(output_prefix, suffix)
    if legacy_output is not None:
        write_tsv_if_changed(df, legacy_output, f"{message_suffix} [legacy compatibility]")


# =============================================================================
# Per-component centrality (layout is handled by graph_plot_style.grid_layout)
# =============================================================================

def compute_component_centrality(
    G_non_iso: nx.Graph,
    *,
    center_top_n: int,
):
    """
    Compute per-node betweenness centrality and degree-based anchor selection.

    Layout is not produced here — call gps.grid_layout separately.

    Returns (components_list_sorted, anchor_nodes,
             betweenness_unnormalised, betweenness_normalised,
             selected_by_unnormalised, selected_by_normalised,
             selected_by_degree).
    """
    components_list = sorted(
        nx.connected_components(G_non_iso),
        key=lambda c: (-len(c), sorted(c)[0]),
    )
    anchor_nodes: list[str] = []
    betweenness_unnormalised: dict[str, float] = {}
    betweenness_normalised: dict[str, float] = {}
    selected_by_unnormalised: dict[int, str | None] = {}
    selected_by_normalised: dict[int, str | None] = {}
    selected_by_degree: dict[int, str | None] = {}

    for i, comp in enumerate(components_list):
        comp_nodes = sorted(comp)
        sub = G_non_iso.subgraph(comp_nodes).copy()
        n = sub.number_of_nodes()

        if n == 1:
            betweenness_unnormalised[comp_nodes[0]] = 0.0
            betweenness_normalised[comp_nodes[0]] = 0.0
            selected_by_unnormalised[i] = None
            selected_by_normalised[i] = None
            selected_by_degree[i] = None
        else:
            bc_unnorm = nx.betweenness_centrality(sub, normalized=False, weight=None)
            bc_norm = nx.betweenness_centrality(sub, normalized=True, weight=None)
            betweenness_unnormalised.update(bc_unnorm)
            betweenness_normalised.update(bc_norm)

            selected_by_unnormalised[i] = max(
                comp_nodes, key=lambda node: (bc_unnorm.get(node, 0.0), node))
            selected_by_normalised[i] = max(
                comp_nodes, key=lambda node: (bc_norm.get(node, 0.0), node))
            selected_by_degree[i] = max(
                comp_nodes, key=lambda node: (sub.degree(node), node))

            if n >= 3:
                ranked = sorted(
                    bc_unnorm.items(), key=lambda x: (x[1], x[0]), reverse=True,
                )
                take = min(center_top_n, max(1, n // 12)) if n >= 12 else 1
                anchor_nodes.extend([node for node, _ in ranked[:take]])

    seen: set[str] = set()
    anchor_unique: list[str] = []
    for node in anchor_nodes:
        if node not in seen:
            anchor_unique.append(node)
            seen.add(node)

    return (
        components_list, anchor_unique,
        betweenness_unnormalised, betweenness_normalised,
        selected_by_unnormalised, selected_by_normalised, selected_by_degree,
    )


def component_band_mid_y_axes(
    ax,
    pos: dict[str, np.ndarray],
    components_list: list[set[str]],
    component_indices: tuple[int, ...] = (0, 1),
) -> tuple[float, str, float, float]:
    """Return the C-band midpoint in axes coordinates after axes limits exist."""
    y_values: list[float] = []
    labels: list[str] = []
    for idx in component_indices:
        if idx >= len(components_list):
            continue
        nodes = [n for n in components_list[idx] if n in pos]
        if not nodes:
            continue
        y_values.extend(float(pos[n][1]) for n in nodes)
        labels.append(f"C{idx + 1}")
    if not y_values:
        return 1.0, "fallback", float("nan"), float("nan")
    y_min = min(y_values)
    y_max = max(y_values)
    y_data = (y_min + y_max) / 2.0
    _, y_axes = ax.transAxes.inverted().transform(
        ax.transData.transform((0.0, y_data))
    )
    return float(y_axes), "/".join(labels), float(y_min), float(y_max)


def multicohort_adaptive_full_graph_layout(
    G_non_iso: nx.Graph,
    isolates: list[str],
    *,
    seed: int,
    isolate_spacing: float,
    isolate_panel_padding_x: float,
    isolate_panel_padding_y: float,
    isolate_label_band_frac: float,
    isolate_max_per_row: int,
    context: str,
) -> tuple[
    dict[str, np.ndarray],
    tuple[float, float, float, float] | None,
    list[dict[str, float | int | str | tuple[float, float, float, float]]],
    dict[str, tuple[float, float, float, float]],
]:
    """
    Lay out the multicohort full resolved graph as a compact thesis panel.

    Connected non-isolate components occupy the first two rows of a 3-column
    grid. Component cell size is scaled from node count plus edge count, with a
    smaller density term, so dense components get more visual radius without
    changing topology. Degree-0 isolates occupy a compact third-row panel.
    """
    components = [sorted(comp) for comp in gps.ordered_non_isolate_components(G_non_iso)]
    if not components:
        return {}, None, [], {}

    n_cols = 3
    n_component_rows = int(math.ceil(len(components) / float(n_cols)))
    max_nodes = max(len(nodes) for nodes in components)
    max_edges = max(
        G_non_iso.subgraph(nodes).number_of_edges()
        for nodes in components
    )
    max_edges = max(max_edges, 1)

    specs: list[dict[str, float | int | str | list[str]]] = []
    readability_expansion = {
        1: (1.75, 1.58),
        2: (1.85, 1.70),
        5: (1.45, 1.38),
    }
    for idx, nodes in enumerate(components, start=1):
        sub = G_non_iso.subgraph(nodes)
        n_nodes = len(nodes)
        n_edges = sub.number_of_edges()
        density = gps.component_edge_density(G_non_iso, nodes)
        node_score = math.sqrt(float(n_nodes) / float(max_nodes))
        edge_score = math.sqrt(float(n_edges) / float(max_edges))
        density_score = min(float(density) / 0.75, 1.0)
        complexity = (
            0.35 * node_score
            + 0.50 * edge_score
            + 0.15 * density_score
        )
        if n_nodes <= 2:
            complexity = min(complexity, 0.34)

        label_pad = min(
            max((len(str(node)) for node in nodes), default=0) * 0.16,
            5.5,
        )
        target_width = 20.0 + 54.0 * complexity
        target_height = 14.0 + 35.0 * complexity
        if n_nodes <= 2:
            target_width = min(target_width, 35.0)
            target_height = min(target_height, 22.0)
        expansion_x, expansion_y = readability_expansion.get(idx, (1.0, 1.0))
        target_width *= expansion_x
        target_height *= expansion_y
        cell_width = target_width + 15.0 + label_pad
        cell_height = target_height + 12.0
        row = (idx - 1) // n_cols
        col = (idx - 1) % n_cols
        specs.append({
            "idx": idx,
            "component": f"C{idx}",
            "nodes": nodes,
            "row": row,
            "col": col,
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            "density": density,
            "complexity": complexity,
            "readability_expand_x": expansion_x,
            "readability_expand_y": expansion_y,
            "target_width": target_width,
            "target_height": target_height,
            "cell_width": cell_width,
            "cell_height": cell_height,
        })

    col_gap = 8.0
    row_gap = 7.0
    col_widths = [
        max(
            float(spec["cell_width"])
            for spec in specs
            if int(spec["col"]) == col
        )
        for col in range(n_cols)
    ]
    row_heights = [
        max(
            float(spec["cell_height"])
            for spec in specs
            if int(spec["row"]) == row
        )
        for row in range(n_component_rows)
    ]
    grid_width = sum(col_widths) + col_gap * float(n_cols - 1)
    grid_height = sum(row_heights) + row_gap * float(max(0, n_component_rows - 1))
    grid_x0 = -grid_width / 2.0
    grid_x1 = grid_width / 2.0
    grid_y1 = 0.0
    grid_y0 = -grid_height

    col_lefts: list[float] = []
    cursor = grid_x0
    for width in col_widths:
        col_lefts.append(cursor)
        cursor += width + col_gap

    row_tops: list[float] = []
    cursor = grid_y1
    for height in row_heights:
        row_tops.append(cursor)
        cursor -= height + row_gap

    pos_non_iso: dict[str, np.ndarray] = {}
    placements: list[dict[str, float | int | str | tuple[float, float, float, float]]] = []
    for spec in specs:
        nodes = list(spec["nodes"])  # type: ignore[arg-type]
        row = int(spec["row"])
        col = int(spec["col"])
        cell_x0 = col_lefts[col]
        cell_x1 = cell_x0 + col_widths[col]
        cell_y1 = row_tops[row]
        cell_y0 = cell_y1 - row_heights[row]
        cx = (cell_x0 + cell_x1) / 2.0
        cy = (cell_y0 + cell_y1) / 2.0

        local_pos = gps._layout_component_for_cell(
            G_non_iso,
            nodes,
            seed=seed,
            target_width=float(spec["target_width"]),
            target_height=float(spec["target_height"]),
            spring_iters=450,
        )
        for node, p in local_pos.items():
            pos_non_iso[node] = np.array(
                [cx + float(p[0]), cy + float(p[1])],
                dtype=float,
            )

        pts = np.array([pos_non_iso[n] for n in nodes], dtype=float)
        bbox = (
            float(pts[:, 0].min()), float(pts[:, 1].min()),
            float(pts[:, 0].max()), float(pts[:, 1].max()),
        )
        cell_bbox = (float(cell_x0), float(cell_y0), float(cell_x1), float(cell_y1))
        placement = {
            "component": str(spec["component"]),
            "component_index": int(spec["idx"]),
            "row": row + 1,
            "col": col + 1,
            "n": int(spec["n_nodes"]),
            "edges": int(spec["n_edges"]),
            "density": float(spec["density"]),
            "complexity": float(spec["complexity"]),
            "readability_expand_x": float(spec["readability_expand_x"]),
            "readability_expand_y": float(spec["readability_expand_y"]),
            "target_width": float(spec["target_width"]),
            "target_height": float(spec["target_height"]),
            "cell_bbox": cell_bbox,
            "bbox": bbox,
        }
        placements.append(placement)
        print(
            f"[QC] {context} adaptive 3x3 placement: "
            f"{placement['component']} row={row + 1} col={col + 1} "
            f"n={placement['n']} edges={placement['edges']} "
            f"density={placement['density']:.3f} "
            f"complexity={placement['complexity']:.3f} "
            f"readability_expand=({placement['readability_expand_x']:.2f},"
            f"{placement['readability_expand_y']:.2f}) "
            f"target=({placement['target_width']:.2f},"
            f"{placement['target_height']:.2f}) "
            f"cell={gps._fmt_bbox(cell_bbox)} bbox={gps._fmt_bbox(bbox)}"
        )

    bottom_gap = 6.0
    bottom_row_height = 18.5
    bottom_y1 = grid_y0 - bottom_gap
    bottom_y0 = bottom_y1 - bottom_row_height
    left_two_width = col_widths[0] + col_gap + col_widths[1]
    left_two_x0 = grid_x0
    left_two_x1 = left_two_x0 + left_two_width
    legend_x0 = col_lefts[2]
    legend_x1 = grid_x1

    iso_bbox = None
    pos_iso: dict[str, np.ndarray] = {}
    if isolates:
        iso_sorted = sorted(isolates)
        n_iso = len(iso_sorted)
        n_iso_rows = int(math.ceil(n_iso / float(max(1, isolate_max_per_row))))
        spacing = min(max(float(isolate_spacing), 14.0), 20.0)
        panel_width = max(
            82.0,
            min(n_iso, max(1, isolate_max_per_row)) * spacing
            + 2.0 * max(float(isolate_panel_padding_x), 8.0),
        )
        panel_width = min(panel_width, left_two_width)
        panel_cx = (left_two_x0 + left_two_x1) / 2.0
        iso_bbox = (
            panel_cx - panel_width / 2.0,
            bottom_y0,
            panel_cx + panel_width / 2.0,
            bottom_y1,
        )
        label_band_frac = float(np.clip(isolate_label_band_frac, 0.26, 0.38))
        label_band_height = (iso_bbox[3] - iso_bbox[1]) * label_band_frac
        node_y0 = iso_bbox[1] + label_band_height + float(isolate_panel_padding_y)
        node_y1 = iso_bbox[3] - float(isolate_panel_padding_y)
        row_step = max((node_y1 - node_y0) / float(n_iso_rows + 1), 1.0)
        for row in range(n_iso_rows):
            row_start = row * max(1, isolate_max_per_row)
            row_nodes = iso_sorted[row_start:row_start + max(1, isolate_max_per_row)]
            if len(row_nodes) == 1:
                xs = np.array([panel_cx], dtype=float)
            else:
                xs = np.linspace(
                    iso_bbox[0] + float(isolate_panel_padding_x),
                    iso_bbox[2] - float(isolate_panel_padding_x),
                    len(row_nodes),
                )
            y_iso = node_y1 - float(row + 1) * row_step
            for node, x in zip(row_nodes, xs):
                pos_iso[node] = np.array([float(x), float(y_iso)], dtype=float)

    legend_region = (
        float(legend_x0),
        float(bottom_y0),
        float(legend_x1),
        float(bottom_y1),
    )
    regions = {
        "title_region": (
            float(grid_x0), 1.0, float(grid_x1), 7.0,
        ),
        "component_grid": (
            float(grid_x0), float(grid_y0), float(grid_x1), float(grid_y1),
        ),
        "legend_region": legend_region,
        "footnote_region": (
            float(grid_x0), float(bottom_y0 - 6.0),
            float(grid_x1), float(bottom_y0 - 1.0),
        ),
    }
    if iso_bbox is not None:
        regions["isolate_box"] = iso_bbox
        regions["isolate_region"] = iso_bbox
    print(
        f"[QC] {context} adaptive 3x3 regions: "
        f"grid={gps._fmt_bbox(regions['component_grid'])} "
        f"isolate={gps._fmt_bbox(iso_bbox) if iso_bbox else 'none'} "
        f"legend={gps._fmt_bbox(legend_region)} "
        f"bottom_row=({bottom_y0:.2f},{bottom_y1:.2f}) "
        f"column_widths={[round(v, 2) for v in col_widths]} "
        f"row_heights={[round(v, 2) for v in row_heights]}"
    )
    gps.report_planned_layout_clearance(context, placements, regions)
    return {**pos_non_iso, **pos_iso}, iso_bbox, placements, regions


def compact_resolved_legend_handles(handles: list) -> list:
    """Shorten verbose resolved-graph legend text for compact full-graph panels."""
    replacements = {
        "Cell line; node fill colour denotes connected component only":
            "Cell line (fill = component)",
        "Bridge-like anchor, within component": "Bridge-like anchor",
        "Most-connected anchor, within component": "Most-connected anchor",
        "Both centrality annotations": "Both anchor markers",
        "Isolate (degree 0)": "Isolate",
        "Edge width denotes support count": "Edge width = support",
        "Single-supported edge (dashed)": "Single-supported edge",
    }
    compact = []
    for item in handles:
        if isinstance(item, tuple) and len(item) == 3:
            label = replacements.get(str(item[2]), str(item[2]))
            compact.append((item[0], item[1], label))
            continue
        label = item.get_label()
        if label.startswith("Support count = "):
            label = label.replace("Support count = ", "Support = ")
        else:
            label = replacements.get(label, label)
        item.set_label(label)
        compact.append(item)
    return compact


def remove_edge_support_legend_handles(handles: list) -> list:
    """Remove support-count/dashed-edge legend entries from resolved graph panels."""
    filtered = []
    drop_labels = {
        "Edge width denotes support count",
        "Edge width = support",
        "Single-supported edge (dashed)",
        "Single-supported edge",
    }
    for item in handles:
        label = item[2] if isinstance(item, tuple) and len(item) == 3 else item.get_label()
        label = str(label)
        if label.startswith("Support count = ") or label.startswith("Support = "):
            continue
        if label in drop_labels:
            continue
        filtered.append(item)
    return filtered


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Plot resolved cell-line neighbourhood graph (strict + "
                    "packed components + §5 visual language).",
    )
    ap.add_argument("resolved_tsv", help="Path to resolved_dsmz_neighbours.tsv")
    ap.add_argument("output_prefix", help="Output prefix (no extension)")
    ap.add_argument("label", help="Plot label (e.g. BRCA, NBL, RBL)")

    ap.add_argument("--center-top-n", type=int, default=3)
    ap.add_argument("--fig-w", type=float, default=22.0)
    ap.add_argument("--fig-h", type=float, default=16.0)
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--seed", type=int, default=gps.LAYOUT_SEED)
    ap.add_argument("--edge-width-min", type=float, default=1.2)
    ap.add_argument("--edge-width-max", type=float, default=5.0)
    ap.add_argument("--font-size", type=int, default=11)
    ap.add_argument("--node-size", type=float, default=gps.NODE_SIZE)
    ap.add_argument("--legend-mode",
                    choices=("keep", "none", "separate"),
                    default="keep")
    ap.add_argument("--legend-out-prefix", default=None)
    ap.add_argument("--use-halo-labels", action="store_true",
                    help="Draw node/component labels with halo strokes and no label boxes.")
    ap.add_argument("--halo-linewidth", type=float, default=3.6,
                    help="White halo stroke width for node labels.")
    ap.add_argument("--halo-color", default="white",
                    help="Halo stroke colour for node and component labels.")
    ap.add_argument("--dense-component-min-size", type=int, default=8)
    ap.add_argument("--dense-component-min-density", type=float, default=0.30)
    ap.add_argument("--dense-component-expand-factor", type=float, default=2.00)
    ap.add_argument("--adaptive-component-layout", action="store_true",
                    help="Expand each non-isolate component after layout using structural burden metrics.")
    ap.add_argument("--adaptive-component-min-nodes", type=int, default=4)
    ap.add_argument("--adaptive-component-min-density", type=float, default=0.15)
    ap.add_argument("--adaptive-component-min-expand", type=float, default=1.0)
    ap.add_argument("--adaptive-component-max-expand", type=float, default=4.0)
    ap.add_argument("--adaptive-component-density-weight", type=float, default=1.8)
    ap.add_argument("--adaptive-component-size-weight", type=float, default=0.14)
    ap.add_argument("--adaptive-component-label-weight", type=float, default=0.30)
    ap.add_argument("--position-fill-x", type=float, default=0.0)
    ap.add_argument("--position-fill-y", type=float, default=0.0)
    ap.add_argument("--axis-margin-frac", type=float, default=0.06)
    ap.add_argument("--label-overlap-avoidance", action="store_true")
    ap.add_argument("--label-padding-factor", type=float, default=1.35)
    ap.add_argument("--label-bbox-margin", type=float, default=0.0)
    ap.add_argument("--max-label-overlap-iterations", type=int, default=1200)
    ap.add_argument("--min-label-separation", type=float, default=0.0)
    ap.add_argument("--avoid-isolate-overlap", action="store_true")
    ap.add_argument("--min-component-isolate-gap", type=float, default=4.0)
    ap.add_argument("--isolate-panel-reserved-height", type=float, default=0.0)
    ap.add_argument("--adaptive-isolate-y-shift", action="store_true")
    ap.add_argument("--component-isolate-gap-factor", type=float, default=1.0)
    ap.add_argument("--component-label-gap", type=float, default=2.0)
    ap.add_argument("--component-grid-ncols", type=int, default=2)
    ap.add_argument("--component-grid-nrows", type=int, default=2)
    ap.add_argument("--component-cell-width", type=float, default=66.0)
    ap.add_argument("--component-cell-height", type=float, default=68.0)
    ap.add_argument("--component-row-gap", type=float, default=14.0)
    ap.add_argument("--component-col-gap", type=float, default=18.0)
    ap.add_argument("--legend-width", type=float, default=62.0)
    ap.add_argument("--isolate-region-height", type=float, default=10.5)
    ap.add_argument("--footnote-region-height", type=float, default=7.0)
    ap.add_argument("--component-cell-fill-x", type=float, default=0.78)
    ap.add_argument("--component-cell-fill-y", type=float, default=0.72)
    ap.add_argument("--dense-component-cell-fill-x", type=float, default=0.90)
    ap.add_argument("--dense-component-cell-fill-y", type=float, default=0.82)
    ap.add_argument("--c2-component-cell-fill-x", type=float, default=0.84)
    ap.add_argument("--c2-component-cell-fill-y", type=float, default=0.78)
    ap.add_argument("--c1-x-expand-factor", type=float, default=3.00)
    ap.add_argument("--c1-y-expand-factor", type=float, default=3.00)
    ap.add_argument("--c2-x-expand-factor", type=float, default=1.35)
    ap.add_argument("--c2-y-expand-factor", type=float, default=1.05)
    ap.add_argument("--isolate-gap", type=float, default=gps.ADAPTIVE_ISOLATE_GAP)
    ap.add_argument("--isolate-spacing", type=float, default=gps.ADAPTIVE_ISOLATE_SPACING)
    ap.add_argument("--isolate-box-height-factor", type=float, default=1.40)
    ap.add_argument("--isolate-label-y-frac", type=float, default=0.010)
    ap.add_argument("--isolate-label-left-pad", type=float, default=0.0)
    ap.add_argument("--compact-isolate-panel", action="store_true")
    ap.add_argument("--isolate-panel-width-mode",
                    choices=("layout", "content"), default="layout")
    ap.add_argument("--isolate-panel-spacing", type=float, default=None)
    ap.add_argument("--min-isolate-label-spacing", type=float, default=0.0)
    ap.add_argument("--isolate-panel-padding-x", type=float, default=12.0)
    ap.add_argument("--isolate-panel-padding-y", type=float, default=1.4)
    ap.add_argument("--isolate-label-band-frac", type=float, default=0.34)
    ap.add_argument("--isolate-max-per-row", type=int, default=6)
    ap.add_argument("--isolate-panel-width-frac", type=float, default=0.62)
    ap.add_argument("--legend-alignment-mode",
                    choices=("default", "component-band"),
                    default="component-band")
    args = ap.parse_args()

    label_upper = args.label.upper()
    is_nbl  = label_upper == "NBL"
    is_brca = label_upper == "BRCA"
    is_rbl  = label_upper == "RBL"
    is_multicohort = label_upper == "MULTICOHORT_CANCER"
    if is_multicohort:
        args.fig_w = 21.0
        args.fig_h = 14.2
        args.font_size = max(args.font_size, 12)
        print(
            "[QC] Multicohort resolved full-graph thesis framing: "
            f"fig_w={args.fig_w:.1f} fig_h={args.fig_h:.1f} "
            f"font_size={args.font_size}"
        )

    # ------------------------------------------------------------------ Load
    try:
        df = pd.read_csv(args.resolved_tsv, sep="\t", dtype=str)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to read TSV: {args.resolved_tsv}\n{e}\n")
        sys.exit(1)
    if "cell_line" not in df.columns:
        sys.stderr.write("[ERROR] TSV must contain column: cell_line\n")
        sys.exit(1)
    if "final_neighbors" in df.columns:
        nb_col = "final_neighbors"
    elif "final_neighbours" in df.columns:
        nb_col = "final_neighbours"
    else:
        sys.stderr.write("[ERROR] TSV must contain column: final_neighbors or final_neighbours\n")
        sys.exit(1)

    df["cell_line"] = df["cell_line"].astype(str).str.strip()
    df[nb_col] = df[nb_col].fillna("").astype(str)
    df["cell_line"] = df["cell_line"].map(canonical_cell_line_id)
    df[nb_col] = df[nb_col].map(canonical_neighbour_string)
    df = df[df["cell_line"] != ""].copy()
    if df["cell_line"].duplicated().any():
        agg_cols = {
            nb_col: lambda values: canonical_neighbour_string(
                ";".join([str(v) for v in values if str(v).strip()])
            )
        }
        for col in ("reps_used", "best_overall_dir", "winner_dir"):
            if col in df.columns:
                agg_cols[col] = lambda values: ";".join(
                    sorted({str(v).strip() for v in values if str(v).strip()})
                )
        df = df.groupby("cell_line", as_index=False).agg(agg_cols)
        df["n_final"] = df[nb_col].map(lambda x: len(parse_neighbours(x)))

    node_set = set(df["cell_line"].unique())
    short_map = build_shortname_map(node_set)
    if len(set(short_map.values())) != len(short_map.values()):
        raise ValueError("[ERROR] Non-unique short names detected.")

    # ------------------------------------------------------------ Build graph
    G = nx.Graph()
    G.add_nodes_from(sorted(node_set))
    self_loops = 0
    ignored_neighbours: set[str] = set()
    for _, row in df.iterrows():
        src = row["cell_line"]
        for tgt in parse_neighbours(row[nb_col]):
            if tgt == src:
                self_loops += 1
                continue
            if tgt not in node_set:
                ignored_neighbours.add(tgt)
                continue
            G.add_edge(src, tgt)

    isolates = sorted(n for n in G.nodes() if G.degree(n) == 0)
    components_all = list(nx.connected_components(G))
    print("=" * 60)
    print(f"Nodes={G.number_of_nodes()} Edges={G.number_of_edges()} "
          f"Components={len(components_all)} Isolates={len(isolates)}")
    if self_loops:
        print(f"Self-loop entries ignored: {self_loops}")
    if ignored_neighbours:
        print(f"Neighbours ignored (not in cell_line): {len(ignored_neighbours)}")

    # ------------------------------------------------------------- Shortnames TSV
    mapping_out = sidecar_output_path(args.output_prefix, "cell_line_display_names.tsv")
    mapping_df = pd.DataFrame(
        [{"long_id": k, "short_id": v} for k, v in sorted(short_map.items())]
    )
    write_tsv_if_changed(mapping_df, mapping_out)
    write_legacy_sidecar_if_needed(
        mapping_df, args.output_prefix, "cell_line_display_names.tsv"
    )

    # -------------------------------------------------------------------- Edges TSV
    edges_output = sidecar_output_path(args.output_prefix, "resolved_cell_line_neighbourhood_graph_edges.tsv")
    output_dir = Path(edges_output).parent
    # Per-direction edge files live at
    #   tumour_neighbourhoods/<direction>/final_consensus/cell_line_similarity_graph_edges_<direction>.tsv
    # output_dir is the …/final_consensus_all/plots/ directory, so we go up
    # two levels to reach tumour_neighbourhoods/ (sibling of final_consensus_all).
    tumour_nh_root = output_dir.parent.parent

    edge_records = []
    if "best_overall_dir" in df.columns and "winner_dir" in df.columns:
        directions_used: set[str] = set()
        for _, row in df.iterrows():
            if pd.notna(row.get("best_overall_dir")):
                directions_used.add(str(row["best_overall_dir"]))
            if pd.notna(row.get("winner_dir")):
                directions_used.add(str(row["winner_dir"]))
        if tumour_nh_root.exists():
            for direction_dir in tumour_nh_root.iterdir():
                if direction_dir.is_dir() and direction_dir.name != "final_consensus_all":
                    direction = direction_dir.name
                    edge_file = (direction_dir / "final_consensus"
                                 / f"cell_line_similarity_graph_edges_{direction}.tsv")
                    if edge_file.exists():
                        directions_used.add(direction)
        for direction in sorted(directions_used):
            edge_file = (tumour_nh_root / direction / "final_consensus"
                         / f"cell_line_similarity_graph_edges_{direction}.tsv")
            if not edge_file.exists():
                continue
            try:
                dir_edges = pd.read_csv(edge_file, sep="\t")
                if "cell_line1" in dir_edges.columns and "cell_line2" in dir_edges.columns:
                    for _, edge_row in dir_edges.iterrows():
                        n1 = str(edge_row["cell_line1"]).strip()
                        n2 = str(edge_row["cell_line2"]).strip()
                        n1 = canonical_cell_line_id(n1)
                        n2 = canonical_cell_line_id(n2)
                        sim = (float(edge_row.get("similarity", 1.0))
                               if "similarity" in dir_edges.columns else 1.0)
                        if (n1 != n2 and n1 in node_set and n2 in node_set
                                and G.has_edge(n1, n2)):
                            if n1 > n2:
                                n1, n2 = n2, n1
                            edge_records.append({
                                "node1": n1, "node2": n2,
                                "direction": direction, "similarity": sim,
                            })
            except Exception as e:
                print(f"[WARNING] Failed to read {edge_file}: {e}", file=sys.stderr)

    edges_source = "direction_files"
    n_edges_written = 0
    edge_support: dict[tuple[str, str], int] = {}
    if edge_records:
        edges_df = pd.DataFrame(edge_records)
        methods_map: dict[str, str] = {}
        for _, row in df.iterrows():
            cell_line = str(row["cell_line"]).strip()
            methods_str = ""
            if pd.notna(row.get("reps_used")):
                methods_str = str(row["reps_used"])
            elif (pd.notna(row.get("best_overall_dir")) and
                  pd.notna(row.get("winner_dir"))):
                methods_str = (f"best_overall={row['best_overall_dir']};"
                               f"winner={row['winner_dir']}")
            methods_map[cell_line] = methods_str
        aggregated = []
        for (n1, n2), group in edges_df.groupby(["node1", "node2"]):
            directions_list = sorted(group["direction"].unique())
            sims = group["similarity"].tolist()
            methods_union = methods_map.get(n1, methods_map.get(n2, ""))
            aggregated.append({
                "node1": n1, "node2": n2,
                "support_directions": len(directions_list),
                "support_weight_mean": np.mean(sims) if sims else 1.0,
                "support_weight_sum": sum(sims),
                "support_weight_max": max(sims) if sims else 1.0,
                "methods_union": methods_union,
            })
            edge_support[(n1, n2)] = len(directions_list)
        rich_edges_df = pd.DataFrame(aggregated)
        rich_edges_df["node1_short"] = rich_edges_df["node1"].map(lambda x: short_map.get(x, x))
        rich_edges_df["node2_short"] = rich_edges_df["node2"].map(lambda x: short_map.get(x, x))
        rich_edges_df = rich_edges_df.sort_values(["node1_short", "node2_short", "node1", "node2"])
        n_edges_written = len(rich_edges_df)
        write_tsv_if_changed(rich_edges_df, edges_output, f" ({len(rich_edges_df)} edges)")
        write_legacy_sidecar_if_needed(
            rich_edges_df,
            args.output_prefix,
            "resolved_cell_line_neighbourhood_graph_edges.tsv",
            f" ({len(rich_edges_df)} edges)",
        )
    else:
        print("[WARN] No direction-specific edge files found; falling back to resolved neighbours.")
        edges_source = "fallback_resolved_neighbors"
        fallback_edges = edges_from_resolved_neighbours(df, nb_col=nb_col, node_set=node_set)
        if fallback_edges:
            fb = pd.DataFrame(fallback_edges)
            fb["node1_short"] = fb["node1"].map(lambda x: short_map.get(x, x))
            fb["node2_short"] = fb["node2"].map(lambda x: short_map.get(x, x))
            n_edges_written = len(fallback_edges)
            write_tsv_if_changed(fb, edges_output, f" ({len(fallback_edges)} edges)")
            write_legacy_sidecar_if_needed(
                fb,
                args.output_prefix,
                "resolved_cell_line_neighbourhood_graph_edges.tsv",
                f" ({len(fallback_edges)} edges)",
            )
        else:
            empty_edges_df = pd.DataFrame(columns=[
                "node1", "node2", "support_directions",
                "support_weight_mean", "support_weight_sum",
                "support_weight_max", "methods_union",
            ])
            write_tsv_if_changed(empty_edges_df, edges_output, " (empty - no edges)")
            write_legacy_sidecar_if_needed(
                empty_edges_df,
                args.output_prefix,
                "resolved_cell_line_neighbourhood_graph_edges.tsv",
                " (empty - no edges)",
            )

    print(f"[QC] Edge export: nodes={G.number_of_nodes()} "
          f"edges_written={n_edges_written} source={edges_source}")

    # ----------------------------------------------------- Layout + centrality
    G_non_iso = G.subgraph([n for n in G.nodes() if n not in set(isolates)]).copy()
    (components_list, anchor_nodes,
     betweenness_unnormalised, betweenness_normalised,
     selected_by_unnormalised, selected_by_normalised, selected_by_degree
     ) = compute_component_centrality(
        G_non_iso, center_top_n=args.center_top_n,
    )
    if is_rbl:
        pos, iso_bbox = gps.rbl_one_component_layout(
            G_non_iso, isolates, seed=args.seed, isolate_gap=args.isolate_gap,
            isolate_spacing=args.isolate_spacing,
            isolate_box_height_factor=args.isolate_box_height_factor,
        )
        planner_regions = {}
        if args.adaptive_component_layout:
            gps.expand_adaptive_components_after_layout(
                G,
                pos,
                isolates,
                label_map=short_map,
                min_nodes=args.adaptive_component_min_nodes,
                min_density=args.adaptive_component_min_density,
                min_expand=args.adaptive_component_min_expand,
                max_expand=args.adaptive_component_max_expand,
                density_weight=args.adaptive_component_density_weight,
                size_weight=args.adaptive_component_size_weight,
                label_weight=args.adaptive_component_label_weight,
                context=f"{args.label} resolved neighbours",
            )
    else:
        if is_brca:
            pos, iso_bbox, _grid_placements, planner_regions = gps.planned_component_grid_layout(
                G_non_iso, isolates, seed=args.seed,
                n_cols=args.component_grid_ncols,
                n_rows=args.component_grid_nrows,
                cell_width=args.component_cell_width,
                cell_height=args.component_cell_height,
                row_gap=args.component_row_gap,
                col_gap=args.component_col_gap,
                legend_width=args.legend_width,
                isolate_gap=args.isolate_gap,
                isolate_spacing=args.isolate_spacing,
                isolate_region_height=args.isolate_region_height,
                isolate_box_height_factor=args.isolate_box_height_factor,
                isolate_label_left_pad=args.isolate_label_left_pad,
                isolate_panel_padding_x=args.isolate_panel_padding_x,
                isolate_panel_padding_y=args.isolate_panel_padding_y,
                isolate_label_band_frac=args.isolate_label_band_frac,
                isolate_max_per_row=args.isolate_max_per_row,
                isolate_panel_width_frac=args.isolate_panel_width_frac,
                footnote_region_height=args.footnote_region_height,
                component_cell_fill_x=args.component_cell_fill_x,
                component_cell_fill_y=args.component_cell_fill_y,
                dense_component_cell_fill_x=args.dense_component_cell_fill_x,
                dense_component_cell_fill_y=args.dense_component_cell_fill_y,
                c2_component_cell_fill_x=args.c2_component_cell_fill_x,
                c2_component_cell_fill_y=args.c2_component_cell_fill_y,
                context=f"{args.label} resolved neighbours",
            )
            print(
                f"[QC] {args.label} resolved neighbours deterministic planner: "
                f"row_major_grid={args.component_grid_ncols}x{args.component_grid_nrows} "
                f"cell_width={args.component_cell_width:.2f} "
                f"cell_height={args.component_cell_height:.2f} "
                f"row_gap={args.component_row_gap:.2f} "
                f"col_gap={args.component_col_gap:.2f} "
                f"legend_width={args.legend_width:.2f} "
                f"isolate_region_height={args.isolate_region_height:.2f} "
                f"footnote_region_height={args.footnote_region_height:.2f}"
            )
            for comp_label, bbox in gps.component_layout_bboxes(G, pos, isolates).items():
                print(
                    f"[QC] {args.label} resolved neighbours "
                    f"{comp_label} final bbox={gps._fmt_bbox(bbox)}"
                )
        elif is_multicohort:
            pos, iso_bbox, _grid_placements, planner_regions = (
                multicohort_adaptive_full_graph_layout(
                    G_non_iso,
                    isolates,
                    seed=args.seed,
                    isolate_spacing=args.isolate_spacing,
                    isolate_panel_padding_x=args.isolate_panel_padding_x,
                    isolate_panel_padding_y=args.isolate_panel_padding_y,
                    isolate_label_band_frac=args.isolate_label_band_frac,
                    isolate_max_per_row=args.isolate_max_per_row,
                    context=f"{args.label} resolved neighbours",
                )
            )
            print(
                f"[QC] {args.label} resolved neighbours adaptive full graph: "
                "layout=3 columns x 3 rows; bottom row contains compact isolates "
                "and legend; component visual space scales with node and edge count"
            )
            for comp_label, bbox in gps.component_layout_bboxes(G, pos, isolates).items():
                print(
                    f"[QC] {args.label} resolved neighbours "
                    f"{comp_label} final bbox={gps._fmt_bbox(bbox)}"
                )
        else:
            pos, iso_bbox = gps.grid_layout(
                G_non_iso, isolates, seed=args.seed,
                isolate_gap=args.isolate_gap,
                isolate_spacing=args.isolate_spacing,
                isolate_box_height_factor=args.isolate_box_height_factor,
            )
            planner_regions = {}
        if args.adaptive_component_layout:
            gps.expand_adaptive_components_after_layout(
                G,
                pos,
                isolates,
                label_map=short_map,
                min_nodes=args.adaptive_component_min_nodes,
                min_density=args.adaptive_component_min_density,
                min_expand=args.adaptive_component_min_expand,
                max_expand=args.adaptive_component_max_expand,
                density_weight=args.adaptive_component_density_weight,
                size_weight=args.adaptive_component_size_weight,
                label_weight=args.adaptive_component_label_weight,
                context=f"{args.label} resolved neighbours",
            )
        elif not (is_brca or is_multicohort):
            gps.expand_dense_components_after_layout(
                G, pos, isolates,
                min_size=args.dense_component_min_size,
                min_density=args.dense_component_min_density,
                expand_factor=args.dense_component_expand_factor,
                context=f"{args.label} resolved neighbours",
            )
        if (
            is_brca and iso_bbox is not None
            and args.isolate_label_left_pad > 0
            and not planner_regions.get("isolate_label_band")
        ):
            iso_bbox = gps.pad_isolate_bbox_for_label(
                iso_bbox,
                left_pad=args.isolate_label_left_pad,
            )
            print(
                f"[QC] {args.label} resolved neighbours isolate label left pad: "
                f"{args.isolate_label_left_pad:.2f}"
            )
        gps.separate_focus_component_overlaps(
            G, pos, isolates,
            context=f"{args.label} resolved neighbours",
            focus_component_index=1,
            min_gap=3.0,
        )
        gps.report_component_clearance(
            G, pos, isolates,
            context=f"{args.label} resolved neighbours",
            focus_component_index=1,
            extra_bboxes={"isolate_box": iso_bbox},
        )
        gps.report_component_clearance(
            G, pos, isolates,
            context=f"{args.label} resolved neighbours",
            focus_component_index=2,
            extra_bboxes={"isolate_box": iso_bbox},
        )
    if iso_bbox is not None and args.compact_isolate_panel:
        isolate_panel_spacing = (
            args.isolate_panel_spacing
            if args.isolate_panel_spacing is not None
            else args.isolate_spacing
        )
        non_isolate_bbox = gps.position_bbox(
            pos,
            nodes=[n for n in G.nodes() if n not in set(isolates)],
        )
        iso_bbox = gps.compact_isolate_panel_layout(
            pos,
            isolates,
            iso_bbox,
            spacing=isolate_panel_spacing,
            padding_x=args.isolate_panel_padding_x,
            width_mode=args.isolate_panel_width_mode,
            align_bbox=non_isolate_bbox,
            label_map=short_map,
            min_label_spacing=args.min_isolate_label_spacing,
            label_padding_factor=args.label_padding_factor,
            context=f"{args.label} resolved neighbours",
        )
        print(
            f"[QC] {args.label} resolved neighbours compact isolate panel: "
            f"enabled={args.compact_isolate_panel} mode={args.isolate_panel_width_mode} "
            f"spacing={float(isolate_panel_spacing):.2f} "
            f"padding_x={float(args.isolate_panel_padding_x):.2f} "
            f"bbox={gps._fmt_bbox(iso_bbox) if iso_bbox is not None else 'none'}"
        )
    if args.avoid_isolate_overlap:
        iso_bbox, _gap_report = gps.enforce_component_isolate_gap(
            G,
            pos,
            isolates,
            iso_bbox,
            min_gap=args.min_component_isolate_gap,
            component_isolate_gap_factor=args.component_isolate_gap_factor,
            bbox_margin=args.label_bbox_margin,
            label_padding_factor=args.label_padding_factor,
            labels_by_node=short_map,
            shift_isolate_panel=args.adaptive_isolate_y_shift,
            context=f"{args.label} resolved neighbours",
        )
    print(
        f"[QC] {args.label} resolved neighbours position fill requested: "
        f"fill_x={float(args.position_fill_x):.3f} fill_y={float(args.position_fill_y):.3f} "
        f"collision_aware={'yes' if args.avoid_isolate_overlap else 'no'}"
    )
    fill_transform = gps.rescale_positions_to_target_fill(
        pos,
        fill_x=args.position_fill_x,
        fill_y=args.position_fill_y,
        figure_aspect=args.fig_w / max(args.fig_h, 1e-9),
    )
    if fill_transform is not None:
        iso_bbox = gps.apply_position_transform_to_bbox(iso_bbox, fill_transform)
        planner_regions = {
            key: (
                gps.apply_position_transform_to_bbox(tuple(value), fill_transform)
                if value is not None and len(value) == 4 else value
            )
            for key, value in planner_regions.items()
        }
        print(
            f"[QC] {args.label} resolved neighbours position fill rescale: "
            f"scale_x={fill_transform['scale_x']:.3f} "
            f"scale_y={fill_transform['scale_y']:.3f} "
            f"current_aspect={fill_transform['current_aspect']:.3f} "
            f"target_aspect={fill_transform['target_aspect']:.3f}"
        )
    else:
        print(
            f"[QC] {args.label} resolved neighbours position fill rescale: "
            "scale_x=1.000 scale_y=1.000 effective_fill_unchanged=yes"
        )
    if iso_bbox is not None:
        iso_height = float(iso_bbox[3] - iso_bbox[1])
        print(
            f"[QC] {args.label} resolved neighbours isolate box height: "
            f"old_unscaled_height=4.40 factor={args.isolate_box_height_factor:.2f} "
            f"height={iso_height:.2f} label_y_frac={args.isolate_label_y_frac:.3f}"
        )
        gps.report_isolate_layout_metrics(pos, isolates, iso_bbox, f"{args.label} resolved neighbours")

    for n in G.nodes():
        if n not in pos:
            betweenness_unnormalised.setdefault(n, 0.0)
            betweenness_normalised.setdefault(n, 0.0)

    component_size_by_node = {
        node: len(comp)
        for comp in components_list
        for node in comp
    }
    component_id_by_node = {
        node: idx
        for idx, comp in enumerate(components_list)
        for node in comp
    }
    bridge_anchor_nodes_all = sorted(
        {node for node in selected_by_unnormalised.values() if node is not None}
    )
    most_connected_nodes_all = sorted(
        {node for node in selected_by_degree.values() if node is not None}
    )
    bridge_anchor_nodes = bridge_anchor_nodes_all
    most_connected_nodes = most_connected_nodes_all
    both_centrality_nodes = sorted(set(bridge_anchor_nodes) & set(most_connected_nodes))
    small_component_centrality_nodes = sorted(
        node for node in (set(bridge_anchor_nodes) | set(most_connected_nodes))
        if component_size_by_node.get(node, 0) < 3
    )
    if small_component_centrality_nodes:
        small_summary = ", ".join(
            f"C{component_id_by_node[node] + 1}:{node}"
            for node in small_component_centrality_nodes
        )
        print(
            "[INFO] Displaying selected visual centrality overlays for "
            f"components with <3 nodes: {small_summary}; "
            "centrality interpretation is limited for two-node components."
        )

    # ----------------------------------------------------- Draw
    if is_brca:
        component_palette = gps.BRCA_NO_GREEN_COMPONENT_PALETTE
    elif is_nbl:
        component_palette = gps.NBL_NO_GREEN_COMPONENT_PALETTE
    else:
        component_palette = gps.OKABE_ITO_PALETTE
    node_colour, _, _ = gps.component_colour_map(
        G, isolates=isolates, palette=component_palette
    )

    fig, ax = plt.subplots(figsize=(args.fig_w, args.fig_h), dpi=args.dpi)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    if iso_bbox is not None:
        if is_brca or is_multicohort:
            gps.draw_isolate_panel(
                ax,
                iso_bbox,
                isolate_nodes=isolates,
                pos=pos,
                labels_by_node={n: short_map.get(n, n) for n in isolates},
                panel_padding_x=args.isolate_panel_padding_x,
                panel_padding_y=args.isolate_panel_padding_y,
                label_band_frac=args.isolate_label_band_frac,
                node_label_fontsize=max(args.font_size - 1, 10),
            )
        else:
            gps.draw_isolate_zone(
                ax,
                iso_bbox,
                label_y_frac=args.isolate_label_y_frac,
            )

    # Edges. Multicohort thesis resolved graph uses uniform solid edges; support
    # counts remain in the sidecar table but do not drive figure styling.
    has_dashed = False
    support_values_for_legend: list[int] = []
    if G.number_of_edges():
        edges_ordered = sorted(G.edges())
        if is_multicohort:
            nx.draw_networkx_edges(
                G, pos, ax=ax,
                edgelist=edges_ordered,
                width=2.0,
                alpha=0.82, edge_color="black", style="solid",
            )
        else:
            def _support(e: tuple[str, str]) -> int | None:
                key = tuple(sorted(e))
                return edge_support.get(key)

            supports = [_support(e) for e in edges_ordered]
            support_values_for_legend = [(s if s is not None else 1) for s in supports]
            widths = gps.edge_widths_from_support(
                support_values_for_legend,
                min_width=args.edge_width_min,
                max_width=args.edge_width_max,
            )

            solid_idx = [i for i, s in enumerate(supports) if (s is not None and s >= 2)]
            dashed_idx = [i for i in range(len(supports)) if i not in solid_idx]

            if solid_idx:
                nx.draw_networkx_edges(
                    G, pos, ax=ax,
                    edgelist=[edges_ordered[i] for i in solid_idx],
                    width=[widths[i] for i in solid_idx],
                    alpha=0.85, edge_color="black", style="solid",
                )
            if dashed_idx:
                has_dashed = True
                nx.draw_networkx_edges(
                    G, pos, ax=ax,
                    edgelist=[edges_ordered[i] for i in dashed_idx],
                    width=[widths[i] for i in dashed_idx],
                    alpha=0.85, edge_color="black", style=(0, (5, 3)),
                )

    non_iso_nodes = [n for n in G.nodes() if n not in set(isolates)]
    if non_iso_nodes:
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            nodelist=non_iso_nodes,
            node_size=args.node_size,
            node_color=[node_colour[n] for n in non_iso_nodes],
            edgecolors=gps.NODE_EDGE_COLOUR,
            linewidths=1.3, alpha=0.95,
        )

    if isolates and not (is_brca or is_multicohort):
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            nodelist=isolates,
            node_size=int(args.node_size * 0.85),
            node_color=gps.ISOLATE_COLOUR,
            node_shape="D",
            edgecolors=gps.NODE_EDGE_COLOUR,
            linewidths=1.3, alpha=0.95,
        )

    # Centrality overlays: ring for max betweenness, square for max degree.
    if bridge_anchor_nodes:
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            nodelist=bridge_anchor_nodes,
            node_size=int(args.node_size * 1.55),
            node_color="none",
            edgecolors=gps.ANCHOR_RING_COLOUR,
            linewidths=3.6, alpha=1.0,
        )
    if most_connected_nodes:
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            nodelist=most_connected_nodes,
            node_size=int(args.node_size * 1.95),
            node_color="none",
            node_shape="s",
            edgecolors=gps.MOST_CONNECTED_COLOUR,
            linewidths=3.0, alpha=1.0,
        )

    labels = {n: short_map.get(n, n) for n in G.nodes()}
    if not is_brca:
        gps.draw_component_labels(
            ax, G, pos, isolates=isolates, min_size=2,
            min_y_offset=0.30 if is_rbl else (0.55 if is_multicohort else 0.40),
            y_offset_scale=0.045 if (is_rbl or is_multicohort) else 0.08,
            same_row_top_n=2 if is_nbl else 0,
            n_cols_for_rows=3 if is_multicohort else 0,
            use_bbox=False,
            font_color="black",
            halo_width=max(1.4, float(args.halo_linewidth) * 0.65),
        )
    label_nodes = [
        n for n in G.nodes()
        if not ((is_brca or is_multicohort) and n in set(isolates))
    ]
    text_items = gps.halo_labels(
        ax, G, pos,
        {n: labels[n] for n in label_nodes},
        font_size=args.font_size,
        font_color="black" if args.use_halo_labels else "#1A1A2E",
        halo_width=float(args.halo_linewidth) if args.use_halo_labels else 0.0,
        halo_color=args.halo_color,
    )
    decollide_iters = int(args.max_label_overlap_iterations)
    decollide_step = max(0.25, float(args.min_label_separation))
    decollide_expand = float(args.label_padding_factor)
    decollide_char_width = 8.5
    decollide_char_height = 11.0
    if is_multicohort:
        decollide_iters = max(decollide_iters, 1200)
        decollide_step = max(decollide_step, 0.42)
        decollide_expand = max(decollide_expand, 1.55)
        decollide_char_width = 10.0
        decollide_char_height = 12.8
    if args.label_overlap_avoidance:
        gps.decollide_labels(
            ax,
            text_items,
            iters=decollide_iters,
            step=decollide_step,
            char_width_pts=decollide_char_width,
            char_height_pts=decollide_char_height,
            expand=decollide_expand,
            bbox_margin=args.label_bbox_margin,
            min_separation=args.min_label_separation,
            halo_width=float(args.halo_linewidth) if args.use_halo_labels else 0.0,
            halo_color=args.halo_color,
            context=f"{args.label} resolved neighbours",
        )
    else:
        gps.decollide_labels(
            ax,
            text_items,
            halo_width=float(args.halo_linewidth) if args.use_halo_labels else 0.0,
            halo_color=args.halo_color,
            context=f"{args.label} resolved neighbours",
        )
    for text in text_items.values():
        if not args.use_halo_labels:
            text.set_path_effects([])

    title = (
        "Resolved multicohort cell-line neighbourhood graph"
        if is_multicohort
        else f"{args.label} resolved neighbours"
    )
    ax.set_title(
        title,
        fontsize=15 if is_multicohort else 14,
        pad=12 if is_multicohort else 20,
    )
    caption = (
        "node fill colour denotes connected component only. "
        "C labels are ordered by component size within this graph; "
        "component colours should not be compared across cohorts or graph types.\n"
        "Resolved edges are final resolved-neighbour relationships. "
        "Bridge-like and most-connected anchors are selected within connected components "
        "and displayed wherever an anchor assignment is available; centrality "
        "interpretation is limited for two-node components."
    )
    if is_rbl:
        caption = (
            "node fill colour denotes connected component only. "
            "C labels are ordered by component size within this graph;\n"
            "component colours should not be compared across cohorts or graph types. "
            "Resolved edges are final resolved-neighbour relationships.\n"
            "Bridge-like and most-connected anchors are selected within connected components "
            "and displayed wherever an anchor assignment is available; centrality "
            "interpretation is limited for two-node components."
        )
    if is_multicohort:
        caption = (
            "Node fill colour denotes connected component only; C labels are ordered "
            "by component size within this graph. Resolved edges are drawn uniformly.\n"
            "Bridge-like and most-connected anchors are selected within connected "
            "components; interpretation is limited for two-node components."
        )
    gps.add_caption(fig, caption, fontsize=7 if is_multicohort else 8)

    extra_bboxes = []
    if iso_bbox is not None:
        extra_bboxes.append(iso_bbox)
    if (is_brca or is_multicohort) and planner_regions.get("legend_region") is not None:
        extra_bboxes.append(planner_regions["legend_region"])
    gps.frame_axes(
        ax, pos, margin=args.axis_margin_frac,
        extra_bboxes=extra_bboxes if extra_bboxes else None,
    )
    if is_multicohort:
        gps.decollide_labels(
            ax,
            text_items,
            iters=900,
            step=0.34,
            char_width_pts=9.8,
            char_height_pts=12.8,
            expand=1.42,
        )
        gps.frame_axes(
            ax,
            pos,
            margin=0.035,
            extra_bboxes=extra_bboxes if extra_bboxes else None,
        )
    if is_brca:
        gps.place_component_labels_safely(
            ax,
            G,
            pos,
            isolates=isolates,
            labels_by_node=labels,
            label_gap=args.component_label_gap,
            context=f"{args.label} resolved neighbours",
            extra_data_bboxes={
                "isolate_box": iso_bbox,
                "legend_region": planner_regions.get("legend_region"),
            },
            fontsize=10,
            use_bbox=False,
            halo_width=max(1.4, float(args.halo_linewidth) * 0.65),
            halo_color=args.halo_color,
        )
        gps.frame_axes(
            ax, pos, margin=args.axis_margin_frac,
            extra_bboxes=extra_bboxes if extra_bboxes else None,
        )
    legend_anchor = (0.72, 0.90) if is_rbl else (1.01, 1.0)
    legend_loc = "upper left"
    legend_bbox_transform = None
    legend_coordinate_system = "axes"
    if not is_rbl and args.legend_alignment_mode == "component-band":
        legend_mid_y, legend_band, band_ymin, band_ymax = component_band_mid_y_axes(
            ax, pos, components_list, component_indices=(0, 1)
        )
        legend_anchor = (1.01, legend_mid_y)
        legend_loc = "center left"
        print(
            f"[QC] Resolved-graph legend y-anchor derived from "
            f"{legend_band} band: data_y_range=({band_ymin:.2f},{band_ymax:.2f}) "
            f"anchor=({legend_anchor[0]:.3f},{legend_anchor[1]:.3f})"
        )
    if is_brca and planner_regions.get("legend_region") is not None:
        legend_anchor = gps.legend_anchor_from_region(planner_regions["legend_region"])
        legend_loc = "center left"
        legend_bbox_transform = ax.transData
        legend_coordinate_system = "data"
        print(
            f"[QC] BRCA resolved neighbours legend placement: "
            f"loc={legend_loc} anchor=({legend_anchor[0]:.2f},{legend_anchor[1]:.2f}) "
            f"coordinate_system=data reserved_region="
            f"{gps._fmt_bbox(planner_regions['legend_region'])}"
        )
    if is_multicohort and planner_regions.get("legend_region") is not None:
        legend_region = planner_regions["legend_region"]
        legend_anchor = (
            float(legend_region[0]) + 1.0,
            float(legend_region[3]) - 1.0,
        )
        legend_loc = "upper left"
        legend_bbox_transform = ax.transData
        legend_coordinate_system = "data"
        print(
            f"[QC] Multicohort resolved neighbours legend placement: "
            f"loc={legend_loc} anchor=({legend_anchor[0]:.2f},"
            f"{legend_anchor[1]:.2f}) coordinate_system=data "
            f"reserved_region={gps._fmt_bbox(legend_region)}"
        )
    legend_handles = gps.legend_handles_resolved(
        has_isolates=bool(isolates),
        has_bridge_anchor=bool(bridge_anchor_nodes),
        has_most_connected=bool(most_connected_nodes),
        has_both=bool(both_centrality_nodes),
        has_dashed=has_dashed,
        support_values=support_values_for_legend,
        edge_width_min=args.edge_width_min,
        edge_width_max=args.edge_width_max,
    )
    if is_multicohort:
        legend_handles = compact_resolved_legend_handles(legend_handles)
        legend_handles = remove_edge_support_legend_handles(legend_handles)
    legend_title = "Encoding" if is_multicohort else "Resolved-graph legend"
    legend_spacious = not is_multicohort
    legend_labelspacing = 1.05 if is_multicohort else 1.5
    if args.legend_mode == "keep":
        gps.place_legend(
            ax,
            legend_handles,
            title=legend_title,
            spacious=legend_spacious,
            labelspacing_override=legend_labelspacing,
            bbox_to_anchor=legend_anchor,
            loc=legend_loc,
            bbox_transform=legend_bbox_transform,
            coordinate_system=legend_coordinate_system,
        )
    elif args.legend_mode == "separate":
        legend_out_prefix = args.legend_out_prefix or f"{args.output_prefix}_legend"
        gps.save_legend_only(
            legend_handles,
            legend_out_prefix,
            title=legend_title,
            dpi=args.dpi,
            spacious=legend_spacious,
            labelspacing_override=legend_labelspacing,
        )

    if args.legend_mode != "keep":
        fig.tight_layout(rect=[0.0, 0.04, 0.98, 0.97])
    else:
        if is_rbl:
            fig.tight_layout(rect=[0.02, 0.09, 0.98, 0.95])
            gps.align_legend_to_component_box(
                fig, ax, pos, non_iso_nodes, "RBL resolved neighbours",
                component_texts=[labels[n] for n in non_iso_nodes],
                iso_bbox=iso_bbox,
                target_gap_frac=gps.RBL_SINGLE_COMPONENT_LEGEND_GAP_FRAC,
            )
        elif is_multicohort:
            fig.tight_layout(rect=[0.01, 0.055, 0.99, 0.94])
        else:
            fig.tight_layout(rect=[0.0, 0.06, 0.98, 0.95] if is_brca else [0.0, 0.04, 0.80, 0.97])
    gps.report_rendered_layout_clearance(
        fig, ax, G, pos, isolates,
        f"{args.label} resolved neighbours",
        labels_by_node=labels,
        iso_bbox=iso_bbox,
    )
    if is_nbl:
        gps.report_c_label_legend_metrics(fig, ax, "NBL resolved neighbours")
    save_kw = dict(facecolor="white")
    if is_multicohort:
        save_kw.update(bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{args.output_prefix}.png", dpi=args.dpi, **save_kw)
    fig.savefig(f"{args.output_prefix}.pdf", **save_kw)
    fig.savefig(f"{args.output_prefix}.svg", **save_kw)
    print(f"[OK] Saved: {args.output_prefix}.{{png,pdf,svg}}")
    plt.close(fig)

    # -------------------------------------------------------------------- Node stats TSV
    node_stats_output = sidecar_output_path(args.output_prefix, "resolved_cell_line_neighbourhood_graph_node_stats.tsv")
    comp_by_node = {node: i for i, comp in enumerate(components_list) for node in comp}
    # Isolates each form their own (unpacked) component; encode them with a
    # negative-indexed pseudo-component so consumers can identify them by
    # `is_isolate=True` and still keep a non-null component id.
    for k, n in enumerate(isolates):
        comp_by_node[n] = -(k + 1)

    anchor_set = set(anchor_nodes)
    node_stats_rows = []
    for node in sorted(G.nodes()):
        cid = comp_by_node.get(node, -1)
        sel_unnorm = selected_by_unnormalised.get(cid)
        sel_norm = selected_by_normalised.get(cid)
        sel_degree = selected_by_degree.get(cid)
        node_stats_rows.append({
            "cell_line": node,
            "cell_line_short": short_map.get(node, node),
            "degree": G.degree(node),
            "betweenness": betweenness_unnormalised.get(node, 0.0),
            "betweenness_unnormalised": betweenness_unnormalised.get(node, 0.0),
            "betweenness_normalised": betweenness_normalised.get(node, 0.0),
            "component": cid,
            "is_isolate": node in set(isolates),
            "is_central": node in anchor_set,
            "selected_by_unnormalised": node == sel_unnorm if sel_unnorm else False,
            "selected_by_normalised": node == sel_norm if sel_norm else False,
            "selected_by_degree": node == sel_degree if sel_degree else False,
            "canonical_selected": node == sel_unnorm if sel_unnorm else False,
            "canonical_bridge_selected": node == sel_unnorm if sel_unnorm else False,
            "most_connected_selected": node == sel_degree if sel_degree else False,
            "centrality_metric": "degree_and_betweenness",
            "centrality_normalised": False,
            "centrality_weighted": False,
            "centrality_scope": "within_component",
            "centrality_tie_break": "highest_value_then_node_id",
        })
    node_stats_df = pd.DataFrame(node_stats_rows)
    write_tsv_if_changed(node_stats_df, node_stats_output)
    write_legacy_sidecar_if_needed(
        node_stats_df,
        args.output_prefix,
        "resolved_cell_line_neighbourhood_graph_node_stats.tsv",
    )

    # ----------------------------------------------------------------- Anchor audit TSV
    anchor_audit_output = sidecar_output_path(args.output_prefix, "resolved_cell_line_neighbourhood_anchor_centrality_audit.tsv")
    audit_df = node_stats_df.rename(columns={
        "cell_line": "node_id",
        "cell_line_short": "display_label",
        "component": "component_id",
    })[[
        "component_id", "node_id", "display_label", "degree",
        "betweenness_normalised", "betweenness_unnormalised",
        "selected_by_normalised", "selected_by_unnormalised",
        "selected_by_degree", "canonical_selected", "canonical_bridge_selected",
        "most_connected_selected",
    ]]
    audit_df.insert(0, "cohort", args.label)
    audit_df["changed_relative_to_legacy"] = (
        audit_df["selected_by_normalised"] != audit_df["selected_by_unnormalised"]
    )
    audit_df["degree_vs_betweenness_node_differs"] = (
        audit_df["selected_by_degree"] != audit_df["selected_by_unnormalised"]
    )
    audit_df["canonical_metric"] = "degree_and_unnormalised_betweenness"
    audit_df["alternative_metric"] = "normalised_betweenness"
    audit_df["centrality_scope"] = "within_component"
    audit_df["centrality_weighted"] = False
    audit_df["centrality_tie_break"] = "highest_value_then_node_id"
    audit_df["degree_anchor_selected"] = audit_df["most_connected_selected"]
    audit_df["bridge_betweenness_selected"] = audit_df["canonical_bridge_selected"]
    audit_df["anchor_selected"] = (
        audit_df["degree_anchor_selected"] | audit_df["bridge_betweenness_selected"]
    )

    def _anchor_reason(row):
        by_deg = bool(row["degree_anchor_selected"])
        by_btw = bool(row["bridge_betweenness_selected"])
        if by_deg and by_btw:
            return "degree;betweenness_unnormalised"
        if by_deg:
            return "degree"
        if by_btw:
            return "betweenness_unnormalised"
        return ""

    audit_df["anchor_selection_reason"] = audit_df.apply(_anchor_reason, axis=1)
    write_tsv_if_changed(audit_df, anchor_audit_output)
    write_legacy_sidecar_if_needed(
        audit_df,
        args.output_prefix,
        "resolved_cell_line_neighbourhood_anchor_centrality_audit.tsv",
    )


if __name__ == "__main__":
    main()
