#!/usr/bin/env python3
"""
plot_consensus_graph.py
-------------------------------------------------
Visualises the undirected support-threshold consensus cell-line similarity
network for one cohort.

Visual language (single source of truth: scripts/graph_plot_style.py):

  - Node fill encodes connected component using a colour-blind-friendly
    qualitative palette with no green hue. Isolates share one neutral grey.
    Component colours are structural identifiers only.
  - Isolates (degree 0) are drawn as diamonds inside a labelled
    "Isolates (degree 0)" strip below the rest of the layout.
  - Edge width is proportional to support count (support_directions).
    The consensus network has no single-supported edges by construction
    (support floor = ceil(n_directions / 2)), so solid/dashed style is
    not used here.
  - Labels are short names with a thin white halo (no boxed background).
  - A legend documents every encoding actually present, placed outside
    the axes so it never occludes nodes.
  - No central / anchor marker on the consensus figure.

Layout: all cohorts use graph_plot_style.grid_layout (adaptive component
shelves, independent per-component spring layout, fixed seed).
NBL support-threshold figures additionally compact the first two components
after layout, without changing graph content.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import networkx as nx
import numpy as np
import pandas as pd

plt.rcParams["svg.fonttype"] = "none"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import graph_plot_style as gps  # noqa: E402


RBL_SUPPORT_RESOLVED_REFERENCE_ISOLATE_GAP = 4.6
NBL_SUPPORT_COMPONENT_GAP_FRACTION = 0.57
NBL_SUPPORT_COMPONENT_MIN_GAP = 9.5
NBL_SUPPORT_LEGEND_GAP_FRAC = 0.025
NBL_SUPPORT_ISOLATE_LABEL_Y_FRAC = 0.02
MULTICOHORT_MAJORITY_GRID_N_COLS = 3
MULTICOHORT_MAJORITY_GRID_N_ROWS = 3
MULTICOHORT_UNION_LAYOUT_GAP = 14.0
MULTICOHORT_UNION_COMPONENT_PALETTE = (
    "#5F8FB7",  # muted blue
    "#C58AA5",  # muted reddish-purple
    "#D6A65F",  # muted amber
    "#8F9BC5",  # muted periwinkle
    "#A9825A",  # muted brown
)


def compact_nbl_support_components(
    G,
    pos,
    isolates,
    *,
    gap_fraction=NBL_SUPPORT_COMPONENT_GAP_FRACTION,
    min_gap=NBL_SUPPORT_COMPONENT_MIN_GAP,
):
    # Move C2 toward C1 after the deterministic component layout is complete.
    components = gps.ordered_non_isolate_components(G, isolates=isolates)
    if len(components) < 2:
        return

    c1, c2 = components[0], components[1]
    c1_pts = np.array([pos[n] for n in c1 if n in pos], dtype=float)
    c2_pts = np.array([pos[n] for n in c2 if n in pos], dtype=float)
    if not c1_pts.size or not c2_pts.size:
        return

    c1_x0, c1_x1 = float(c1_pts[:, 0].min()), float(c1_pts[:, 0].max())
    c2_x0, c2_x1 = float(c2_pts[:, 0].min()), float(c2_pts[:, 0].max())
    c1_mid = (c1_x0 + c1_x1) / 2.0
    c2_mid = (c2_x0 + c2_x1) / 2.0

    if c2_mid >= c1_mid:
        current_gap = c2_x0 - c1_x1
        direction = 1.0
    else:
        current_gap = c1_x0 - c2_x1
        direction = -1.0
    if current_gap <= 0:
        return

    target_gap = max(
        current_gap * float(gap_fraction),
        float(min_gap),
    )
    if target_gap >= current_gap:
        return

    delta = direction * (target_gap - current_gap)
    for node in c2:
        if node in pos:
            pos[node] = np.array(
                [float(pos[node][0]) + delta, float(pos[node][1])],
                dtype=float,
            )

    print(
        "[QC] NBL support-threshold C1/C2 gap compaction: "
        f"old_gap={current_gap:.2f} new_gap={target_gap:.2f} "
        f"reduction={(1.0 - target_gap / current_gap) * 100.0:.1f}%"
    )


def component_band_mid_y_axes(ax, pos, nodes):
    pts = np.array([pos[n] for n in nodes if n in pos], dtype=float)
    if not pts.size:
        return 0.5, float("nan"), float("nan")
    y_min = float(pts[:, 1].min())
    y_max = float(pts[:, 1].max())
    mid_y = (y_min + y_max) / 2.0
    _, axes_y = ax.transAxes.inverted().transform(
        ax.transData.transform((0.0, mid_y))
    )
    return float(np.clip(axes_y, 0.15, 0.85)), y_min, y_max


def widen_nbl_support_target_isolates(G, pos, iso_bbox, min_gap=13.0):
    # Widen the NBL support isolate row around long labels without changing topology.
    targets = ('CCLF_PEDS_0051_T', 'GI_ME_N', 'LAN2')
    if not all(node in G and node in pos for node in targets):
        return iso_bbox
    if any(G.degree(node) != 0 for node in targets):
        return iso_bbox

    isolate_nodes = sorted(n for n in G.nodes() if G.degree(n) == 0 and n in pos)
    if not isolate_nodes:
        return iso_bbox

    original_x = [float(pos[node][0]) for node in isolate_nodes]
    centre = (min(original_x) + max(original_x)) / 2.0
    start_x = centre - (len(isolate_nodes) - 1) * min_gap / 2.0

    for i, node in enumerate(isolate_nodes):
        pos[node] = np.array([start_x + i * min_gap, float(pos[node][1])], dtype=float)

    if iso_bbox is not None:
        _, y0, _, y1 = iso_bbox
        pad = 6.0
        row_x = [float(pos[node][0]) for node in isolate_nodes]
        iso_bbox = (min(row_x) - pad, float(y0), max(row_x) + pad, float(y1))

    gaps = {
        f'{left}->{right}': abs(float(pos[right][0]) - float(pos[left][0]))
        for left, right in zip(targets, targets[1:])
    }
    print(
        '[INFO] NBL support-threshold isolate row spacing: '
        f'ordered isolates={isolate_nodes}; '
        f'target_gaps=' + ','.join(f'{k}={v:.2f}' for k, v in gaps.items())
    )
    return iso_bbox


def multicohort_majority_adaptive_grid_layout(
    G_non_iso,
    isolates,
    *,
    seed,
    context="MULTICOHORT_CANCER majority-threshold consensus",
):
    """
    Place multicohort majority components in a structured 3 x 3 grid.

    Component cell size is driven by both node count and edge count. The helper
    changes only plotting coordinates; graph nodes, edges, support values, and
    component membership are untouched.
    """
    components = [sorted(comp) for comp in gps.ordered_non_isolate_components(G_non_iso)]
    if not components:
        return {}, None, [], {}

    n_cols = MULTICOHORT_MAJORITY_GRID_N_COLS
    n_rows = max(
        MULTICOHORT_MAJORITY_GRID_N_ROWS,
        int(np.ceil(len(components) / float(n_cols))),
    )
    col_gap = 15.0
    row_gap = 13.0
    legend_min_width = 62.0
    isolate_min_height = 22.0
    bottom_row_min_height = 30.0

    specs = []
    for idx, nodes in enumerate(components, start=1):
        n_nodes = len(nodes)
        n_edges = G_non_iso.subgraph(nodes).number_of_edges()
        density = gps.component_edge_density(G_non_iso, nodes)
        edge_term = np.sqrt(max(float(n_edges), 0.0))
        node_term = np.sqrt(max(float(n_nodes), 1.0))
        density_boost = 1.0 + min(float(density), 0.85) * 0.24

        cell_width = (22.0 + 6.8 * node_term + 5.1 * edge_term) * density_boost
        cell_height = (16.0 + 4.8 * node_term + 3.7 * edge_term) * density_boost
        cell_width = float(np.clip(cell_width, 36.0, 84.0))
        cell_height = float(np.clip(cell_height, 26.0, 62.0))
        complexity = float(n_nodes + 0.55 * n_edges)
        specs.append({
            "idx": idx,
            "nodes": nodes,
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            "density": density,
            "complexity": complexity,
            "cell_width": cell_width,
            "cell_height": cell_height,
        })

    complexities = [float(s["complexity"]) for s in specs]
    cmin, cmax = min(complexities), max(complexities)
    denom = max(cmax - cmin, 1e-9)

    col_widths = [0.0] * n_cols
    row_heights = [0.0] * n_rows
    for spec in specs:
        row = (int(spec["idx"]) - 1) // n_cols
        col = (int(spec["idx"]) - 1) % n_cols
        col_widths[col] = max(col_widths[col], float(spec["cell_width"]))
        row_heights[row] = max(row_heights[row], float(spec["cell_height"]))
    if len(components) <= 7 and n_rows >= 3:
        col_widths[2] = max(col_widths[2], legend_min_width)
        row_heights[2] = max(row_heights[2], bottom_row_min_height, isolate_min_height)
    col_widths = [max(w, 34.0) for w in col_widths]
    row_heights = [max(h, 24.0) for h in row_heights]

    total_width = sum(col_widths) + col_gap * (n_cols - 1)
    total_height = sum(row_heights) + row_gap * (n_rows - 1)
    grid_x0 = -total_width / 2.0
    grid_y1 = 0.0
    grid_y0 = grid_y1 - total_height

    x_starts = []
    cursor = grid_x0
    for width in col_widths:
        x_starts.append(cursor)
        cursor += width + col_gap
    y_tops = []
    cursor = grid_y1
    for height in row_heights:
        y_tops.append(cursor)
        cursor -= height + row_gap

    pos_non_iso = {}
    placements = []
    for spec in specs:
        idx = int(spec["idx"])
        row = (idx - 1) // n_cols
        col = (idx - 1) % n_cols
        cell_x0 = x_starts[col]
        cell_x1 = cell_x0 + col_widths[col]
        cell_y1 = y_tops[row]
        cell_y0 = cell_y1 - row_heights[row]
        cx = (cell_x0 + cell_x1) / 2.0
        cy = (cell_y0 + cell_y1) / 2.0

        norm = (float(spec["complexity"]) - cmin) / denom
        fill_x = 0.78 + 0.18 * norm
        fill_y = 0.72 + 0.20 * norm
        if idx <= 3:
            fill_x = max(fill_x, 0.91)
            fill_y = max(fill_y, 0.87)
        if int(spec["n_nodes"]) <= 2:
            fill_x = min(fill_x, 0.70)
            fill_y = min(fill_y, 0.52)
        target_width = col_widths[col] * float(np.clip(fill_x, 0.50, 0.94))
        target_height = row_heights[row] * float(np.clip(fill_y, 0.42, 0.92))
        local_pos = gps._layout_component_for_cell(
            G_non_iso,
            list(spec["nodes"]),
            seed=seed,
            target_width=target_width,
            target_height=target_height,
            spring_iters=700,
            spring_k=0.82 if int(spec["n_nodes"]) >= 8 else None,
        )
        for node, p in local_pos.items():
            pos_non_iso[node] = np.array(
                [cx + float(p[0]), cy + float(p[1])],
                dtype=float,
            )

        pts = np.array([pos_non_iso[n] for n in spec["nodes"] if n in pos_non_iso], dtype=float)
        bbox = (
            float(pts[:, 0].min()), float(pts[:, 1].min()),
            float(pts[:, 0].max()), float(pts[:, 1].max()),
        )
        cell_bbox = (float(cell_x0), float(cell_y0), float(cell_x1), float(cell_y1))
        placement = {
            "component": f"C{idx}",
            "component_index": idx,
            "row": row + 1,
            "col": col + 1,
            "n": int(spec["n_nodes"]),
            "edges": int(spec["n_edges"]),
            "density": float(spec["density"]),
            "complexity": float(spec["complexity"]),
            "cell_bbox": cell_bbox,
            "bbox": bbox,
            "fill_x": float(fill_x),
            "fill_y": float(fill_y),
        }
        placements.append(placement)
        print(
            f"[QC] {context} adaptive placement: C{idx} "
            f"row={row + 1} col={col + 1} n={int(spec['n_nodes'])} "
            f"edges={int(spec['n_edges'])} density={float(spec['density']):.3f} "
            f"complexity={float(spec['complexity']):.2f} "
            f"cell={gps._fmt_bbox(cell_bbox)} bbox={gps._fmt_bbox(bbox)} "
            f"fill=({float(fill_x):.2f},{float(fill_y):.2f})"
        )

    regions = {
        "component_grid": (grid_x0, grid_y0, grid_x0 + total_width, grid_y1),
        "title_region": (grid_x0, grid_y1 + 2.0, grid_x0 + total_width, grid_y1 + 8.0),
        "footnote_region": (grid_x0, grid_y0 - 8.0, grid_x0 + total_width, grid_y0 - 2.0),
    }

    legend_region = None
    if n_rows >= 3:
        legend_region = (
            float(x_starts[2]),
            float(y_tops[2] - row_heights[2]),
            float(x_starts[2] + col_widths[2]),
            float(y_tops[2]),
        )
        regions["legend_region"] = legend_region

    iso_bbox = None
    pos_iso = {}
    if isolates:
        iso_sorted = sorted(isolates)
        if n_rows >= 3 and len(components) <= 7:
            iso_bbox = (
                float(x_starts[1]),
                float(y_tops[2] - row_heights[2]),
                float(x_starts[1] + col_widths[1]),
                float(y_tops[2]),
            )
        else:
            panel_width = min(max(total_width * 0.45, len(iso_sorted) * 12.0 + 18.0), total_width)
            iso_y1 = grid_y0 - 5.0
            iso_bbox = (
                -panel_width / 2.0,
                iso_y1 - isolate_min_height,
                panel_width / 2.0,
                iso_y1,
            )
        label_band_frac = 0.33
        panel_padding_x = 8.0
        panel_padding_y = 1.2
        label_band_h = (float(iso_bbox[3]) - float(iso_bbox[1])) * label_band_frac
        node_y0 = float(iso_bbox[1]) + label_band_h + panel_padding_y
        node_y1 = float(iso_bbox[3]) - panel_padding_y
        y_iso = (node_y0 + node_y1) / 2.0
        if len(iso_sorted) == 1:
            xs = np.array([(float(iso_bbox[0]) + float(iso_bbox[2])) / 2.0])
        else:
            xs = np.linspace(
                float(iso_bbox[0]) + panel_padding_x,
                float(iso_bbox[2]) - panel_padding_x,
                len(iso_sorted),
            )
        pos_iso = {
            node: np.array([float(x), y_iso], dtype=float)
            for node, x in zip(iso_sorted, xs)
        }
        regions["isolate_box"] = iso_bbox
        regions["isolate_region"] = iso_bbox
        print(
            f"[QC] {context} isolate panel: bbox={gps._fmt_bbox(iso_bbox)} "
            f"isolates={','.join(iso_sorted)}"
        )

    gps.report_planned_layout_clearance(context, placements, regions)
    return {**pos_non_iso, **pos_iso}, iso_bbox, placements, regions


def layout_union_large_component_by_modules(
    G_non_iso,
    nodes,
    *,
    seed,
    target_width,
    target_height,
    context,
):
    """
    Spread a large union-supported component by internal graph modules.

    This is still one connected component in the plotted graph; only the node
    coordinates are arranged module-by-module to make provenance plots readable.
    """
    ordered_nodes = sorted(nodes)
    sub = G_non_iso.subgraph(ordered_nodes).copy()
    communities = [
        sorted(comm)
        for comm in nx.algorithms.community.greedy_modularity_communities(sub)
    ]
    communities = sorted(communities, key=lambda c: (-len(c), c[0]))
    if len(communities) < 2 or len(communities) > 6:
        return gps._layout_component_for_cell(
            G_non_iso,
            ordered_nodes,
            seed=seed,
            target_width=target_width,
            target_height=target_height,
            spring_iters=950,
            spring_k=0.70,
        )

    target_width = float(target_width)
    target_height = float(target_height)
    if len(communities) == 3:
        region_specs = [
            (-0.29, -0.08, 0.54, 0.45),
            (0.30, -0.34, 0.45, 0.43),
            (-0.36, 0.43, 0.38, 0.30),
        ]
    else:
        n_cols = 2 if len(communities) <= 4 else 3
        n_rows = int(np.ceil(len(communities) / float(n_cols)))
        region_specs = []
        for i in range(len(communities)):
            row = i // n_cols
            col = i % n_cols
            x_frac = (col + 0.5) / n_cols - 0.5
            y_frac = 0.5 - (row + 0.5) / n_rows
            region_specs.append((
                x_frac,
                y_frac,
                0.82 / n_cols,
                0.78 / n_rows,
            ))

    local = {}
    for idx, comm_nodes in enumerate(communities, start=1):
        x_frac, y_frac, w_frac, h_frac = region_specs[min(idx - 1, len(region_specs) - 1)]
        comm_w = target_width * float(w_frac)
        comm_h = target_height * float(h_frac)
        comm_cx = target_width * float(x_frac)
        comm_cy = target_height * float(y_frac)
        comm_pos = gps._layout_component_for_cell(
            sub,
            comm_nodes,
            seed=seed + idx,
            target_width=comm_w,
            target_height=comm_h,
            spring_iters=500,
            spring_k=None,
        )
        for node, p in comm_pos.items():
            local[node] = np.array(
                [comm_cx + float(p[0]), comm_cy + float(p[1])],
                dtype=float,
            )
        print(
            f"[QC] {context} C1 module layout: Q{idx} "
            f"n={len(comm_nodes)} edges={sub.subgraph(comm_nodes).number_of_edges()} "
            f"center=({comm_cx:.2f},{comm_cy:.2f}) "
            f"target=({comm_w:.2f},{comm_h:.2f}) "
            f"nodes={';'.join(comm_nodes)}"
        )
    return local


def multicohort_union_supported_adaptive_layout(
    G_non_iso,
    isolates,
    *,
    seed,
    spacing_scale=1.0,
    context="MULTICOHORT_CANCER union-supported-edges network",
):
    """
    Place multicohort union-supported components in a two-zone provenance layout.

    The union graph has one large, dense provenance component plus a smaller
    component and isolates. The large component receives most of the page, while
    small structures and the legend are reserved in a compact right column.
    """
    components = [sorted(comp) for comp in gps.ordered_non_isolate_components(G_non_iso)]
    if not components:
        return {}, None, [], {}
    spacing_scale = max(1.0, float(spacing_scale))

    specs = []
    for idx, nodes in enumerate(components, start=1):
        n_nodes = len(nodes)
        n_edges = G_non_iso.subgraph(nodes).number_of_edges()
        density = gps.component_edge_density(G_non_iso, nodes)
        complexity = float(n_nodes + 0.45 * n_edges)
        specs.append({
            "idx": idx,
            "nodes": nodes,
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            "density": density,
            "complexity": complexity,
        })

    max_complexity = max(float(s["complexity"]) for s in specs)
    largest = specs[0]
    full_label_boost = min(spacing_scale, 1.45)
    large_width = float(np.clip(
        82.0
        + 5.7 * np.sqrt(max(float(largest["n_nodes"]), 1.0))
        + 3.4 * np.sqrt(max(float(largest["n_edges"]), 0.0))
        + 8.0,
        130.0,
        185.0,
    )) * full_label_boost
    large_height = float(np.clip(
        62.0
        + 4.5 * np.sqrt(max(float(largest["n_nodes"]), 1.0))
        + 2.7 * np.sqrt(max(float(largest["n_edges"]), 0.0))
        + 10.0,
        118.0,
        160.0,
    )) * full_label_boost

    right_width = 70.0 * min(spacing_scale, 1.35)
    top_small_height = 50.0 * min(spacing_scale, 1.35)
    isolate_height = 14.0 * min(spacing_scale, 1.25) if isolates else 0.0
    legend_height = 56.0 * min(spacing_scale, 1.30)
    gap = MULTICOHORT_UNION_LAYOUT_GAP * spacing_scale
    print(
        f"[QC] {context} adaptive spacing: scale={spacing_scale:.2f} "
        f"large_region=({large_width:.2f},{large_height:.2f}) "
        f"right_width={right_width:.2f} gap={gap:.2f}"
    )
    right_height = top_small_height + 8.0 + isolate_height + 10.0 + legend_height
    overall_height = max(large_height, right_height)
    total_width = large_width + gap + right_width
    x0 = -total_width / 2.0
    y1 = 0.0

    large_region = (x0, y1 - overall_height, x0 + large_width, y1)
    right_x0 = x0 + large_width + gap
    right_x1 = right_x0 + right_width
    c2_region = (right_x0, y1 - top_small_height, right_x1, y1)
    isolate_region = None
    if isolates:
        iso_y1 = c2_region[1] - 8.0
        iso_width = min(42.0 * min(spacing_scale, 1.45), right_width * 0.66)
        iso_cx = (right_x0 + right_x1) / 2.0
        isolate_region = (
            iso_cx - iso_width / 2.0,
            iso_y1 - isolate_height,
            iso_cx + iso_width / 2.0,
            iso_y1,
        )
        legend_y1 = isolate_region[1] - 10.0
    else:
        legend_y1 = c2_region[1] - 10.0
    legend_region = (right_x0, legend_y1 - legend_height, right_x1, legend_y1)

    pos_non_iso = {}
    placements = []
    fallback_regions = []
    if len(specs) > 2:
        spare_bottom = min(large_region[1] + 34.0, legend_region[1] - 8.0)
        for extra_i, spec in enumerate(specs[2:], start=0):
            fallback_regions.append((
                right_x0,
                spare_bottom - extra_i * 28.0,
                right_x1,
                spare_bottom + 22.0 - extra_i * 28.0,
            ))

    for spec in specs:
        idx = int(spec["idx"])
        if idx == 1:
            region = large_region
            fill_x, fill_y = 0.95, 0.93
            spring_iters = 1200
            spring_k = 0.86 if spacing_scale > 1.0 else 0.70
        elif idx == 2:
            region = c2_region
            fill_x, fill_y = 0.90, 0.88
            spring_iters = 650
            spring_k = 0.75 if spacing_scale > 1.0 else None
        else:
            region = fallback_regions[idx - 3] if idx - 3 < len(fallback_regions) else c2_region
            fill_x, fill_y = 0.72, 0.64
            spring_iters = 400
            spring_k = None

        rx0, ry0, rx1, ry1 = region
        region_w = float(rx1 - rx0)
        region_h = float(ry1 - ry0)
        complexity_ratio = float(spec["complexity"]) / max(max_complexity, 1e-9)
        target_width = region_w * min(0.95, fill_x + 0.03 * complexity_ratio)
        target_height = region_h * min(0.93, fill_y + 0.03 * complexity_ratio)
        if idx == 1 and int(spec["n_nodes"]) >= 24:
            local_pos = layout_union_large_component_by_modules(
                G_non_iso,
                list(spec["nodes"]),
                seed=seed,
                target_width=target_width,
                target_height=target_height,
                context=context,
            )
        else:
            local_pos = gps._layout_component_for_cell(
                G_non_iso,
                list(spec["nodes"]),
                seed=seed,
                target_width=target_width,
                target_height=target_height,
                spring_iters=spring_iters,
                spring_k=spring_k,
            )
        cx = (float(rx0) + float(rx1)) / 2.0
        cy = (float(ry0) + float(ry1)) / 2.0
        for node, p in local_pos.items():
            pos_non_iso[node] = np.array(
                [cx + float(p[0]), cy + float(p[1])],
                dtype=float,
            )

        pts = np.array([pos_non_iso[n] for n in spec["nodes"] if n in pos_non_iso], dtype=float)
        bbox = (
            float(pts[:, 0].min()), float(pts[:, 1].min()),
            float(pts[:, 0].max()), float(pts[:, 1].max()),
        )
        placement = {
            "component": f"C{idx}",
            "component_index": idx,
            "row": 1 if idx == 1 else idx,
            "col": 1 if idx == 1 else 2,
            "n": int(spec["n_nodes"]),
            "edges": int(spec["n_edges"]),
            "density": float(spec["density"]),
            "complexity": float(spec["complexity"]),
            "cell_bbox": tuple(float(v) for v in region),
            "bbox": bbox,
            "fill_x": float(fill_x),
            "fill_y": float(fill_y),
        }
        placements.append(placement)
        print(
            f"[QC] {context} adaptive placement: C{idx} "
            f"n={int(spec['n_nodes'])} edges={int(spec['n_edges'])} "
            f"density={float(spec['density']):.3f} "
            f"complexity={float(spec['complexity']):.2f} "
            f"region={gps._fmt_bbox(region)} bbox={gps._fmt_bbox(bbox)} "
            f"target=({target_width:.2f},{target_height:.2f})"
        )

    pos_iso = {}
    iso_bbox = None
    if isolates and isolate_region is not None:
        iso_sorted = sorted(isolates)
        iso_bbox = tuple(float(v) for v in isolate_region)
        iy0, iy1 = float(iso_bbox[1]), float(iso_bbox[3])
        ix0, ix1 = float(iso_bbox[0]), float(iso_bbox[2])
        label_band_h = (iy1 - iy0) * 0.34
        y_iso = (iy0 + label_band_h + 1.0 + iy1 - 1.0) / 2.0
        if len(iso_sorted) == 1:
            xs = np.array([(ix0 + ix1) / 2.0])
        else:
            iso_x_pad = 8.0 * min(spacing_scale, 1.35)
            xs = np.linspace(ix0 + iso_x_pad, ix1 - iso_x_pad, len(iso_sorted))
        pos_iso = {
            node: np.array([float(x), y_iso], dtype=float)
            for node, x in zip(iso_sorted, xs)
        }
        print(
            f"[QC] {context} isolate panel: bbox={gps._fmt_bbox(iso_bbox)} "
            f"isolates={','.join(iso_sorted)}"
        )

    regions = {
        "component_grid": (x0, y1 - overall_height, x0 + total_width, y1),
        "title_region": (x0, y1 + 2.0, x0 + total_width, y1 + 8.0),
        "legend_region": tuple(float(v) for v in legend_region),
        "footnote_region": (x0, y1 - overall_height - 8.0, x0 + total_width, y1 - overall_height - 2.0),
    }
    if iso_bbox is not None:
        regions["isolate_box"] = iso_bbox
        regions["isolate_region"] = iso_bbox
    gps.report_planned_layout_clearance(context, placements, regions)
    return {**pos_non_iso, **pos_iso}, iso_bbox, placements, regions


def union_thesis_label_nodes(
    G,
    isolates,
    *,
    pos=None,
    max_component_labels=14,
    min_component_nodes=16,
    min_component_density=0.25,
    min_label_distance=0.0,
    context_label="union-supported",
):
    """
    Select thesis-readable labels for dense union-supported graphs.

    All nodes remain drawn. The largest component is labelled by stable
    within-component degree/module representatives; smaller or sparse
    components and isolates keep full labels.
    """
    isolate_set = set(isolates)
    components = gps.ordered_non_isolate_components(G, isolates=isolates)
    if not components:
        return sorted(G.nodes())

    selected: set[str] = set(isolate_set)
    for idx, comp in enumerate(components, start=1):
        nodes = sorted(comp)
        sub = G.subgraph(nodes).copy()
        density = gps.component_edge_density(G, nodes)
        if (
            len(nodes) <= int(max_component_labels)
            or len(nodes) < int(min_component_nodes)
            or density < float(min_component_density)
        ):
            selected.update(nodes)
            print(
                f"[QC] {context_label} union-supported thesis labels: "
                f"C{idx} full_labels={len(nodes)}/{len(nodes)} "
                f"density={density:.3f}"
            )
            continue

        communities = [
            sorted(comm)
            for comm in nx.algorithms.community.greedy_modularity_communities(sub)
        ]
        communities = sorted(communities, key=lambda c: (-len(c), c[0]))

        candidate_order: list[str] = []
        for comm in communities:
            comm_ranked = sorted(
                comm,
                key=lambda n: (-sub.degree(n), n),
            )
            candidate_order.extend(comm_ranked[:4])

        degree_ranked = [n for n, _ in sorted(sub.degree, key=lambda x: (-x[1], x[0]))]
        candidate_order.extend(degree_ranked[:8])

        if sub.number_of_nodes() > 2:
            between = nx.betweenness_centrality(sub, normalized=False, weight=None)
            bridge_ranked = sorted(nodes, key=lambda n: (-between.get(n, 0.0), n))
            candidate_order.extend(bridge_ranked[:4])

        chosen_for_component: list[str] = []
        for node in candidate_order:
            if node in chosen_for_component:
                continue
            if pos is not None and float(min_label_distance) > 0.0:
                node_xy = np.array(pos.get(node, [np.nan, np.nan]), dtype=float)
                if np.isfinite(node_xy).all():
                    too_close = False
                    for previous in chosen_for_component:
                        previous_xy = np.array(pos.get(previous, [np.nan, np.nan]), dtype=float)
                        if (
                            np.isfinite(previous_xy).all()
                            and float(np.linalg.norm(node_xy - previous_xy)) < float(min_label_distance)
                        ):
                            too_close = True
                            break
                    if too_close:
                        continue
            selected.add(node)
            chosen_for_component.append(node)
            if len(chosen_for_component) >= int(max_component_labels):
                break
        if not chosen_for_component and candidate_order:
            selected.add(candidate_order[0])
            chosen_for_component.append(candidate_order[0])
        shown = sorted(selected.intersection(nodes))
        print(
            f"[QC] {context_label} union-supported thesis labels: "
            f"C{idx} selected_labels={len(shown)}/{len(nodes)} "
            f"density={density:.3f} min_label_distance={float(min_label_distance):.2f} "
            f"labels={';'.join(shown)}"
        )

    label_nodes = sorted(n for n in selected if n in G)
    print(
        f"[QC] {context_label} union-supported thesis labels: "
        f"total_labels={len(label_nodes)}/{G.number_of_nodes()} "
        f"isolates={';'.join(sorted(isolate_set))}"
    )
    return label_nodes


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Plot the support-threshold consensus cell-line similarity "
            "network. Width-encodes support_directions; component-coloured "
            "nodes; halo labels; isolate strip; legend."
        )
    )
    ap.add_argument("--edges", required=True,
                    help="Consensus edge TSV (node1, node2 required; "
                         "support_directions optional but recommended).")
    ap.add_argument("--out", required=True,
                    help="Output path prefix; .png/.pdf/.svg are appended.")
    ap.add_argument("--label", default="BRCA",
                    help="Dataset label shown in the figure title.")
    ap.add_argument("--graph-mode",
                    choices=("majority_threshold", "union_supported_edges"),
                    default="majority_threshold",
                    help="Graph product being plotted.")
    ap.add_argument("--graph-title", default=None,
                    help="Optional explicit figure title suffix, overriding the graph-mode title.")
    ap.add_argument("--node-stats-out", default=None,
                    help="Optional node-stat TSV written from the plotted graph.")
    ap.add_argument("--components-out", default=None,
                    help="Optional component TSV written from the plotted graph.")
    ap.add_argument("--node-labels-out", default=None,
                    help="Optional TSV documenting available node labels and whether each is shown.")
    ap.add_argument("--node-label-mode",
                    choices=("auto", "selected", "full"),
                    default="auto",
                    help="Node labels to draw. Auto keeps selected labels for dense union graphs and full labels elsewhere.")
    ap.add_argument("--selected-label-max-component-labels", type=int, default=14,
                    help="Maximum labels to show for one dense component when selected union labelling is active.")
    ap.add_argument("--selected-label-min-component-nodes", type=int, default=16,
                    help="Minimum component size before selected union labelling can hide labels.")
    ap.add_argument("--selected-label-min-component-density", type=float, default=0.25,
                    help="Minimum component density before selected union labelling can hide labels.")
    ap.add_argument("--selected-label-min-distance", type=float, default=0.0,
                    help="Minimum data-coordinate distance between selected labels in a dense union component.")

    ap.add_argument("--fig-w",  type=float, default=22.0)
    ap.add_argument("--fig-h",  type=float, default=16.0)
    ap.add_argument("--dpi",    type=int,   default=300)
    ap.add_argument("--use-halo-labels", action="store_true",
                    help="Keep white path-effect halos on node labels.")
    ap.add_argument("--halo-linewidth", type=float, default=3.6,
                    help="White halo stroke width for node labels.")
    ap.add_argument("--halo-color", default="white",
                    help="Halo stroke colour for node and component labels.")
    ap.add_argument("--match-resolved-style", action="store_true",
                    help="Use resolved-graph label and framing defaults where enabled by config.")
    ap.add_argument("--position-fill-x", type=float, default=0.0,
                    help="Target horizontal fill fraction for final node extent; <=0 disables.")
    ap.add_argument("--position-fill-y", type=float, default=0.0,
                    help="Target vertical fill fraction for final node extent; <=0 disables.")
    ap.add_argument("--axis-margin-frac", type=float, default=None,
                    help="Final axis margin fraction; default keeps graph-type behaviour.")

    ap.add_argument("--nodes", default=None,
                    help="Optional TSV listing all expected nodes so isolates appear.")
    ap.add_argument("--nodes-col", default="short_id",
                    help="Column in --nodes TSV to use as node IDs.")
    ap.add_argument("--font-size",  type=int,   default=10)
    ap.add_argument("--node-size", type=float, default=gps.NODE_SIZE)

    ap.add_argument("--edge-width-min", type=float, default=1.2)
    ap.add_argument("--edge-width-max", type=float, default=5.0)
    ap.add_argument("--edge-alpha",     type=float, default=0.85)

    ap.add_argument("--component-label-gap", type=float, default=0.40,
                    help="Data-coordinate gap between component labels and component tops.")
    ap.add_argument("--nbl-component-gap-fraction", type=float,
                    default=NBL_SUPPORT_COMPONENT_GAP_FRACTION,
                    help="NBL-only fraction of the original C1/C2 horizontal gap to retain.")
    ap.add_argument("--nbl-component-min-gap", type=float,
                    default=NBL_SUPPORT_COMPONENT_MIN_GAP,
                    help="NBL-only minimum retained C1/C2 horizontal gap in layout units.")
    ap.add_argument("--legend-anchor-x", type=float, default=None,
                    help="Optional legend x anchor in axes coordinates.")
    ap.add_argument("--legend-anchor-y", type=float, default=None,
                    help="Optional legend y anchor in axes coordinates.")
    ap.add_argument("--legend-align-mode",
                    choices=("default", "component-band"),
                    default=None,
                    help="Legend vertical alignment mode; NBL defaults to component-band.")
    ap.add_argument("--legend-gap-frac", type=float, default=None,
                    help="Rendered legend gap from the component box as a fraction of figure width.")
    ap.add_argument("--legend-mode",
                    choices=("keep", "none", "separate"),
                    default="keep")
    ap.add_argument("--legend-out-prefix", default=None)
    ap.add_argument("--isolate-label-y-frac", type=float, default=None,
                    help="Isolate-panel label y position as a fraction of isolate-box height.")
    ap.add_argument("--isolate-box-height-factor", type=float, default=2.00,
                    help="Multiplier for isolate-box vertical height.")
    ap.add_argument("--isolate-gap", type=float, default=gps.ADAPTIVE_ISOLATE_GAP,
                    help="Vertical gap between component grid and isolate panel.")
    ap.add_argument("--isolate-label-left-pad", type=float, default=0.0,
                    help="Left padding added to isolate box for its label.")
    ap.add_argument("--isolate-spacing", type=float, default=gps.ADAPTIVE_ISOLATE_SPACING,
                    help="Horizontal spacing between isolate diamonds.")
    ap.add_argument("--compact-isolate-panel", action="store_true",
                    help="Size the isolate strip from isolate content rather than full layout width.")
    ap.add_argument("--isolate-panel-width-mode",
                    choices=("layout", "content"), default="layout")
    ap.add_argument("--isolate-panel-spacing", type=float, default=None,
                    help="Horizontal spacing for compact isolate-panel placement.")
    ap.add_argument("--min-isolate-label-spacing", type=float, default=0.0,
                    help="Minimum data-coordinate spacing between isolate labels in compact panels.")
    ap.add_argument("--isolate-panel-padding-x", type=float, default=12.0)
    ap.add_argument("--isolate-panel-padding-y", type=float, default=1.4)
    ap.add_argument("--isolate-label-band-frac", type=float, default=0.34)
    ap.add_argument("--isolate-max-per-row", type=int, default=6)
    ap.add_argument("--isolate-panel-width-frac", type=float, default=0.62)
    ap.add_argument("--dense-component-min-size", type=int, default=8,
                    help="Minimum component size for post-layout dense-component expansion.")
    ap.add_argument("--dense-component-min-density", type=float, default=0.30,
                    help="Minimum component edge density for post-layout expansion.")
    ap.add_argument("--dense-component-expand-factor", type=float, default=1.85,
                    help="Centroid expansion factor for dense components after layout.")
    ap.add_argument("--adaptive-component-layout", action="store_true",
                    help="Expand each non-isolate component after layout using structural burden metrics.")
    ap.add_argument("--adaptive-component-min-nodes", type=int, default=4,
                    help="Minimum component node count eligible for adaptive expansion.")
    ap.add_argument("--adaptive-component-min-density", type=float, default=0.15,
                    help="Minimum component edge density eligible for adaptive expansion.")
    ap.add_argument("--adaptive-component-min-expand", type=float, default=1.0,
                    help="Lower bound for qualifying adaptive component expansion factors.")
    ap.add_argument("--adaptive-component-max-expand", type=float, default=4.0,
                    help="Upper bound for adaptive component expansion factors.")
    ap.add_argument("--adaptive-component-density-weight", type=float, default=1.8,
                    help="Weight applied to density above the adaptive density threshold.")
    ap.add_argument("--adaptive-component-size-weight", type=float, default=0.14,
                    help="Weight applied per node above the adaptive node threshold.")
    ap.add_argument("--adaptive-component-label-weight", type=float, default=0.30,
                    help="Weight applied to max label length above eight characters.")
    ap.add_argument("--label-overlap-avoidance", action="store_true",
                    help="Run label de-collision and report before/after overlap estimates.")
    ap.add_argument("--label-padding-factor", type=float, default=1.35,
                    help="Estimated label bbox expansion factor used for label collision avoidance.")
    ap.add_argument("--label-bbox-margin", type=float, default=0.0,
                    help="Extra data-coordinate label bbox margin for collision estimates.")
    ap.add_argument("--max-label-overlap-iterations", type=int, default=1200,
                    help="Maximum label de-collision iterations.")
    ap.add_argument("--min-label-separation", type=float, default=0.0,
                    help="Minimum data-coordinate label separation in overlap estimates.")
    ap.add_argument("--avoid-isolate-overlap", action="store_true",
                    help="Shift isolate panels down when graph components are too close.")
    ap.add_argument("--min-component-isolate-gap", type=float, default=4.0,
                    help="Minimum data-coordinate vertical gap from components to isolate panel.")
    ap.add_argument("--isolate-panel-reserved-height", type=float, default=0.0,
                    help="Reserved isolate-panel height for reports/config traceability.")
    ap.add_argument("--adaptive-isolate-y-shift", action="store_true",
                    help="Allow automatic downward isolate-panel shifts when enforcing gaps.")
    ap.add_argument("--component-isolate-gap-factor", type=float, default=1.0,
                    help="Multiplier applied to the configured component-isolate gap.")
    ap.add_argument("--component-grid-ncols", type=int, default=3,
                    help="BRCA fixed-grid component columns for row-major packing.")
    ap.add_argument("--component-grid-nrows", type=int, default=2,
                    help="BRCA fixed-grid component rows for row-major packing.")
    ap.add_argument("--component-cell-width", type=float, default=54.0,
                    help="BRCA fixed-grid component cell width in data units.")
    ap.add_argument("--component-cell-height", type=float, default=27.0,
                    help="BRCA fixed-grid component cell height in data units.")
    ap.add_argument("--component-row-gap", type=float, default=14.0)
    ap.add_argument("--component-col-gap", type=float, default=12.0)
    ap.add_argument("--legend-width", type=float, default=54.0)
    ap.add_argument("--isolate-region-height", type=float, default=10.5)
    ap.add_argument("--footnote-region-height", type=float, default=7.0)
    ap.add_argument("--component-cell-fill-x", type=float, default=0.78)
    ap.add_argument("--component-cell-fill-y", type=float, default=0.72)
    ap.add_argument("--dense-component-cell-fill-x", type=float, default=0.90)
    ap.add_argument("--dense-component-cell-fill-y", type=float, default=0.82)
    ap.add_argument("--c2-component-cell-fill-x", type=float, default=0.86)
    ap.add_argument("--c2-component-cell-fill-y", type=float, default=0.80)
    ap.add_argument("--c1-x-expand-factor", type=float, default=4.00)
    ap.add_argument("--c1-y-expand-factor", type=float, default=4.00)
    ap.add_argument("--c2-x-expand-factor", type=float, default=2.00)
    ap.add_argument("--c2-y-expand-factor", type=float, default=2.00)

    ap.add_argument("--seed",  type=int,   default=gps.LAYOUT_SEED)

    args = ap.parse_args()

    label_upper = args.label.upper()
    is_nbl  = label_upper == "NBL"
    is_brca = label_upper == "BRCA"
    is_rbl  = label_upper == "RBL"
    is_multicohort_majority = (
        label_upper == "MULTICOHORT_CANCER"
        and args.graph_mode == "majority_threshold"
    )
    is_multicohort_union = (
        label_upper == "MULTICOHORT_CANCER"
        and args.graph_mode == "union_supported_edges"
    )
    is_multicohort_multirep_planned = (
        is_multicohort_majority or is_multicohort_union
    )
    node_label_mode = args.node_label_mode
    if node_label_mode == "auto":
        node_label_mode = "selected" if args.graph_mode == "union_supported_edges" else "full"
    if node_label_mode == "selected" and args.graph_mode != "union_supported_edges":
        print(
            "[QC] Selected-label mode is only specialised for union-supported graphs; "
            "using full labels."
        )
        node_label_mode = "full"
    is_multicohort_union_full_labels = (
        is_multicohort_union and node_label_mode == "full"
    )
    use_halo_labels = bool(
        args.use_halo_labels
        or args.match_resolved_style
        or is_multicohort_union
        or is_multicohort_majority
    )
    graph_title = args.graph_title or (
        "Majority-threshold consensus network"
        if args.graph_mode == "majority_threshold"
        else "Union-of-supported-edges network"
    )
    context_label = f"{args.label} {graph_title}"
    if is_multicohort_union:
        context_label = "Multicohort union-supported cell-line similarity network"
        if is_multicohort_union_full_labels:
            context_label += " (full node labels)"
    if is_multicohort_majority:
        args.fig_w = 21.5
        args.fig_h = 14.8
        args.font_size = max(args.font_size, 12)
        print(
            "[QC] MULTICOHORT_CANCER majority-threshold compact framing: "
            f"fig_w={args.fig_w:.1f} fig_h={args.fig_h:.1f} "
            f"font_size={args.font_size}"
        )
    if is_multicohort_union:
        if is_multicohort_union_full_labels:
            args.fig_w = 32.0
            args.fig_h = 23.5
            args.font_size = 9
        else:
            args.fig_w = 22.0
            args.fig_h = 15.2
            args.font_size = 9
        print(
            "[QC] MULTICOHORT_CANCER union-supported framing: "
            f"label_mode={node_label_mode} "
            f"fig_w={args.fig_w:.1f} fig_h={args.fig_h:.1f} "
            f"font_size={args.font_size}"
        )

    # ---------- Load edges + node universe ----------
    df = pd.read_csv(args.edges, sep="\t")
    for col in ("node1", "node2"):
        if col not in df.columns:
            raise SystemExit(
                f"Missing column '{col}' in {args.edges}. "
                f"Columns present: {list(df.columns)}"
            )
    has_support = "support_directions" in df.columns

    G = nx.Graph()
    if args.nodes:
        nodes_df = pd.read_csv(args.nodes, sep="\t")
        if args.nodes_col not in nodes_df.columns:
            raise SystemExit(
                f"--nodes-col '{args.nodes_col}' not found in {args.nodes}. "
                f"Available columns: {list(nodes_df.columns)}"
            )
        for n in nodes_df[args.nodes_col].dropna().astype(str).str.strip().unique():
            if n:
                G.add_node(n)

    edge_support: dict[tuple[str, str], float] = {}
    for _, r in df.iterrows():
        u = str(r["node1"]).strip()
        v = str(r["node2"]).strip()
        if not u or not v or u == v:
            continue
        G.add_edge(u, v)
        if has_support:
            try:
                edge_support[tuple(sorted((u, v)))] = float(r["support_directions"])
            except (TypeError, ValueError):
                pass

    isolates = sorted(n for n in G.nodes() if G.degree(n) == 0)
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    n_comps_total = nx.number_connected_components(G)
    print(f"[INFO] Graph: {n_nodes} nodes, {n_edges} edges, "
          f"{n_comps_total} components ({len(isolates)} isolates)")

    # ---------- Layout: adaptive component shelves via shared module ----------
    G_non_iso = G.subgraph([n for n in G.nodes() if n not in set(isolates)]).copy()
    if is_rbl:
        pos, iso_bbox = gps.rbl_one_component_layout(
            G_non_iso, isolates, seed=args.seed,
            isolate_gap=RBL_SUPPORT_RESOLVED_REFERENCE_ISOLATE_GAP,
            isolate_spacing=args.isolate_spacing,
            isolate_box_height_factor=args.isolate_box_height_factor,
        )
        planner_regions = {}
        if args.adaptive_component_layout:
            gps.expand_adaptive_components_after_layout(
                G,
                pos,
                isolates,
                min_nodes=args.adaptive_component_min_nodes,
                min_density=args.adaptive_component_min_density,
                min_expand=args.adaptive_component_min_expand,
                max_expand=args.adaptive_component_max_expand,
                density_weight=args.adaptive_component_density_weight,
                size_weight=args.adaptive_component_size_weight,
                label_weight=args.adaptive_component_label_weight,
                context=f"{args.label} support-threshold consensus",
            )
    else:
        if is_multicohort_union:
            pos, iso_bbox, _grid_placements, planner_regions = (
                multicohort_union_supported_adaptive_layout(
                    G_non_iso,
                    isolates,
                    seed=args.seed,
                    spacing_scale=1.42 if is_multicohort_union_full_labels else 1.0,
                    context=f"{args.label} union-supported consensus",
                )
            )
            print(
                f"[QC] {args.label} union-supported consensus adaptive planner: "
                "large_component_region plus compact right column; "
                "component_region_size=nodes+edges"
            )
        elif is_multicohort_majority:
            pos, iso_bbox, _grid_placements, planner_regions = (
                multicohort_majority_adaptive_grid_layout(
                    G_non_iso,
                    isolates,
                    seed=args.seed,
                    context=f"{args.label} majority-threshold consensus",
                )
            )
            print(
                f"[QC] {args.label} majority-threshold consensus adaptive planner: "
                f"grid={MULTICOHORT_MAJORITY_GRID_N_COLS}x"
                f"{MULTICOHORT_MAJORITY_GRID_N_ROWS} "
                "component_cell_size=nodes+edges"
            )
        elif is_brca:
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
                context=f"{args.label} support-threshold consensus",
            )
            print(
                f"[QC] {args.label} support-threshold consensus deterministic planner: "
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
                    f"[QC] {args.label} support-threshold consensus "
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
                f"[QC] {args.label} support-threshold consensus isolate label left pad: "
                f"{args.isolate_label_left_pad:.2f}"
            )
        if args.adaptive_component_layout:
            gps.expand_adaptive_components_after_layout(
                G,
                pos,
                isolates,
                min_nodes=args.adaptive_component_min_nodes,
                min_density=args.adaptive_component_min_density,
                min_expand=args.adaptive_component_min_expand,
                max_expand=args.adaptive_component_max_expand,
                density_weight=args.adaptive_component_density_weight,
                size_weight=args.adaptive_component_size_weight,
                label_weight=args.adaptive_component_label_weight,
                context=f"{args.label} support-threshold consensus",
            )
        elif not (is_brca or is_multicohort_union or is_multicohort_majority):
            gps.expand_dense_components_after_layout(
                G, pos, isolates,
                min_size=args.dense_component_min_size,
                min_density=args.dense_component_min_density,
                expand_factor=args.dense_component_expand_factor,
                context=f"{args.label} support-threshold consensus",
            )
        gps.separate_focus_component_overlaps(
            G, pos, isolates,
            context=f"{args.label} support-threshold consensus",
            focus_component_index=1,
            min_gap=3.0,
        )
        gps.separate_component_overlaps(
            G, pos, isolates,
            context=f"{args.label} support-threshold consensus",
            min_gap=3.0,
        )
        gps.report_component_clearance(
            G, pos, isolates,
            context=f"{args.label} support-threshold consensus",
            focus_component_index=1,
            extra_bboxes={"isolate_box": iso_bbox},
        )
        gps.report_component_clearance(
            G, pos, isolates,
            context=f"{args.label} support-threshold consensus",
            focus_component_index=2,
            extra_bboxes={"isolate_box": iso_bbox},
        )
    if is_nbl:
        compact_nbl_support_components(
            G, pos, isolates,
            gap_fraction=args.nbl_component_gap_fraction,
            min_gap=args.nbl_component_min_gap,
        )
        if not args.compact_isolate_panel:
            iso_bbox = widen_nbl_support_target_isolates(G, pos, iso_bbox)
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
            label_map={n: n for n in G.nodes()},
            min_label_spacing=args.min_isolate_label_spacing,
            label_padding_factor=args.label_padding_factor,
            context=f"{args.label} support-threshold consensus",
        )
        print(
            f"[QC] {args.label} support-threshold consensus compact isolate panel: "
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
            labels_by_node={n: n for n in G.nodes()},
            shift_isolate_panel=args.adaptive_isolate_y_shift,
            context=f"{args.label} support-threshold consensus",
        )
    print(
        f"[QC] {args.label} support-threshold consensus position fill requested: "
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
            f"[QC] {args.label} support-threshold consensus position fill rescale: "
            f"scale_x={fill_transform['scale_x']:.3f} "
            f"scale_y={fill_transform['scale_y']:.3f} "
            f"current_aspect={fill_transform['current_aspect']:.3f} "
            f"target_aspect={fill_transform['target_aspect']:.3f}"
        )
    else:
        print(
            f"[QC] {args.label} support-threshold consensus position fill rescale: "
            "scale_x=1.000 scale_y=1.000 effective_fill_unchanged=yes"
        )
    if is_rbl:
        gps.report_isolate_layout_metrics(pos, isolates, iso_bbox, "RBL support-threshold consensus")
    if iso_bbox is not None:
        iso_height = float(iso_bbox[3] - iso_bbox[1])
        print(
            f"[QC] {args.label} support-threshold consensus isolate box height: "
            f"old_unscaled_height=4.40 factor={args.isolate_box_height_factor:.2f} "
            f"height={iso_height:.2f} "
            f"label_y_frac={(args.isolate_label_y_frac if args.isolate_label_y_frac is not None else (NBL_SUPPORT_ISOLATE_LABEL_Y_FRAC if is_nbl else 0.010)):.3f}"
        )
        gps.report_isolate_layout_metrics(pos, isolates, iso_bbox, f"{args.label} support-threshold consensus")

    # ---------- Colour map ----------
    if is_brca:
        component_palette = gps.BRCA_NO_GREEN_COMPONENT_PALETTE
    elif is_nbl:
        component_palette = gps.NBL_NO_GREEN_COMPONENT_PALETTE
    elif is_multicohort_union:
        component_palette = MULTICOHORT_UNION_COMPONENT_PALETTE
    else:
        component_palette = gps.OKABE_ITO_PALETTE
    node_colour, _, _ = gps.component_colour_map(
        G, isolates=isolates, palette=component_palette
    )

    # ---------- Draw ----------
    fig, ax = plt.subplots(figsize=(args.fig_w, args.fig_h), dpi=args.dpi)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    if iso_bbox is not None:
        isolate_label_y_frac = (
            args.isolate_label_y_frac
            if args.isolate_label_y_frac is not None
            else (NBL_SUPPORT_ISOLATE_LABEL_Y_FRAC if is_nbl else 0.010)
        )
        if is_brca:
            gps.draw_isolate_panel(
                ax,
                iso_bbox,
                isolate_nodes=isolates,
                pos=pos,
                labels_by_node={n: n for n in isolates},
                panel_padding_x=args.isolate_panel_padding_x,
                panel_padding_y=args.isolate_panel_padding_y,
                label_band_frac=args.isolate_label_band_frac,
                node_label_fontsize=args.font_size,
            )
        else:
            gps.draw_isolate_zone(
                ax,
                iso_bbox,
                label_y_frac=isolate_label_y_frac,
            )
        if is_nbl:
            print(
                "[QC] NBL support-threshold isolate label y fraction: "
                f"{float(isolate_label_y_frac):.3f}"
            )

    support_values_for_legend: list[float] = []
    if G.number_of_edges():
        edges_ordered = sorted(G.edges())
        if has_support:
            support_per_edge = [
                edge_support.get(tuple(sorted(e)), 0.0) for e in edges_ordered
            ]
            support_values_for_legend = support_per_edge
            widths = gps.edge_widths_from_support(
                support_per_edge,
                min_width=args.edge_width_min,
                max_width=args.edge_width_max,
            )
        else:
            widths = [args.edge_width_min] * len(edges_ordered)
        edge_style_by_key: dict[tuple[str, str], str] = {}
        if "edge_style" in df.columns:
            for _, r in df.iterrows():
                u = str(r["node1"]).strip()
                v = str(r["node2"]).strip()
                if u and v:
                    edge_style_by_key[tuple(sorted((u, v)))] = str(r["edge_style"]).strip()
        edge_index = {tuple(sorted(e)): i for i, e in enumerate(edges_ordered)}
        for style_name, mpl_style in (("solid", "solid"), ("dashed", (0, (5, 3)))):
            styled_edges = [
                edge for edge in edges_ordered
                if edge_style_by_key.get(tuple(sorted(edge)), "solid") == style_name
            ]
            if not styled_edges:
                continue
            styled_widths = [
                widths[edge_index[tuple(sorted(edge))]] for edge in styled_edges
            ]
            nx.draw_networkx_edges(
                G, pos, ax=ax, edgelist=styled_edges,
                width=styled_widths, alpha=args.edge_alpha,
                edge_color="black", style=mpl_style,
            )

    non_iso_nodes = [n for n in G.nodes() if n not in set(isolates)]
    if non_iso_nodes:
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            nodelist=non_iso_nodes,
            node_size=args.node_size,
            node_color=[node_colour[n] for n in non_iso_nodes],
            edgecolors=gps.NODE_EDGE_COLOUR,
            linewidths=1.4, alpha=0.95,
        )

    if isolates:
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            nodelist=isolates,
            node_size=int(args.node_size * 0.85),
            node_color=gps.ISOLATE_COLOUR,
            node_shape="D",
            edgecolors=gps.NODE_EDGE_COLOUR,
            linewidths=1.4, alpha=0.95,
        )

    if is_multicohort_union:
        gps.draw_component_labels(
            ax, G, pos, isolates=isolates, min_size=2,
            min_y_offset=2.0,
            y_offset_scale=0.025,
            same_row_top_n=0,
            n_cols_for_rows=0,
            use_bbox=False,
            font_color="black" if use_halo_labels else "#1A1A2E",
            halo_width=max(1.4, float(args.halo_linewidth) * 0.65),
        )
    elif not (is_brca or is_multicohort_multirep_planned):
        gps.draw_component_labels(
            ax, G, pos, isolates=isolates, min_size=2,
            min_y_offset=args.component_label_gap,
            same_row_top_n=2 if is_nbl else 0,
            n_cols_for_rows=0,
            use_bbox=False,
            font_color="black" if use_halo_labels else "#1A1A2E",
            halo_width=max(1.4, float(args.halo_linewidth) * 0.65),
        )
    if args.graph_mode == "union_supported_edges" and node_label_mode == "selected":
        label_nodes = union_thesis_label_nodes(
            G,
            isolates,
            pos=pos,
            max_component_labels=args.selected_label_max_component_labels,
            min_component_nodes=args.selected_label_min_component_nodes,
            min_component_density=args.selected_label_min_component_density,
            min_label_distance=args.selected_label_min_distance,
            context_label=args.label,
        )
    else:
        label_nodes = list(G.nodes())
    label_nodes = sorted(label_nodes)
    all_nodes_labelled = set(label_nodes) == set(G.nodes())
    print(
        f"[QC] {args.label} {args.graph_mode} node label mode: "
        f"mode={node_label_mode} shown={len(label_nodes)}/{G.number_of_nodes()} "
        f"all_nodes_labelled={'yes' if all_nodes_labelled else 'no'}"
    )
    text_items = gps.halo_labels(
        ax, G, pos, {n: n for n in label_nodes}, font_size=args.font_size,
        font_color="black" if use_halo_labels else "#1A1A2E",
        halo_width=float(args.halo_linewidth) if use_halo_labels else 0.0,
        halo_color=args.halo_color,
    )
    decollide_iters = int(args.max_label_overlap_iterations)
    decollide_step = max(0.25, float(args.min_label_separation))
    decollide_expand = float(args.label_padding_factor)
    decollide_char_width = 8.5
    decollide_char_height = 11.0
    if is_multicohort_union_full_labels:
        decollide_iters = max(decollide_iters, 5200)
        decollide_step = max(decollide_step, 0.72)
        decollide_expand = max(decollide_expand, 1.95)
        decollide_char_width = 10.8
        decollide_char_height = 13.2
    elif is_multicohort_union:
        decollide_iters = max(decollide_iters, 1800)
        decollide_step = max(decollide_step, 0.50)
        decollide_expand = max(decollide_expand, 1.55)
        decollide_char_width = 9.4
        decollide_char_height = 11.8
    elif is_multicohort_majority:
        decollide_iters = max(decollide_iters, 1800)
        decollide_step = max(decollide_step, 0.46)
        decollide_expand = max(decollide_expand, 1.62)
        decollide_char_width = 10.2
        decollide_char_height = 13.0
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
            halo_width=float(args.halo_linewidth) if use_halo_labels else 0.0,
            halo_color=args.halo_color,
            context=f"{args.label} {args.graph_mode} consensus",
        )
    else:
        gps.decollide_labels(
            ax,
            text_items,
            halo_width=float(args.halo_linewidth) if use_halo_labels else 0.0,
            halo_color=args.halo_color,
            context=f"{args.label} {args.graph_mode} consensus",
        )
    for text in text_items.values():
        if not use_halo_labels:
            text.set_path_effects([])

    ax.set_title(
        context_label,
        fontsize=args.font_size + 5, pad=20, color="#1A1A2E",
    )
    caption = (
        "node fill colour denotes connected component only. "
        "C labels are ordered by component size within this graph; "
        "component colours should not be compared across cohorts or graph types. "
        "Edges aggregate support across feature--distance representations."
    )
    if is_multicohort_union:
        if is_multicohort_union_full_labels:
            caption = (
                "node colour denotes connected component within this graph only; "
                "component colours should not be compared across cohorts or graph types. "
                "All node labels are shown in this full-label supplementary version; "
                "all nodes and edges are drawn. Edge width and line style denote representation support."
            )
        else:
            caption = (
                "node colour denotes connected component within this graph only; "
                "component colours should not be compared across cohorts or graph types. "
                "Selected node labels are shown for readability; all node labels are available in "
                "the full-label supplementary version and sidecar table. "
                "All nodes and edges are drawn. Edge width and line style denote representation support."
            )
    if is_rbl:
        caption = (
            "node fill colour denotes connected component only. "
            "C labels are ordered by component size within this graph;\n"
            "component colours should not be compared across cohorts or graph types. "
            "Edges aggregate support across feature--distance representations."
        )
    gps.add_caption(fig, caption, fontsize=8)

    extra_bboxes = []
    if iso_bbox is not None:
        extra_bboxes.append(iso_bbox)
    if (is_brca or is_multicohort_multirep_planned) and planner_regions.get("legend_region") is not None:
        extra_bboxes.append(planner_regions["legend_region"])
    axis_margin = (
        args.axis_margin_frac
        if args.axis_margin_frac is not None
        else (0.03 if is_multicohort_union else (0.04 if is_rbl else 0.06))
    )
    gps.frame_axes(
        ax, pos, margin=axis_margin,
        extra_bboxes=extra_bboxes if extra_bboxes else None,
        fill_x=args.position_fill_x if args.position_fill_x > 0.0 else None,
        fill_y=args.position_fill_y if args.position_fill_y > 0.0 else None,
    )
    if is_brca or is_multicohort_majority:
        gps.place_component_labels_safely(
            ax,
            G,
            pos,
            isolates=isolates,
            labels_by_node={n: n for n in G.nodes()},
            label_gap=args.component_label_gap,
            context=f"{args.label} {args.graph_mode} consensus",
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
            ax, pos, margin=axis_margin,
            extra_bboxes=extra_bboxes if extra_bboxes else None,
            fill_x=args.position_fill_x if args.position_fill_x > 0.0 else None,
            fill_y=args.position_fill_y if args.position_fill_y > 0.0 else None,
        )
    legend_align_mode = args.legend_align_mode or (
        "default" if (is_rbl or is_brca or is_multicohort_multirep_planned) else "component-band"
    )
    default_legend_x = 0.81 if is_nbl else (0.72 if is_rbl else 1.01)
    default_legend_y = 1.30 if is_brca else (0.50 if is_multicohort_multirep_planned else (0.90 if is_rbl else 1.0))
    legend_x = args.legend_anchor_x if args.legend_anchor_x is not None else default_legend_x
    if args.legend_anchor_y is not None:
        legend_y = args.legend_anchor_y
    elif legend_align_mode == "component-band":
        legend_y, band_ymin, band_ymax = component_band_mid_y_axes(ax, pos, non_iso_nodes)
        print(
            f"[QC] Support-threshold legend y-anchor derived from component band: "
            f"data_y_range=({band_ymin:.2f},{band_ymax:.2f}) "
            f"anchor=({legend_x:.3f},{legend_y:.3f})"
        )
    else:
        legend_y = default_legend_y
    legend_bbox_to_anchor = (legend_x, legend_y)
    legend_loc = "center left" if legend_align_mode == "component-band" else "upper left"
    legend_bbox_transform = None
    legend_coordinate_system = "axes"
    if (is_brca or is_multicohort_multirep_planned) and planner_regions.get("legend_region") is not None:
        requested_anchor = (float(legend_x), float(legend_y))
        requested_clips = (
            requested_anchor[0] < 0.0 or requested_anchor[0] > 1.0
            or requested_anchor[1] < 0.0 or requested_anchor[1] > 1.0
        )
        if requested_clips:
            print(
                f"[QC] {args.label} {args.graph_mode} consensus requested legend anchor "
                f"would clip outside axes: requested_axes_anchor="
                f"({requested_anchor[0]:.3f},{requested_anchor[1]:.3f}); "
                "fallback=reserved_legend_region"
            )
        legend_bbox_to_anchor = gps.legend_anchor_from_region(planner_regions["legend_region"])
        legend_loc = "center left"
        legend_bbox_transform = ax.transData
        legend_coordinate_system = "data"
    if is_brca or is_multicohort_multirep_planned:
        print(
            f"[QC] {args.label} {args.graph_mode} consensus legend placement: "
            f"mode={legend_align_mode} loc={legend_loc} "
            f"anchor=({legend_bbox_to_anchor[0]:.2f},{legend_bbox_to_anchor[1]:.2f}) "
            f"coordinate_system={legend_coordinate_system} "
            f"reserved_region="
            f"{gps._fmt_bbox(planner_regions['legend_region']) if planner_regions.get('legend_region') is not None else 'NA'}"
        )
    legend_handles = gps.legend_handles_consensus(
        has_isolates=bool(isolates),
        support_values=support_values_for_legend,
        edge_width_min=args.edge_width_min,
        edge_width_max=args.edge_width_max,
    )
    if is_multicohort_union and legend_handles:
        for handle in legend_handles:
            if getattr(handle, "get_label", lambda: "")().startswith("Cell line;"):
                handle.set_label(
                    "Cell line; node colour denotes connected component within this graph only"
                )
                if hasattr(handle, "set_markerfacecolor"):
                    handle.set_markerfacecolor(MULTICOHORT_UNION_COMPONENT_PALETTE[0])
    has_dashed_edges = (
        "edge_style" in df.columns
        and any(str(v).strip() == "dashed" for v in df["edge_style"])
    )
    if has_dashed_edges:
        legend_handles.append(mlines.Line2D(
            [], [], color="black", linewidth=1.5, linestyle=(0, (5, 3)),
            label="Support = 1; one feature--distance representation",
        ))
        legend_handles.append(mlines.Line2D(
            [], [], color="black", linewidth=2.5, linestyle="solid",
            label="Support >= 2; multiple feature--distance representations",
        ))
    legend_title = "Encoding" if is_multicohort_union else f"{graph_title} legend"
    legend_labelspacing = (
        2.15 if is_multicohort_union_full_labels
        else (2.00 if is_multicohort_union else (1.95 if is_multicohort_multirep_planned else 1.5))
    )
    if args.legend_mode == "keep":
        gps.place_legend(
            ax,
            legend_handles,
            title=legend_title,
            spacious=True,
            labelspacing_override=legend_labelspacing,
            bbox_to_anchor=legend_bbox_to_anchor,
            loc=legend_loc,
            bbox_transform=legend_bbox_transform,
            coordinate_system=legend_coordinate_system,
        )
    elif args.legend_mode == "separate":
        legend_out_prefix = args.legend_out_prefix or f"{args.out}_legend"
        gps.save_legend_only(
            legend_handles,
            legend_out_prefix,
            title=legend_title,
            dpi=args.dpi,
            spacious=True,
            labelspacing_override=legend_labelspacing,
        )

    # ---------- Save ----------
    # tight_layout pre-commits the geometry so savefig does not re-derive
    # the bounding box (which is non-deterministic when the legend is outside
    # the axes and the tight-bbox calculation involves text extents).
    if args.legend_mode != "keep":
        fig.tight_layout(rect=[0.0, 0.04, 0.98, 0.97])
    else:
        if is_rbl:
            fig.tight_layout(rect=[0.02, 0.09, 0.98, 0.95])
            gps.align_legend_to_component_box(
                fig, ax, pos, non_iso_nodes, "RBL support-threshold consensus",
                component_texts=non_iso_nodes,
                iso_bbox=iso_bbox,
                target_gap_frac=gps.RBL_SINGLE_COMPONENT_LEGEND_GAP_FRAC,
            )
        else:
            if is_multicohort_union:
                fig.tight_layout(rect=[0.005, 0.050, 0.995, 0.955])
            elif is_multicohort_multirep_planned:
                fig.tight_layout(rect=[0.01, 0.055, 0.99, 0.95])
            else:
                fig.tight_layout(rect=[0.0, 0.06, 0.98, 0.95] if is_brca else [0.0, 0.04, 0.80, 0.97])
            if is_nbl and legend_align_mode == "component-band":
                gps.align_legend_to_component_box(
                    fig, ax, pos, non_iso_nodes, "NBL support-threshold consensus",
                    component_texts=[*non_iso_nodes, "C1", "C2"],
                    iso_bbox=iso_bbox,
                    target_gap_frac=(
                        args.legend_gap_frac
                        if args.legend_gap_frac is not None
                        else NBL_SUPPORT_LEGEND_GAP_FRAC
                    ),
                )
    gps.report_rendered_layout_clearance(
        fig, ax, G, pos, isolates,
        f"{args.label} {args.graph_mode} consensus",
        labels_by_node={n: n for n in G.nodes()},
        iso_bbox=iso_bbox,
    )
    if is_nbl:
        gps.report_c_label_legend_metrics(fig, ax, "NBL support-threshold consensus")
    save_kw = dict(facecolor="white", transparent=False)
    if is_rbl:
        save_kw.update(bbox_inches="tight", pad_inches=0.12)
    fig.savefig(f"{args.out}.png", dpi=args.dpi, **save_kw)
    fig.savefig(f"{args.out}.pdf", **save_kw)
    fig.savefig(f"{args.out}.svg", **save_kw)
    component_by_node: dict[str, str] = {}
    component_rows = []
    for idx, comp in enumerate(gps.ordered_non_isolate_components(G, isolates=isolates), start=1):
        component = f"C{idx}"
        comp_nodes = sorted(comp)
        subgraph = G.subgraph(comp_nodes)
        for node in comp_nodes:
            component_by_node[node] = component
        component_rows.append({
            "cohort": args.label,
            "graph_product": args.graph_mode,
            "component": component,
            "is_isolate_component": False,
            "node_count": len(comp_nodes),
            "edge_count": subgraph.number_of_edges(),
            "nodes": ";".join(comp_nodes),
        })
    for idx, node in enumerate(isolates, start=1):
        component = f"I{idx}"
        component_by_node[node] = component
        component_rows.append({
            "cohort": args.label,
            "graph_product": args.graph_mode,
            "component": component,
            "is_isolate_component": True,
            "node_count": 1,
            "edge_count": 0,
            "nodes": node,
        })
    if args.node_labels_out:
        node_labels_path = Path(args.node_labels_out)
        node_labels_path.parent.mkdir(parents=True, exist_ok=True)
        visible_label_nodes = set(label_nodes)
        label_rows = [
            {
                "cohort": args.label,
                "graph_product": args.graph_mode,
                "node_label_mode": node_label_mode,
                "node": node,
                "label": node,
                "shown_in_figure": node in visible_label_nodes,
                "degree": int(G.degree(node)),
                "component": component_by_node.get(node, ""),
                "is_isolate": G.degree(node) == 0,
            }
            for node in sorted(G.nodes())
        ]
        pd.DataFrame(label_rows).to_csv(node_labels_path, sep="\t", index=False)
        print(
            f"[OK] Saved node label sidecar: {node_labels_path} "
            f"shown={sum(row['shown_in_figure'] for row in label_rows)}/"
            f"{len(label_rows)} mode={node_label_mode}"
        )
    if args.node_stats_out:
        node_stats_path = Path(args.node_stats_out)
        node_stats_path.parent.mkdir(parents=True, exist_ok=True)
        node_rows = [
            {
                "cohort": args.label,
                "graph_product": args.graph_mode,
                "node": node,
                "degree": int(G.degree(node)),
                "component": component_by_node.get(node, ""),
                "is_isolate": G.degree(node) == 0,
            }
            for node in sorted(G.nodes())
        ]
        pd.DataFrame(node_rows).to_csv(node_stats_path, sep="\t", index=False)
        print(f"[OK] Saved node stats: {node_stats_path}")
    if args.components_out:
        components_path = Path(args.components_out)
        components_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(component_rows).to_csv(components_path, sep="\t", index=False)
        print(f"[OK] Saved components: {components_path}")
    plt.close(fig)
    print(f"[OK] Saved: {args.out}.png / .pdf / .svg")


if __name__ == "__main__":
    main()
