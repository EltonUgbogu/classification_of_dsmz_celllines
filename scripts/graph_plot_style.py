"""
graph_plot_style.py
-------------------
Shared visual language for the cell-line similarity / neighbourhood graph
figures (plot_consensus_graph.py and visualize_resolved_dsmz_graph.py).

The visual language is documented in the plotting brief, section 5. Anchors
referenced in comments below:

  5.3  Isolates drawn as a diamond (neutral grey) in a labelled "Isolates
       (degree 0)" strip, never scattered as strays along the frame.
  5.4  Node fill encodes connected component using a colour-blind-friendly
       qualitative palette with no green hue. Isolates share one neutral grey.
       Component colours are structural identifiers only.
  5.6  Bridge-like anchors are drawn as a non-fill ring. Most-connected nodes
       are drawn as non-fill square outlines so both annotations can be
       overlaid without changing node fill.
  5.7  Edge width is proportional to support count. Where solid/dashed is
       meaningful (resolved graph with possible single-supported edges),
       dashed = single-supported, solid = multi-supported.
  5.8  Labels are short names rendered with a thin white halo for legibility
       (matplotlib.patheffects.withStroke), identical style for every node.
  5.9  Legend documents every encoding actually present, placed outside the
       axes so it cannot occlude nodes; swatches match the rendered glyphs
       (rings are rings in the legend, diamonds are diamonds, etc).

All public helpers are deterministic given an input. Colour assignment is
stable across runs because components are sorted (largest first, ties by
lexicographically smallest node id) before being assigned palette indices.
"""

from __future__ import annotations

import math
from typing import Iterable, Mapping, Sequence

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.legend_handler import HandlerTuple
import networkx as nx
import numpy as np


OKABE_ITO_PALETTE: tuple[str, ...] = (
    "#0072B2",  # blue
    "#56B4E9",  # sky-blue
    "#CC79A7",  # reddish-purple
    "#E69F00",  # orange
    "#F0E442",  # yellow
    "#8DA0CB",  # periwinkle
    "#A6761D",  # brown
)
NBL_NO_GREEN_COMPONENT_PALETTE: tuple[str, ...] = (
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#CC79A7",  # reddish-purple
    "#8DA0CB",  # periwinkle
    "#A6761D",  # brown
    "#D55E00",  # vermillion
    "#9467BD",  # violet
    "#56B4E9",  # sky-blue
)
BRCA_NO_GREEN_COMPONENT_PALETTE: tuple[str, ...] = NBL_NO_GREEN_COMPONENT_PALETTE
ISOLATE_COLOUR: str = "#7F7F7F"
ANCHOR_RING_COLOUR: str = "#111111"
MOST_CONNECTED_COLOUR: str = "#B00020"
NODE_EDGE_COLOUR: str = "#1A1A2E"

# Shared layout constants — identical for every cohort and both graph types.
LAYOUT_SEED: int = 42
NODE_SIZE: int = 1400  # scatter area (pts²); constant node radius across all figures


def component_colour_map(
    G: nx.Graph,
    isolates: Iterable[str] | None = None,
    palette: tuple[str, ...] = OKABE_ITO_PALETTE,
    isolate_colour: str = ISOLATE_COLOUR,
) -> tuple[dict[str, str], dict[int, str], list[int]]:
    """
    Assign one palette colour per connected component (non-isolate).

    Sort order is largest component first, then lexicographic node id, so
    colour assignment is stable across runs. Isolates (passed explicitly to
    avoid duplicating the caller's degree-0 detection) share `isolate_colour`.

    Returns (per_node_colour, per_component_colour, component_indices_in_order).
    """
    isolate_set = set(isolates) if isolates is not None else {n for n in G.nodes() if G.degree(n) == 0}
    comps_non_iso = ordered_non_isolate_components(G, isolates=isolate_set)

    per_node: dict[str, str] = {}
    per_comp: dict[int, str] = {}
    comp_order: list[int] = []
    for i, comp in enumerate(comps_non_iso):
        colour = palette[i % len(palette)]
        per_comp[i] = colour
        comp_order.append(i)
        for node in comp:
            per_node[node] = colour

    for node in isolate_set:
        per_node[node] = isolate_colour

    return per_node, per_comp, comp_order


def ordered_non_isolate_components(
    G: nx.Graph,
    isolates: Iterable[str] | None = None,
) -> list[set[str]]:
    """Return connected components ordered by size, excluding degree-0 isolates."""
    isolate_set = set(isolates) if isolates is not None else {n for n in G.nodes() if G.degree(n) == 0}
    return [
        set(c) for c in sorted(
            nx.connected_components(G),
            key=lambda c: (-len(c), sorted(c)[0]),
        )
        if not (len(c) == 1 and next(iter(c)) in isolate_set)
    ]


def component_edge_density(G: nx.Graph, nodes: Sequence[str]) -> float:
    """Return simple undirected edge density for a component node list."""
    n = len(nodes)
    if n < 2:
        return 0.0
    sub = G.subgraph(nodes)
    return float(2 * sub.number_of_edges()) / float(n * (n - 1))


def component_layout_metrics(
    G: nx.Graph,
    isolates: Iterable[str] | None = None,
    label_map: Mapping[str, str] | None = None,
) -> list[dict[str, object]]:
    """
    Return per-component layout burden metrics for non-isolate components.

    Component labels (C1, C2, ...) are display ranks derived from the current
    graph structure. They are reported for QC only and are not used as stable
    biological identifiers.
    """
    metrics: list[dict[str, object]] = []
    for idx, comp in enumerate(ordered_non_isolate_components(G, isolates=isolates), start=1):
        nodes = sorted(str(n) for n in comp)
        n_nodes = len(nodes)
        sub = G.subgraph(nodes)
        n_edges = int(sub.number_of_edges())
        possible_edges = float(n_nodes * (n_nodes - 1) / 2.0)
        density = float(n_edges / possible_edges) if possible_edges > 0 else 0.0
        avg_degree = float(2.0 * n_edges / n_nodes) if n_nodes > 0 else 0.0

        labels: list[str] = []
        for node in nodes:
            raw_label = label_map.get(node, node) if label_map is not None else node
            labels.append(str(raw_label if raw_label is not None else node))
        label_lengths = [len(label) for label in labels]
        max_label_len = max(label_lengths) if label_lengths else 0
        mean_label_len = float(np.mean(label_lengths)) if label_lengths else 0.0

        metrics.append({
            "component": f"C{idx}",
            "component_index": idx,
            "nodes": nodes,
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            "possible_edges": possible_edges,
            "density": density,
            "avg_degree": avg_degree,
            "max_label_len": max_label_len,
            "mean_label_len": mean_label_len,
        })
    return metrics


def adaptive_component_expand_factors(
    G: nx.Graph,
    isolates: Iterable[str] | None = None,
    label_map: Mapping[str, str] | None = None,
    *,
    min_nodes: int = 4,
    min_density: float = 0.15,
    min_expand: float = 1.0,
    max_expand: float = 4.0,
    density_weight: float = 1.8,
    size_weight: float = 0.14,
    label_weight: float = 0.30,
    context: str = "graph",
) -> tuple[dict[int, tuple[float, float]], list[dict[str, object]]]:
    """
    Compute bounded per-component expansion factors from structural burden.

    Dense multi-node components with more nodes, higher density, and longer
    labels receive larger spacing. Degree-0 isolates are excluded before
    metrics are computed.
    """
    min_nodes = int(min_nodes)
    min_density = float(min_density)
    min_expand = max(1.0, float(min_expand))
    max_expand = max(min_expand, float(max_expand))
    density_weight = float(density_weight)
    size_weight = float(size_weight)
    label_weight = float(label_weight)

    component_factors: dict[int, tuple[float, float]] = {}
    rows = component_layout_metrics(G, isolates=isolates, label_map=label_map)
    for row in rows:
        component_index = int(row["component_index"])
        n_nodes = int(row["n_nodes"])
        density = float(row["density"])
        max_label_len = int(row["max_label_len"])

        qualifies = n_nodes >= min_nodes and density >= min_density
        size_term = 0.0
        density_term = 0.0
        label_term = 0.0
        expand_factor = 1.0
        if qualifies:
            size_term = max(0.0, float(n_nodes - min_nodes)) * size_weight
            density_term = max(0.0, density - min_density) * density_weight
            label_term = max(0.0, float(max_label_len - 8)) / 20.0 * label_weight
            expand_factor = min_expand + size_term + density_term + label_term
            expand_factor = max(min_expand, min(max_expand, expand_factor))
            if expand_factor > 1.0:
                component_factors[component_index] = (expand_factor, expand_factor)

        row["qualifies_for_adaptive_expansion"] = qualifies
        row["size_term"] = size_term
        row["density_term"] = density_term
        row["label_term"] = label_term
        row["expand_factor"] = expand_factor
        print(
            f"[QC] {context} adaptive component layout metrics: {row['component']} "
            f"n_nodes={n_nodes} n_edges={int(row['n_edges'])} "
            f"possible_edges={float(row['possible_edges']):.1f} "
            f"density={density:.3f} avg_degree={float(row['avg_degree']):.2f} "
            f"max_label_len={max_label_len} "
            f"mean_label_len={float(row['mean_label_len']):.1f} "
            f"size_term={size_term:.3f} density_term={density_term:.3f} "
            f"label_term={label_term:.3f} expand_factor={expand_factor:.2f} "
            f"expanded={'yes' if expand_factor > 1.0 else 'no'}"
        )
    return component_factors, rows


def expand_adaptive_components_after_layout(
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Iterable[str] | None = None,
    label_map: Mapping[str, str] | None = None,
    *,
    min_nodes: int = 4,
    min_density: float = 0.15,
    min_expand: float = 1.0,
    max_expand: float = 4.0,
    density_weight: float = 1.8,
    size_weight: float = 0.14,
    label_weight: float = 0.30,
    context: str = "graph",
) -> list[dict[str, object]]:
    """
    Expand qualifying components around their centroids after initial layout.

    Isolates are deliberately excluded from the component ranking and are not
    moved by this operation.
    """
    component_factors, rows = adaptive_component_expand_factors(
        G,
        isolates=isolates,
        label_map=label_map,
        min_nodes=min_nodes,
        min_density=min_density,
        min_expand=min_expand,
        max_expand=max_expand,
        density_weight=density_weight,
        size_weight=size_weight,
        label_weight=label_weight,
        context=context,
    )
    expand_components_after_layout(
        G,
        pos,
        isolates=isolates,
        component_factors=component_factors,
        context=f"{context} adaptive",
    )
    return rows


def component_layout_bboxes(
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Iterable[str] | None = None,
) -> dict[str, tuple[float, float, float, float]]:
    """Return data-coordinate bboxes for non-isolate components C1, C2, ..."""
    bboxes: dict[str, tuple[float, float, float, float]] = {}
    for idx, comp in enumerate(ordered_non_isolate_components(G, isolates=isolates), start=1):
        pts = np.array([pos[n] for n in comp if n in pos], dtype=float)
        if not pts.size:
            continue
        bboxes[f"C{idx}"] = (
            float(pts[:, 0].min()), float(pts[:, 1].min()),
            float(pts[:, 0].max()), float(pts[:, 1].max()),
        )
    return bboxes


def _fmt_bbox(bbox: tuple[float, float, float, float]) -> str:
    return "(" + ",".join(f"{float(v):.2f}" for v in bbox) + ")"


def expand_components_after_layout(
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Iterable[str] | None = None,
    *,
    component_factors: Mapping[int | str, tuple[float, float]],
    context: str = "graph",
) -> list[dict[str, float | int | str | tuple[float, float, float, float]]]:
    """
    Expand selected components anisotropically around their own centroid.

    The component centroid is preserved, so this changes visual spacing only;
    graph topology, component membership, edge endpoints by node identity, and
    centrality selections are untouched.
    """
    expanded: list[dict[str, float | int | str | tuple[float, float, float, float]]] = []
    components = ordered_non_isolate_components(G, isolates=isolates)
    for idx, comp in enumerate(components, start=1):
        factors = component_factors.get(idx, component_factors.get(f"C{idx}"))
        if factors is None:
            continue
        x_factor, y_factor = float(factors[0]), float(factors[1])
        if x_factor <= 1.0 and y_factor <= 1.0:
            continue

        nodes = sorted(n for n in comp if n in pos)
        if len(nodes) < 2:
            continue
        coords = np.vstack([pos[n] for n in nodes])
        centroid = coords.mean(axis=0)
        before = (
            float(coords[:, 0].min()), float(coords[:, 1].min()),
            float(coords[:, 0].max()), float(coords[:, 1].max()),
        )
        for node in nodes:
            delta = np.array(pos[node], dtype=float) - centroid
            pos[node] = np.array([
                float(centroid[0]) + x_factor * float(delta[0]),
                float(centroid[1]) + y_factor * float(delta[1]),
            ], dtype=float)
        after_pts = np.vstack([pos[n] for n in nodes])
        after = (
            float(after_pts[:, 0].min()), float(after_pts[:, 1].min()),
            float(after_pts[:, 0].max()), float(after_pts[:, 1].max()),
        )
        density = component_edge_density(G, nodes)
        expanded.append({
            "component": f"C{idx}",
            "component_index": idx,
            "n": len(nodes),
            "density": density,
            "x_factor": x_factor,
            "y_factor": y_factor,
            "bbox_before": before,
            "bbox_after": after,
        })
        print(
            f"[QC] {context} component anisotropic layout expansion: C{idx} "
            f"n={len(nodes)} density={density:.3f} "
            f"x_factor={x_factor:.2f} y_factor={y_factor:.2f} "
            f"bbox_before={_fmt_bbox(before)} bbox_after={_fmt_bbox(after)}"
        )
    return expanded


def expand_dense_components_after_layout(
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Iterable[str] | None = None,
    *,
    min_size: int = 8,
    min_density: float = 0.30,
    expand_factor: float = 1.0,
    x_expand_factor: float | None = None,
    y_expand_factor: float | None = None,
    context: str = "graph",
) -> list[dict[str, float | int | str]]:
    """
    Expand dense connected components around their own centroid after layout.

    This preserves each component's relative topology and centroid while adding
    internal spacing only to large, dense components where node labels and
    overlays are likely to be crowded.
    """
    x_factor = float(expand_factor if x_expand_factor is None else x_expand_factor)
    y_factor = float(expand_factor if y_expand_factor is None else y_expand_factor)
    if x_factor <= 1.0 and y_factor <= 1.0:
        return []

    expanded: list[dict[str, float | int | str]] = []
    components = ordered_non_isolate_components(G, isolates=isolates)
    for idx, comp in enumerate(components, start=1):
        nodes = sorted(n for n in comp if n in pos)
        if len(nodes) < min_size:
            continue
        density = component_edge_density(G, nodes)
        if density < min_density:
            continue
        coords = np.vstack([pos[n] for n in nodes])
        centroid = coords.mean(axis=0)
        before = (
            float(coords[:, 0].min()), float(coords[:, 1].min()),
            float(coords[:, 0].max()), float(coords[:, 1].max()),
        )
        for node in nodes:
            delta = np.array(pos[node], dtype=float) - centroid
            pos[node] = np.array([
                float(centroid[0]) + x_factor * float(delta[0]),
                float(centroid[1]) + y_factor * float(delta[1]),
            ], dtype=float)
        after_pts = np.vstack([pos[n] for n in nodes])
        after = (
            float(after_pts[:, 0].min()), float(after_pts[:, 1].min()),
            float(after_pts[:, 0].max()), float(after_pts[:, 1].max()),
        )
        expanded.append({
            "component": f"C{idx}",
            "component_index": idx,
            "n": len(nodes),
            "density": density,
            "factor": float(expand_factor),
            "x_factor": x_factor,
            "y_factor": y_factor,
            "bbox_before": before,
            "bbox_after": after,
        })
        print(
            f"[QC] {context} dense component layout expansion: C{idx} "
            f"n={len(nodes)} density={density:.3f} "
            f"x_factor={x_factor:.2f} y_factor={y_factor:.2f} "
            f"bbox_before={_fmt_bbox(before)} bbox_after={_fmt_bbox(after)}"
        )
    return expanded


def _bbox_clearance(
    a: tuple[float, float, float, float],
    b: tuple[float, float, float, float],
) -> tuple[float, bool]:
    """Return axis-aligned clearance between two bboxes and whether they overlap."""
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    x_gap = max(bx0 - ax1, ax0 - bx1, 0.0)
    y_gap = max(by0 - ay1, ay0 - by1, 0.0)
    overlap = x_gap == 0.0 and y_gap == 0.0
    if overlap:
        return 0.0, True
    if x_gap == 0.0:
        return float(y_gap), False
    if y_gap == 0.0:
        return float(x_gap), False
    return float(math.hypot(x_gap, y_gap)), False


def separate_focus_component_overlaps(
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Iterable[str] | None = None,
    *,
    context: str = "graph",
    focus_component_index: int = 1,
    min_gap: float = 3.0,
) -> list[dict[str, float | str]]:
    """
    Translate non-focus components away if their layout bbox overlaps C1.

    This is a layout safeguard after centroid expansion; it never changes graph
    topology or within-component node coordinates.
    """
    components = ordered_non_isolate_components(G, isolates=isolates)
    if focus_component_index < 1 or focus_component_index > len(components):
        return []

    moved: list[dict[str, float | str]] = []
    focus_nodes = [n for n in components[focus_component_index - 1] if n in pos]
    if not focus_nodes:
        return moved

    def _bbox(nodes: Sequence[str]) -> tuple[float, float, float, float] | None:
        pts = np.array([pos[n] for n in nodes if n in pos], dtype=float)
        if not pts.size:
            return None
        return (
            float(pts[:, 0].min()), float(pts[:, 1].min()),
            float(pts[:, 0].max()), float(pts[:, 1].max()),
        )

    focus_label = f"C{focus_component_index}"
    for idx, comp in enumerate(components, start=1):
        if idx == focus_component_index:
            continue
        focus_bbox = _bbox(focus_nodes)
        other_nodes = [n for n in comp if n in pos]
        other_bbox = _bbox(other_nodes)
        if focus_bbox is None or other_bbox is None:
            continue
        _, overlap = _bbox_clearance(focus_bbox, other_bbox)
        if not overlap:
            continue

        fx0, _, fx1, _ = focus_bbox
        ox0, _, ox1, _ = other_bbox
        focus_cx = (fx0 + fx1) / 2.0
        other_cx = (ox0 + ox1) / 2.0
        if other_cx >= focus_cx:
            dx = (fx1 - ox0) + float(min_gap)
        else:
            dx = -((ox1 - fx0) + float(min_gap))
        for node in other_nodes:
            pos[node] = np.array(
                [float(pos[node][0]) + dx, float(pos[node][1])],
                dtype=float,
            )
        moved.append({"component": f"C{idx}", "dx": float(dx), "dy": 0.0})
        print(
            f"[QC] {context} layout overlap safeguard: moved C{idx} "
            f"away from {focus_label} by dx={dx:.2f}"
        )
    return moved


def separate_component_overlaps(
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Iterable[str] | None = None,
    *,
    context: str = "graph",
    min_gap: float = 3.0,
    max_passes: int = 4,
) -> list[dict[str, float | str]]:
    """
    Translate later ordered components until component bboxes no longer overlap.

    This post-layout safeguard preserves graph topology and the relative
    coordinates inside each connected component. It only moves whole components
    when dense-component expansion or compact packing causes component bboxes to
    overlap after the initial layout.
    """
    components = ordered_non_isolate_components(G, isolates=isolates)
    moved: list[dict[str, float | str]] = []

    def _bbox(nodes: Sequence[str]) -> tuple[float, float, float, float] | None:
        pts = np.array([pos[n] for n in nodes if n in pos], dtype=float)
        if not pts.size:
            return None
        return (
            float(pts[:, 0].min()), float(pts[:, 1].min()),
            float(pts[:, 0].max()), float(pts[:, 1].max()),
        )

    for pass_idx in range(1, int(max_passes) + 1):
        moved_this_pass = False
        for i, comp_i in enumerate(components):
            nodes_i = [n for n in comp_i if n in pos]
            bbox_i = _bbox(nodes_i)
            if bbox_i is None:
                continue
            ix0, iy0, ix1, iy1 = bbox_i
            icx = (ix0 + ix1) / 2.0
            icy = (iy0 + iy1) / 2.0
            for j in range(i + 1, len(components)):
                comp_j = components[j]
                nodes_j = [n for n in comp_j if n in pos]
                bbox_j = _bbox(nodes_j)
                if bbox_j is None:
                    continue
                clearance, overlap = _bbox_clearance(bbox_i, bbox_j)
                if not overlap:
                    continue

                jx0, jy0, jx1, jy1 = bbox_j
                jcx = (jx0 + jx1) / 2.0
                jcy = (jy0 + jy1) / 2.0
                x_overlap = min(ix1, jx1) - max(ix0, jx0)
                y_overlap = min(iy1, jy1) - max(iy0, jy0)

                if x_overlap <= y_overlap:
                    direction = 1.0 if jcx >= icx else -1.0
                    if jcx == icx:
                        direction = 1.0
                    dx = direction * (float(x_overlap) + float(min_gap))
                    dy = 0.0
                else:
                    direction = 1.0 if jcy >= icy else -1.0
                    if jcy == icy:
                        direction = -1.0
                    dx = 0.0
                    dy = direction * (float(y_overlap) + float(min_gap))

                for node in nodes_j:
                    pos[node] = np.array(
                        [float(pos[node][0]) + dx, float(pos[node][1]) + dy],
                        dtype=float,
                    )
                moved_this_pass = True
                moved.append({
                    "component": f"C{j + 1}",
                    "against": f"C{i + 1}",
                    "pass": float(pass_idx),
                    "dx": float(dx),
                    "dy": float(dy),
                    "previous_clearance": float(clearance),
                })
                print(
                    f"[QC] {context} general layout overlap safeguard: "
                    f"moved C{j + 1} away from C{i + 1} on pass {pass_idx} "
                    f"by dx={dx:.2f} dy={dy:.2f}"
                )
        if not moved_this_pass:
            break
    return moved


def report_component_clearance(
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Iterable[str] | None = None,
    *,
    context: str = "graph",
    focus_component_index: int = 1,
    extra_bboxes: dict[str, tuple[float, float, float, float] | None] | None = None,
) -> dict[str, float | str | bool]:
    """
    Report layout bbox clearance from a focus component to other plot objects.

    The calculation is intentionally conservative and data-coordinate based; it
    catches component/isolate-box overlaps after post-layout transformations.
    """
    components = ordered_non_isolate_components(G, isolates=isolates)
    bboxes: dict[str, tuple[float, float, float, float]] = {}
    for idx, comp in enumerate(components, start=1):
        pts = np.array([pos[n] for n in comp if n in pos], dtype=float)
        if pts.size:
            bboxes[f"C{idx}"] = (
                float(pts[:, 0].min()), float(pts[:, 1].min()),
                float(pts[:, 0].max()), float(pts[:, 1].max()),
            )
    if extra_bboxes:
        for label, bbox in extra_bboxes.items():
            if bbox is not None:
                bboxes[label] = tuple(float(v) for v in bbox)

    focus = f"C{focus_component_index}"
    if focus not in bboxes:
        return {}
    focus_bbox = bboxes[focus]
    reports = []
    min_label = ""
    min_clearance = math.inf
    any_overlap = False
    for label, bbox in bboxes.items():
        if label == focus:
            continue
        clearance, overlap = _bbox_clearance(focus_bbox, bbox)
        any_overlap = any_overlap or overlap
        if clearance < min_clearance:
            min_clearance = clearance
            min_label = label
        reports.append(f"{label}={clearance:.2f}{'*' if overlap else ''}")

    if not reports:
        return {}
    print(
        f"[QC] {context} {focus} layout clearance: "
        f"min_target={min_label} min_clearance={min_clearance:.2f} "
        f"overlap={'yes' if any_overlap else 'no'}; "
        + ", ".join(reports)
    )
    return {
        "focus": focus,
        "min_target": min_label,
        "min_clearance": float(min_clearance),
        "overlap": bool(any_overlap),
    }


def component_label_by_node(
    G: nx.Graph,
    isolates: Iterable[str] | None = None,
) -> dict[str, str]:
    """Map each non-isolate node to C1, C2, etc. ordered by component size."""
    out: dict[str, str] = {}
    for idx, comp in enumerate(ordered_non_isolate_components(G, isolates=isolates), start=1):
        label = f"C{idx}"
        for node in comp:
            out[node] = label
    return out

def edge_widths_from_support(
    support_values: list[float],
    min_width: float = 0.9,
    max_width: float = 4.5,
) -> list[float]:
    """
    Map a per-edge support count to a per-edge line width.

    Width range is clipped to [min_width, max_width]. If all support values
    are equal (or the input is empty), all widths are min_width — flat support
    should not produce arbitrary contrast.
    """
    if not support_values:
        return []
    finite = [s for s in support_values if s is not None and np.isfinite(s)]
    if not finite:
        return [min_width] * len(support_values)
    lo, hi = min(finite), max(finite)
    if hi <= lo:
        return [min_width] * len(support_values)
    return [
        min_width + (max_width - min_width) * ((s - lo) / (hi - lo))
        if (s is not None and np.isfinite(s))
        else min_width
        for s in support_values
    ]


def representative_support_values(support_values: Sequence[float | int | None]) -> list[float]:
    """Pick min/middle/max observed support counts for a compact legend."""
    finite = sorted({
        float(s) for s in support_values
        if s is not None and np.isfinite(s)
    })
    if not finite:
        return []
    if len(finite) <= 3:
        return finite
    return [finite[0], finite[len(finite) // 2], finite[-1]]


def edge_width_for_support(
    value: float,
    support_values: Sequence[float | int | None],
    min_width: float,
    max_width: float,
) -> float:
    """Return the same width mapping used by edge_widths_from_support for one value."""
    finite = [
        float(s) for s in support_values
        if s is not None and np.isfinite(s)
    ]
    if not finite:
        return min_width
    lo, hi = min(finite), max(finite)
    if hi <= lo:
        return min_width
    return min_width + (max_width - min_width) * ((float(value) - lo) / (hi - lo))


def support_count_legend_handles(
    support_values: Sequence[float | int | None],
    min_width: float,
    max_width: float,
) -> list:
    """Build representative edge-width legend handles using observed support counts."""
    handles: list = []
    for value in representative_support_values(support_values):
        width = edge_width_for_support(value, support_values, min_width, max_width)
        label_value = int(value) if float(value).is_integer() else value
        handles.append(mlines.Line2D(
            [], [], color="black", linewidth=width,
            label=f"Support count = {label_value}",
        ))
    return handles


def halo_labels(
    ax,
    G,
    pos,
    labels: dict[str, str],
    font_size: int = 9,
    font_color: str = "#1A1A2E",
    halo_width: float = 3.3,
    halo_color: str = "white",
) -> dict:
    """
    Draw node labels with a thin white halo (path-effect stroke).

    Returns the Text-artist dict so the caller can pass it to decollide_labels.
    """
    text_items = nx.draw_networkx_labels(
        G, pos, labels=labels,
        font_size=font_size, font_weight="bold",
        font_color=font_color, ax=ax,
    )
    for t in text_items.values():
        t.set_path_effects([pe.withStroke(linewidth=halo_width, foreground=halo_color)])
    return text_items


def decollide_labels(
    ax,
    text_items: dict,
    *,
    iters: int = 300,
    step: float = 0.25,
    char_width_pts: float = 8.5,
    char_height_pts: float = 11.0,
    expand: float = 1.20,
    bbox_margin: float = 0.0,
    min_separation: float = 0.0,
    halo_width: float = 3.3,
    halo_color: str = "white",
    context: str | None = None,
) -> dict[str, float | int]:
    """
    Iteratively push overlapping node labels apart in data coordinates.

    Uses a purely data-coordinate approach (no renderer calls) so the result
    is fully deterministic across process invocations. Label widths/heights are
    estimated from character count and nominal font metrics; then the repel
    runs in data space. Path-effect halos are reapplied after all iterations.
    """
    texts = [t for t in text_items.values() if t.get_text()]
    if not texts:
        return {"labels": 0, "overlaps_before": 0, "overlaps_after": 0, "iterations": 0}

    # Convert nominal character dimensions from points to data coordinates.
    # ax.transData maps data→display; its inverse maps display→data.
    # We compute the scale factor at the centre of the axes.
    try:
        inv = ax.transData.inverted()
        cx, cy = ax.get_xlim()[0], ax.get_ylim()[0]
        p0 = ax.transData.transform((cx, cy))
        p1 = ax.transData.transform((cx + 1.0, cy + 1.0))
        pts_per_data_x = abs(p1[0] - p0[0])  # display pts per data unit (x)
        pts_per_data_y = abs(p1[1] - p0[1])  # display pts per data unit (y)
    except Exception:
        return {"labels": len(texts), "overlaps_before": 0, "overlaps_after": 0, "iterations": 0}

    if pts_per_data_x < 1e-9 or pts_per_data_y < 1e-9:
        return {"labels": len(texts), "overlaps_before": 0, "overlaps_after": 0, "iterations": 0}

    bbox_margin = max(0.0, float(bbox_margin))
    min_separation = max(0.0, float(min_separation))

    def _half_size(t):
        n = len(t.get_text())
        hw = (n * char_width_pts * expand * 0.5) / pts_per_data_x + bbox_margin
        hh = (char_height_pts * expand * 0.5) / pts_per_data_y + bbox_margin
        return hw, hh

    def _count_overlaps() -> int:
        overlaps = 0
        for i, ti in enumerate(texts):
            xi, yi = ti.get_position()
            hwi, hhi = _half_size(ti)
            for j in range(i + 1, len(texts)):
                tj = texts[j]
                xj, yj = tj.get_position()
                hwj, hhj = _half_size(tj)
                if (
                    abs(xi - xj) < hwi + hwj + min_separation
                    and abs(yi - yj) < hhi + hhj + min_separation
                ):
                    overlaps += 1
        return overlaps

    overlaps_before = _count_overlaps()
    iterations_used = 0
    for iter_idx in range(int(iters)):
        moved = False
        iterations_used = iter_idx + 1
        for i, ti in enumerate(texts):
            xi, yi = ti.get_position()
            hwi, hhi = _half_size(ti)
            for j in range(i + 1, len(texts)):
                tj = texts[j]
                xj, yj = tj.get_position()
                hwj, hhj = _half_size(tj)
                # AABB overlap
                if (
                    abs(xi - xj) >= hwi + hwj + min_separation
                    or abs(yi - yj) >= hhi + hhj + min_separation
                ):
                    continue
                dx, dy = xi - xj, yi - yj
                d = math.hypot(dx, dy)
                if d < 1e-9:
                    dx, dy, d = step, 0.0, step
                s = step / d
                ti.set_position((xi + dx * s, yi + dy * s))
                tj.set_position((xj - dx * s, yj - dy * s))
                xi, yi = ti.get_position()
                moved = True
        if not moved:
            break

    overlaps_after = _count_overlaps()
    for t in texts:
        t.set_path_effects([pe.withStroke(linewidth=float(halo_width), foreground=halo_color)])
    if context is not None:
        print(
            f"[QC] {context} label overlap avoidance: "
            f"labels={len(texts)} overlaps_before={overlaps_before} "
            f"overlaps_after={overlaps_after} iterations={iterations_used} "
            f"padding_factor={float(expand):.2f} bbox_margin={bbox_margin:.2f} "
            f"min_label_separation={min_separation:.2f}"
        )
    return {
        "labels": len(texts),
        "overlaps_before": overlaps_before,
        "overlaps_after": overlaps_after,
        "iterations": iterations_used,
    }


def draw_component_labels(
    ax,
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Iterable[str] | None = None,
    min_size: int = 2,
    min_y_offset: float = 0.35,
    y_offset_scale: float = 0.08,
    same_row_top_n: int = 0,
    n_cols_for_rows: int = 0,
    use_bbox: bool = False,
    font_color: str = "#1A1A2E",
    halo_width: float = 1.4,
) -> None:
    """
    Draw C1, C2, ... labels for non-isolate components ordered by size.

    When n_cols_for_rows > 0 (grid layout), all C-labels within the same grid
    row are placed at a common y (the maximum node-top in that row), ensuring
    row-mates share a common top edge.
    """
    components = ordered_non_isolate_components(G, isolates=isolates)

    if n_cols_for_rows > 0:
        # Compute the highest label-y per row across all components in that row.
        row_y: dict[int, float] = {}
        for idx, comp in enumerate(components):
            if len(comp) < min_size:
                continue
            pts = np.array([pos[n] for n in comp if n in pos], dtype=float)
            if not pts.size:
                continue
            row = idx // n_cols_for_rows
            span = np.ptp(pts, axis=0)
            candidate = float(pts[:, 1].max()) + max(min_y_offset, float(span[1]) * y_offset_scale)
            row_y[row] = max(row_y.get(row, -np.inf), candidate)

        for idx, comp in enumerate(components, start=1):
            if len(comp) < min_size:
                continue
            pts = np.array([pos[n] for n in comp if n in pos], dtype=float)
            if not pts.size:
                continue
            x = float(pts[:, 0].mean())
            row = (idx - 1) // n_cols_for_rows
            label_y = row_y.get(row, float(pts[:, 1].max()) + min_y_offset)
            text = ax.text(
                x, label_y, f"C{idx}",
                fontsize=10, fontweight="bold", color=font_color,
                ha="center", va="bottom", zorder=5,
                bbox=(
                    dict(boxstyle="round,pad=0.22", facecolor="white",
                         edgecolor="#555555", linewidth=0.8, alpha=0.92)
                    if use_bbox else None
                ),
            )
            text.set_path_effects([pe.withStroke(linewidth=halo_width, foreground="white")])
        return

    # Legacy path (same_row_top_n or no grid)
    same_row_y = None
    if same_row_top_n > 0 and components:
        selected = components[:same_row_top_n]
        selected_pts = np.array([
            pos[n] for comp in selected for n in comp if n in pos
        ], dtype=float)
        if selected_pts.size:
            selected_span = np.ptp(selected_pts, axis=0)
            same_row_y = float(selected_pts[:, 1].max()) + max(
                min_y_offset, float(selected_span[1]) * y_offset_scale,
            )

    for idx, comp in enumerate(components, start=1):
        if len(comp) < min_size:
            continue
        pts = np.array([pos[n] for n in comp if n in pos], dtype=float)
        if not pts.size:
            continue
        x = float(pts[:, 0].mean())
        y = float(pts[:, 1].max())
        span = np.ptp(pts, axis=0)
        y_offset = max(min_y_offset, float(span[1]) * y_offset_scale)
        label_y = same_row_y if same_row_y is not None and idx <= same_row_top_n else y + y_offset
        text = ax.text(
            x, label_y, f"C{idx}",
            fontsize=10, fontweight="bold", color=font_color,
            ha="center", va="bottom", zorder=5,
            bbox=(
                dict(
                    boxstyle="round,pad=0.22",
                    facecolor="white",
                    edgecolor="#555555",
                    linewidth=0.8,
                    alpha=0.92,
                )
                if use_bbox else None
            ),
        )
        text.set_path_effects([pe.withStroke(linewidth=halo_width, foreground="white")])


def _component_display_bbox_for_safe_label(
    fig,
    ax,
    pos: dict[str, np.ndarray],
    nodes: Sequence[str],
    *,
    labels_by_node: Mapping[str, str] | None = None,
    node_size_factor: float = 1.95,
    safety_px: float = 10.0,
):
    """Rendered component bbox including node glyphs, node labels, and safety padding."""
    from matplotlib.transforms import Bbox

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    node_radius_px = math.sqrt(NODE_SIZE * float(node_size_factor) / math.pi) * fig.dpi / 72.0
    label_map = labels_by_node or {node: node for node in nodes}
    wanted_text = {str(label_map.get(node, node)) for node in nodes}

    bboxes = []
    for node in nodes:
        if node not in pos:
            continue
        x, y = ax.transData.transform(pos[node])
        bboxes.append(Bbox.from_extents(
            x - node_radius_px, y - node_radius_px,
            x + node_radius_px, y + node_radius_px,
        ))
    for text in ax.texts:
        if text.get_text() in wanted_text:
            bboxes.append(text.get_window_extent(renderer=renderer))
    bbox = _union_display_bboxes(bboxes)
    if bbox is None:
        return None
    pad = max(0.0, float(safety_px))
    return Bbox.from_extents(
        float(bbox.x0) - pad,
        float(bbox.y0) - pad,
        float(bbox.x1) + pad,
        float(bbox.y1) + pad,
    )


def _data_bbox_from_display_bbox(ax, bbox) -> tuple[float, float, float, float]:
    """Convert a display-coordinate bbox to a data-coordinate bbox."""
    inv = ax.transData.inverted()
    p0 = inv.transform((float(bbox.x0), float(bbox.y0)))
    p1 = inv.transform((float(bbox.x1), float(bbox.y1)))
    return (
        min(float(p0[0]), float(p1[0])),
        min(float(p0[1]), float(p1[1])),
        max(float(p0[0]), float(p1[0])),
        max(float(p0[1]), float(p1[1])),
    )


def _component_label_candidate(
    data_bbox: tuple[float, float, float, float],
    candidate: str,
    gap: float,
) -> tuple[float, float, str, str]:
    """Return x, y, ha, va for one component-label candidate."""
    x0, y0, x1, y1 = data_bbox
    cx = 0.5 * (x0 + x1)
    cy = 0.5 * (y0 + y1)
    gap = float(gap)
    if candidate == "above_center":
        return cx, y1 + gap, "center", "bottom"
    if candidate == "above_left":
        return x0, y1 + gap, "left", "bottom"
    if candidate == "above_right":
        return x1, y1 + gap, "right", "bottom"
    if candidate == "below_center":
        return cx, y0 - gap, "center", "top"
    if candidate == "left_center":
        return x0 - gap, cy, "right", "center"
    if candidate == "right_center":
        return x1 + gap, cy, "left", "center"
    return cx, y1 + gap, "center", "bottom"


def place_component_labels_safely(
    ax,
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    *,
    isolates: Iterable[str] | None = None,
    labels_by_node: Mapping[str, str] | None = None,
    min_size: int = 2,
    label_gap: float = 2.0,
    context: str = "graph",
    extra_data_bboxes: Mapping[str, tuple[float, float, float, float] | None] | None = None,
    candidate_positions: Sequence[str] = (
        "above_center", "above_left", "above_right",
        "below_center", "left_center", "right_center",
    ),
    min_clearance_px: float = 8.0,
    boundary_clearance_px: float = 2.0,
    fontsize: int = 10,
    use_bbox: bool = False,
    halo_width: float = 1.4,
    halo_color: str = "white",
) -> dict[str, dict[str, object]]:
    """
    Place C1, C2, ... labels outside expanded rendered component bboxes.

    Components are still ordered by size/label. The helper is display-bbox
    aware: it sees node glyphs, node labels, resolved anchor overlays via an
    enlarged node radius, title, isolate/legend boxes passed by the caller, and
    already placed component labels.
    """
    from matplotlib.transforms import Bbox

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes_bbox = ax.get_window_extent(renderer=renderer)

    components = ordered_non_isolate_components(G, isolates=isolates)
    component_bboxes: dict[str, object] = {}
    component_data_bboxes: dict[str, tuple[float, float, float, float]] = {}
    for idx, comp in enumerate(components, start=1):
        if len(comp) < min_size:
            continue
        nodes = sorted(n for n in comp if n in pos)
        if not nodes:
            continue
        bbox = _component_display_bbox_for_safe_label(
            fig, ax, pos, nodes,
            labels_by_node=labels_by_node,
            node_size_factor=1.95,
            safety_px=12.0,
        )
        if bbox is None:
            continue
        label = f"C{idx}"
        component_bboxes[label] = bbox
        component_data_bboxes[label] = _data_bbox_from_display_bbox(ax, bbox)

    occupied: dict[str, object] = dict(component_bboxes)
    if extra_data_bboxes:
        for label, bbox in extra_data_bboxes.items():
            if bbox is not None:
                occupied[label] = display_bbox_from_data_bbox(ax, bbox)
    if ax.title is not None and ax.title.get_text():
        occupied["title"] = ax.title.get_window_extent(renderer=renderer)
    legend = ax.get_legend()
    if legend is not None:
        occupied["legend"] = legend.get_window_extent(renderer=renderer)

    results: dict[str, dict[str, object]] = {}
    placed_bboxes: dict[str, object] = {}
    for label in sorted(component_data_bboxes, key=lambda x: int(x[1:])):
        data_bbox = component_data_bboxes[label]
        h = max(data_bbox[3] - data_bbox[1], 1.0)
        gap = max(float(label_gap), 0.08 * h)
        chosen = None
        best_reason = ""
        for candidate in candidate_positions:
            x, y, ha, va = _component_label_candidate(data_bbox, candidate, gap)
            text = ax.text(
                x, y, label,
                fontsize=fontsize, fontweight="bold", color="#1A1A2E",
                ha=ha, va=va, zorder=6,
                bbox=(
                    dict(boxstyle="round,pad=0.22", facecolor="white",
                         edgecolor="#555555", linewidth=0.8, alpha=0.94)
                    if use_bbox else None
                ),
            )
            text.set_path_effects([pe.withStroke(linewidth=float(halo_width), foreground=halo_color)])
            fig.canvas.draw()
            bbox = text.get_window_extent(renderer=fig.canvas.get_renderer())
            boundary = min(
                float(bbox.x0) - float(axes_bbox.x0),
                float(axes_bbox.x1) - float(bbox.x1),
                float(bbox.y0) - float(axes_bbox.y0),
                float(axes_bbox.y1) - float(bbox.y1),
            )
            conflicts = []
            min_clearance = math.inf
            for occ_label, occ_bbox in {**occupied, **placed_bboxes}.items():
                if occ_label == label:
                    occ_name = f"{label}_component"
                else:
                    occ_name = occ_label
                clearance, overlap = _display_clearance(bbox, occ_bbox)
                min_clearance = min(min_clearance, clearance)
                if overlap or clearance < float(min_clearance_px):
                    conflicts.append(f"{occ_name}:{clearance:.1f}{'*' if overlap else ''}")
            if boundary < float(boundary_clearance_px):
                conflicts.append(f"boundary:{boundary:.1f}")
            if not conflicts:
                chosen = (candidate, text, bbox, x, y, min_clearance, boundary)
                break
            best_reason = ",".join(conflicts[:4])
            text.remove()

        if chosen is None:
            candidate = "above_center_forced"
            x, y, ha, va = _component_label_candidate(data_bbox, "above_center", gap * 1.5)
            text = ax.text(
                x, y, label,
                fontsize=fontsize, fontweight="bold", color="#1A1A2E",
                ha=ha, va=va, zorder=6,
                bbox=(
                    dict(boxstyle="round,pad=0.22", facecolor="white",
                         edgecolor="#555555", linewidth=0.8, alpha=0.94)
                    if use_bbox else None
                ),
            )
            text.set_path_effects([pe.withStroke(linewidth=float(halo_width), foreground=halo_color)])
            fig.canvas.draw()
            bbox = text.get_window_extent(renderer=fig.canvas.get_renderer())
            min_clearance = 0.0
            boundary = 0.0
        else:
            candidate, text, bbox, x, y, min_clearance, boundary = chosen
        placed_bboxes[f"{label}_label"] = bbox
        occupied[label] = component_bboxes[label]
        results[label] = {
            "x": float(x),
            "y": float(y),
            "candidate": candidate,
            "min_clearance_px": float(min_clearance),
            "boundary_px": float(boundary),
            "fallback_reason": best_reason,
        }
        print(
            f"[QC] {context} component label placement: {label} "
            f"candidate={candidate} xy=({float(x):.2f},{float(y):.2f}) "
            f"gap={gap:.2f} min_clearance_px={float(min_clearance):.2f} "
            f"boundary_px={float(boundary):.2f} "
            f"fallback_reason={best_reason or 'none'}"
        )
    return results


def add_caption(fig, text: str, fontsize: int = 8) -> None:
    """Place a short explanatory caption below the axes."""
    fig.text(
        0.5, 0.018, text,
        ha="center", va="bottom",
        fontsize=fontsize, color="#242424",
        linespacing=1.25,
    )


# Adaptive component layout geometry. The legacy GRID_* names are retained for
# backward-compatible callers, but grid_layout no longer uses a fixed 27-unit
# row-major step.
GRID_N_COLS: int = 3
GRID_CELL_HALFBOX: float = 9.0
GRID_CELL_PAD: float = 2.5
GRID_GUTTER: float = 8.0
GRID_CELL_SIDE: float = 2.0 * GRID_CELL_HALFBOX + 2.0 * GRID_CELL_PAD
GRID_STEP: float = GRID_CELL_SIDE + GRID_GUTTER

ADAPTIVE_H_GUTTER: float = 8.0
ADAPTIVE_V_GUTTER: float = 5.0
ADAPTIVE_ISOLATE_GAP: float = 4.0
ADAPTIVE_ISOLATE_SPACING: float = 7.0
RBL_SINGLE_COMPONENT_ISOLATE_GAP: float = 2.6
RBL_SINGLE_COMPONENT_LEGEND_GAP_FRAC: float = 0.05


def _component_base_box(n_nodes: int) -> tuple[float, float]:
    """Recommended layout box size by connected-component node count."""
    if n_nodes <= 2:
        return 14.0, 8.0
    if n_nodes <= 4:
        return 18.0, 10.0
    if n_nodes <= 8:
        return 24.0, 14.0
    if n_nodes <= 12:
        return 30.0, 18.0
    if n_nodes <= 18:
        return 36.0, 22.0
    return 42.0, 26.0


def _component_box_spec(nodes: Sequence[str]) -> dict[str, float | int]:
    """
    Return layout and packing dimensions for one non-isolate component.

    Node coordinates are scaled into the base box. The packed component box is
    widened for labels and made slightly taller because labels are drawn at the
    node centres and C-labels are placed above components.
    """
    n_nodes = len(nodes)
    layout_w, layout_h = _component_base_box(n_nodes)
    max_label_len = max((len(str(node)) for node in nodes), default=0)
    label_w = min(0.30 * float(max_label_len), 6.0)
    return {
        "n_nodes": n_nodes,
        "layout_w": layout_w,
        "layout_h": layout_h,
        "box_w": layout_w + label_w,
        "box_h": layout_h + 2.0,
    }


def _rescale_to_box(local_pos: dict, width: float, height: float) -> dict:
    """Rescale spring-layout positions into a centred rectangular box."""
    pts = np.array(list(local_pos.values()), dtype=float)
    if not pts.size:
        return {}
    if len(pts) == 1:
        return {n: np.array([0.0, 0.0], dtype=float) for n in local_pos}

    mn, mx = pts.min(axis=0), pts.max(axis=0)
    span = mx - mn
    centre = (mn + mx) / 2.0
    usable_w = max(float(width), 1e-9)
    usable_h = max(float(height), 1e-9)
    finite_spans = [span[axis] for axis in (0, 1) if span[axis] > 1e-9]
    if not finite_spans:
        scale = 1.0
    elif span[0] <= 1e-9:
        scale = usable_h / float(span[1])
    elif span[1] <= 1e-9:
        scale = usable_w / float(span[0])
    else:
        scale = min(usable_w / float(span[0]), usable_h / float(span[1]))

    return {
        n: (np.array(p, dtype=float) - centre) * scale
        for n, p in local_pos.items()
    }


def _layout_component(
    G_non_iso: nx.Graph,
    nodes: Sequence[str],
    *,
    seed: int,
    spring_iters: int,
    layout_w: float,
    layout_h: float,
) -> dict[str, np.ndarray]:
    """Deterministically lay out one component with sorted node insertion."""
    if len(nodes) == 1:
        return {nodes[0]: np.array([0.0, 0.0], dtype=float)}

    sub_stable = nx.Graph()
    sub_stable.add_nodes_from(nodes)
    sub_stable.add_edges_from(sorted(G_non_iso.subgraph(nodes).edges()))
    raw = nx.spring_layout(sub_stable, seed=seed, iterations=spring_iters)
    return _rescale_to_box(raw, layout_w, layout_h)


def _target_shelf_width(specs: Sequence[dict], h_gutter: float) -> float:
    """Choose a compact row width that still lets small components share rows."""
    if not specs:
        return 0.0
    widths = [float(spec["box_w"]) for spec in specs]
    if len(widths) == 1:
        return widths[0]

    total_width = sum(widths) + h_gutter * (len(widths) - 1)
    largest = max(widths)
    small_row_n = min(4, len(widths))
    small_row_width = (
        sum(sorted(widths)[:small_row_n]) + h_gutter * (small_row_n - 1)
    )
    area = sum(float(spec["box_w"]) * float(spec["box_h"]) for spec in specs)
    area_width = math.sqrt(area) * 1.45
    target = max(largest * 1.35, small_row_width, area_width, 64.0)
    return min(total_width, max(largest, min(target, 104.0)))


def _pack_component_specs(
    specs: Sequence[dict],
    *,
    h_gutter: float,
    v_gutter: float,
) -> tuple[list[dict], float, float, float]:
    """
    Pack component boxes into deterministic shelves sorted by component size.

    Returns the updated specs and the component-grid x extent plus bottom y.
    """
    if not specs:
        return [], 0.0, 0.0, 0.0

    target_width = _target_shelf_width(specs, h_gutter)
    rows: list[list[dict]] = []
    current: list[dict] = []
    current_w = 0.0
    for i, spec in enumerate(specs):
        box_w = float(spec["box_w"])
        remaining = specs[i:]
        keep_tail_pairs_together = (
            int(spec.get("n_nodes", 0)) == 2
            and current
            and any(int(s.get("n_nodes", 0)) != 2 for s in current)
            and len(remaining) >= 2
            and all(int(s.get("n_nodes", 0)) == 2 for s in remaining)
        )
        if keep_tail_pairs_together:
            rows.append(current)
            current = [spec]
            current_w = box_w
            continue
        next_w = box_w if not current else current_w + h_gutter + box_w
        if current and next_w > target_width:
            rows.append(current)
            current = [spec]
            current_w = box_w
        else:
            current.append(spec)
            current_w = next_w
    if current:
        rows.append(current)

    y_top = 0.0
    packed: list[dict] = []
    grid_x0 = math.inf
    grid_x1 = -math.inf
    grid_y0 = math.inf
    for row in rows:
        row_w = sum(float(spec["box_w"]) for spec in row) + h_gutter * (len(row) - 1)
        row_h = max(float(spec["box_h"]) for spec in row)
        x_left = -row_w / 2.0
        cursor = x_left
        for spec in row:
            box_w = float(spec["box_w"])
            box_h = float(spec["box_h"])
            placed = dict(spec)
            placed["cx"] = cursor + box_w / 2.0
            placed["cy"] = y_top - box_h / 2.0
            placed["bbox"] = (cursor, y_top - box_h, cursor + box_w, y_top)
            packed.append(placed)
            grid_x0 = min(grid_x0, cursor)
            grid_x1 = max(grid_x1, cursor + box_w)
            grid_y0 = min(grid_y0, y_top - box_h)
            cursor += box_w + h_gutter
        y_top -= row_h + v_gutter

    return packed, grid_x0, grid_y0, grid_x1 if packed else 0.0


def grid_layout(
    G_non_iso: nx.Graph,
    isolates: list[str],
    *,
    seed: int = LAYOUT_SEED,
    n_cols: int = GRID_N_COLS,
    cell_halfbox: float = GRID_CELL_HALFBOX,
    cell_pad: float = GRID_CELL_PAD,
    gutter: float = ADAPTIVE_H_GUTTER,
    isolate_gap: float = ADAPTIVE_ISOLATE_GAP,
    isolate_spacing: float = ADAPTIVE_ISOLATE_SPACING,
    isolate_box_height_factor: float = 1.0,
    spring_iters: int = 300,
) -> tuple[dict[str, np.ndarray], tuple[float, float, float, float] | None]:
    """
    Lay out non-isolate components in adaptive shelves, isolates in a band below.

    Algorithm
    ---------
    1. Sort connected components by decreasing node count, then node id.
    2. Rebuild each component with sorted node insertion and run a deterministic
       spring_layout(component, seed=seed).
    3. Rescale each component into a node-count-specific rectangular layout box.
    4. Pack component boxes into shelf rows with compact gutters, keeping
       2-node components on shared rows where possible.
    5. Place degree-0 isolates in a separate deterministic band below the grid.

    Returns
    -------
    pos : dict[str, np.ndarray]
        Positions for every node (non-isolates + isolates).
    iso_bbox : (x0, y0, x1, y1) or None
        Bounding box of the isolate band for draw_isolate_zone.
    """
    del n_cols, cell_halfbox, cell_pad  # retained in the signature for compatibility

    components = [
        sorted(comp)
        for comp in sorted(
            nx.connected_components(G_non_iso),
            key=lambda c: (-len(c), sorted(c)[0]),
        )
    ]
    specs: list[dict] = []
    for idx, nodes in enumerate(components, start=1):
        spec = _component_box_spec(nodes)
        spec["idx"] = idx
        spec["nodes"] = nodes
        specs.append(spec)

    packed_specs, grid_x0, grid_y0, grid_x1 = _pack_component_specs(
        specs, h_gutter=max(float(gutter), ADAPTIVE_H_GUTTER),
        v_gutter=ADAPTIVE_V_GUTTER,
    )

    pos_non_iso: dict[str, np.ndarray] = {}
    for spec in packed_specs:
        nodes = list(spec["nodes"])
        local_pos = _layout_component(
            G_non_iso, nodes, seed=seed, spring_iters=spring_iters,
            layout_w=float(spec["layout_w"]), layout_h=float(spec["layout_h"]),
        )
        cx = float(spec["cx"])
        cy = float(spec["cy"])
        for node, p in local_pos.items():
            pos_non_iso[node] = np.array(
                [cx + float(p[0]), cy + float(p[1])], dtype=float
            )

    if not isolates:
        return pos_non_iso, None

    iso_sorted = sorted(isolates)
    n_iso = len(iso_sorted)
    if packed_specs:
        mid_x = (grid_x0 + grid_x1) / 2.0
        component_bottom = grid_y0
    else:
        mid_x = 0.0
        component_bottom = 0.0

    spacing = max(float(isolate_spacing), ADAPTIVE_ISOLATE_SPACING)
    gap = min(max(float(isolate_gap), 4.0), 12.0)
    old_band_pad_y = 2.2
    band_pad_y = old_band_pad_y * max(1.0, float(isolate_box_height_factor))
    y_iso = component_bottom - gap - band_pad_y
    if n_iso == 1:
        iso_xs = np.array([mid_x])
    else:
        total_w = (n_iso - 1) * spacing
        iso_xs = np.linspace(mid_x - total_w / 2.0, mid_x + total_w / 2.0, n_iso)

    iso_pos = {
        n: np.array([float(x), y_iso], dtype=float)
        for n, x in zip(iso_sorted, iso_xs)
    }

    px = max(spacing * 0.55, 3.5)
    iso_bbox = (
        float(iso_xs.min() - px), float(y_iso - band_pad_y),
        float(iso_xs.max() + px), float(y_iso + band_pad_y),
    )
    return {**pos_non_iso, **iso_pos}, iso_bbox


def component_grid_layout(
    G_non_iso: nx.Graph,
    isolates: list[str],
    *,
    seed: int = LAYOUT_SEED,
    n_cols: int = 3,
    n_rows: int | None = None,
    cell_width: float = 54.0,
    cell_height: float = 27.0,
    isolate_gap: float = ADAPTIVE_ISOLATE_GAP,
    isolate_spacing: float = ADAPTIVE_ISOLATE_SPACING,
    isolate_box_height_factor: float = 1.0,
    spring_iters: int = 300,
) -> tuple[
    dict[str, np.ndarray],
    tuple[float, float, float, float] | None,
    list[dict[str, float | int | str | tuple[float, float, float, float]]],
]:
    """
    Lay out components in deterministic row-major grid cells.

    Components are already ordered by size/label via ordered_non_isolate_components:
    C1, C2, C3 occupy the first row when n_cols=3; C4, C5, C6 occupy
    the second row, with the same row-major fallback for other component counts.
    """
    components = [sorted(comp) for comp in ordered_non_isolate_components(G_non_iso)]
    if not components:
        return {}, None, []

    n_cols = max(1, int(n_cols))
    required_rows = int(math.ceil(len(components) / float(n_cols)))
    n_rows_eff = max(required_rows, int(n_rows) if n_rows is not None else required_rows)
    cell_width = float(cell_width)
    cell_height = float(cell_height)

    pos_non_iso: dict[str, np.ndarray] = {}
    placements: list[dict[str, float | int | str | tuple[float, float, float, float]]] = []
    for idx, nodes in enumerate(components, start=1):
        row = (idx - 1) // n_cols
        col = (idx - 1) % n_cols
        cx = (float(col) - (float(n_cols) - 1.0) / 2.0) * cell_width
        cy = -float(row) * cell_height

        spec = _component_box_spec(nodes)
        local_pos = _layout_component(
            G_non_iso, nodes, seed=seed, spring_iters=spring_iters,
            layout_w=float(spec["layout_w"]), layout_h=float(spec["layout_h"]),
        )
        for node, p in local_pos.items():
            pos_non_iso[node] = np.array(
                [cx + float(p[0]), cy + float(p[1])], dtype=float
            )

        pts = np.array([pos_non_iso[n] for n in nodes if n in pos_non_iso], dtype=float)
        bbox = (
            float(pts[:, 0].min()), float(pts[:, 1].min()),
            float(pts[:, 0].max()), float(pts[:, 1].max()),
        )
        placement = {
            "component": f"C{idx}",
            "component_index": idx,
            "row": row + 1,
            "col": col + 1,
            "center_x": cx,
            "center_y": cy,
            "bbox": bbox,
        }
        placements.append(placement)
        print(
            f"[QC] fixed component grid placement: C{idx} "
            f"row={row + 1} col={col + 1} center=({cx:.2f},{cy:.2f}) "
            f"bbox={_fmt_bbox(bbox)}"
        )

    if not isolates:
        return pos_non_iso, None, placements

    iso_sorted = sorted(isolates)
    spacing = max(float(isolate_spacing), ADAPTIVE_ISOLATE_SPACING)
    grid_x0 = -((float(n_cols) - 1.0) * cell_width) / 2.0 - cell_width / 2.0
    grid_x1 = ((float(n_cols) - 1.0) * cell_width) / 2.0 + cell_width / 2.0
    mid_x = (grid_x0 + grid_x1) / 2.0
    pts_non_iso = np.array(list(pos_non_iso.values()), dtype=float)
    component_bottom = float(pts_non_iso[:, 1].min()) if pts_non_iso.size else -float(n_rows_eff - 1) * cell_height

    gap = min(max(float(isolate_gap), 4.0), 12.0)
    band_pad_y = 2.2 * max(1.0, float(isolate_box_height_factor))
    y_iso = component_bottom - gap - band_pad_y
    if len(iso_sorted) == 1:
        iso_xs = np.array([mid_x])
    else:
        total_w = (len(iso_sorted) - 1) * spacing
        iso_xs = np.linspace(mid_x - total_w / 2.0, mid_x + total_w / 2.0, len(iso_sorted))

    iso_pos = {
        n: np.array([float(x), y_iso], dtype=float)
        for n, x in zip(iso_sorted, iso_xs)
    }
    px = max(spacing * 0.55, 3.5)
    iso_bbox = (
        float(iso_xs.min() - px), float(y_iso - band_pad_y),
        float(iso_xs.max() + px), float(y_iso + band_pad_y),
    )
    return {**pos_non_iso, **iso_pos}, iso_bbox, placements


def _clip_fill(value: float, *, default: float) -> float:
    """Keep requested cell-fill fractions inside a useful plotting range."""
    try:
        out = float(value)
    except (TypeError, ValueError):
        out = float(default)
    return float(np.clip(out, 0.20, 0.96))


def _stable_component_graph(G_non_iso: nx.Graph, nodes: Sequence[str]) -> nx.Graph:
    """Build a component graph with deterministic insertion order."""
    sub = nx.Graph()
    ordered_nodes = sorted(nodes)
    sub.add_nodes_from(ordered_nodes)
    sub.add_edges_from(sorted(G_non_iso.subgraph(ordered_nodes).edges()))
    return sub


def _layout_component_for_cell(
    G_non_iso: nx.Graph,
    nodes: Sequence[str],
    *,
    seed: int,
    target_width: float,
    target_height: float,
    spring_iters: int = 300,
    spring_k: float | None = None,
) -> dict[str, np.ndarray]:
    """
    Deterministically lay out one component and fit it to a target cell area.

    Small components use Kamada-Kawai for readability; larger components use a
    seeded spring layout. The final scaling is anisotropic by design so dense
    label-heavy components use the requested cell fill without spilling out of
    their assigned grid cell.
    """
    ordered_nodes = sorted(nodes)
    if not ordered_nodes:
        return {}
    if len(ordered_nodes) == 1:
        return {ordered_nodes[0]: np.array([0.0, 0.0], dtype=float)}

    sub = _stable_component_graph(G_non_iso, ordered_nodes)
    if len(ordered_nodes) <= 20:
        raw = nx.kamada_kawai_layout(sub, weight=None, scale=1.0)
        layout_name = "kamada_kawai"
    else:
        k = spring_k if spring_k is not None else max(0.35, 2.0 / math.sqrt(len(ordered_nodes)))
        raw = nx.spring_layout(
            sub, seed=seed, k=k, iterations=int(spring_iters), scale=1.0
        )
        layout_name = "spring"

    pts = np.array([raw[n] for n in ordered_nodes], dtype=float)
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    centre = (mins + maxs) / 2.0
    spans = maxs - mins
    eps = 1e-9
    scale_x = float(target_width) / float(spans[0]) if spans[0] > eps else 1.0
    scale_y = float(target_height) / float(spans[1]) if spans[1] > eps else 1.0

    fitted = {
        node: np.array([
            (float(raw[node][0]) - float(centre[0])) * scale_x,
            (float(raw[node][1]) - float(centre[1])) * scale_y,
        ], dtype=float)
        for node in ordered_nodes
    }
    print(
        f"[QC] planned local component layout: n={len(ordered_nodes)} "
        f"method={layout_name} target=({float(target_width):.2f},"
        f"{float(target_height):.2f}) raw_span=({float(spans[0]):.3f},"
        f"{float(spans[1]):.3f})"
    )
    return fitted


def planned_component_grid_layout(
    G_non_iso: nx.Graph,
    isolates: list[str],
    *,
    seed: int = LAYOUT_SEED,
    n_cols: int = 3,
    n_rows: int | None = None,
    cell_width: float = 54.0,
    cell_height: float = 42.0,
    row_gap: float = 12.0,
    col_gap: float = 12.0,
    legend_width: float = 50.0,
    isolate_gap: float = ADAPTIVE_ISOLATE_GAP,
    isolate_spacing: float = ADAPTIVE_ISOLATE_SPACING,
    isolate_region_height: float = 9.0,
    isolate_box_height_factor: float | None = None,
    isolate_label_left_pad: float = 0.0,
    isolate_panel_padding_x: float = 12.0,
    isolate_panel_padding_y: float = 1.4,
    isolate_label_band_frac: float = 0.34,
    isolate_max_per_row: int = 6,
    isolate_panel_width_frac: float = 0.62,
    footnote_region_height: float = 6.0,
    title_region_height: float = 6.0,
    component_cell_fill_x: float = 0.78,
    component_cell_fill_y: float = 0.72,
    dense_component_cell_fill_x: float = 0.90,
    dense_component_cell_fill_y: float = 0.82,
    c2_component_cell_fill_x: float = 0.84,
    c2_component_cell_fill_y: float = 0.78,
    dense_component_indices: Sequence[int] = (1,),
    secondary_component_indices: Sequence[int] = (2,),
    spring_iters: int = 300,
    spring_k: float | None = None,
    context: str = "planned component grid",
) -> tuple[
    dict[str, np.ndarray],
    tuple[float, float, float, float] | None,
    list[dict[str, float | int | str | tuple[float, float, float, float]]],
    dict[str, tuple[float, float, float, float]],
]:
    """
    Deterministic region-aware grid planner for thesis graph figures.

    The planner separates the figure into named data-space regions: component
    grid, right-side legend column, isolate strip, title clearance, and
    footnote clearance. Components are laid out locally, then scaled to a
    target fill of their assigned row-major cell. No component is expanded
    after placement.
    """
    components = [sorted(comp) for comp in ordered_non_isolate_components(G_non_iso)]
    if not components:
        return {}, None, [], {}

    n_cols = max(1, int(n_cols))
    required_rows = int(math.ceil(len(components) / float(n_cols)))
    n_rows_eff = max(required_rows, int(n_rows) if n_rows is not None else required_rows)
    cell_width = float(cell_width)
    cell_height = float(cell_height)
    row_gap = max(0.0, float(row_gap))
    col_gap = max(0.0, float(col_gap))
    legend_width = max(0.0, float(legend_width))
    isolate_region_height = max(4.4, float(isolate_region_height))
    if isolate_box_height_factor is not None:
        isolate_region_height = max(
            isolate_region_height,
            4.4 * max(1.0, float(isolate_box_height_factor)),
        )
    isolate_panel_padding_x = max(1.0, float(isolate_panel_padding_x))
    isolate_panel_padding_y = max(0.4, float(isolate_panel_padding_y))
    isolate_label_band_frac = float(np.clip(isolate_label_band_frac, 0.20, 0.45))
    isolate_max_per_row = max(1, int(isolate_max_per_row))
    isolate_panel_width_frac = float(np.clip(isolate_panel_width_frac, 0.40, 0.95))
    footnote_region_height = max(0.0, float(footnote_region_height))
    title_region_height = max(0.0, float(title_region_height))

    grid_width = n_cols * cell_width + (n_cols - 1) * col_gap
    grid_height = n_rows_eff * cell_height + (n_rows_eff - 1) * row_gap
    grid_x0 = -grid_width / 2.0
    grid_x1 = grid_width / 2.0
    grid_y1 = 0.0
    grid_y0 = -grid_height

    legend_x0 = grid_x1 + col_gap
    legend_bbox = (
        legend_x0,
        grid_y1 - cell_height,
        legend_x0 + legend_width,
        grid_y1,
    )
    isolate_y1 = grid_y0 - max(float(isolate_gap), 0.0)
    isolate_region_bbox = (
        grid_x0,
        isolate_y1 - isolate_region_height,
        grid_x1,
        isolate_y1,
    )
    regions = {
        "title_region": (
            grid_x0,
            grid_y1 + 2.0,
            legend_bbox[2],
            grid_y1 + 2.0 + title_region_height,
        ),
        "component_grid": (grid_x0, grid_y0, grid_x1, grid_y1),
        "legend_region": legend_bbox,
        "isolate_region": isolate_region_bbox,
        "footnote_region": (
            grid_x0,
            isolate_region_bbox[1] - 2.0 - footnote_region_height,
            legend_bbox[2],
            isolate_region_bbox[1] - 2.0,
        ),
    }
    print(
        f"[QC] {context} planned regions: "
        f"title={_fmt_bbox(regions['title_region'])} "
        f"grid={_fmt_bbox(regions['component_grid'])} "
        f"legend={_fmt_bbox(regions['legend_region'])} "
        f"isolate={_fmt_bbox(regions['isolate_region'])} "
        f"footnote={_fmt_bbox(regions['footnote_region'])}"
    )

    default_fill_x = _clip_fill(component_cell_fill_x, default=0.78)
    default_fill_y = _clip_fill(component_cell_fill_y, default=0.72)
    dense_fill_x = _clip_fill(dense_component_cell_fill_x, default=0.90)
    dense_fill_y = _clip_fill(dense_component_cell_fill_y, default=0.82)
    c2_fill_x = _clip_fill(c2_component_cell_fill_x, default=0.84)
    c2_fill_y = _clip_fill(c2_component_cell_fill_y, default=0.78)
    dense_set = {int(i) for i in dense_component_indices}
    secondary_set = {int(i) for i in secondary_component_indices}

    pos_non_iso: dict[str, np.ndarray] = {}
    placements: list[dict[str, float | int | str | tuple[float, float, float, float]]] = []
    for idx, nodes in enumerate(components, start=1):
        row = (idx - 1) // n_cols
        col = (idx - 1) % n_cols
        cell_x0 = grid_x0 + col * (cell_width + col_gap)
        cell_x1 = cell_x0 + cell_width
        cell_y1 = grid_y1 - row * (cell_height + row_gap)
        cell_y0 = cell_y1 - cell_height
        cx = (cell_x0 + cell_x1) / 2.0
        cy = (cell_y0 + cell_y1) / 2.0

        if idx in dense_set:
            fill_x, fill_y = dense_fill_x, dense_fill_y
            fill_role = "dense"
        elif idx in secondary_set:
            fill_x, fill_y = c2_fill_x, c2_fill_y
            fill_role = "secondary"
        else:
            fill_x, fill_y = default_fill_x, default_fill_y
            fill_role = "default"

        target_width = cell_width * fill_x
        target_height = cell_height * fill_y
        local_pos = _layout_component_for_cell(
            G_non_iso, nodes, seed=seed,
            target_width=target_width, target_height=target_height,
            spring_iters=spring_iters, spring_k=spring_k,
        )
        for node, p in local_pos.items():
            pos_non_iso[node] = np.array(
                [cx + float(p[0]), cy + float(p[1])], dtype=float
            )

        pts = np.array([pos_non_iso[n] for n in nodes if n in pos_non_iso], dtype=float)
        bbox = (
            float(pts[:, 0].min()), float(pts[:, 1].min()),
            float(pts[:, 0].max()), float(pts[:, 1].max()),
        )
        cell_bbox = (float(cell_x0), float(cell_y0), float(cell_x1), float(cell_y1))
        density = component_edge_density(G_non_iso, nodes)
        placement = {
            "component": f"C{idx}",
            "component_index": idx,
            "row": row + 1,
            "col": col + 1,
            "center_x": cx,
            "center_y": cy,
            "cell_bbox": cell_bbox,
            "bbox": bbox,
            "fill_x": fill_x,
            "fill_y": fill_y,
            "density": density,
            "n": len(nodes),
            "fill_role": fill_role,
        }
        placements.append(placement)
        cell_margins = (
            bbox[0] - cell_bbox[0],
            bbox[1] - cell_bbox[1],
            cell_bbox[2] - bbox[2],
            cell_bbox[3] - bbox[3],
        )
        outside_cell = any(float(m) < -1e-6 for m in cell_margins)
        print(
            f"[QC] {context} placement: C{idx} row={row + 1} col={col + 1} "
            f"n={len(nodes)} density={density:.3f} fill_role={fill_role} "
            f"fill=({fill_x:.2f},{fill_y:.2f}) "
            f"cell={_fmt_bbox(cell_bbox)} bbox={_fmt_bbox(bbox)} "
            f"cell_margins=({cell_margins[0]:.2f},{cell_margins[1]:.2f},"
            f"{cell_margins[2]:.2f},{cell_margins[3]:.2f}) "
            f"outside_cell={'yes' if outside_cell else 'no'}"
        )

    iso_bbox = None
    pos_iso: dict[str, np.ndarray] = {}
    if isolates:
        iso_sorted = sorted(isolates)
        n_iso = len(iso_sorted)
        n_iso_rows = int(math.ceil(n_iso / float(isolate_max_per_row)))
        marker_diameter = max(4.0, float(isolate_spacing) * 0.25)
        label_text_height = 1.8
        node_label_height = 2.2
        label_band_height = max(
            1.2 * label_text_height,
            isolate_region_height * isolate_label_band_frac,
        )
        node_row_height = max(
            2.2 * marker_diameter + node_label_height,
            6.0,
        )
        computed_panel_height = (
            label_band_height
            + n_iso_rows * node_row_height
            + 2.0 * isolate_panel_padding_y
        )
        if computed_panel_height > isolate_region_height:
            isolate_region_height = computed_panel_height
            isolate_y1 = grid_y0 - max(float(isolate_gap), 0.0)
            isolate_region_bbox = (
                grid_x0,
                isolate_y1 - isolate_region_height,
                grid_x1,
                isolate_y1,
            )
            regions["isolate_region"] = isolate_region_bbox
            regions["footnote_region"] = (
                grid_x0,
                isolate_region_bbox[1] - 2.0 - footnote_region_height,
                legend_bbox[2],
                isolate_region_bbox[1] - 2.0,
            )

        spacing = max(float(isolate_spacing), ADAPTIVE_ISOLATE_SPACING)
        iso_y0, iso_y1 = isolate_region_bbox[1], isolate_region_bbox[3]
        panel_width = max(
            grid_width * isolate_panel_width_frac,
            min(n_iso, isolate_max_per_row) * spacing + 2.0 * isolate_panel_padding_x,
            72.0,
        )
        iso_bbox = (
            -panel_width / 2.0,
            float(iso_y0),
            panel_width / 2.0,
            float(iso_y1),
        )
        label_band_height = (float(iso_bbox[3]) - float(iso_bbox[1])) * isolate_label_band_frac
        label_band_y0 = float(iso_bbox[1])
        label_band_y1 = label_band_y0 + label_band_height
        node_band_y0 = label_band_y1 + isolate_panel_padding_y
        node_band_y1 = float(iso_bbox[3]) - isolate_panel_padding_y
        if node_band_y1 <= node_band_y0:
            node_band_y1 = node_band_y0 + 1.0
        node_band_height = node_band_y1 - node_band_y0
        row_step = node_band_height / float(n_iso_rows + 1)

        pos_iso = {}
        for row in range(n_iso_rows):
            row_start = row * isolate_max_per_row
            row_nodes = iso_sorted[row_start:row_start + isolate_max_per_row]
            row_count = len(row_nodes)
            if row_count == 1:
                xs = np.array([0.0])
            else:
                xs = np.linspace(
                    float(iso_bbox[0]) + isolate_panel_padding_x,
                    float(iso_bbox[2]) - isolate_panel_padding_x,
                    row_count,
                )
            y_iso = node_band_y1 - float(row + 1) * row_step
            for node, x in zip(row_nodes, xs):
                pos_iso[node] = np.array([float(x), float(y_iso)], dtype=float)

        regions["isolate_label_band"] = (
            float(iso_bbox[0]), label_band_y0,
            float(iso_bbox[2]), label_band_y1,
        )
        regions["isolate_node_band"] = (
            float(iso_bbox[0]), node_band_y0,
            float(iso_bbox[2]), node_band_y1,
        )
        print(
            f"[QC] {context} isolate region: bbox={_fmt_bbox(iso_bbox)} "
            f"label_band={_fmt_bbox(regions['isolate_label_band'])} "
            f"node_band={_fmt_bbox(regions['isolate_node_band'])} "
            f"rows={n_iso_rows} max_per_row={isolate_max_per_row} "
            f"spacing_min={spacing:.2f} padding=({isolate_panel_padding_x:.2f},"
            f"{isolate_panel_padding_y:.2f})"
        )
    else:
        print(f"[QC] {context} isolate region: no isolates")

    if iso_bbox is not None:
        regions["isolate_box"] = iso_bbox
    report_planned_layout_clearance(context, placements, regions)
    return {**pos_non_iso, **pos_iso}, iso_bbox, placements, regions


def legend_anchor_from_region(
    region_bbox: tuple[float, float, float, float],
    *,
    x_pad_frac: float = 0.05,
    y_frac: float = 0.50,
) -> tuple[float, float]:
    """Return a data-coordinate legend anchor inside a reserved region."""
    x0, y0, x1, y1 = region_bbox
    x = float(x0) + (float(x1) - float(x0)) * float(x_pad_frac)
    y = float(y0) + (float(y1) - float(y0)) * float(y_frac)
    return x, y


def report_planned_layout_clearance(
    context: str,
    placements: Sequence[Mapping[str, object]],
    regions: Mapping[str, tuple[float, float, float, float]],
) -> dict[str, object]:
    """Log data-coordinate clearance between planned major layout regions."""
    bboxes: dict[str, tuple[float, float, float, float]] = {
        str(p["component"]): tuple(float(v) for v in p["bbox"])  # type: ignore[index]
        for p in placements
        if "component" in p and "bbox" in p
    }
    for key in ("legend_region", "isolate_box", "isolate_region", "title_region", "footnote_region"):
        if key in regions:
            bboxes[key] = tuple(float(v) for v in regions[key])

    pair_specs: list[tuple[str, str]] = []
    comp_labels = sorted(
        [k for k in bboxes if k.startswith("C")],
        key=lambda x: int(x[1:]) if x[1:].isdigit() else 999,
    )
    for i, a in enumerate(comp_labels):
        for b in comp_labels[i + 1:]:
            pair_specs.append((a, b))
        for b in ("legend_region", "isolate_box", "title_region"):
            if b in bboxes:
                pair_specs.append((a, b))
    for a, b in (
        ("isolate_box", "footnote_region"),
        ("isolate_region", "footnote_region"),
        ("legend_region", "title_region"),
    ):
        if a in bboxes and b in bboxes:
            pair_specs.append((a, b))

    min_pair = None
    any_overlap = False
    for a, b in pair_specs:
        clearance, overlap = _bbox_clearance(bboxes[a], bboxes[b])
        any_overlap = any_overlap or overlap
        if min_pair is None or clearance < float(min_pair[2]):
            min_pair = (a, b, clearance, overlap)
        print(
            f"[QC] {context} planned clearance {a}_vs_{b}: "
            f"clearance={clearance:.2f} overlap={'yes' if overlap else 'no'}"
        )
    if min_pair is not None:
        print(
            f"[QC] {context} planned clearance summary: "
            f"min_pair={min_pair[0]}_vs_{min_pair[1]} "
            f"min_clearance={float(min_pair[2]):.2f} "
            f"overlap={'yes' if any_overlap else 'no'}"
        )
    return {"any_overlap": any_overlap, "min_pair": min_pair}


def rbl_one_component_layout(
    G_non_iso: nx.Graph,
    isolates: list[str],
    *,
    seed: int = LAYOUT_SEED,
    isolate_gap: float = RBL_SINGLE_COMPONENT_ISOLATE_GAP,
    isolate_spacing: float = ADAPTIVE_ISOLATE_SPACING,
    isolate_box_height_factor: float = 1.0,
    spring_iters: int = 300,
) -> tuple[dict[str, np.ndarray], tuple[float, float, float, float] | None]:
    """
    Compact RBL-only layout for sparse resolved/support-threshold graphs.

    RBL has a single non-isolate component in these thesis figures, so this
    bypasses the multi-component shelf reservation and treats C1, the legend,
    and the isolate band as the only layout boxes.
    """
    components = [
        sorted(comp)
        for comp in sorted(
            nx.connected_components(G_non_iso),
            key=lambda c: (-len(c), sorted(c)[0]),
        )
    ]
    if len(components) != 1:
        return grid_layout(
            G_non_iso, isolates, seed=seed, isolate_gap=isolate_gap,
            isolate_spacing=isolate_spacing,
            isolate_box_height_factor=isolate_box_height_factor,
            spring_iters=spring_iters,
        )

    nodes = components[0]
    spec = _component_box_spec(nodes)
    local_pos = _layout_component(
        G_non_iso, nodes, seed=seed, spring_iters=spring_iters,
        layout_w=float(spec["layout_w"]), layout_h=float(spec["layout_h"]),
    )
    pos_non_iso = {
        node: np.array([float(p[0]), float(p[1])], dtype=float)
        for node, p in local_pos.items()
    }
    if not isolates:
        return pos_non_iso, None

    iso_sorted = sorted(isolates)
    spacing = max(float(isolate_spacing), ADAPTIVE_ISOLATE_SPACING)
    band_pad_y = 3.0 * max(1.0, float(isolate_box_height_factor))
    component_bottom = -float(spec["box_h"]) / 2.0
    iso_top = component_bottom - max(float(isolate_gap), 1.8)
    y_iso = iso_top - band_pad_y
    if len(iso_sorted) == 1:
        iso_xs = np.array([0.0])
    else:
        total_w = (len(iso_sorted) - 1) * spacing
        iso_xs = np.linspace(-total_w / 2.0, total_w / 2.0, len(iso_sorted))

    iso_pos = {
        node: np.array([float(x), y_iso], dtype=float)
        for node, x in zip(iso_sorted, iso_xs)
    }
    px = max(spacing * 0.55, 3.5)
    iso_bbox = (
        float(iso_xs.min() - px), float(y_iso - band_pad_y),
        float(iso_xs.max() + px), float(iso_top),
    )
    return {**pos_non_iso, **iso_pos}, iso_bbox


def place_isolates_strip(
    pos: dict[str, np.ndarray],
    isolates: list[str],
    pad_below: float = 1.5,
    spacing: float = 1.6,
    max_per_row: int = 6,
    compact: bool = False,
    single_row: bool = False,
) -> tuple[dict[str, np.ndarray], tuple[float, float, float, float] | None]:
    """
    Lay out isolate nodes in a compact labelled strip below the rest of the
    layout. Returns (positions_dict, strip_bbox_for_zone_label).

    By default this preserves the legacy full-width strip used by consensus
    figures. Resolved plots pass compact=True to keep isolates in a small panel
    whose width is determined by isolate count, not by the full graph span.
    """
    if not isolates:
        return {}, None

    non_iso_pts = np.array(
        [pos[n] for n in pos if n not in set(isolates)], dtype=float
    ) if pos else np.empty((0, 2))

    n_iso = len(isolates)
    if not compact:
        if non_iso_pts.size:
            xmin, ymin = non_iso_pts.min(axis=0)
            xmax, _ = non_iso_pts.max(axis=0)
            mid_x = 0.5 * (xmin + xmax)
            y_strip = ymin - pad_below
            x_span = max(xmax - xmin, len(isolates) * spacing * 0.5)
        else:
            mid_x = 0.0
            y_strip = 0.0
            x_span = max(1.0, len(isolates) * spacing * 0.5)

        total_width = max(x_span * 0.85, (n_iso - 1) * spacing)
        if n_iso == 1:
            xs = np.array([mid_x])
        else:
            xs = np.linspace(mid_x - total_width / 2, mid_x + total_width / 2, n_iso)

        out = {n: np.array([x, y_strip], dtype=float) for n, x in zip(sorted(isolates), xs)}
        pad_x = max(spacing * 0.9, total_width * 0.06)
        pad_y = spacing * 0.9
        bbox = (
            float(xs.min() - pad_x),
            float(y_strip - pad_y),
            float(xs.max() + pad_x),
            float(y_strip + pad_y),
        )
        return out, bbox

    if single_row:
        n_cols = n_iso
    elif n_iso <= 3:
        n_cols = n_iso
    else:
        n_cols = max(1, min(max_per_row, int(np.ceil(np.sqrt(n_iso)))))
    n_rows = int(np.ceil(n_iso / n_cols))

    if non_iso_pts.size:
        xmin, ymin = non_iso_pts.min(axis=0)
        xmax, _ = non_iso_pts.max(axis=0)
        mid_x = 0.5 * (xmin + xmax)
        y_top = ymin - pad_below
    else:
        mid_x = 0.0
        y_top = 0.0

    grid_width = (n_cols - 1) * spacing
    row_gap = spacing * 1.05
    out: dict[str, np.ndarray] = {}
    for idx, node in enumerate(sorted(isolates)):
        row = idx // n_cols
        col = idx % n_cols
        row_count = min(n_cols, n_iso - row * n_cols)
        row_width = (row_count - 1) * spacing
        x0 = mid_x - row_width / 2.0
        out[node] = np.array([x0 + col * spacing, y_top - row * row_gap], dtype=float)

    pts = np.array(list(out.values()), dtype=float)
    pad_x = max(0.75, spacing * 0.52, grid_width * 0.08)
    pad_y = max(0.75, spacing * 0.58)
    bbox = (
        float(pts[:, 0].min() - pad_x),
        float(pts[:, 1].min() - pad_y),
        float(pts[:, 0].max() + pad_x),
        float(pts[:, 1].max() + pad_y),
    )
    return out, bbox


def draw_isolate_panel(
    ax,
    bbox: tuple[float, float, float, float],
    *,
    isolate_nodes: Sequence[str] | None = None,
    pos: Mapping[str, np.ndarray] | None = None,
    labels_by_node: Mapping[str, str] | None = None,
    label_text: str = "Isolates (degree 0)",
    label_fontsize: int = 8,
    node_label_fontsize: int = 9,
    panel_padding_x: float = 12.0,
    panel_padding_y: float = 1.4,
    label_band_frac: float = 0.34,
    marker_size: int | None = None,
    draw_nodes: bool = True,
) -> tuple[float, float, float, float]:
    """
    Draw a structured isolate panel with a lower-right label band and node band.

    The node coordinates are supplied by the layout planner; this helper owns
    the shared visual treatment for the dashed panel, the bottom-right caption, and
    optionally the isolate diamonds plus their labels.
    """
    x0, y0, x1, y1 = (float(v) for v in bbox)
    panel_w = x1 - x0
    panel_h = y1 - y0
    panel_padding_x = max(0.5, float(panel_padding_x))
    panel_padding_y = max(0.3, float(panel_padding_y))
    label_band_frac = float(np.clip(label_band_frac, 0.20, 0.45))
    label_band_h = panel_h * label_band_frac
    label_band_y0 = y0
    label_band_y1 = y0 + label_band_h
    node_band_y0 = label_band_y1 + panel_padding_y
    node_band_y1 = y1 - panel_padding_y

    ax.add_patch(mpatches.Rectangle(
        (x0, y0), panel_w, panel_h,
        fill=False, edgecolor="#777777", linewidth=1.0,
        linestyle=(0, (4, 3)), zorder=0.5,
    ))
    header = ax.text(
        x1 - panel_padding_x,
        y0 + panel_padding_y,
        label_text,
        fontsize=label_fontsize, color="#444444", style="italic",
        ha="right", va="bottom", zorder=0.8,
        bbox=dict(
            boxstyle="round,pad=0.10",
            facecolor="white", edgecolor="none", alpha=0.82,
        ),
    )
    header.set_path_effects([pe.withStroke(linewidth=1.1, foreground="white")])

    nodes = [n for n in (isolate_nodes or []) if pos is not None and n in pos]
    if draw_nodes and nodes:
        xs = [float(pos[n][0]) for n in nodes]
        ys = [float(pos[n][1]) for n in nodes]
        ax.scatter(
            xs, ys,
            s=int(NODE_SIZE * 0.85) if marker_size is None else int(marker_size),
            marker="D", facecolors=ISOLATE_COLOUR, edgecolors=NODE_EDGE_COLOUR,
            linewidths=1.4, alpha=0.95, zorder=2.4,
        )
        for node in nodes:
            label = str(labels_by_node.get(node, node) if labels_by_node else node)
            text = ax.text(
                float(pos[node][0]), float(pos[node][1]), label,
                fontsize=node_label_fontsize, fontweight="bold",
                color="#1A1A2E", ha="center", va="center",
                zorder=4.2,
                bbox=dict(
                    boxstyle="round,pad=0.08",
                    facecolor="white", edgecolor="none", alpha=0.72,
                ),
            )
            text.set_path_effects([pe.withStroke(linewidth=2.6, foreground="white")])

        label_clearance = min(float(pos[n][1]) for n in nodes) - label_band_y1
        border_clearance = min(
            min(float(pos[n][0]) - x0, x1 - float(pos[n][0])) for n in nodes
        )
        header_to_nodes = min(float(pos[n][1]) - label_band_y1 for n in nodes)
        print(
            "[QC] isolate panel layered draw: "
            f"bbox={_fmt_bbox((x0, y0, x1, y1))} "
            f"label_xy=({x1 - panel_padding_x:.2f},{y0 + panel_padding_y:.2f}) "
            f"label_band=({x0:.2f},{label_band_y0:.2f},{x1:.2f},{label_band_y1:.2f}) "
            f"node_band=({x0:.2f},{node_band_y0:.2f},{x1:.2f},{node_band_y1:.2f}) "
            f"n_nodes={len(nodes)} label_to_node_min={header_to_nodes:.2f} "
            f"node_to_border_min={border_clearance:.2f} "
            f"node_to_label_band_min={label_clearance:.2f}"
        )
    else:
        print(
            "[QC] isolate panel layered draw: "
            f"bbox={_fmt_bbox((x0, y0, x1, y1))} "
            f"label_xy=({x1 - panel_padding_x:.2f},{y0 + panel_padding_y:.2f}) "
            "n_nodes=0"
        )

    return (x0, y0, x1, y1)


def position_bbox(
    pos: Mapping[str, np.ndarray],
    nodes: Iterable[str] | None = None,
) -> tuple[float, float, float, float] | None:
    """Return the data-coordinate bounding box for selected node positions."""
    selected = list(nodes) if nodes is not None else list(pos.keys())
    pts = np.array([pos[n] for n in selected if n in pos], dtype=float)
    if pts.size == 0:
        return None
    return (
        float(pts[:, 0].min()),
        float(pts[:, 1].min()),
        float(pts[:, 0].max()),
        float(pts[:, 1].max()),
    )


def _transform_bbox(
    bbox: tuple[float, float, float, float] | None,
    *,
    center_x: float,
    center_y: float,
    scale_x: float,
    scale_y: float,
) -> tuple[float, float, float, float] | None:
    if bbox is None:
        return None
    x0, y0, x1, y1 = (float(v) for v in bbox)
    pts = np.array([[x0, y0], [x0, y1], [x1, y0], [x1, y1]], dtype=float)
    pts[:, 0] = center_x + (pts[:, 0] - center_x) * scale_x
    pts[:, 1] = center_y + (pts[:, 1] - center_y) * scale_y
    return (
        float(pts[:, 0].min()),
        float(pts[:, 1].min()),
        float(pts[:, 0].max()),
        float(pts[:, 1].max()),
    )


def apply_position_transform_to_bbox(
    bbox: tuple[float, float, float, float] | None,
    transform: Mapping[str, float] | None,
) -> tuple[float, float, float, float] | None:
    """Apply a position rescaling transform to an auxiliary layout bbox."""
    if bbox is None or not transform:
        return bbox
    return _transform_bbox(
        bbox,
        center_x=float(transform["center_x"]),
        center_y=float(transform["center_y"]),
        scale_x=float(transform["scale_x"]),
        scale_y=float(transform["scale_y"]),
    )


def rescale_positions_to_target_fill(
    pos: dict[str, np.ndarray],
    *,
    fill_x: float | None = None,
    fill_y: float | None = None,
    figure_aspect: float | None = None,
    nodes: Iterable[str] | None = None,
    max_scale: float = 3.0,
) -> dict[str, float] | None:
    """
    Rescale final positions so the node extent better matches the figure aspect.

    This keeps the fixed thesis figure size unchanged and adjusts only final
    data-coordinate positions before tight axis framing.
    """
    if not pos:
        return None
    fill_x = float(fill_x) if fill_x is not None else 0.0
    fill_y = float(fill_y) if fill_y is not None else 0.0
    if fill_x <= 0.0 and fill_y <= 0.0:
        return None
    if figure_aspect is None or float(figure_aspect) <= 0.0:
        return None

    bbox = position_bbox(pos, nodes=nodes)
    if bbox is None:
        return None
    x0, y0, x1, y1 = bbox
    width = max(x1 - x0, 1e-9)
    height = max(y1 - y0, 1e-9)
    current_aspect = width / height
    target_aspect = float(figure_aspect)
    if fill_x > 0.0 and fill_y > 0.0:
        target_aspect *= fill_x / fill_y
    elif fill_x > 0.0:
        target_aspect *= fill_x
    elif fill_y > 0.0:
        target_aspect /= fill_y
    if target_aspect <= 0.0:
        return None

    scale_x = 1.0
    scale_y = 1.0
    if current_aspect < target_aspect:
        scale_x = min(max_scale, target_aspect / current_aspect)
    elif current_aspect > target_aspect:
        scale_y = min(max_scale, current_aspect / target_aspect)
    if abs(scale_x - 1.0) < 1e-6 and abs(scale_y - 1.0) < 1e-6:
        return None

    center_x = (x0 + x1) / 2.0
    center_y = (y0 + y1) / 2.0
    for node, point in pos.items():
        arr = np.array(point, dtype=float)
        arr[0] = center_x + (arr[0] - center_x) * scale_x
        arr[1] = center_y + (arr[1] - center_y) * scale_y
        pos[node] = arr

    return {
        "center_x": center_x,
        "center_y": center_y,
        "scale_x": scale_x,
        "scale_y": scale_y,
        "target_aspect": target_aspect,
        "current_aspect": current_aspect,
    }


def compact_isolate_panel_layout(
    pos: dict[str, np.ndarray],
    isolates: Sequence[str],
    bbox: tuple[float, float, float, float] | None,
    *,
    spacing: float,
    padding_x: float,
    width_mode: str = "layout",
    align_bbox: tuple[float, float, float, float] | None = None,
    label_map: Mapping[str, str] | None = None,
    min_label_spacing: float = 0.0,
    label_padding_factor: float = 0.55,
    context: str = "graph",
) -> tuple[float, float, float, float] | None:
    """Compact the isolate strip around its content without changing membership."""
    if bbox is None or not isolates or width_mode != "content":
        return bbox

    x0, y0, x1, y1 = (float(v) for v in bbox)
    nodes = [n for n in sorted(isolates) if n in pos]
    if not nodes:
        return bbox
    center_x = (x0 + x1) / 2.0
    if align_bbox is not None:
        center_x = (float(align_bbox[0]) + float(align_bbox[2])) / 2.0

    spacing_requested = max(0.1, float(spacing))
    padding_x = max(0.1, float(padding_x))
    label_lengths = []
    for node in nodes:
        raw_label = label_map.get(node, node) if label_map is not None else node
        label_lengths.append(len(str(raw_label if raw_label is not None else node)))
    max_label_len = max(label_lengths) if label_lengths else 0
    # Isolate labels sit on a single horizontal strip; they need less spacing
    # than the two-dimensional node-label de-collision envelope.
    label_extent = float(max_label_len) * float(label_padding_factor) * 0.55
    label_spacing = max(
        0.0,
        float(min_label_spacing),
        label_extent,
    )
    spacing = max(spacing_requested, label_spacing)
    label_half_extent = 0.5 * max(0.0, label_extent)
    node_span = spacing * max(len(nodes) - 1, 0)
    panel_width = max(
        node_span + 2.0 * (padding_x + label_half_extent),
        spacing + 2.0 * (padding_x + label_half_extent),
    )
    y_center = float(np.mean([pos[n][1] for n in nodes]))

    if len(nodes) == 1:
        pos[nodes[0]] = np.array([center_x, y_center], dtype=float)
    else:
        start_x = center_x - node_span / 2.0
        for idx, node in enumerate(nodes):
            pos[node] = np.array([start_x + idx * spacing, y_center], dtype=float)

    new_bbox = (
        center_x - panel_width / 2.0,
        y0,
        center_x + panel_width / 2.0,
        y1,
    )
    print(
        f"[QC] {context} compact isolate spacing: "
        f"requested_spacing={spacing_requested:.2f} effective_spacing={spacing:.2f} "
        f"min_label_spacing={float(min_label_spacing):.2f} "
        f"max_label_len={max_label_len} label_padding_factor={float(label_padding_factor):.2f} "
        f"label_half_extent={label_half_extent:.2f} "
        f"bbox={_fmt_bbox(new_bbox)}"
    )
    return new_bbox


def enforce_component_isolate_gap(
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Sequence[str],
    iso_bbox: tuple[float, float, float, float] | None,
    *,
    min_gap: float = 4.0,
    component_isolate_gap_factor: float = 1.0,
    bbox_margin: float = 0.0,
    label_padding_factor: float = 0.0,
    labels_by_node: Mapping[str, str] | None = None,
    shift_isolate_panel: bool = True,
    context: str = "graph",
) -> tuple[tuple[float, float, float, float] | None, dict[str, float | bool | str]]:
    """
    Ensure the isolate panel has a minimum vertical data-coordinate gap.

    The operation is purely visual: if the non-isolate component bbox is too
    close to or overlaps the isolate panel, isolates and the isolate bbox are
    shifted downward as a group. Non-isolate component coordinates are not
    altered.
    """
    if iso_bbox is None or not isolates:
        return iso_bbox, {
            "enabled": True,
            "available": False,
            "overlap_or_near_overlap": False,
            "shift_y": 0.0,
        }
    non_iso_nodes = [n for n in G.nodes() if n not in set(isolates) and n in pos]
    comp_bbox = position_bbox(pos, nodes=non_iso_nodes)
    if comp_bbox is None:
        return iso_bbox, {
            "enabled": True,
            "available": False,
            "overlap_or_near_overlap": False,
            "shift_y": 0.0,
        }

    label_map = labels_by_node or {node: node for node in non_iso_nodes}
    max_label_len = max((len(str(label_map.get(node, node))) for node in non_iso_nodes), default=0)
    label_margin = max(0.0, float(max_label_len - 8)) * 0.04 * max(0.0, float(label_padding_factor))
    margin = max(0.0, float(bbox_margin)) + label_margin
    effective_min_gap = max(0.0, float(min_gap)) * max(0.0, float(component_isolate_gap_factor))

    comp_bottom = float(comp_bbox[1]) - margin
    iso_top = float(iso_bbox[3])
    gap_before = comp_bottom - iso_top
    overlap_or_near = gap_before < effective_min_gap
    shift_y = effective_min_gap - gap_before if overlap_or_near else 0.0

    new_bbox = iso_bbox
    if overlap_or_near and shift_isolate_panel and shift_y > 0.0:
        for node in isolates:
            if node in pos:
                pos[node] = np.array(
                    [float(pos[node][0]), float(pos[node][1]) - shift_y],
                    dtype=float,
                )
        new_bbox = (
            float(iso_bbox[0]),
            float(iso_bbox[1]) - shift_y,
            float(iso_bbox[2]),
            float(iso_bbox[3]) - shift_y,
        )
    gap_after = (float(comp_bbox[1]) - margin) - float(new_bbox[3]) if new_bbox is not None else gap_before
    print(
        f"[QC] {context} component-isolate gap: "
        f"component_bbox={_fmt_bbox(comp_bbox)} isolate_bbox_before={_fmt_bbox(iso_bbox)} "
        f"bbox_margin={margin:.2f} min_gap={effective_min_gap:.2f} "
        f"gap_before={gap_before:.2f} overlap_or_near={'yes' if overlap_or_near else 'no'} "
        f"shift_isolate_y={shift_y:.2f} isolate_bbox_after={_fmt_bbox(new_bbox)} "
        f"gap_after={gap_after:.2f}"
    )
    return new_bbox, {
        "enabled": True,
        "available": True,
        "overlap_or_near_overlap": bool(overlap_or_near),
        "shift_y": float(shift_y),
        "gap_before": float(gap_before),
        "gap_after": float(gap_after),
        "min_gap": float(effective_min_gap),
        "component_bbox": _fmt_bbox(comp_bbox),
        "isolate_bbox": _fmt_bbox(new_bbox),
    }


def draw_isolate_zone(
    ax,
    bbox: tuple[float, float, float, float],
    label: str = "Isolates (degree 0)",
    label_y_frac: float = 0.010,
) -> None:
    """
    Draw a soft labelled rectangle around the isolate strip so degree-0 nodes
    cannot be misread as strays. No fill — just a thin dashed border with the
    label aligned to its bottom-left.
    """
    x0, y0, x1, y1 = bbox
    ax.add_patch(mpatches.Rectangle(
        (x0, y0), x1 - x0, y1 - y0,
        fill=False, edgecolor="#777777", linewidth=1.0,
        linestyle=(0, (4, 3)), zorder=0.5,
    ))
    label_y_frac = float(np.clip(label_y_frac, 0.005, 0.40))
    ax.text(
        x0 + (x1 - x0) * 0.02, y0 + (y1 - y0) * label_y_frac,
        label,
        fontsize=8, color="#444444", style="italic",
        ha="left", va="bottom", zorder=0.6,
    )


def frame_axes(
    ax,
    pos: dict[str, np.ndarray],
    margin: float = 0.06,
    extra_bboxes: Sequence[tuple[float, float, float, float]] | None = None,
    include_text: bool = True,
    fill_x: float | None = None,
    fill_y: float | None = None,
) -> None:
    """
    Crop axes to the layout bounding box with equal physical padding on all sides.

    Uses max(span_x, span_y) * margin as padding so the margin is the same number
    of data units horizontally and vertically, which with set_aspect("equal") produces
    visually equal margins. Optional bboxes and existing data-coordinate text
    artists are included so isolate bands and labels do not drive fixed padding.
    """
    base_points = [np.array(p, dtype=float) for p in pos.values()]
    if extra_bboxes:
        for x0, y0, x1, y1 in extra_bboxes:
            base_points.extend([
                np.array([x0, y0], dtype=float),
                np.array([x0, y1], dtype=float),
                np.array([x1, y0], dtype=float),
                np.array([x1, y1], dtype=float),
            ])
    if not base_points:
        return

    points = list(base_points)
    for _ in range(2 if include_text else 1):
        pts = np.array(points, dtype=float)
        mn, mx = pts.min(axis=0), pts.max(axis=0)
        span = mx - mn
        span = np.where(span < 1e-9, 1.0, span)
        use_fill = (
            fill_x is not None and float(fill_x) > 0.0
        ) or (
            fill_y is not None and float(fill_y) > 0.0
        )
        if use_fill:
            x_fill = float(fill_x) if fill_x is not None and float(fill_x) > 0.0 else 1.0 / (1.0 + 2.0 * margin)
            y_fill = float(fill_y) if fill_y is not None and float(fill_y) > 0.0 else 1.0 / (1.0 + 2.0 * margin)
            x_fill = float(np.clip(x_fill, 0.05, 0.98))
            y_fill = float(np.clip(y_fill, 0.05, 0.98))
            x_span = float(span[0]) / x_fill
            y_span = float(span[1]) / y_fill
            try:
                fig_w, fig_h = ax.figure.get_size_inches()
                axes_box = ax.get_position()
                axes_aspect = (fig_w * axes_box.width) / max(fig_h * axes_box.height, 1e-9)
                if axes_aspect > 0.0:
                    if x_span / max(y_span, 1e-9) < axes_aspect:
                        x_span = y_span * axes_aspect
                    else:
                        y_span = x_span / axes_aspect
            except Exception:
                pass
            center = (mn + mx) / 2.0
            ax.set_xlim(center[0] - x_span / 2.0, center[0] + x_span / 2.0)
            ax.set_ylim(center[1] - y_span / 2.0, center[1] + y_span / 2.0)
        else:
            pad = float(span.max()) * margin
            ax.set_xlim(mn[0] - pad, mx[0] + pad)
            ax.set_ylim(mn[1] - pad, mx[1] + pad)
        ax.set_aspect("equal")

        if not include_text:
            break
        try:
            ax.figure.canvas.draw()
            renderer = ax.figure.canvas.get_renderer()
            inv = ax.transData.inverted()
            text_points = []
            for text in ax.texts:
                bbox = text.get_window_extent(renderer=renderer).transformed(inv)
                text_points.extend([
                    np.array([bbox.x0, bbox.y0], dtype=float),
                    np.array([bbox.x0, bbox.y1], dtype=float),
                    np.array([bbox.x1, bbox.y0], dtype=float),
                    np.array([bbox.x1, bbox.y1], dtype=float),
                ])
            points = base_points + text_points
        except Exception:
            break
    ax.axis("off")


def legend_handles_resolved(
    has_isolates: bool,
    has_bridge_anchor: bool,
    has_most_connected: bool,
    has_both: bool,
    has_dashed: bool,
    support_values: Sequence[float | int | None] | None = None,
    edge_width_min: float = 0.9,
    edge_width_max: float = 4.5,
) -> list:
    """
    Build legend handles for the resolved-graph figure. Only the encodings
    actually present in the rendered figure are documented.
    """
    handles: list = [
        mlines.Line2D(
            [], [], marker="o", color="w",
            markerfacecolor=OKABE_ITO_PALETTE[0],
            markeredgecolor=NODE_EDGE_COLOUR, markeredgewidth=1.2,
            markersize=11, linestyle="None",
            label="Cell line; node fill colour denotes connected component only",
        ),
    ]
    if has_bridge_anchor:
        handles.append(mlines.Line2D(
            [], [], marker="o", color="w",
            markerfacecolor="white",
            markeredgecolor=ANCHOR_RING_COLOUR, markeredgewidth=3.6,
            markersize=12, linestyle="None",
            label="Bridge-like anchor, within component",
        ))
    if has_most_connected:
        handles.append(mlines.Line2D(
            [], [], marker="s", color="w",
            markerfacecolor="white",
            markeredgecolor=MOST_CONNECTED_COLOUR, markeredgewidth=3.0,
            markersize=12, linestyle="None",
            label="Most-connected anchor, within component",
        ))
    if has_both:
        both_square = mlines.Line2D(
            [], [], marker="s", color="w",
            markerfacecolor="none",
            markeredgecolor=MOST_CONNECTED_COLOUR, markeredgewidth=3.0,
            markersize=13, linestyle="None",
        )
        both_circle = mlines.Line2D(
            [], [], marker="o", color="w",
            markerfacecolor="none",
            markeredgecolor=ANCHOR_RING_COLOUR, markeredgewidth=3.4,
            markersize=11, linestyle="None",
        )
        handles.append((both_square, both_circle, "Both centrality annotations"))
    if has_isolates:
        handles.append(mlines.Line2D(
            [], [], marker="D", color="w",
            markerfacecolor=ISOLATE_COLOUR,
            markeredgecolor=NODE_EDGE_COLOUR, markeredgewidth=1.2,
            markersize=10, linestyle="None",
            label="Isolate (degree 0)",
        ))
    support_handles = support_count_legend_handles(
        support_values or [], edge_width_min, edge_width_max,
    )
    if support_handles:
        handles.extend(support_handles)
    else:
        handles.append(mlines.Line2D(
            [], [], color="black", linewidth=2.6,
            label="Edge width denotes support count",
        ))
    if has_dashed:
        handles.append(mlines.Line2D(
            [], [], color="black", linewidth=1.2, linestyle=(0, (5, 3)),
            label="Single-supported edge (dashed)",
        ))
    return handles


def legend_handles_consensus(
    has_isolates: bool,
    support_values: Sequence[float | int | None] | None = None,
    edge_width_min: float = 0.9,
    edge_width_max: float = 4.5,
) -> list:
    """
    Build legend handles for the consensus-graph figure. The consensus
    network has no single-supported edges by construction (support floor =
    ceil(n_directions / 2)), so dashed encoding is omitted; no anchor marker
    either (the pipeline computes only degree for the consensus network).
    """
    handles: list = [
        mlines.Line2D(
            [], [], marker="o", color="w",
            markerfacecolor=OKABE_ITO_PALETTE[0],
            markeredgecolor=NODE_EDGE_COLOUR, markeredgewidth=1.2,
            markersize=11, linestyle="None",
            label="Cell line; node fill colour denotes connected component only",
        ),
    ]
    if has_isolates:
        handles.append(mlines.Line2D(
            [], [], marker="D", color="w",
            markerfacecolor=ISOLATE_COLOUR,
            markeredgecolor=NODE_EDGE_COLOUR, markeredgewidth=1.2,
            markersize=10, linestyle="None",
            label="Isolate (degree 0)",
        ))
    support_handles = support_count_legend_handles(
        support_values or [], edge_width_min, edge_width_max,
    )
    if support_handles:
        handles.extend(support_handles)
    else:
        handles.append(mlines.Line2D(
            [], [], color="black", linewidth=2.6,
            label="Edge width denotes support count",
        ))
    return handles


def place_legend(
    ax,
    handles: list,
    title: str = "Visual encoding",
    spacious: bool = False,
    labelspacing_override: float | None = None,
    bbox_to_anchor: tuple[float, float] = (1.01, 1.0),
    loc: str = "upper left",
    bbox_transform=None,
    coordinate_system: str = "axes",
):
    """Place the legend with an explicit anchor."""
    labels = [
        item[2] if isinstance(item, tuple) and len(item) == 3 else None
        for item in handles
    ]
    legend_handles = [
        item[:2] if isinstance(item, tuple) and len(item) == 3 else item
        for item in handles
    ]
    if any(label is not None for label in labels):
        resolved_labels = [
            label if label is not None else handle.get_label()
            for handle, label in zip(legend_handles, labels)
        ]
    else:
        resolved_labels = None

    legend_spacing = (
        labelspacing_override
        if labelspacing_override is not None
        else (1.5 if spacious else 0.55)
    )
    handle_text_pad = 1.05 if spacious else 0.8
    border_pad = 0.9 if spacious else 0.4
    handle_height = 1.5 if spacious else 0.7
    legend = ax.legend(
        handles=legend_handles, labels=resolved_labels, loc=loc,
        bbox_to_anchor=bbox_to_anchor, borderaxespad=0.0,
        bbox_transform=bbox_transform,
        fontsize=8, frameon=True, framealpha=0.95,
        title=title, title_fontsize=9,
        labelspacing=legend_spacing, handletextpad=handle_text_pad,
        borderpad=border_pad, handleheight=handle_height,
        handler_map={tuple: HandlerTuple(ndivide=1)},
    )
    legend._codex_bbox_to_anchor_axes = bbox_to_anchor
    legend._codex_bbox_coordinate_system = coordinate_system
    return legend


def save_legend_only(
    handles: list,
    output_prefix: str,
    *,
    title: str = "Visual encoding",
    dpi: int = 300,
    spacious: bool = False,
    labelspacing_override: float | None = None,
    fig_width: float = 7.5,
) -> None:
    """Save a standalone legend using the same handles as the main plot."""
    if not output_prefix:
        return

    import matplotlib.pyplot as plt

    n_items = max(1, len(handles))
    fig_height = max(1.6, 0.42 * n_items + 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)
    ax.axis("off")
    place_legend(
        ax,
        handles,
        title=title,
        spacious=spacious,
        labelspacing_override=labelspacing_override,
        bbox_to_anchor=(0.0, 0.5),
        loc="center left",
    )
    fig.canvas.draw()
    save_kw = dict(
        bbox_inches="tight",
        pad_inches=0.08,
        facecolor="white",
        transparent=False,
    )
    fig.savefig(f"{output_prefix}.png", dpi=dpi, **save_kw)
    fig.savefig(f"{output_prefix}.pdf", **save_kw)
    fig.savefig(f"{output_prefix}.svg", **save_kw)
    plt.close(fig)
    print(f"[OK] Saved legend: {output_prefix}.{{png,pdf,svg}}")


def report_c_label_legend_metrics(fig, ax, context: str) -> dict[str, float]:
    """Report rendered C1/C2 alignment and C2-label-to-legend gap."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labels = {text.get_text(): text for text in ax.texts if text.get_text() in {"C1", "C2"}}
    legend = ax.get_legend()
    if "C1" not in labels or "C2" not in labels or legend is None:
        print(f"[QC] {context} C-label/legend metrics unavailable")
        return {}

    c1_bbox = labels["C1"].get_window_extent(renderer=renderer)
    c2_bbox = labels["C2"].get_window_extent(renderer=renderer)
    legend_bbox = legend.get_window_extent(renderer=renderer)
    delta_top_px = abs(float(c1_bbox.y1) - float(c2_bbox.y1))
    gap_px = float(legend_bbox.x0) - float(c2_bbox.x1)
    px_to_mm = 25.4 / float(fig.dpi)
    metrics = {
        "c1_c2_top_delta_px": delta_top_px,
        "c1_c2_top_delta_mm": delta_top_px * px_to_mm,
        "c2_to_legend_gap_px": gap_px,
        "c2_to_legend_gap_mm": gap_px * px_to_mm,
    }
    print(
        f"[QC] {context} C1/C2 top alignment: "
        f"delta_px={metrics['c1_c2_top_delta_px']:.2f} "
        f"delta_mm={metrics['c1_c2_top_delta_mm']:.2f}"
    )
    print(
        f"[QC] {context} C2_to_legend_gap: "
        f"gap_px={metrics['c2_to_legend_gap_px']:.2f} "
        f"gap_mm={metrics['c2_to_legend_gap_mm']:.2f}"
    )
    return metrics


def report_component_label_legend_gap(
    fig, ax, context: str, component_label: str = "C1"
) -> dict[str, float]:
    """Report rendered component-label-to-legend gap for sparse layouts."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_text = None
    for text in ax.texts:
        if text.get_text() == component_label:
            label_text = text
            break
    legend = ax.get_legend()
    if label_text is None or legend is None:
        print(f"[QC] {context} {component_label}_to_legend_gap unavailable")
        return {}

    label_bbox = label_text.get_window_extent(renderer=renderer)
    legend_bbox = legend.get_window_extent(renderer=renderer)
    gap_px = float(legend_bbox.x0) - float(label_bbox.x1)
    px_to_mm = 25.4 / float(fig.dpi)
    metrics = {
        f"{component_label.lower()}_to_legend_gap_px": gap_px,
        f"{component_label.lower()}_to_legend_gap_mm": gap_px * px_to_mm,
    }
    print(
        f"[QC] {context} {component_label}_to_legend_gap: "
        f"gap_px={gap_px:.2f} gap_mm={gap_px * px_to_mm:.2f}"
    )
    return metrics


def report_isolate_layout_metrics(
    pos: dict[str, np.ndarray],
    isolates: list[str],
    iso_bbox,
    context: str,
) -> dict[str, object]:
    """Report isolate ordering, spacing, and band placement for layout QC."""
    ordered = [node for node in sorted(isolates) if node in pos]
    if not ordered:
        print(f"[QC] {context} isolate layout: no isolates")
        return {"ordered_isolates": []}
    xs = [float(pos[node][0]) for node in ordered]
    ys = [float(pos[node][1]) for node in ordered]
    gaps = [xs[i + 1] - xs[i] for i in range(len(xs) - 1)]
    bbox_text = None
    if iso_bbox is not None:
        bbox_text = tuple(round(float(v), 2) for v in iso_bbox)
    print(
        f"[QC] {context} isolate layout: "
        f"ordered_isolates={ordered}; "
        f"x_gaps={[round(g, 2) for g in gaps]}; "
        f"y_values={[round(y, 2) for y in ys]}; "
        f"iso_bbox={bbox_text}"
    )
    return {
        "ordered_isolates": ordered,
        "x_gaps": gaps,
        "y_values": ys,
        "iso_bbox": iso_bbox,
    }


def align_legend_midline_to_component_label(
    fig,
    ax,
    context: str,
    component_label: str = "C1",
    target_gap_mm: float | None = None,
    max_iter: int = 5,
) -> dict[str, float]:
    # Align the rendered legend midpoint, and optionally its gap, to a component label.
    legend = ax.get_legend()
    if legend is None:
        print(f"[QC] {context} {component_label}_legend_midline unavailable")
        return {}

    x_anchor, y_anchor = getattr(legend, "_codex_bbox_to_anchor_axes", (1.01, 1.0))
    px_to_mm = 25.4 / float(fig.dpi)
    target_gap_px = None if target_gap_mm is None else target_gap_mm / px_to_mm

    for _ in range(max_iter):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        label_text = next((text for text in ax.texts if text.get_text() == component_label), None)
        if label_text is None:
            print(f"[QC] {context} {component_label}_legend_midline unavailable")
            return {}

        label_bbox = label_text.get_window_extent(renderer=renderer)
        legend_bbox = legend.get_window_extent(renderer=renderer)
        axes_bbox = ax.get_window_extent(renderer=renderer)
        label_mid = (float(label_bbox.y0) + float(label_bbox.y1)) / 2.0
        legend_mid = (float(legend_bbox.y0) + float(legend_bbox.y1)) / 2.0
        delta_y_px = label_mid - legend_mid
        delta_x_px = 0.0
        if target_gap_px is not None:
            gap_px = float(legend_bbox.x0) - float(label_bbox.x1)
            delta_x_px = target_gap_px - gap_px

        if abs(delta_y_px) <= 0.5 and abs(delta_x_px) <= 0.5:
            break
        x_anchor += delta_x_px / float(axes_bbox.width)
        y_anchor += delta_y_px / float(axes_bbox.height)
        legend.set_bbox_to_anchor((x_anchor, y_anchor), transform=ax.transAxes)
        legend._codex_bbox_to_anchor_axes = (x_anchor, y_anchor)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_text = next((text for text in ax.texts if text.get_text() == component_label), None)
    if label_text is None:
        print(f"[QC] {context} {component_label}_legend_midline unavailable")
        return {}
    label_bbox = label_text.get_window_extent(renderer=renderer)
    legend_bbox = legend.get_window_extent(renderer=renderer)
    label_mid = (float(label_bbox.y0) + float(label_bbox.y1)) / 2.0
    legend_mid = (float(legend_bbox.y0) + float(legend_bbox.y1)) / 2.0
    mid_delta_px = abs(label_mid - legend_mid)
    gap_px = float(legend_bbox.x0) - float(label_bbox.x1)
    print(
        f"[QC] {context} {component_label}_legend_midline: "
        f"delta_px={mid_delta_px:.2f} delta_mm={mid_delta_px * px_to_mm:.2f}"
    )
    print(
        f"[QC] {context} {component_label}_to_legend_gap: "
        f"gap_px={gap_px:.2f} gap_mm={gap_px * px_to_mm:.2f}"
    )
    return {
        f"{component_label.lower()}_legend_midline_delta_px": mid_delta_px,
        f"{component_label.lower()}_to_legend_gap_px": gap_px,
    }


def _union_display_bboxes(bboxes):
    """Return one display-coordinate bbox covering all non-empty bboxes."""
    from matplotlib.transforms import Bbox

    bboxes = [bbox for bbox in bboxes if bbox is not None]
    if not bboxes:
        return None
    return Bbox.from_extents(
        min(float(b.x0) for b in bboxes),
        min(float(b.y0) for b in bboxes),
        max(float(b.x1) for b in bboxes),
        max(float(b.y1) for b in bboxes),
    )


def component_render_bbox(
    fig,
    ax,
    pos: dict[str, np.ndarray],
    component_nodes: Sequence[str],
    component_texts: Sequence[str] | None = None,
):
    """
    Rendered display bbox for one component, including node glyphs and labels.

    This is used for RBL spacing QC/alignment because its sparse plots need the
    legend positioned relative to the actual C1 content, not a shelf-grid cell.
    """
    from matplotlib.transforms import Bbox

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    node_radius_px = math.sqrt(NODE_SIZE * 1.95 / math.pi) * fig.dpi / 72.0
    wanted_text = set(component_texts or []) | {"C1"}
    bboxes = []
    for node in component_nodes:
        if node not in pos:
            continue
        x, y = ax.transData.transform(pos[node])
        bboxes.append(Bbox.from_extents(
            x - node_radius_px, y - node_radius_px,
            x + node_radius_px, y + node_radius_px,
        ))
    for text in ax.texts:
        if text.get_text() in wanted_text:
            bboxes.append(text.get_window_extent(renderer=renderer))
    return _union_display_bboxes(bboxes)


def display_bbox_from_data_bbox(ax, data_bbox: tuple[float, float, float, float]):
    """Convert a data-coordinate bbox to display coordinates."""
    from matplotlib.transforms import Bbox

    x0, y0, x1, y1 = data_bbox
    p0 = ax.transData.transform((x0, y0))
    p1 = ax.transData.transform((x1, y1))
    return Bbox.from_extents(
        min(float(p0[0]), float(p1[0])),
        min(float(p0[1]), float(p1[1])),
        max(float(p0[0]), float(p1[0])),
        max(float(p0[1]), float(p1[1])),
    )


def pad_isolate_bbox_for_label(
    bbox: tuple[float, float, float, float] | None,
    *,
    left_pad: float = 0.0,
) -> tuple[float, float, float, float] | None:
    """Expand the isolate box leftward so its caption clears the first isolate."""
    if bbox is None:
        return None
    x0, y0, x1, y1 = bbox
    return (float(x0) - max(0.0, float(left_pad)), float(y0), float(x1), float(y1))


def _display_clearance(a, b) -> tuple[float, bool]:
    ax0, ay0, ax1, ay1 = float(a.x0), float(a.y0), float(a.x1), float(a.y1)
    bx0, by0, bx1, by1 = float(b.x0), float(b.y0), float(b.x1), float(b.y1)
    x_gap = max(bx0 - ax1, ax0 - bx1, 0.0)
    y_gap = max(by0 - ay1, ay0 - by1, 0.0)
    overlap = x_gap == 0.0 and y_gap == 0.0
    if overlap:
        return 0.0, True
    if x_gap == 0.0:
        return float(y_gap), False
    if y_gap == 0.0:
        return float(x_gap), False
    return float(math.hypot(x_gap, y_gap)), False


def _node_label_render_bbox(
    fig,
    ax,
    pos: dict[str, np.ndarray],
    nodes: Sequence[str],
    labels_by_node: Mapping[str, str] | None = None,
    extra_texts: Sequence[str] | None = None,
):
    """Rendered bbox covering node glyphs and their labels."""
    from matplotlib.transforms import Bbox

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    node_radius_px = math.sqrt(NODE_SIZE * 1.95 / math.pi) * fig.dpi / 72.0
    label_map = labels_by_node or {node: node for node in nodes}
    wanted_text = {str(label_map.get(node, node)) for node in nodes}
    wanted_text.update(str(text) for text in (extra_texts or []))

    bboxes = []
    for node in nodes:
        if node not in pos:
            continue
        x, y = ax.transData.transform(pos[node])
        bboxes.append(Bbox.from_extents(
            x - node_radius_px, y - node_radius_px,
            x + node_radius_px, y + node_radius_px,
        ))
    for text in ax.texts:
        if text.get_text() in wanted_text:
            bboxes.append(text.get_window_extent(renderer=renderer))
    return _union_display_bboxes(bboxes)


def report_rendered_layout_clearance(
    fig,
    ax,
    G: nx.Graph,
    pos: dict[str, np.ndarray],
    isolates: Iterable[str] | None,
    context: str,
    *,
    labels_by_node: Mapping[str, str] | None = None,
    iso_bbox: tuple[float, float, float, float] | None = None,
) -> dict[str, object]:
    """
    Report rendered display-coordinate clearances for final figure layout QC.

    This is deliberately inspection-only: it does not move artists or change
    graph data. Units are display pixels plus millimetres at the figure DPI.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_to_mm = 25.4 / float(fig.dpi)

    bboxes: dict[str, object] = {}
    components = ordered_non_isolate_components(G, isolates=isolates)
    for idx, comp in enumerate(components, start=1):
        nodes = sorted(n for n in comp if n in pos)
        bbox = _node_label_render_bbox(
            fig, ax, pos, nodes,
            labels_by_node=labels_by_node,
            extra_texts=[f"C{idx}"],
        )
        if bbox is not None:
            bboxes[f"C{idx}"] = bbox

    if isolates:
        iso_nodes = sorted(n for n in isolates if n in pos)
        bbox = _node_label_render_bbox(
            fig, ax, pos, iso_nodes,
            labels_by_node=labels_by_node,
        )
        if bbox is not None:
            bboxes["isolate_nodes"] = bbox
    if iso_bbox is not None:
        bboxes["isolate_box"] = display_bbox_from_data_bbox(ax, iso_bbox)
    isolate_label = next((text for text in ax.texts if text.get_text() == "Isolates (degree 0)"), None)
    if isolate_label is not None:
        bboxes["isolate_label"] = isolate_label.get_window_extent(renderer=renderer)
    legend = ax.get_legend()
    if legend is not None:
        bboxes["legend"] = legend.get_window_extent(renderer=renderer)
    if ax.title is not None and ax.title.get_text():
        bboxes["title"] = ax.title.get_window_extent(renderer=renderer)
    if fig.texts:
        caption_bbox = _union_display_bboxes([
            text.get_window_extent(renderer=renderer)
            for text in fig.texts
        ])
        if caption_bbox is not None:
            bboxes["footnote"] = caption_bbox

    pair_specs: list[tuple[str, str]] = []
    for a, b in (("C1", "C2"), ("C1", "C3"), ("C1", "C4"), ("C2", "legend"),
                 ("legend", "title"), ("isolate_box", "footnote"),
                 ("isolate_label", "isolate_nodes")):
        if a in bboxes and b in bboxes:
            pair_specs.append((a, b))
    for comp_label in sorted(k for k in bboxes if k.startswith("C")):
        for target in ("isolate_box", "title"):
            if target in bboxes:
                pair_specs.append((comp_label, target))

    min_pair = None
    any_overlap = False
    seen: set[tuple[str, str]] = set()
    for a, b in pair_specs:
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        clearance, overlap = _display_clearance(bboxes[a], bboxes[b])
        any_overlap = any_overlap or overlap
        if min_pair is None or clearance < float(min_pair[2]):
            min_pair = (a, b, clearance, overlap)
        print(
            f"[QC] {context} rendered clearance {a}_vs_{b}: "
            f"clearance_px={clearance:.2f} clearance_mm={clearance * px_to_mm:.2f} "
            f"overlap={'yes' if overlap else 'no'}"
        )

    fig_bbox = fig.bbox
    min_boundary = math.inf
    min_boundary_label = ""
    any_outside = False
    for label, bbox in bboxes.items():
        boundary = min(
            float(bbox.x0) - float(fig_bbox.x0),
            float(fig_bbox.x1) - float(bbox.x1),
            float(bbox.y0) - float(fig_bbox.y0),
            float(fig_bbox.y1) - float(bbox.y1),
        )
        min_boundary = min(min_boundary, boundary)
        if boundary == min_boundary:
            min_boundary_label = label
        any_outside = any_outside or boundary < 0.0
        print(
            f"[QC] {context} rendered boundary {label}: "
            f"min_px={boundary:.2f} min_mm={boundary * px_to_mm:.2f} "
            f"outside={'yes' if boundary < 0.0 else 'no'}"
        )

    if min_pair is not None:
        print(
            f"[QC] {context} rendered clearance summary: "
            f"min_pair={min_pair[0]}_vs_{min_pair[1]} "
            f"min_clearance_px={float(min_pair[2]):.2f} "
            f"min_clearance_mm={float(min_pair[2]) * px_to_mm:.2f} "
            f"overlap={'yes' if any_overlap else 'no'} "
            f"min_boundary={min_boundary_label} "
            f"min_boundary_px={min_boundary:.2f} "
            f"outside={'yes' if any_outside else 'no'}"
        )
    return {
        "bboxes": bboxes,
        "any_overlap": any_overlap,
        "any_outside": any_outside,
        "min_boundary_label": min_boundary_label,
        "min_boundary_px": min_boundary,
    }


def align_legend_to_component_box(
    fig,
    ax,
    pos: dict[str, np.ndarray],
    component_nodes: Sequence[str],
    context: str,
    *,
    component_texts: Sequence[str] | None = None,
    iso_bbox: tuple[float, float, float, float] | None = None,
    target_gap_frac: float = RBL_SINGLE_COMPONENT_LEGEND_GAP_FRAC,
    max_iter: int = 6,
) -> dict[str, float]:
    """
    Align legend beside the rendered C1 component box for RBL sparse figures.

    The target horizontal gap is a fraction of the current figure width. The
    vertical midpoints of C1 and the legend are aligned so they share a row.
    """
    legend = ax.get_legend()
    if legend is None:
        print(f"[QC] {context} component_to_legend unavailable")
        return {}

    x_anchor, y_anchor = getattr(legend, "_codex_bbox_to_anchor_axes", (1.01, 1.0))
    target_gap_px = float(fig.bbox.width) * float(target_gap_frac)

    for _ in range(max_iter):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        component_bbox = component_render_bbox(
            fig, ax, pos, component_nodes, component_texts=component_texts,
        )
        if component_bbox is None:
            print(f"[QC] {context} component_to_legend unavailable")
            return {}
        legend_bbox = legend.get_window_extent(renderer=renderer)
        axes_bbox = ax.get_window_extent(renderer=renderer)
        gap_px = float(legend_bbox.x0) - float(component_bbox.x1)
        component_mid = (float(component_bbox.y0) + float(component_bbox.y1)) / 2.0
        legend_mid = (float(legend_bbox.y0) + float(legend_bbox.y1)) / 2.0
        delta_x_px = target_gap_px - gap_px
        delta_y_px = component_mid - legend_mid
        if abs(delta_x_px) <= 0.5 and abs(delta_y_px) <= 0.5:
            break
        x_anchor += delta_x_px / float(axes_bbox.width)
        y_anchor += delta_y_px / float(axes_bbox.height)
        legend.set_bbox_to_anchor((x_anchor, y_anchor), transform=ax.transAxes)
        legend._codex_bbox_to_anchor_axes = (x_anchor, y_anchor)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    component_bbox = component_render_bbox(
        fig, ax, pos, component_nodes, component_texts=component_texts,
    )
    legend_bbox = legend.get_window_extent(renderer=renderer)
    gap_px = float(legend_bbox.x0) - float(component_bbox.x1)
    mid_delta_px = abs(
        ((float(component_bbox.y0) + float(component_bbox.y1)) / 2.0)
        - ((float(legend_bbox.y0) + float(legend_bbox.y1)) / 2.0)
    )
    px_to_mm = 25.4 / float(fig.dpi)
    metrics = {
        "component_to_legend_gap_px": gap_px,
        "component_to_legend_gap_mm": gap_px * px_to_mm,
        "component_to_legend_gap_frac": gap_px / float(fig.bbox.width),
        "component_legend_midline_delta_px": mid_delta_px,
        "component_legend_midline_delta_mm": mid_delta_px * px_to_mm,
    }
    if iso_bbox is not None:
        iso_display = display_bbox_from_data_bbox(ax, iso_bbox)
        iso_gap_px = float(component_bbox.y0) - float(iso_display.y1)
        metrics["component_to_isolate_gap_px"] = iso_gap_px
        metrics["component_to_isolate_gap_mm"] = iso_gap_px * px_to_mm
        metrics["component_to_isolate_gap_frac_height"] = iso_gap_px / float(fig.bbox.height)
    print(
        f"[QC] {context} component_to_legend_gap: "
        f"gap_px={metrics['component_to_legend_gap_px']:.2f} "
        f"gap_mm={metrics['component_to_legend_gap_mm']:.2f} "
        f"gap_frac={metrics['component_to_legend_gap_frac']:.3f}"
    )
    print(
        f"[QC] {context} component_legend_midline: "
        f"delta_px={metrics['component_legend_midline_delta_px']:.2f} "
        f"delta_mm={metrics['component_legend_midline_delta_mm']:.2f}"
    )
    if iso_bbox is not None:
        print(
            f"[QC] {context} component_to_isolate_gap: "
            f"gap_px={metrics['component_to_isolate_gap_px']:.2f} "
            f"gap_mm={metrics['component_to_isolate_gap_mm']:.2f} "
            f"gap_frac_height={metrics['component_to_isolate_gap_frac_height']:.3f}"
        )
    return metrics
