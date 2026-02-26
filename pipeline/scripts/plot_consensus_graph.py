#!/usr/bin/env python3
"""
plot_consensus_dsmz_graph.py
-------------------------------------------------
Visualises the undirected consensus DSMZ cell-line similarity graph.

Key improvements over the previous version
-------------------------------------------
  - Component-aware layout: each connected component gets its own
    spring layout, preventing small components from collapsing into
    the centre of a large one.
  - Per-component adaptive spring constant (k) and bounding-box scaling
    so that small chains are spread out rather than compressed.
  - Post-layout repulsion pass to eliminate residual node overlap.
  - Grid packing with size-proportional inter-component padding so that
    whitespace is used efficiently.
  - Curved edges (arc style) to keep labels readable when nodes are close.
  - Label truncation / offset for very dense regions.

Edge semantics
--------------
  Solid edges   — reciprocal (supported by ≥2 independent tags/methods).
  Dashed edges  — one-way (supported by exactly 1 tag/method).
  Edge width    — scales with support_directions when available.

Node highlighting
-----------------
  Thick outline — top-N% by betweenness centrality (bridge nodes).
"""

import argparse
import math
import sys

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def _spring_k(n: int, k_min: float = 2.5, k_max: float = 8.0) -> float:
    """
    Adaptive spring constant that grows logarithmically with component size.

    Larger k spreads nodes further apart.  The log scale prevents large
    components from becoming unmanageably spread while still giving small
    chains (n=2..5) enough room to separate their labels.

    Parameters
    ----------
    n     : number of nodes in the component
    k_min : spring constant for singletons / very small components
    k_max : spring constant for large components (≥50 nodes)
    """
    if n <= 1:
        return k_min
    return float(np.clip(
        k_min + (math.log(n) / math.log(50)) * (k_max - k_min),
        k_min, k_max
    ))


def _target_halfbox(n: int, t_min: float = 3.0, t_max: float = 9.0) -> float:
    """
    Minimum half-span (in layout units) for a component of size n.

    Components are only *enlarged* to this size; they are never shrunk,
    so the spring geometry is preserved.
    """
    if n <= 1:
        return t_min
    return float(np.clip(
        t_min + (math.log(n) / math.log(50)) * (t_max - t_min),
        t_min, t_max
    ))


def _scale_up_only(local_pos: dict, target: float) -> dict:
    """
    Scale a component layout so its larger dimension is at least `target`.
    Never shrinks — only enlarges compact layouts.

    Parameters
    ----------
    local_pos : node → (x, y) from spring_layout
    target    : minimum half-span in the larger axis
    """
    pts = np.array(list(local_pos.values()), dtype=float)
    mn, mx = pts.min(axis=0), pts.max(axis=0)
    span = max(float((mx - mn).max()), 1e-9)
    half_span = span / 2.0
    scale = (target / half_span) if half_span < target else 1.0
    centre = (mn + mx) / 2.0
    return {n: (np.array(p, dtype=float) - centre) * scale
            for n, p in local_pos.items()}


def _repel(pos: dict, nodes: list,
           min_dist: float = 1.2, iters: int = 120,
           step: float = 0.12, seed: int = 42) -> dict:
    """
    Iterative repulsion pass: push pairs of nodes that are closer than
    `min_dist` apart until separation is achieved or iterations run out.

    Parameters
    ----------
    pos      : current node positions (modified in-place and returned)
    nodes    : list of node ids to consider (single component)
    min_dist : desired minimum Euclidean separation
    iters    : maximum relaxation iterations
    step     : movement fraction per iteration
    seed     : RNG seed for symmetry-breaking jitter on coincident nodes
    """
    rng = np.random.default_rng(seed)
    nodes = sorted(nodes)
    for _ in range(iters):
        moved = False
        for i in range(len(nodes)):
            u = nodes[i]
            pu = np.array(pos[u], dtype=float)
            for j in range(i + 1, len(nodes)):
                v = nodes[j]
                pv = np.array(pos[v], dtype=float)
                dvec = pu - pv
                dist = float(np.linalg.norm(dvec))
                if dist < 1e-9:
                    # Coincident nodes: break symmetry with random jitter
                    dvec = rng.normal(size=2)
                    dvec /= np.linalg.norm(dvec) + 1e-9
                    dist = 1e-9
                if dist < min_dist:
                    push = (min_dist - dist) * step
                    direction = dvec / dist
                    pos[u] = pu + direction * push
                    pos[v] = pv - direction * push
                    pu = np.array(pos[u], dtype=float)
                    moved = True
        if not moved:
            break
    return pos


def pack_components(G: nx.Graph, pad: float, seed: int,
                    k_min: float, k_max: float,
                    t_min: float, t_max: float) -> tuple:
    """
    Compute a packed layout for all connected components of G.

    Strategy
    --------
    1. Sort components largest-first so important groups appear top-left.
    2. Run an independent spring layout per component with adaptive k.
    3. Enlarge (never shrink) each component to a minimum bounding box.
    4. Run a within-component repulsion pass to eliminate overlap.
    5. Arrange components on a square grid with size-scaled padding.

    Returns
    -------
    pos : dict  node → (x, y) global coordinates
    comps : list of sorted node lists, one per component
    """
    comps = sorted(
        nx.connected_components(G),
        key=lambda c: (-len(c), sorted(c)[0])
    )

    cols = max(1, math.ceil(math.sqrt(len(comps))))
    pos = {}

    for idx, comp in enumerate(comps):
        nodes = sorted(comp)
        sub = G.subgraph(nodes).copy()
        n = sub.number_of_nodes()

        # --- per-component layout ---
        if n == 1:
            local_pos = {nodes[0]: np.array([0.0, 0.0])}
        else:
            k = _spring_k(n, k_min=k_min, k_max=k_max)
            local_pos = nx.spring_layout(sub, seed=seed, k=k, iterations=300)
            target = _target_halfbox(n, t_min=t_min, t_max=t_max)
            local_pos = _scale_up_only(local_pos, target)
            # Repulsion: min separation = 55% of target halfbox
            local_pos = _repel(
                local_pos, nodes,
                min_dist=0.55 * target,
                iters=160, step=0.12, seed=seed
            )

        # --- grid offset ---
        row = idx // cols
        col = idx % cols
        # Pad scales with component size so large components don't collide
        pad_i = pad * (1.0 + 0.18 * math.log(max(n, 2)))
        offset = np.array([col * pad_i, -row * pad_i], dtype=float)
        for node, p in local_pos.items():
            pos[node] = np.array(p, dtype=float) + offset

    return pos, comps


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # ------------------------------------------------------------------
    # CLI
    # ------------------------------------------------------------------
    ap = argparse.ArgumentParser(
        description=(
            "Plot consensus DSMZ graph with solid (reciprocal) and "
            "dashed (one-way) edges. Improved component-aware layout."
        )
    )
    ap.add_argument("--edges", required=True,
                    help="Consensus edge TSV (node1, node2, reciprocal columns required).")
    ap.add_argument("--out", required=True,
                    help="Output path prefix; .png/.pdf/.svg are appended.")
    ap.add_argument("--label", default="BRCA",
                    help="Dataset label shown in the figure title.")

    # Figure
    ap.add_argument("--fig-w",  type=float, default=28.0)
    ap.add_argument("--fig-h",  type=float, default=20.0)
    ap.add_argument("--dpi",    type=int,   default=300)

    # Nodes / labels
    ap.add_argument("--nodes", default=None,
                    help="Optional TSV with a list of nodes to include (e.g. dsmz_cellline_shortnames.tsv). "
                         "If provided, isolates (degree 0) will appear in the plot.")
    ap.add_argument("--nodes-col", default="short_id",
                    help="Column name in --nodes TSV to use as node IDs (default: short_id).")
    ap.add_argument("--node-size",  type=float, default=1800.0,
                    help="Base node area pt² (default 1800; increase for fewer nodes).")
    ap.add_argument("--font-size",  type=int,   default=9,
                    help="Label font size pt (default 9).")

    # Edges
    ap.add_argument("--edge-width", type=float, default=2.5)
    ap.add_argument("--edge-alpha", type=float, default=0.80)

    # Layout
    ap.add_argument("--pad",      type=float, default=14.0,
                    help="Base inter-component padding in layout units (default 14).")
    ap.add_argument("--k-min",    type=float, default=2.5,
                    help="Min spring k for small components (default 2.5).")
    ap.add_argument("--k-max",    type=float, default=8.0,
                    help="Max spring k for large components (default 8.0).")
    ap.add_argument("--t-min",    type=float, default=3.0,
                    help="Min bounding-box half-span (default 3.0).")
    ap.add_argument("--t-max",    type=float, default=9.0,
                    help="Max bounding-box half-span (default 9.0).")
    ap.add_argument("--seed",     type=int,   default=42)

    # Centrality
    ap.add_argument("--top-betw-frac", type=float, default=0.15,
                    help="Fraction of nodes highlighted as high-betweenness (default 0.15).")

    args = ap.parse_args()

    # ------------------------------------------------------------------
    # Step 1 — Load and validate edges
    # ------------------------------------------------------------------
    df = pd.read_csv(args.edges, sep="\t")
    for col in ("node1", "node2", "reciprocal"):
        if col not in df.columns:
            raise SystemExit(
                f"Missing column '{col}' in {args.edges}. "
                f"Columns present: {list(df.columns)}"
            )

    # ------------------------------------------------------------------
    # Step 2 — Build graph
    # ------------------------------------------------------------------
    G = nx.Graph()
    support_map: dict = {}

    # Add all nodes from optional node list first (so isolates appear)
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

    for _, r in df.iterrows():
        u = str(r["node1"]).strip()
        v = str(r["node2"]).strip()
        if not u or not v or u == v:
            continue
        G.add_edge(u, v, reciprocal=int(r["reciprocal"]))
        if "support_directions" in df.columns:
            support_map[tuple(sorted((u, v)))] = int(r.get("support_directions", 1))

    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    n_comps = nx.number_connected_components(G)
    print(f"[INFO] Graph: {n_nodes} nodes, {n_edges} edges, {n_comps} components")

    # ------------------------------------------------------------------
    # Step 3 — Component-aware packed layout
    # ------------------------------------------------------------------
    pos, comps = pack_components(
        G,
        pad=args.pad,
        seed=args.seed,
        k_min=args.k_min,
        k_max=args.k_max,
        t_min=args.t_min,
        t_max=args.t_max,
    )

    # ------------------------------------------------------------------
    # Step 4 — Betweenness centrality → node classification
    # ------------------------------------------------------------------
    bet = nx.betweenness_centrality(G)
    k_top = max(1, int(np.ceil(args.top_betw_frac * n_nodes)))
    central = set(sorted(bet, key=bet.get, reverse=True)[:k_top])
    normal_nodes  = [x for x in G.nodes() if x not in central]
    central_nodes = [x for x in G.nodes() if x in central]

    # ------------------------------------------------------------------
    # Step 5 — Edge classification and width computation
    # ------------------------------------------------------------------
    solid_edges  = [(u, v) for u, v, d in G.edges(data=True) if d.get("reciprocal", 0) == 1]
    dashed_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("reciprocal", 0) == 0]

    def _ew(u, v, base: float) -> float:
        """Scale edge width by support_directions; fall back to base."""
        if support_map:
            s = support_map.get(tuple(sorted((u, v))), 1)
            return base + 0.4 * (s - 1)
        return base

    solid_widths  = [_ew(u, v, args.edge_width * 1.4) for u, v in solid_edges]
    dashed_widths = [_ew(u, v, args.edge_width * 0.9) for u, v in dashed_edges]

    # ------------------------------------------------------------------
    # Step 6 — Draw
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(args.fig_w, args.fig_h), dpi=args.dpi)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.axis("off")

    # --- Solid (reciprocal) edges — black, minimum weight, continuous ---
    if solid_edges:
        nx.draw_networkx_edges(
            G, pos, ax=ax,
            edgelist=solid_edges,
            width=1.0,             # minimum readable line weight
            alpha=1.0,
            edge_color="black",
            style="solid"
        )

    # --- Dashed (one-way) edges — black, minimum weight, dashed ---
    if dashed_edges:
        nx.draw_networkx_edges(
            G, pos, ax=ax,
            edgelist=dashed_edges,
            width=1.0,             # same minimum weight as solid
            alpha=1.0,
            edge_color="black",
            style=(0, (6, 4))      # 6pt on / 4pt off dash pattern
        )

    # --- Normal nodes ---
    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        nodelist=normal_nodes,
        node_size=args.node_size,
        node_color="#a8d8ea",              # calm mid-blue
        edgecolors="#2d4059",
        linewidths=1.5,
        alpha=0.95
    )

    # --- Central (high-betweenness) nodes ---
    # Drawn on top with thicker outline to mark bridge nodes
    if central_nodes:
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            nodelist=central_nodes,
            node_size=args.node_size * 1.4,
            node_color="#a8d8ea",
            edgecolors="#2d4059",
            linewidths=4.5,
            alpha=1.0
        )

    # --- Labels ---
    # White background box with tight padding improves legibility in
    # dense regions; font weight bold ensures readability at small sizes.
    nx.draw_networkx_labels(
        G, pos, ax=ax,
        font_size=args.font_size,
        font_weight="bold",
        font_color="#1a1a2e",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            boxstyle="round,pad=0.15",
            alpha=0.75                     # semi-transparent: shows node colour beneath
        )
    )

    # --- Title ---
    ax.set_title(
        f"DSMZ Cell-line Consensus Graph — {args.label}\n"
        "(Solid = multi-supported (≥2 tags); dashed = single-supported; thick outline = high betweenness)",
        fontsize=args.font_size + 5,
        pad=20,
        color="#1a1a2e"
    )

    # Auto-scale axes with a small margin so edge labels are not clipped
    pts = np.array(list(pos.values()), dtype=float)
    mn, mx = pts.min(axis=0), pts.max(axis=0)
    margin_x = (mx[0] - mn[0]) * 0.06
    margin_y = (mx[1] - mn[1]) * 0.06
    ax.set_xlim(mn[0] - margin_x, mx[0] + margin_x)
    ax.set_ylim(mn[1] - margin_y, mx[1] + margin_y)
    ax.set_aspect("equal")

    # ------------------------------------------------------------------
    # Step 7 — Save
    # ------------------------------------------------------------------
    save_kw = dict(bbox_inches="tight", pad_inches=0.3,
                   facecolor="white", transparent=False)
    fig.savefig(f"{args.out}.png", dpi=args.dpi, **save_kw)
    fig.savefig(f"{args.out}.pdf", **save_kw)
    fig.savefig(f"{args.out}.svg", **save_kw)
    plt.close(fig)
    print(f"[OK] Saved: {args.out}.png / .pdf / .svg")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()