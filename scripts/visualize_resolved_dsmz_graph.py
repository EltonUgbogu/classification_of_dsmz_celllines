#!/usr/bin/env python3
"""
=============================================================================
visualise_resolved_dsmz_graph.py
=============================================================================

DESCRIPTION
-----------
This script generates visualisations of cell line
similarity networks derived from tumour neighbourhood analysis. The
visualisation employs graph-theoretic principles to reveal the structure
of relationships between DSMZ cell lines based on their transcriptomic
similarity to tumour samples.

BACKGROUND: CELL LINE SIMILARITY NETWORKS
-----------------------------------------
In cancer genomics research, cell lines serve as in vitro models for
studying tumour biology. However, not all cell lines faithfully recapitulate
the molecular characteristics of their purported tissue of origin. By
constructing similarity networks based on gene expression profiles, one
can identify:

    (1) Cell lines that cluster together, suggesting shared molecular
        features and potentially similar biological behaviour.

    (2) Isolated cell lines that lack strong similarity to any other
        cell line, which may represent unique molecular subtypes or
        indicate quality/identity concerns.

    (3) Central ``hub'' cell lines that bridge multiple groups, which
        may represent transitional phenotypes or broadly representative
        models.

GRAPH-THEORETIC CONCEPTS
------------------------
The script employs several fundamental concepts from graph theory and
network analysis:

A connected component is a maximal subgraph in which every pair of nodes
is connected by a path. In the context of cell line networks:

    - Each component represents a group of cell lines with mutual
      similarity relationships (direct or transitive).
    - Isolated nodes (components of size 1) represent cell lines with
      no detected similarity to others under the applied thresholds.
    - The number and size distribution of components reveals the overall
      fragmentation or cohesion of the cell line population.

Betweenness centrality quantifies the importance of a node as a bridge
between other nodes. For a node $v$, it is defined as:

In biological networks, high betweenness centrality identifies:
    - Cell lines that connect otherwise separate clusters
    - Potential ``bridge'' phenotypes between molecular subtypes
    - Key nodes whose removal would fragment the network

INPUT REQUIREMENTS
------------------
The script requires a TSV file (resolved_dsmz_neighbors.tsv) containing:

    - cell_line:        Cell line identifier (defines the node universe)
    - final_neighbours:  Semicolon-separated list of neighbour cell lines

STRICT INPUT RULES:
    - Nodes are derived ONLY from the cell_line column
    - Edges are derived ONLY from the final_neighbours column
    - Neighbours not present in the cell_line column are ignored
    - Self-loops (cell line listed as its own neighbour) are ignored

This strict approach ensures that the visualisation accurately represents
the resolved neighbour relationships without introducing spurious nodes
from incomplete or inconsistent data.

OUTPUT FILES
------------
The script produces:

    1. {output_prefix}.png / .pdf
       Network visualisation with:
       - Packed component layout (one spring layout per component)
       - Adaptive spacing based on component size
       - Central nodes highlighted with thicker outlines

    2. dsmz_cellline_graph_edges.tsv
       Edge list with aggregated support metrics:
       - node1, node2: Edge endpoints (alphabetically ordered)
       - support_directions: Number of feature-distance directions
         supporting this edge
       - support_weight_mean/sum/max: Aggregated similarity scores

    3. dsmz_cellline_graph_node_stats.tsv
       Per-node statistics:
       - degree: Number of edges incident to the node
       - betweenness: Betweenness centrality score
       - component: Component membership index
       - is_isolate: Boolean flag for isolated nodes
       - is_central: Boolean flag for highlighted central nodes

USAGE
-----
    python visualise_resolved_dsmz_graph.py \\
        resolved_dsmz_neighbors.tsv \\
        output/dsmz_graph \\
        BRCA \\
        --pad 12.0 \\
        --k-min 1.8 --k-max 5.0 \\
        --t-min 3.4 --t-max 8.0

DEPENDENCIES
------------
Python packages: pandas, networkx, matplotlib, numpy

=============================================================================
"""

import argparse
import datetime
import sys
import math
import re
from collections import defaultdict
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np

# Ensure underscores render literally (avoid mathtext subscripts)
mpl.rcParams["text.usetex"] = False
mpl.rcParams["mathtext.default"] = "regular"

# Colour blind safe Wong palette; last slot reserved for isolates
CB_PALETTE = [
    "#4477AA", "#EE6677", "#228833", "#CCBB44",
    "#66CCEE", "#AA3377", "#BB5566", "#DDAA33",
]
ISOLATE_COLOUR = "#BBBBBB"


# =============================================================================
# Utility Functions
# =============================================================================

def parse_neighbours(s: str):
    """
    Parses a semicolon-separated neighbour string into a list of identifiers.

    In the resolved neighbours file, each cell line's neighbours are stored
    as a single string with semicolon delimiters (e.g., "MCF7;T47D;BT474").
    This function handles edge cases including:

        - None or NaN values (empty neighbour lists)
        - Empty strings (isolated cell lines)
        - Whitespace around identifiers

    Parameters
    ----------
    s : str or None
        Semicolon-separated string of neighbour identifiers.

    Returns
    -------
    list of str
        List of cleaned neighbour identifiers (empty list if no neighbours).
    """
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return []

    s = str(s).strip()
    if not s:
        return []

    return [x.strip() for x in s.split(";") if x.strip()]


def build_shortname_map(node_ids):
    """
    Naming convention:
      - base: RBL_<num> extracted from long_id (or WERI_RB1, Y_79)
      - if base occurs once -> short_id = base
      - if base occurs multiple times -> short_id = base_<lastToken>
        where lastToken is the trailing underscore token of the long_id (must be digits)
    Special cases:
      - WERI_RB1 -> WERI_RB1
      - Y_79     -> Y_79
    """
    node_ids = sorted(set(node_ids))

    base_for = {}
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

    short_map = {}

    for base, longs in by_base.items():
        longs = sorted(longs)

        if len(longs) == 1:
            short_map[longs[0]] = base
            continue

        for long_id in longs:
            last = str(long_id).split("_")[-1]
            if not last.isdigit():
                raise ValueError(
                    f"[shortname] Duplicate base {base} but long_id has no numeric trailing token: {long_id}"
                )
            short_map[long_id] = f"{base}_{last}"

    return short_map


def edges_from_resolved_neighbors(df, nb_col, node_set):
    """
    Derive edges directly from resolved neighbors column.
    Enforces strict node universe filtering (only nodes in node_set).
    """
    edges = []
    seen = set()
    for _, row in df.iterrows():
        src = str(row.get("cell_line", "")).strip()
        if not src:
            continue

        neigh_str = str(row.get(nb_col, "") or "")
        for tgt in parse_neighbours(neigh_str):
            if not tgt or tgt.upper() in {"NA", "NAN", "."}:
                continue
            if tgt == src:
                continue
            if tgt not in node_set:
                continue

            node1, node2 = sorted([src, tgt])
            if (node1, node2) in seen:
                continue
            seen.add((node1, node2))
            edges.append({
                "node1": node1,
                "node2": node2,
                "support_directions": "",
                "support_weight_mean": "",
                "support_weight_sum": "",
                "support_weight_max": "",
                "methods_union": "fallback=resolved_neighbors"
            })
    return edges


def dynamic_spring_k(n: int, k_min: float, k_max: float) -> float:
    """
    Computes an adaptive spring constant based on component size.

    The spring constant $k$ in force-directed layouts controls the natural
    length of edges. For larger components, a larger $k$ produces more
    spread-out layouts, preventing node overlap and improving readability.

    The scaling follows a logarithmic function to provide smooth growth:

    \begin{equation}
        k = k_{\min} + \frac{\log(n)}{\log(50)} \cdot (k_{\max} - k_{\min})
    \end{equation}

    The logarithmic scaling ensures that:
        - Small components (n < 10) remain compact
        - Medium components (10 < n < 50) receive moderate spacing
        - Large components (n > 50) approach maximum spacing

    Parameters
    ----------
    n : int
        Number of nodes in the component.
    k_min : float
        Minimum spring constant (for smallest components).
    k_max : float
        Maximum spring constant (for largest components).

    Returns
    -------
    float
        Adaptive spring constant for the component.
    """
    if n <= 1:
        return k_min

    k = k_min + (math.log(n) / math.log(50)) * (k_max - k_min)
    return max(k_min, min(k, k_max))


def dynamic_target_box(n: int, t_min: float, t_max: float) -> float:
    """
    Computes an adaptive bounding box size for component normalisation.

    After computing the spring layout, node positions are normalised to
    fit within a standardised bounding box. Larger components receive
    larger boxes to accommodate more nodes and labels without overlap.

    The scaling follows the same logarithmic function as the spring
    constant to maintain proportional spacing:

    \begin{equation}
        t = t_{\min} + \frac{\log(n)}{\log(50)} \cdot (t_{\max} - t_{\min})
    \end{equation}

    Parameters
    ----------
    n : int
        Number of nodes in the component.
    t_min : float
        Minimum half-box size (component spans $[-t, t]$ in each dimension).
    t_max : float
        Maximum half-box size.

    Returns
    -------
    float
        Adaptive half-box size for the component.
    """
    if n <= 1:
        return t_min

    t = t_min + (math.log(n) / math.log(50)) * (t_max - t_min)
    return max(t_min, min(t, t_max))


def normalise_layout(local_pos, target: float):
    """
    Normalises component layout coordinates to a standardised bounding box.

    Force-directed layouts produce coordinates in arbitrary ranges depending
    on the initial positions and convergence behaviour. This function
    rescales coordinates to fit within $[-\text{target}, \text{target}]$
    in both dimensions, ensuring consistent visual sizing across components.

    The normalisation process:
        1. Computes the bounding box of current positions
        2. Rescales coordinates to $[0, 1]$
        3. Shifts and scales to $[-\text{target}, \text{target}]$

    Parameters
    ----------
    local_pos : dict
        Dictionary mapping node identifiers to (x, y) coordinate arrays.
    target : float
        Half-size of the target bounding box.

    Returns
    -------
    dict
        Dictionary mapping node identifiers to normalised coordinates.
    """
    pts = np.array(list(local_pos.values()), dtype=float)

    minxy = pts.min(axis=0)
    maxxy = pts.max(axis=0)

    span = np.maximum(maxxy - minxy, 1e-9)

    out = {}
    for node, p in local_pos.items():
        q = (np.array(p) - minxy) / span
        q = (q - 0.5) * (2 * target)
        out[node] = q

    return out


def scale_layout_to_min_box(local_pos, target: float):
    """
    Ensure component is at least target-sized, but NEVER shrink it.
    Preserves spring_layout geometry; scales up only when too compact.

    Unlike normalise_layout(), this preserves the natural spacing from
    spring_layout and only enlarges components that are too small.

    Parameters
    ----------
    local_pos : dict
        Dictionary mapping node identifiers to (x, y) coordinate arrays.
    target : float
        Minimum half-span (larger dimension) for the component.

    Returns
    -------
    dict
        Dictionary mapping node identifiers to scaled coordinates (centered).
    """
    pts = np.array(list(local_pos.values()), dtype=float)
    minxy = pts.min(axis=0)
    maxxy = pts.max(axis=0)
    span = np.maximum(maxxy - minxy, 1e-9)

    half_span = 0.5 * max(span[0], span[1])

    s = (target / half_span) if half_span < target else 1.0

    center = (minxy + maxxy) / 2.0
    return {node: (np.array(p, dtype=float) - center) * s
            for node, p in local_pos.items()}


def normalize_pos_to_square(pos):
    """
    Normalises global layout coordinates to a unit square.

    This rescales packed component coordinates so the overall layout
    fills a square canvas without distorting topology.
    """
    pts = np.array(list(pos.values()), dtype=float)

    minxy = pts.min(axis=0)
    maxxy = pts.max(axis=0)
    span = max(maxxy[0] - minxy[0], maxxy[1] - minxy[1], 1e-9)

    out = {}
    for node, p in pos.items():
        q = (np.array(p) - minxy) / span
        out[node] = q

    return out


# =============================================================================
# Automatic Node Separation (Collision Avoidance)
# =============================================================================

def repel_close_nodes(pos, G, min_dist=0.55, iters=80, step=0.06, seed=42):
    """
    Generic post-layout repulsion to reduce node/label collisions.
    Works per connected component. Deterministic given seed.

    Parameters
    ----------
    pos : dict
        Dictionary mapping node identifiers to (x, y) coordinates.
    G : networkx.Graph
        Input graph (used to identify connected components).
    min_dist : float
        Desired minimum Euclidean distance between nodes (in layout units).
    iters : int
        Number of relaxation iterations.
    step : float
        Movement step size per iteration.
    seed : int
        Random seed for deterministic jitter when nodes have identical coordinates.

    Returns
    -------
    dict
        Updated position dictionary with nodes pushed apart.
    """
    rng = np.random.default_rng(seed)

    for comp in nx.connected_components(G):
        nodes = sorted(comp)
        if len(nodes) <= 1:
            continue

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
                        jitter = rng.normal(0, 1, size=2)
                        jitter = jitter / (np.linalg.norm(jitter) + 1e-9)
                        dvec = jitter
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


# =============================================================================
# Component Layout and Packing
# =============================================================================

def pack_components(G: nx.Graph, pad: float, seed: int,
                    k_min: float, k_max: float,
                    t_min: float, t_max: float,
                    center_top_n: int):
    """
    Computes a packed layout for all connected components of the graph.

    This function implements a multi-step layout strategy optimised for
    networks with multiple disconnected components of varying sizes:

    \paragraph{Step 1: Component Identification and Sorting}
    Connected components are identified using depth-first search and
    sorted by size (descending). This ensures that larger, more important
    components are placed first in the grid layout.

    \paragraph{Step 2: Per-Component Spring Layout}
    Each component receives its own force-directed layout with adaptive
    parameters:
        - Spring constant $k$ scales with component size
        - Target bounding box scales with component size
        - Singleton components are placed at the origin

    \paragraph{Step 3: Centrality Computation}
    Betweenness centrality is computed within each component to identify
    structurally important nodes. The top-ranked nodes (by centrality)
    are flagged for visual highlighting:
        - Small components (n < 12): highlight top 1 node
        - Large components (n >= 12): highlight up to center_top_n nodes

    \paragraph{Step 4: Grid Packing}
    Components are arranged on a grid with configurable padding. The
    number of columns is set to $\lceil\sqrt{n_{\text{components}}}\rceil$
    to produce an approximately square arrangement.

    Parameters
    ----------
    G : networkx.Graph
        Input graph (may be disconnected).
    pad : float
        Distance between component centres in the grid layout.
    seed : int
        Random seed for reproducible spring layouts.
    k_min, k_max : float
        Range for adaptive spring constant.
    t_min, t_max : float
        Range for adaptive bounding box size.
    center_top_n : int
        Maximum number of central nodes to highlight per large component.

    Returns
    -------
    pos : dict
        Dictionary mapping node identifiers to (x, y) coordinates.
    components_list : list of set
        List of node sets, one per connected component.
    central_nodes : list
        List of node identifiers flagged as central (high betweenness).
    """
    components_list = sorted(
        nx.connected_components(G),
        key=lambda c: (-len(c), sorted(c)[0]),
        reverse=False
    )

    cols = max(1, math.ceil(math.sqrt(len(components_list))))

    pos = {}
    central_nodes = []

    for i, comp in enumerate(components_list):
        comp_nodes = sorted(comp)
        sub = G.subgraph(comp_nodes).copy()
        n = sub.number_of_nodes()

        if n == 1:
            local_pos = {comp_nodes[0]: np.array([0.0, 0.0])}
            central_nodes.append(comp_nodes[0])
        else:
            k = dynamic_spring_k(n, k_min=k_min, k_max=k_max)
            local_pos = nx.spring_layout(sub, seed=seed, k=k, iterations=220)

            target = dynamic_target_box(n, t_min=t_min, t_max=t_max)
            local_pos = scale_layout_to_min_box(local_pos, target=target)

            local_pos = repel_close_nodes(
                local_pos, sub,
                min_dist=0.55 * target,
                iters=140,
                step=0.10,
                seed=seed
            )

            bc = nx.betweenness_centrality(sub, normalized=True)
            ranked = sorted(bc.items(), key=lambda x: (-x[1], x[0]))

            if n >= 12:
                take = min(center_top_n, max(1, n // 12))
            else:
                take = 1
            central_nodes.extend([node for node, _ in ranked[:take]])

        row = i // cols
        col = i % cols
        pad_i = pad * (1.0 + 0.15 * math.log(max(n, 2)))
        offset = np.array([col * pad_i, -row * pad_i], dtype=float)

        for node, p in local_pos.items():
            pos[node] = np.array(p, dtype=float) + offset

    seen = set()
    central_nodes_unique = []
    for n in central_nodes:
        if n not in seen:
            central_nodes_unique.append(n)
            seen.add(n)

    return pos, components_list, central_nodes_unique


# =============================================================================
# Main Execution
# =============================================================================

def main():
    """
    Main entry point for the visualisation script.

    This function orchestrates the complete visualisation pipeline:

        1. Parse command-line arguments
        2. Load and validate input data
        3. Construct the similarity graph
        4. Compute packed component layout
        5. Generate visualisation
        6. Export supplementary data files

    The script enforces strict input rules to ensure reproducibility:
        - Only cell lines listed in the cell_line column become nodes
        - Only relationships in final_neighbours become edges
        - Invalid references (missing cell lines, self-loops) are logged
          and excluded
    """

    # -------------------------------------------------------------------------
    # Argument Parsing
    # -------------------------------------------------------------------------
    ap = argparse.ArgumentParser(
        description="Plot resolved DSMZ graph from resolved_dsmz_neighbors.tsv "
                    "(strict + packed components + centrality highlights)."
    )

    ap.add_argument("resolved_tsv",
                    help="Path to resolved_dsmz_neighbors.tsv")
    ap.add_argument("output_prefix",
                    help="Output prefix (no extension)")
    ap.add_argument("label",
                    help="Plot label (e.g., BRCA, NBL, RBL)")

    ap.add_argument("--pad", type=float, default=6.5,
                    help="Distance between packed components [default: 6.5]")

    ap.add_argument("--k-min", type=float, default=1.2,
                    help="Min spring k for small components [default: 1.2]")
    ap.add_argument("--k-max", type=float, default=3.0,
                    help="Max spring k for large components [default: 3.0]")
    ap.add_argument("--t-min", type=float, default=2.6,
                    help="Min normalised half-box size [default: 2.6]")
    ap.add_argument("--t-max", type=float, default=6.0,
                    help="Max normalised half-box size [default: 6.0]")

    ap.add_argument("--center-top-n", type=int, default=3,
                    help="Max number of central nodes to highlight per large "
                         "component [default: 3]")

    ap.add_argument("--fig-w", type=float, default=24.0,
                    help="Figure width (inches) [default: 24]")
    ap.add_argument("--fig-h", type=float, default=16.0,
                    help="Figure height (inches) [default: 16]")
    ap.add_argument("--dpi", type=int, default=600,
                    help="Raster output DPI [default: 600]")
    ap.add_argument("--seed", type=int, default=42,
                    help="Layout seed [default: 42]")

    # Optional supplementary inputs
    ap.add_argument("--node-stats", type=str, default=None,
                    help="Path to node stats TSV with columns cell_line, component, "
                         "is_isolate, degree, betweenness, community_louv, community_leid")
    ap.add_argument("--edge-stats", type=str, default=None,
                    help="Path to edge stats TSV with support_directions column")
    ap.add_argument("--component-annotations", type=str, default=None,
                    help="Path to component annotations TSV with columns "
                         "resolved_component_id, anchor_cell_line")
    ap.add_argument("--out-layout-coords", type=str, default=None,
                    help="Path to write layout coordinates TSV")
    ap.add_argument("--out-provenance", type=str, default=None,
                    help="Path to write figure provenance TSV")

    args = ap.parse_args()

    # -------------------------------------------------------------------------
    # Data Loading and Validation
    # -------------------------------------------------------------------------
    try:
        df = pd.read_csv(args.resolved_tsv, sep="\t", dtype=str)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to read TSV: {args.resolved_tsv}\n{e}\n")
        sys.exit(1)

    if "cell_line" not in df.columns:
        sys.stderr.write("[ERROR] TSV must contain column: cell_line\n")
        sys.stderr.write(f"Found columns: {', '.join(df.columns)}\n")
        sys.exit(1)

    if "final_neighbors" in df.columns:
        nb_col = "final_neighbors"
    elif "final_neighbours" in df.columns:
        nb_col = "final_neighbours"
    else:
        sys.stderr.write("[ERROR] TSV must contain column: final_neighbors or final_neighbours\n")
        sys.stderr.write(f"Found columns: {', '.join(df.columns)}\n")
        sys.exit(1)

    df["cell_line"] = df["cell_line"].astype(str).str.strip()
    df[nb_col] = df[nb_col].fillna("").astype(str)
    df = df[df["cell_line"] != ""].copy()

    # -------------------------------------------------------------------------
    # Optional node stats — override component and is_isolate if provided
    # -------------------------------------------------------------------------
    node_stats_ext = None
    if args.node_stats:
        try:
            node_stats_ext = pd.read_csv(args.node_stats, sep="\t", dtype=str)
            node_stats_ext["cell_line"] = node_stats_ext["cell_line"].astype(str).str.strip()
            print(f"[INFO] Loaded external node stats: {args.node_stats} "
                  f"({len(node_stats_ext)} rows)")
        except Exception as e:
            sys.stderr.write(f"[WARN] Could not load --node-stats: {e}\n")
            node_stats_ext = None

    # -------------------------------------------------------------------------
    # Optional edge stats — for support-weighted drawing
    # -------------------------------------------------------------------------
    edge_stats_ext = None
    edge_support_map = {}
    if args.edge_stats:
        try:
            edge_stats_ext = pd.read_csv(args.edge_stats, sep="\t")
            if "support_directions" in edge_stats_ext.columns:
                for _, er in edge_stats_ext.iterrows():
                    n1 = str(er.get("node1", "")).strip()
                    n2 = str(er.get("node2", "")).strip()
                    if n1 and n2:
                        key = tuple(sorted((n1, n2)))
                        edge_support_map[key] = int(er["support_directions"])
            print(f"[INFO] Loaded edge stats: {args.edge_stats} "
                  f"({len(edge_stats_ext)} rows, "
                  f"{len(edge_support_map)} edge support entries)")
        except Exception as e:
            sys.stderr.write(f"[WARN] Could not load --edge-stats: {e}\n")

    # -------------------------------------------------------------------------
    # Optional component annotations — anchor nodes
    # -------------------------------------------------------------------------
    anchor_nodes = set()
    if args.component_annotations:
        try:
            ann_df = pd.read_csv(args.component_annotations, sep="\t", dtype=str)
            if "anchor_cell_line" in ann_df.columns:
                anchor_nodes = set(
                    ann_df["anchor_cell_line"].dropna().astype(str).str.strip().tolist()
                )
                print(f"[INFO] Loaded {len(anchor_nodes)} anchor nodes from "
                      f"{args.component_annotations}")
        except Exception as e:
            sys.stderr.write(f"[WARN] Could not load --component-annotations: {e}\n")

    # -------------------------------------------------------------------------
    # Graph Construction
    # -------------------------------------------------------------------------
    node_set = set(df["cell_line"].unique())

    short_map = build_shortname_map(node_set)

    if len(set(short_map.values())) != len(short_map.values()):
        raise ValueError("[ERROR] Non-unique short names detected. Check build_shortname_map() rules.")

    mapping_out = (f"{args.output_prefix.rsplit('/', 1)[0]}/dsmz_cellline_shortnames.tsv"
                   if '/' in args.output_prefix
                   else "dsmz_cellline_shortnames.tsv")
    pd.DataFrame(
        [{"long_id": k, "short_id": v} for k, v in sorted(short_map.items())]
    ).to_csv(mapping_out, sep="\t", index=False)
    print(f"[OK] Saved: {mapping_out}")

    G = nx.Graph()
    G.add_nodes_from(sorted(node_set))

    ignored_neighbours = set()
    self_loops = 0

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

    # -------------------------------------------------------------------------
    # Determine component and isolate info
    # -------------------------------------------------------------------------
    isolates_raw = [n for n in G.nodes() if G.degree(n) == 0]
    comps = list(nx.connected_components(G))

    # Build per-node component index from graph structure (integer label)
    node_to_comp = {}
    for ci, comp in enumerate(comps):
        for nd in comp:
            node_to_comp[nd] = ci

    # Override with external node stats if available
    node_to_comp_ext = {}
    node_is_isolate_ext = {}
    node_degree_ext = {}
    node_betweenness_ext = {}
    node_community_louv_ext = {}
    node_community_leid_ext = {}
    node_display_ext = {}

    if node_stats_ext is not None:
        for _, nr in node_stats_ext.iterrows():
            cl = str(nr["cell_line"]).strip()
            if "component" in node_stats_ext.columns:
                try:
                    node_to_comp_ext[cl] = int(nr["component"])
                except (ValueError, TypeError):
                    pass
            if "is_isolate" in node_stats_ext.columns:
                val = str(nr["is_isolate"]).strip().upper()
                node_is_isolate_ext[cl] = val in {"TRUE", "1", "YES"}
            if "degree" in node_stats_ext.columns:
                try:
                    node_degree_ext[cl] = int(nr["degree"])
                except (ValueError, TypeError):
                    pass
            if "betweenness" in node_stats_ext.columns:
                try:
                    node_betweenness_ext[cl] = float(nr["betweenness"])
                except (ValueError, TypeError):
                    pass
            if "community_louv" in node_stats_ext.columns:
                node_community_louv_ext[cl] = str(nr.get("community_louv", ""))
            if "community_leid" in node_stats_ext.columns:
                node_community_leid_ext[cl] = str(nr.get("community_leid", ""))
            if "cell_line_display" in node_stats_ext.columns:
                node_display_ext[cl] = str(nr.get("cell_line_display", cl))

    def get_component(node):
        if node_to_comp_ext:
            return node_to_comp_ext.get(node, node_to_comp.get(node, -1))
        return node_to_comp.get(node, -1)

    def is_isolate_node(node):
        if node_is_isolate_ext:
            return node_is_isolate_ext.get(node, node in isolates_raw)
        return node in isolates_raw

    def get_degree(node):
        if node_degree_ext:
            return node_degree_ext.get(node, G.degree(node))
        return G.degree(node)

    isolates = [n for n in G.nodes() if is_isolate_node(n)]

    # -------------------------------------------------------------------------
    # Graph Statistics Summary
    # -------------------------------------------------------------------------
    print("=" * 60)
    print("GRAPH STATISTICS (STRICT: cell_line + final_neighbors only)")
    print("=" * 60)
    print(f"Nodes: {G.number_of_nodes()}")
    print(f"Edges: {G.number_of_edges()}")
    print(f"Connected components: {len(comps)}")
    print(f"Isolated nodes: {len(isolates)}")
    if self_loops:
        print(f"Self-loop entries ignored: {self_loops}")
    if ignored_neighbours:
        print(f"Neighbours ignored (not in cell_line): {len(ignored_neighbours)}")
    print()

    # -------------------------------------------------------------------------
    # Layout Computation
    # -------------------------------------------------------------------------
    pos, components_list, central_nodes = pack_components(
        G,
        pad=args.pad,
        seed=args.seed,
        k_min=args.k_min,
        k_max=args.k_max,
        t_min=args.t_min,
        t_max=args.t_max,
        center_top_n=args.center_top_n
    )

    pos = repel_close_nodes(
        pos, G,
        min_dist=0.18 * args.t_min,
        iters=60,
        step=0.04,
        seed=args.seed
    )

    # -------------------------------------------------------------------------
    # Colour assignment per component
    # -------------------------------------------------------------------------
    # Sort unique component IDs so colour assignment is deterministic
    unique_comp_ids = sorted(set(get_component(n) for n in G.nodes()))
    comp_colour_map = {}
    palette_idx = 0
    for cid in unique_comp_ids:
        # Isolate component: check whether all members of this component are isolates
        members = [n for n in G.nodes() if get_component(n) == cid]
        all_isolates = all(is_isolate_node(m) for m in members)
        if all_isolates:
            comp_colour_map[cid] = ISOLATE_COLOUR
        else:
            comp_colour_map[cid] = CB_PALETTE[palette_idx % len(CB_PALETTE)]
            palette_idx += 1

    def node_colour(node):
        if is_isolate_node(node):
            return ISOLATE_COLOUR
        return comp_colour_map.get(get_component(node), CB_PALETTE[0])

    def node_size(node):
        deg = get_degree(node)
        return max(300, 400 + 150 * deg)

    # -------------------------------------------------------------------------
    # Edge classification
    # -------------------------------------------------------------------------
    all_edges = list(G.edges())

    multi_supported_edges = []
    single_supported_edges = []
    no_stats_edges = []

    for u, v in all_edges:
        key = tuple(sorted((u, v)))
        sd = edge_support_map.get(key, None)
        if sd is None:
            no_stats_edges.append((u, v))
        elif sd >= 2:
            multi_supported_edges.append((u, v, sd))
        else:
            single_supported_edges.append((u, v, sd))

    # -------------------------------------------------------------------------
    # Visualisation
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(args.fig_w, args.fig_h), dpi=args.dpi)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # --- Edges: multi supported (solid) ---
    if multi_supported_edges:
        for u, v, sd in multi_supported_edges:
            width = min(6.0, 1.2 + 0.6 * (sd - 1))
            alpha = min(0.9, 0.3 + 0.1 * sd)
            nx.draw_networkx_edges(
                G, pos, ax=ax,
                edgelist=[(u, v)],
                width=width, alpha=alpha,
                edge_color="black", style="solid"
            )

    # --- Edges: single supported (dashed) ---
    if single_supported_edges:
        for u, v, sd in single_supported_edges:
            nx.draw_networkx_edges(
                G, pos, ax=ax,
                edgelist=[(u, v)],
                width=1.2, alpha=0.5,
                edge_color="#444444", style=(0, (5, 4))
            )

    # --- Edges: no edge stats (fallback solid) ---
    if no_stats_edges:
        nx.draw_networkx_edges(
            G, pos, ax=ax,
            edgelist=no_stats_edges,
            width=2.0, alpha=0.7,
            edge_color="black", style="solid"
        )

    # --- Nodes: regular non-isolates (circle) ---
    regular_nodes = [n for n in G.nodes()
                     if not is_isolate_node(n) and n not in anchor_nodes]
    if regular_nodes:
        nx.draw_networkx_nodes(
            G, pos,
            nodelist=regular_nodes,
            node_size=[node_size(n) for n in regular_nodes],
            node_color=[node_colour(n) for n in regular_nodes],
            edgecolors="black",
            linewidths=1.4,
            alpha=0.95,
            ax=ax
        )

    # --- Nodes: anchor nodes (circle, thick border) ---
    anchor_list = [n for n in anchor_nodes if n in G.nodes() and not is_isolate_node(n)]
    if anchor_list:
        nx.draw_networkx_nodes(
            G, pos,
            nodelist=anchor_list,
            node_size=[node_size(n) for n in anchor_list],
            node_color=[node_colour(n) for n in anchor_list],
            edgecolors="black",
            linewidths=4.5,
            alpha=1.0,
            ax=ax
        )

    # --- Nodes: isolates (diamond, separate scatter call) ---
    if isolates:
        ix = [pos[n][0] for n in isolates]
        iy = [pos[n][1] for n in isolates]
        isizes = [node_size(n) for n in isolates]
        # networkx uses pt², scatter uses pt for s; use consistent scale
        ax.scatter(ix, iy,
                   s=isizes,
                   c=ISOLATE_COLOUR,
                   marker="D",
                   edgecolors="black",
                   linewidths=1.4,
                   alpha=0.92,
                   zorder=4)

    # --- Labels ---
    display_labels = {}
    for n in G.nodes():
        if n in node_display_ext:
            display_labels[n] = node_display_ext[n].replace("_", r"\_")
        else:
            display_labels[n] = short_map.get(n, n).replace("_", r"\_")

    nx.draw_networkx_labels(
        G, pos,
        labels=display_labels,
        font_size=8,
        font_weight="bold",
        bbox=dict(facecolor="white", edgecolor="none",
                  boxstyle="round,pad=0.18", alpha=0.80),
        ax=ax
    )

    # --- Title ---
    ax.set_title(
        f"Resolved cell line similarity network — {args.label}",
        fontsize=14, pad=18
    )
    ax.set_aspect("equal")

    pts = np.array(list(pos.values()), dtype=float)
    minxy = pts.min(axis=0)
    maxxy = pts.max(axis=0)
    mx = (maxxy[0] - minxy[0]) * 0.08
    my = (maxxy[1] - minxy[1]) * 0.08
    ax.set_xlim(minxy[0] - mx, maxxy[0] + mx)
    ax.set_ylim(minxy[1] - my, maxxy[1] + my)
    ax.axis("off")

    # --- Legend ---
    legend_handles = []

    # One coloured patch per non-isolate component
    for cid in unique_comp_ids:
        members = [n for n in G.nodes() if get_component(n) == cid]
        all_iso = all(is_isolate_node(m) for m in members)
        if all_iso:
            continue
        colour = comp_colour_map[cid]
        n_members = len(members)
        legend_handles.append(
            mpatches.Patch(facecolor=colour, edgecolor="black", linewidth=1.0,
                           label=f"Component {cid} (n={n_members})")
        )

    # Isolates patch
    if isolates:
        legend_handles.append(
            mpatches.Patch(facecolor=ISOLATE_COLOUR, edgecolor="black", linewidth=1.0,
                           label=f"Isolate (n={len(isolates)})")
        )
        legend_handles.append(
            mlines.Line2D([], [], color=ISOLATE_COLOUR, marker="D",
                          markersize=8, markeredgecolor="black",
                          linewidth=0, label="Isolate (diamond shape)")
        )

    # Anchor nodes
    if anchor_list:
        legend_handles.append(
            mpatches.Patch(facecolor="white", edgecolor="black", linewidth=4.5,
                           label="Anchor cell line (thick border)")
        )

    # Edge styles
    has_multi = bool(multi_supported_edges)
    has_single = bool(single_supported_edges)
    has_fallback = bool(no_stats_edges)

    if has_multi:
        legend_handles.append(
            mlines.Line2D([], [], color="black", linewidth=2.0, linestyle="solid",
                          label="Multi supported edge (≥2 directions)")
        )
    if has_single:
        legend_handles.append(
            mlines.Line2D([], [], color="#444444", linewidth=1.5, linestyle=(0, (5, 4)),
                          label="Single supported edge (1 direction)")
        )
    if has_fallback and not (has_multi or has_single):
        legend_handles.append(
            mlines.Line2D([], [], color="black", linewidth=2.0, linestyle="solid",
                          label="Edge (no support data)")
        )

    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc="lower left",
            fontsize=7,
            framealpha=0.85,
            frameon=True,
            title="Visual encoding",
            title_fontsize=8,
        )

    # -------------------------------------------------------------------------
    # Save figure
    # -------------------------------------------------------------------------
    fig.savefig(f"{args.output_prefix}.png", dpi=args.dpi, bbox_inches="tight", pad_inches=0.25)
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", pad_inches=0.25)
    fig.savefig(f"{args.output_prefix}.svg", bbox_inches="tight", pad_inches=0.25)
    print(f"[OK] Saved: {args.output_prefix}.png, "
          f"{args.output_prefix}.pdf, {args.output_prefix}.svg")

    plt.close(fig)

    # -------------------------------------------------------------------------
    # Edge Statistics Export
    # -------------------------------------------------------------------------
    import os
    from pathlib import Path

    edges_output = (f"{args.output_prefix.rsplit('/', 1)[0]}/"
                    "dsmz_cellline_graph_edges.tsv"
                    if '/' in args.output_prefix
                    else "dsmz_cellline_graph_edges.tsv")

    output_dir = Path(edges_output).parent
    tumour_nh_root = output_dir.parent

    edge_records = []

    if "best_overall_dir" in df.columns and "winner_dir" in df.columns:
        directions_used = set()

        for _, row in df.iterrows():
            if pd.notna(row.get("best_overall_dir")):
                directions_used.add(str(row["best_overall_dir"]))
            if pd.notna(row.get("winner_dir")):
                directions_used.add(str(row["winner_dir"]))

        if tumour_nh_root.exists():
            for direction_dir in tumour_nh_root.iterdir():
                if direction_dir.is_dir() and direction_dir.name != "final_consensus_all":
                    direction = direction_dir.name
                    edge_file = (direction_dir / "final_consensus" /
                                 f"cell_line_similarity_graph_edges_{direction}.tsv")
                    if edge_file.exists():
                        directions_used.add(direction)

        for direction in sorted(directions_used):
            direction_dir = tumour_nh_root / direction / "final_consensus"
            edge_file = direction_dir / f"cell_line_similarity_graph_edges_{direction}.tsv"

            if edge_file.exists():
                try:
                    dir_edges = pd.read_csv(edge_file, sep="\t")

                    if ("cell_line1" in dir_edges.columns and
                            "cell_line2" in dir_edges.columns):
                        for _, edge_row in dir_edges.iterrows():
                            node1 = str(edge_row["cell_line1"]).strip()
                            node2 = str(edge_row["cell_line2"]).strip()

                            similarity = (float(edge_row.get("similarity", 1.0))
                                          if "similarity" in dir_edges.columns
                                          else 1.0)

                            if (node1 in node_set and node2 in node_set and
                                    G.has_edge(node1, node2)):
                                if node1 > node2:
                                    node1, node2 = node2, node1
                                edge_records.append({
                                    "node1": node1,
                                    "node2": node2,
                                    "direction": direction,
                                    "similarity": similarity
                                })
                except Exception as e:
                    print(f"[WARNING] Failed to read {edge_file}: {e}",
                          file=sys.stderr)

    edges_source = "direction_files"
    n_edges_written = 0
    if edge_records:
        edges_df = pd.DataFrame(edge_records)

        methods_map = {}
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
        for (node1, node2), group in edges_df.groupby(["node1", "node2"]):
            directions_list = sorted(group["direction"].unique())
            similarities = group["similarity"].tolist()

            methods_union = methods_map.get(node1, methods_map.get(node2, ""))

            aggregated.append({
                "node1": node1,
                "node2": node2,
                "support_directions": len(directions_list),
                "support_weight_mean": (np.mean(similarities)
                                        if similarities else 1.0),
                "support_weight_sum": sum(similarities),
                "support_weight_max": max(similarities) if similarities else 1.0,
                "methods_union": methods_union
            })

        rich_edges_df = pd.DataFrame(aggregated)
        rich_edges_df["node1_short"] = rich_edges_df["node1"].map(lambda x: short_map.get(x, x))
        rich_edges_df["node2_short"] = rich_edges_df["node2"].map(lambda x: short_map.get(x, x))
        rich_edges_df = rich_edges_df.sort_values(["node1_short", "node2_short", "node1", "node2"])
        rich_edges_df.to_csv(edges_output, sep="\t", index=False)
        n_edges_written = len(rich_edges_df)
        print(f"[OK] Saved: {edges_output} ({len(rich_edges_df)} edges)")
    else:
        print("[WARN] No direction-specific edge files found; "
              "falling back to resolved neighbors.")
        edges_source = "fallback_resolved_neighbors"
        fallback_edges = edges_from_resolved_neighbors(df, nb_col=nb_col, node_set=node_set)
        if fallback_edges:
            fb = pd.DataFrame(fallback_edges)
            fb["node1_short"] = fb["node1"].map(lambda x: short_map.get(x, x))
            fb["node2_short"] = fb["node2"].map(lambda x: short_map.get(x, x))
            fb.to_csv(edges_output, sep="\t", index=False)
            n_edges_written = len(fallback_edges)
            print(f"[OK] Saved: {edges_output} ({len(fallback_edges)} edges)")
        else:
            pd.DataFrame(columns=["node1", "node2", "support_directions",
                                  "support_weight_mean", "support_weight_sum",
                                  "support_weight_max", "methods_union"]).to_csv(
                edges_output, sep="\t", index=False)
            print(f"[OK] Saved: {edges_output} (empty - no edges)")

    print("[QC] Edge export summary:",
          f"nodes={G.number_of_nodes()}",
          f"edges_written={n_edges_written}",
          f"source={edges_source}")

    # -------------------------------------------------------------------------
    # Node Statistics Export
    # -------------------------------------------------------------------------
    node_stats_output = (f"{args.output_prefix.rsplit('/', 1)[0]}/"
                         "dsmz_cellline_graph_node_stats.tsv"
                         if '/' in args.output_prefix
                         else "dsmz_cellline_graph_node_stats.tsv")

    bc = nx.betweenness_centrality(G, normalized=True)

    node_stats = []
    for node in sorted(G.nodes()):
        comp_id = get_component(node)
        is_iso = is_isolate_node(node)
        component_label = "Isolate" if is_iso else f"Component {comp_id}"
        anchor_status = node in anchor_nodes

        # Use external betweenness if available, else computed
        bet_val = node_betweenness_ext.get(node, bc.get(node, 0.0))

        ns_row = {
            "cell_line": node,
            "sample_id": node,
            "cell_line_display": node_display_ext.get(node, short_map.get(node, node)),
            "component": comp_id,
            "is_isolate": str(is_iso).upper(),
            "degree": get_degree(node),
            "betweenness": bet_val,
            "community_louv": node_community_louv_ext.get(node, ""),
            "community_leid": node_community_leid_ext.get(node, ""),
            # New alias and annotation columns
            "resolved_component_id": comp_id,
            "component_label": component_label,
            "anchor_status": str(anchor_status).upper(),
        }
        node_stats.append(ns_row)

    node_stats_df = pd.DataFrame(node_stats)
    node_stats_df.to_csv(node_stats_output, sep="\t", index=False)
    print(f"[OK] Saved: {node_stats_output}")

    # -------------------------------------------------------------------------
    # Optional: layout coordinates TSV
    # -------------------------------------------------------------------------
    if args.out_layout_coords:
        Path(args.out_layout_coords).parent.mkdir(parents=True, exist_ok=True)
        coord_rows = []
        for node in sorted(G.nodes()):
            x, y = float(pos[node][0]), float(pos[node][1])
            coord_rows.append({
                "cohort": args.label,
                "graph_type": "resolved",
                "node_id": node,
                "display_label": node_display_ext.get(node, short_map.get(node, node)),
                "layout_method": "spring_weighted_fixed_seed",
                "layout_seed": args.seed,
                "x": x,
                "y": y,
            })
        pd.DataFrame(coord_rows).to_csv(args.out_layout_coords, sep="\t", index=False)
        print(f"[OK] Wrote layout coords: {args.out_layout_coords}")

    # -------------------------------------------------------------------------
    # Optional: figure provenance TSV
    # -------------------------------------------------------------------------
    if args.out_provenance:
        Path(args.out_provenance).parent.mkdir(parents=True, exist_ok=True)
        enc = (
            "node_colour=component_id; node_shape=isolate_status(circle/diamond); "
            "node_size=degree; border_width=anchor_status; "
            "edge_style=support_type(solid/dashed); edge_width=support_directions"
        )
        prov_row = {
            "figure_file": f"{args.output_prefix}.png",
            "script": __file__,
            "input_files": args.resolved_tsv,
            "date_time": datetime.datetime.now().isoformat(),
            "layout_seed": args.seed,
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "visual_encodings": enc,
            "output_width": args.fig_w,
            "output_height": args.fig_h,
        }
        pd.DataFrame([prov_row]).to_csv(args.out_provenance, sep="\t", index=False)
        print(f"[OK] Wrote provenance: {args.out_provenance}")

    plt.close()


# =============================================================================
# Script Entry Point
# =============================================================================

if __name__ == "__main__":
    main()
