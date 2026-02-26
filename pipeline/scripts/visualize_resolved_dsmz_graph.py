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
import sys
import math
import re
from collections import defaultdict
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# Ensure underscores render literally (avoid mathtext subscripts)
mpl.rcParams["text.usetex"] = False
mpl.rcParams["mathtext.default"] = "regular"


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
    # Handle missing or NaN values
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return []

    s = str(s).strip()
    if not s:
        return []

    # Split on semicolons and clean each identifier
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

    # 1) compute base for each long id
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
        # fallback (should be rare in RBL)
        base_for[s] = s

    # 2) group by base
    by_base = defaultdict(list)
    for long_id, base in base_for.items():
        by_base[base].append(long_id)

    short_map = {}

    for base, longs in by_base.items():
        longs = sorted(longs)

        # singleton base -> no suffix
        if len(longs) == 1:
            short_map[longs[0]] = base
            continue

        # duplicates -> MUST use trailing token
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

    # Logarithmic scaling with base 50 as reference point
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

    # Compute bounding box
    minxy = pts.min(axis=0)
    maxxy = pts.max(axis=0)

    # Avoid division by zero for single-node or collinear layouts
    span = np.maximum(maxxy - minxy, 1e-9)

    out = {}
    for node, p in local_pos.items():
        # Rescale to [0, 1]
        q = (np.array(p) - minxy) / span
        # Shift and scale to [-target, target]
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

    # Half-span in the larger dimension
    half_span = 0.5 * max(span[0], span[1])

    # Scale up only (never shrink)
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

    # Work component-wise so we don't destroy component packing
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
                        # Identical coordinates; break symmetry deterministically
                        jitter = rng.normal(0, 1, size=2)
                        jitter = jitter / (np.linalg.norm(jitter) + 1e-9)
                        dvec = jitter
                        dist = 1e-9

                    if dist < min_dist:
                        # Push apart equally
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
    # Identify and sort components by size (largest first)
    # Tie-break by smallest node name for deterministic ordering across runs
    components_list = sorted(
        nx.connected_components(G),
        key=lambda c: (-len(c), sorted(c)[0]),
        reverse=False
    )

    # Determine grid dimensions for packing
    cols = max(1, math.ceil(math.sqrt(len(components_list))))

    pos = {}
    central_nodes = []

    for i, comp in enumerate(components_list):
        # Sort nodes deterministically for reproducible spring_layout input
        comp_nodes = sorted(comp)
        sub = G.subgraph(comp_nodes).copy()
        n = sub.number_of_nodes()

        # --- Per-component layout ---
        if n == 1:
            # Singleton: place at origin (will be offset by grid packing)
            local_pos = {comp_nodes[0]: np.array([0.0, 0.0])}
            central_nodes.append(comp_nodes[0])
        else:
            # Multi-node component: spring layout with adaptive parameters
            k = dynamic_spring_k(n, k_min=k_min, k_max=k_max)
            local_pos = nx.spring_layout(sub, seed=seed, k=k, iterations=220)

            # Do NOT squash into a box; only scale up if component is too tight
            target = dynamic_target_box(n, t_min=t_min, t_max=t_max)
            local_pos = scale_layout_to_min_box(local_pos, target=target)

            # Improve within-component separation BEFORE packing (component-local)
            local_pos = repel_close_nodes(
                local_pos, sub,
                min_dist=0.55 * target,   # 0.45–0.75 * target (bigger = more spread)
                iters=140,
                step=0.10,
                seed=seed
            )

            # --- Centrality analysis ---
            # Betweenness centrality identifies bridge nodes within the component
            bc = nx.betweenness_centrality(sub, normalized=True)
            ranked = sorted(bc.items(), key=lambda x: x[1], reverse=True)

            # Select top central nodes for highlighting
            # Larger components may have multiple important bridge nodes
            if n >= 12:
                take = min(center_top_n, max(1, n // 12))
            else:
                take = 1
            central_nodes.extend([node for node, _ in ranked[:take]])

        # --- Grid packing offset ---
        row = i // cols
        col = i % cols
        # Scale pad by component size so big components don't collide visually
        pad_i = pad * (1.0 + 0.15 * math.log(max(n, 2)))
        offset = np.array([col * pad_i, -row * pad_i], dtype=float)

        for node, p in local_pos.items():
            pos[node] = np.array(p, dtype=float) + offset

    # Remove duplicate central nodes while preserving order
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

    # Required positional arguments
    ap.add_argument("resolved_tsv",
                    help="Path to resolved_dsmz_neighbors.tsv")
    ap.add_argument("output_prefix",
                    help="Output prefix (no extension)")
    ap.add_argument("label",
                    help="Plot label (e.g., BRCA, NBL, RBL)")

    # Component packing controls
    ap.add_argument("--pad", type=float, default=6.5,
                    help="Distance between packed components [default: 6.5]")

    # Adaptive spring layout parameters
    ap.add_argument("--k-min", type=float, default=1.2,
                    help="Min spring k for small components [default: 1.2]")
    ap.add_argument("--k-max", type=float, default=3.0,
                    help="Max spring k for large components [default: 3.0]")
    ap.add_argument("--t-min", type=float, default=2.6,
                    help="Min normalised half-box size [default: 2.6]")
    ap.add_argument("--t-max", type=float, default=6.0,
                    help="Max normalised half-box size [default: 6.0]")

    # Centrality highlighting
    ap.add_argument("--center-top-n", type=int, default=3,
                    help="Max number of central nodes to highlight per large "
                         "component [default: 3]")

    # Figure aesthetics
    ap.add_argument("--fig-w", type=float, default=24.0,
                    help="Figure width (inches) [default: 24]")
    ap.add_argument("--fig-h", type=float, default=16.0,
                    help="Figure height (inches) [default: 16]")
    ap.add_argument("--dpi", type=int, default=600,
                    help="Raster output DPI [default: 600]")
    ap.add_argument("--seed", type=int, default=42,
                    help="Layout seed [default: 42]")

    args = ap.parse_args()

    # -------------------------------------------------------------------------
    # Data Loading and Validation
    # -------------------------------------------------------------------------
    # The script performs strict validation to ensure data integrity.
    # Malformed input is caught early with informative error messages.

    try:
        df = pd.read_csv(args.resolved_tsv, sep="\t", dtype=str)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to read TSV: {args.resolved_tsv}\n{e}\n")
        sys.exit(1)

    # Validate required columns
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

    # Clean identifiers
    df["cell_line"] = df["cell_line"].astype(str).str.strip()
    df[nb_col] = df[nb_col].fillna("").astype(str)
    df = df[df["cell_line"] != ""].copy()

    # -------------------------------------------------------------------------
    # Graph Construction
    # -------------------------------------------------------------------------
    # The node universe is defined strictly by the cell_line column.
    # This ensures that the graph accurately represents the resolved
    # relationships without introducing spurious nodes from incomplete data.

    node_set = set(df["cell_line"].unique())

    # Build short-name mapping for cleaner labels and exports
    short_map = build_shortname_map(node_set)

    # Sanity check: ensure all short names are unique
    if len(set(short_map.values())) != len(short_map.values()):
        raise ValueError("[ERROR] Non-unique short names detected. Check build_shortname_map() rules.")

    # Optional: write mapping table for reference (highly recommended)
    mapping_out = (f"{args.output_prefix.rsplit('/', 1)[0]}/dsmz_cellline_shortnames.tsv"
                   if '/' in args.output_prefix
                   else "dsmz_cellline_shortnames.tsv")
    pd.DataFrame(
        [{"long_id": k, "short_id": v} for k, v in sorted(short_map.items())]
    ).to_csv(mapping_out, sep="\t", index=False)
    print(f"[OK] Saved: {mapping_out}")

    # Initialise graph with all nodes (including potential isolates)
    G = nx.Graph()
    G.add_nodes_from(sorted(node_set))

    # Track data quality issues for logging
    ignored_neighbours = set()
    self_loops = 0

    # Add edges from the final_neighbors column
    for _, row in df.iterrows():
        src = row["cell_line"]
        for tgt in parse_neighbours(row[nb_col]):
            # Exclude self-loops (cell line listed as its own neighbour)
            if tgt == src:
                self_loops += 1
                continue
            # Exclude references to cell lines not in the node universe
            if tgt not in node_set:
                ignored_neighbours.add(tgt)
                continue
            G.add_edge(src, tgt)

    # -------------------------------------------------------------------------
    # Graph Statistics Summary
    # -------------------------------------------------------------------------
    # Summary statistics provide insight into the network structure and
    # help identify potential data quality issues.

    isolates = [n for n in G.nodes() if G.degree(n) == 0]
    comps = list(nx.connected_components(G))

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
    # The packed component layout ensures that each connected component
    # is visualised optimally, with adaptive spacing based on size and
    # central nodes highlighted for structural interpretation.

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

    # Keep packing coordinates (don't squeeze to 0-1, which removes visual separation)
    # pos = normalize_pos_to_square(pos)  # DISABLED: preserves component separation

    # Mild global repulsion (optional - most separation handled within components)
    # Only needed for edge cases where components are very close after packing
    pos = repel_close_nodes(
        pos, G,
        min_dist=0.18 * args.t_min,  # Very mild - just prevents extreme edge cases
        iters=60,                     # Fewer iterations since main work done per-component
        step=0.04,                    # Smaller step to avoid inflating packed layout
        seed=args.seed
    )

    # -------------------------------------------------------------------------
    # Visualisation
    # -------------------------------------------------------------------------
    # The visualisation employs a layered approach:
    #   1. Edges (bottom layer)
    #   2. Base nodes (middle layer)
    #   3. Central nodes with enhanced styling (top layer)
    #   4. Labels (topmost layer)

    fig, ax = plt.subplots(figsize=(args.fig_w, args.fig_h), dpi=args.dpi)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Layer 1: Edges
    # Thicker, darker edges for better visibility in slides/PDFs
    nx.draw_networkx_edges(G, pos, ax=ax, width=2.0, alpha=0.70, edge_color="black")

    # Layer 2: Base nodes (normal nodes)
    # Light blue fill with thin black outline
    base_node_size = 1400  # good for slides; use 1000 for thesis pdf if too big
    normal_nodes = [n for n in G.nodes() if n not in central_nodes]
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=normal_nodes,
        node_size=base_node_size,
        node_color="lightblue",
        edgecolors="black",
        linewidths=1.4,
        alpha=0.95,
        ax=ax
    )

    # Layer 3: Central nodes (highlighted)
    # Thicker outline and larger size indicate high betweenness
    if central_nodes:
        nx.draw_networkx_nodes(
            G, pos,
            nodelist=central_nodes,
            node_size=int(base_node_size * 1.35),
            node_color="lightblue",
            edgecolors="black",
            linewidths=4.0,
            alpha=1.0,
            ax=ax
        )

    # Layer 4: Labels
    # Larger font with background box for readability
    # Use short names and escape underscores to prevent mathtext interpretation
    safe_labels = {n: short_map.get(n, n).replace("_", r"\_") for n in G.nodes()}
    
    nx.draw_networkx_labels(
        G, pos,
        labels=safe_labels,
        font_size=11,  # slides: 11–14; thesis: 9–11
        font_weight="bold",
        bbox=dict(facecolor="lightblue", edgecolor="none",
                  boxstyle="round,pad=0.18"),
        ax=ax
    )

    # Title with explanatory subtitle
    ax.set_title(
        f"Resolved DSMZ Neighbour Graph - {args.label}\n"
        "(Packed components; adaptive spacing; "
        "thick outline = high betweenness centrality)",
        fontsize=14, pad=18
    )
    ax.set_aspect("equal")
    
    # Autoscale to fit content with margin (preserves component separation)
    pts = np.array(list(pos.values()), dtype=float)
    minxy = pts.min(axis=0)
    maxxy = pts.max(axis=0)
    
    # margin in data-coordinates (8% of span)
    mx = (maxxy[0] - minxy[0]) * 0.08
    my = (maxxy[1] - minxy[1]) * 0.08
    
    ax.set_xlim(minxy[0] - mx, maxxy[0] + mx)
    ax.set_ylim(minxy[1] - my, maxxy[1] + my)
    ax.axis("off")

    # Save in both raster (PNG) and vector (PDF) formats with tight bounding
    # pad_inches prevents edge labels from being clipped
    fig.savefig(f"{args.output_prefix}.png", dpi=args.dpi, bbox_inches="tight", pad_inches=0.25)
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", pad_inches=0.25)
    fig.savefig(f"{args.output_prefix}.svg", bbox_inches="tight", pad_inches=0.25)
    print(f"[OK] Saved: {args.output_prefix}.png, "
          f"{args.output_prefix}.pdf, {args.output_prefix}.svg")

    plt.close(fig)

    # -------------------------------------------------------------------------
    # Edge Statistics Export
    # -------------------------------------------------------------------------
    # The rich edges file aggregates support information across different
    # feature-distance directions, providing a comprehensive view of edge
    # reliability.

    import os
    from pathlib import Path

    edges_output = (f"{args.output_prefix.rsplit('/', 1)[0]}/"
                    "dsmz_cellline_graph_edges.tsv"
                    if '/' in args.output_prefix
                    else "dsmz_cellline_graph_edges.tsv")

    output_dir = Path(edges_output).parent
    tumour_nh_root = output_dir.parent

    # Collect edges from direction-specific files
    edge_records = []

    if "best_overall_dir" in df.columns and "winner_dir" in df.columns:
        directions_used = set()

        # Extract directions from resolved neighbours metadata
        for _, row in df.iterrows():
            if pd.notna(row.get("best_overall_dir")):
                directions_used.add(str(row["best_overall_dir"]))
            if pd.notna(row.get("winner_dir")):
                directions_used.add(str(row["winner_dir"]))

        # Discover additional directions from directory structure
        if tumour_nh_root.exists():
            for direction_dir in tumour_nh_root.iterdir():
                if direction_dir.is_dir() and direction_dir.name != "final_consensus_all":
                    direction = direction_dir.name
                    edge_file = (direction_dir / "final_consensus" /
                                 f"DSMZ_DSMZ_graph_edges_{direction}.tsv")
                    if edge_file.exists():
                        directions_used.add(direction)

        # Read and aggregate edge files
        for direction in sorted(directions_used):
            direction_dir = tumour_nh_root / direction / "final_consensus"
            edge_file = direction_dir / f"DSMZ_DSMZ_graph_edges_{direction}.tsv"

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

                            # Only include edges present in the resolved graph
                            if (node1 in node_set and node2 in node_set and
                                    G.has_edge(node1, node2)):
                                # Normalise edge order (alphabetically)
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

    # Aggregate edges across directions
    edges_source = "direction_files"
    n_edges_written = 0
    if edge_records:
        edges_df = pd.DataFrame(edge_records)

        # Build methods map from resolved neighbours
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

        # Aggregate by edge pair
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
        # Add short name columns
        rich_edges_df["node1_short"] = rich_edges_df["node1"].map(lambda x: short_map.get(x, x))
        rich_edges_df["node2_short"] = rich_edges_df["node2"].map(lambda x: short_map.get(x, x))
        # Sort by short names for nicer viewing (keep long names for reference)
        rich_edges_df = rich_edges_df.sort_values(["node1_short", "node2_short", "node1", "node2"])
        rich_edges_df.to_csv(edges_output, sep="\t", index=False)
        n_edges_written = len(rich_edges_df)
        print(f"[OK] Saved: {edges_output} ({len(rich_edges_df)} edges)")
    else:
        # Fallback: derive edges directly from resolved neighbors
        print("[WARN] No direction-specific edge files found; "
              "falling back to resolved neighbors.")
        edges_source = "fallback_resolved_neighbors"
        fallback_edges = edges_from_resolved_neighbors(df, nb_col=nb_col, node_set=node_set)
        if fallback_edges:
            fb = pd.DataFrame(fallback_edges)
            # Add short name columns
            fb["node1_short"] = fb["node1"].map(lambda x: short_map.get(x, x))
            fb["node2_short"] = fb["node2"].map(lambda x: short_map.get(x, x))
            fb.to_csv(edges_output, sep="\t", index=False)
            n_edges_written = len(fallback_edges)
            print(f"[OK] Saved: {edges_output} ({len(fallback_edges)} edges)")
        else:
            # Create empty file with headers for consistency
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
    # Per-node statistics facilitate downstream analysis and provide
    # quantitative support for visual observations.

    node_stats_output = (f"{args.output_prefix.rsplit('/', 1)[0]}/"
                         "dsmz_cellline_graph_node_stats.tsv"
                         if '/' in args.output_prefix
                         else "dsmz_cellline_graph_node_stats.tsv")

    # Compute betweenness centrality for all nodes
    bc = nx.betweenness_centrality(G, normalized=True)

    node_stats = []
    for node in sorted(G.nodes()):
        node_stats.append({
            "cell_line": node,
            "cell_line_short": short_map.get(node, node),
            "degree": G.degree(node),
            "betweenness": bc.get(node, 0.0),
            "component": next((i for i, comp in enumerate(comps)
                               if node in comp), -1),
            "is_isolate": node in isolates,
            "is_central": node in central_nodes
        })

    node_stats_df = pd.DataFrame(node_stats)
    node_stats_df.to_csv(node_stats_output, sep="\t", index=False)
    print(f"[OK] Saved: {node_stats_output}")

    plt.close()


# =============================================================================
# Script Entry Point
# =============================================================================

if __name__ == "__main__":
    main()


