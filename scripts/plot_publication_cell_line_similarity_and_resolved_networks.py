#!/usr/bin/env python3
"""
plot_publication_cell_line_similarity_and_resolved_networks.py

Produces the two main thesis figures for DSMZ cell line networks across the
BRCA, NBL, and RBL cohorts.

The separate final_consensus_all/community_stability/ outputs are retained as
supplementary graph stability validation only. They are not inputs to these
main figures and must not be interpreted as final consensus Leiden/Louvain
communities.

Main Figure 1 — clinical patient sample referenced cell line similarity
    Fig_cell_line_similarity_networks_clinical_patient_referenced_combined.{pdf,svg,png}
    Panels:
        A. BRCA clinical patient sample referenced cell line similarity network
        B. NBL  clinical patient sample referenced cell line similarity network
        C. RBL  clinical patient sample referenced cell line similarity network
        D. Similarity network construction summary
    Encodings:
        nodes  = DSMZ cell lines
        edges  = Pearson similarity between clinical patient sample
                 neighbourhood consensus profiles (winning direction per cohort)
        node colour = Leiden community (per-cohort palette; Louvain shown only
                       in the supplementary node TSV for comparison)
        edge width  = scaled by similarity value
        diamond markers = isolates (degree 0)

Main Figure 2 — graph based consensus resolved DSMZ neighbourhoods
    Fig_consensus_resolved_DSMZ_neighbourhoods_combined.{pdf,svg,png}
    Panels:
        A. BRCA graph based consensus resolved DSMZ cell line neighbourhoods
        B. NBL  graph based consensus resolved DSMZ cell line neighbourhoods
        C. RBL  graph based consensus resolved DSMZ cell line neighbourhoods
        D. Resolved component summary
    Encodings:
        nodes = DSMZ cell lines
        edges = final stable neighbour relationships
        node colour = final connected component
        component central node = highest unnormalised betweenness centrality
            node within each non-isolate component (highlighted with a thick
            red border). Normalised betweenness is also reported for audit.

Supporting TSVs:
    cell_line_similarity_network_edges.tsv
    cell_line_similarity_network_nodes.tsv
    similarity_network_construction_summary.tsv
    resolved_DSMZ_neighbourhood_node_summary.tsv
    resolved_DSMZ_component_annotations.tsv
    supp_resolved_DSMZ_component_annotations.tsv
    resolved_component_central_nodes.tsv
    anchor_centrality_audit.tsv
    cell_line_network_layout_coordinates.tsv
    network_results_summary.txt
    network_input_validation.tsv
    figure_provenance.tsv
    plot_publication_cell_line_similarity_and_resolved_networks.log

Strict mode (default) requires all three cohorts. --allow-missing-cohorts is
for diagnostic runs only.
"""

from __future__ import annotations

import argparse
import datetime
import logging
import math
import re
import subprocess
import sys
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

mpl.rcParams["text.usetex"] = False
mpl.rcParams["mathtext.default"] = "regular"


# =============================================================================
# Constants
# =============================================================================

CB_PALETTE = [
    "#4477AA", "#EE6677", "#228833", "#CCBB44",
    "#66CCEE", "#AA3377", "#BB5566", "#DDAA33",
]

# Extended categorical palette for the similarity figure: enough distinct
# colours that no two Leiden communities share a colour by coincidence even
# when the upstream Leiden detector returns singleton-dominated communities
# (BRCA = 29 communities at PCA_corr resolution).
LEIDEN_PALETTE = (
    CB_PALETTE
    + ["#882255", "#117733", "#999933", "#332288", "#44AA99", "#88CCEE",
       "#AA4499", "#DDDDDD"]
    + ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
       "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
       "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94",
       "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5"]
)

ISOLATE_COLOUR = "#AAAAAA"
CENTRAL_BORDER_COLOUR = "#B22222"

COHORTS = ["BRCA", "NBL", "RBL"]
COHORT_PANEL_COLOURS = {
    "BRCA": "#4477AA",
    "NBL":  "#EE6677",
    "RBL":  "#228833",
}

# RBL has no winning_direction.txt; this cohort default is used unless
# the user supplies --rbl-direction.
DEFAULT_RBL_DIRECTION = "hvg_euc"

SIM_EDGE_FILE_TEMPLATES = [
    "cell_line_similarity_graph_edges_{dir}.tsv",
    "DSMZ_DSMZ_graph_edges_{dir}.tsv",
]
SIM_NODE_ANN_TEMPLATES = [
    "cell_line_similarity_graph_node_annotations_{dir}.tsv",
    "DSMZ_DSMZ_graph_node_annotations_{dir}.tsv",
]
SIM_ISOLATES_TEMPLATES = [
    "cell_line_similarity_graph_isolates_{dir}.tsv",
    "DSMZ_DSMZ_graph_isolates_{dir}.tsv",
]

CITATION_FIG1 = (
    "Source TSVs: cell_line_similarity_network_edges.tsv, "
    "cell_line_similarity_network_nodes.tsv, "
    "similarity_network_construction_summary.tsv."
)
CITATION_FIG2 = (
    "Source TSVs: resolved_DSMZ_neighbourhood_node_summary.tsv, "
    "resolved_DSMZ_component_annotations.tsv, "
    "resolved_component_central_nodes.tsv."
)

log = logging.getLogger(__name__)


# =============================================================================
# Layout helpers (verbatim from plot_publication_cell_line_component_networks)
# =============================================================================

def _parse_neighbours(s):
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return []
    s = str(s).strip()
    if not s or s.lower() in ("nan", "none"):
        return []
    return [x.strip() for x in re.split(r"[;,]", s) if x.strip()]


def _dynamic_k(n, k_min=1.2, k_max=3.0):
    if n <= 1:
        return k_min
    return max(k_min, min(k_max,
                          k_min + (math.log(n) / math.log(50)) * (k_max - k_min)))


def _dynamic_box(n, t_min=2.6, t_max=6.0):
    if n <= 1:
        return t_min
    return max(t_min, min(t_max,
                          t_min + (math.log(n) / math.log(50)) * (t_max - t_min)))


def _scale_to_box(local_pos, target):
    pts = np.array(list(local_pos.values()), dtype=float)
    span = np.maximum(pts.max(axis=0) - pts.min(axis=0), 1e-9)
    half = 0.5 * max(span)
    s = (target / half) if half < target else 1.0
    centre = (pts.min(axis=0) + pts.max(axis=0)) / 2.0
    return {n: (np.array(p, float) - centre) * s for n, p in local_pos.items()}


def _repel(pos, G, min_d=0.55, iters=80, step=0.06, seed=42):
    rng = np.random.default_rng(seed)
    for comp in nx.connected_components(G):
        nodes = sorted(comp)
        if len(nodes) <= 1:
            continue
        for _ in range(iters):
            moved = False
            for i in range(len(nodes)):
                u = nodes[i]
                pu = np.array(pos[u], float)
                for j in range(i + 1, len(nodes)):
                    v = nodes[j]
                    pv = np.array(pos[v], float)
                    dvec = pu - pv
                    dist = float(np.linalg.norm(dvec))
                    if dist < 1e-9:
                        dvec = rng.normal(0, 1, 2)
                        dist = 1e-9
                    if dist < min_d:
                        push = (min_d - dist) * step * dvec / dist
                        pos[u] = pu + push
                        pos[v] = pv - push
                        pu = np.array(pos[u], float)
                        moved = True
            if not moved:
                break
    return pos


def pack_components(G, pad=8.0, seed=42, k_min=1.2, k_max=3.0,
                    t_min=2.6, t_max=6.0):
    comps = sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)[0]))
    cols = max(1, math.ceil(math.sqrt(len(comps))))
    pos = {}
    for i, comp in enumerate(comps):
        nodes = sorted(comp)
        sub = G.subgraph(nodes).copy()
        n = sub.number_of_nodes()
        if n == 1:
            lp = {nodes[0]: np.zeros(2)}
        else:
            k = _dynamic_k(n, k_min, k_max)
            lp = nx.spring_layout(sub, seed=seed, k=k, iterations=220)
            target = _dynamic_box(n, t_min, t_max)
            lp = _scale_to_box(lp, target)
            lp = _repel(lp, sub, min_d=0.55 * target,
                        iters=140, step=0.10, seed=seed)
        row_, col_ = i // cols, i % cols
        pad_i = pad * (1.0 + 0.15 * math.log(max(n, 2)))
        offset = np.array([col_ * pad_i, -row_ * pad_i], float)
        for nd, p in lp.items():
            pos[nd] = np.array(p, float) + offset
    return pos, comps


def make_display_label(cell_line_id):
    m = re.match(r"^NG-\d+_(.+)_lib\d+_\d+_(\d+)$", str(cell_line_id))
    if m:
        return f"{m.group(1)}_r{m.group(2)}"
    return cell_line_id


# =============================================================================
# Generic utilities
# =============================================================================

def get_git_commit():
    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
            cwd=Path(__file__).parent,
        )
        if r.returncode == 0:
            return r.stdout.strip()
    except Exception:
        pass
    return "unavailable"


def load_tsv(path, dtype=str):
    if path is None:
        return None, "no path provided"
    p = Path(path)
    if not p.exists():
        return None, f"file not found: {p}"
    try:
        return pd.read_csv(p, sep="\t", dtype=dtype), None
    except Exception as exc:
        return None, f"read error: {exc}"


def first_existing(paths):
    for p in paths:
        if p is not None and Path(p).exists():
            return Path(p)
    return None


def assign_palette_colours(group_ids, isolate_ids=None, palette=None):
    """Map a sequence of group IDs to palette colours; isolates → grey."""
    isolate_ids = isolate_ids or set()
    palette = palette or CB_PALETTE
    seen, colours = {}, {}
    palette_idx = 0
    for gid in group_ids:
        if gid in isolate_ids:
            colours[gid] = ISOLATE_COLOUR
            continue
        if gid not in seen:
            seen[gid] = palette[palette_idx % len(palette)]
            palette_idx += 1
        colours[gid] = seen[gid]
    return colours


# =============================================================================
# Direction discovery and similarity-graph loading
# =============================================================================

def detect_winning_direction(tn_dir, override=None, fallback=None,
                              required_cells=None):
    """Pick the winning direction for a cohort.

    Order of resolution:
      1. CLI override (--<cohort>-direction).
      2. Contents of `final_consensus_all/winning_direction.txt`.
      3. `fallback` parameter (used for RBL where the file is absent).
      4. Top-ranked direction in `p_consensus_direction_summary.tsv` that
         contains all cells in `required_cells`, if supplied.
    """
    if override:
        return override.strip(), "CLI --<cohort>-direction override"

    base = Path(tn_dir) / "final_consensus_all"
    wf = base / "winning_direction.txt"
    if wf.exists():
        try:
            text = wf.read_text().strip().splitlines()
            if text:
                d = text[0].strip()
                if d:
                    return d, f"winning_direction.txt = {d}"
        except Exception:
            pass

    summary = base / "p_consensus_direction_summary.tsv"
    if summary.exists():
        try:
            df = pd.read_csv(summary, sep="\t")
            df["frac_ge_thr"] = pd.to_numeric(
                df["frac_ge_thr"], errors="coerce")
            df["n_cell_lines"] = pd.to_numeric(
                df["n_cell_lines"], errors="coerce")
            # Filter to directions that cover all required resolved cells, then
            # rank by consensus quality (frac_ge_thr). This avoids picking a
            # direction with many extra non-DSMZ cells but poor consensus.
            df = df.sort_values("frac_ge_thr", ascending=False)
            for _, row in df.iterrows():
                d = str(row["direction"]).strip()
                if not d:
                    continue
                if (required_cells is not None and
                        not _direction_covers(tn_dir, d, required_cells)):
                    continue
                return d, ("p_consensus_direction_summary.tsv top by "
                          f"frac_ge_thr (covering all resolved cells) = {d}")
        except Exception:
            pass

    if fallback:
        return fallback, f"hard-coded fallback = {fallback}"
    return None, "no winning direction could be determined"


def _direction_covers(tn_dir, direction, required_cells):
    base = Path(tn_dir) / direction / "final_consensus"
    nodes_fp = first_existing(
        base / t.format(dir=direction) for t in SIM_NODE_ANN_TEMPLATES)
    if nodes_fp is None:
        return False
    try:
        df = pd.read_csv(nodes_fp, sep="\t", dtype=str)
    except Exception:
        return False
    if "cell_line" not in df.columns:
        return False
    have = set(df["cell_line"].dropna().astype(str).str.strip())
    return required_cells.issubset(have)


def discover_similarity_files(tn_dir, direction):
    base = Path(tn_dir) / direction / "final_consensus"
    edges = first_existing(
        base / t.format(dir=direction) for t in SIM_EDGE_FILE_TEMPLATES)
    nodes = first_existing(
        base / t.format(dir=direction) for t in SIM_NODE_ANN_TEMPLATES)
    iso   = first_existing(
        base / t.format(dir=direction) for t in SIM_ISOLATES_TEMPLATES)
    return edges, nodes, iso


def discover_resolved_files(tn_dir):
    base = Path(tn_dir) / "final_consensus_all"
    resolved = first_existing(
        base / fn for fn in
        ("resolved_dsmz_neighbours.tsv", "resolved_dsmz_neighbors.tsv"))
    node_stats = first_existing([
        base / "patient_referenced_aggregated_cell_line_similarity_graph_node_stats.tsv",
        base / "dsmz_cellline_graph_node_stats.tsv",  # legacy fallback (pre-rename)
    ])
    edges = first_existing([
        base / "plots" / "patient_referenced_resolved_cell_line_neighbourhood_graph_edges.tsv",
        base / "plots" / "resolved_cell_line_neighbourhood_graph_edges.tsv",  # legacy fallback (pre-prefix)
        base / "plots" / "dsmz_cellline_graph_edges.tsv",                     # legacy fallback (pre-rename)
        base / "dsmz_cellline_graph_edges.tsv",
    ])
    anchors = first_existing([base / "anchor_components.tsv"])
    return resolved, node_stats, edges, anchors


# =============================================================================
# Similarity network construction (Main Figure 1)
# =============================================================================

def build_similarity_graph(edges_df, nodes_df, isolates_df=None):
    """Return G with similarity-weighted edges; isolates included as
    nodes with degree 0; node attributes carry community_louv/leid + degree.
    """
    G = nx.Graph()
    if (edges_df is not None
            and {"cell_line1", "cell_line2"}.issubset(edges_df.columns)):
        for _, row in edges_df.iterrows():
            u = str(row["cell_line1"]).strip()
            v = str(row["cell_line2"]).strip()
            if not u or not v or u == "nan" or v == "nan":
                continue
            try:
                w = float(row.get("similarity", "nan"))
            except (TypeError, ValueError):
                w = float("nan")
            G.add_edge(u, v, similarity=w)

    if nodes_df is not None and "cell_line" in nodes_df.columns:
        for _, row in nodes_df.iterrows():
            cl = str(row["cell_line"]).strip()
            if not cl or cl == "nan":
                continue
            G.add_node(cl)
            for col in ("community_louv", "community_leid",
                        "degree", "betweenness", "component",
                        "mean_edge_sim", "max_edge_sim",
                        "cell_line_display", "is_outlier"):
                if col in nodes_df.columns:
                    val = row.get(col, "")
                    G.nodes[cl][col] = val

    if isolates_df is not None and "cell_line" in isolates_df.columns:
        for _, row in isolates_df.iterrows():
            cl = str(row["cell_line"]).strip()
            if not cl or cl == "nan":
                continue
            G.add_node(cl)
            for col in ("cell_line_display", "is_outlier",
                        "mean_edge_sim", "max_edge_sim"):
                if col in isolates_df.columns:
                    G.nodes[cl].setdefault(col, row.get(col, ""))

    return G


def _leiden_id(node_attrs, default="NA"):
    v = node_attrs.get("community_leid", default)
    s = "" if v is None else str(v).strip()
    if not s or s.lower() == "nan":
        return default
    return s


def _louvain_id(node_attrs, default="NA"):
    v = node_attrs.get("community_louv", default)
    s = "" if v is None else str(v).strip()
    if not s or s.lower() == "nan":
        return default
    return s


def draw_similarity_panel(ax, cohort, direction, G, pos, node_colours,
                          isolates, panel_letter, font_size=8):
    isolate_set = set(isolates)

    sims = [d.get("similarity", float("nan")) for _, _, d in G.edges(data=True)]
    sims_clean = [s for s in sims if not (s is None or
                                          (isinstance(s, float) and math.isnan(s)))]
    if sims_clean:
        smin = min(sims_clean)
        smax = max(sims_clean)
        srange = smax - smin if smax > smin else 1e-9
    else:
        smin = smax = 0.0
        srange = 1.0

    for u, v, data in G.edges(data=True):
        s = data.get("similarity", float("nan"))
        if isinstance(s, float) and math.isnan(s):
            w = 1.5
            alpha = 0.5
        else:
            w = 0.8 + 4.5 * ((float(s) - smin) / srange)
            alpha = 0.35 + 0.55 * ((float(s) - smin) / srange)
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)],
            width=w, alpha=alpha,
            edge_color="#555555", ax=ax,
        )

    def _size(n):
        try:
            d = float(G.nodes[n].get("degree", G.degree(n)))
        except (TypeError, ValueError):
            d = G.degree(n)
        return max(300, 400 + 120 * d)

    regular = [n for n in G.nodes() if n not in isolate_set]
    if regular:
        nx.draw_networkx_nodes(
            G, pos, nodelist=regular,
            node_color=[node_colours.get(n, ISOLATE_COLOUR) for n in regular],
            node_size=[_size(n) for n in regular],
            linewidths=1.0, edgecolors="#333333", ax=ax,
        )

    iso_in_pos = [n for n in isolate_set if n in pos]
    if iso_in_pos:
        ax.scatter(
            [pos[n][0] for n in iso_in_pos],
            [pos[n][1] for n in iso_in_pos],
            s=[_size(n) * 0.55 for n in iso_in_pos],
            c=[node_colours.get(n, ISOLATE_COLOUR) for n in iso_in_pos],
            marker="D", edgecolors="#333333", linewidths=1.0, zorder=4,
        )

    labels = {n: (str(G.nodes[n].get("cell_line_display", "")).strip()
                  or make_display_label(n))
              for n in G.nodes() if n in pos}
    nx.draw_networkx_labels(
        G, pos, labels=labels,
        font_size=max(5, font_size - 2), font_color="#1a1a1a", ax=ax,
    )

    ax.set_title(
        f"Panel {panel_letter} — {cohort} similarity network ({direction})",
        fontsize=font_size + 2, pad=8,
    )
    ax.axis("off")

    handles = [
        mpatches.Patch(facecolor=ISOLATE_COLOUR, edgecolor="#333333",
                       linewidth=1.0, label="Isolate (degree 0)"),
        mlines.Line2D([], [], marker="o", color="none",
                      markerfacecolor=CB_PALETTE[0],
                      markeredgecolor="#333333", markersize=8,
                      label="Leiden community (per-cohort palette)"),
        mlines.Line2D([], [], color="#555555", linewidth=4.5,
                      label="High Pearson similarity"),
        mlines.Line2D([], [], color="#555555", linewidth=1.0,
                      label="Low Pearson similarity"),
    ]
    ax.legend(handles=handles, loc="lower left",
              fontsize=max(5, font_size - 2),
              framealpha=0.85, ncol=1, handlelength=1.8)


# =============================================================================
# Resolved network construction (Main Figure 2)
# =============================================================================

def _edge_support_map(edges_df):
    m = {}
    if (edges_df is not None and
            {"node1", "node2", "support_directions"}.issubset(edges_df.columns)):
        for _, er in edges_df.iterrows():
            key = tuple(sorted((str(er["node1"]).strip(),
                                str(er["node2"]).strip())))
            try:
                m[key] = float(er["support_directions"])
            except (ValueError, TypeError):
                pass
    return m


def build_resolved_graph(resolved_df, node_stats_df=None):
    nb_col = next((c for c in ("final_neighbours", "final_neighbors")
                   if c in resolved_df.columns), None)
    if nb_col is None or "cell_line" not in resolved_df.columns:
        return None, {}, [], None

    G = nx.Graph()
    for _, row in resolved_df.iterrows():
        cl = str(row["cell_line"]).strip()
        if not cl or cl == "nan":
            continue
        G.add_node(cl)
        for nb in _parse_neighbours(row.get(nb_col, "")):
            if nb and nb != cl:
                G.add_edge(cl, nb)

    if node_stats_df is not None and "cell_line" in node_stats_df.columns:
        for cl in node_stats_df["cell_line"].dropna().astype(str).str.strip():
            if cl and cl != "nan":
                G.add_node(cl)

    node_to_comp = {}
    for cid, comp in enumerate(
        sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)[0]))
    ):
        for node in comp:
            node_to_comp[node] = cid

    isolates = [n for n in G.nodes() if G.degree(n) == 0]
    return G, node_to_comp, isolates, nb_col


def assign_component_colours(G, node_to_comp, isolates):
    isolate_set = set(isolates)
    comp_members = defaultdict(list)
    for n in G.nodes():
        comp_members[node_to_comp[n]].append(n)

    colour_map = {}
    palette_idx = 0
    for cid, members in sorted(comp_members.items()):
        if all(m in isolate_set for m in members) or len(members) == 1:
            colour_map[cid] = ISOLATE_COLOUR
        else:
            colour_map[cid] = CB_PALETTE[palette_idx % len(CB_PALETTE)]
            palette_idx += 1
    return colour_map


def central_node_per_component(G, node_to_comp, isolate_set):
    """Return component descriptors for degree and betweenness centrality.

    For each non-isolate connected component, the highest-degree node is
    reported as the most connected cell line and the highest unnormalised
    betweenness node is reported as the bridge/anchor-like cell line.
    Normalised betweenness is computed in parallel for audit/provenance only.
    """
    comp_members = defaultdict(list)
    for n, c in node_to_comp.items():
        comp_members[c].append(n)
    central = {}
    btw_unnorm = {n: 0.0 for n in G.nodes()}
    btw_norm = {n: 0.0 for n in G.nodes()}
    selected_by_unnorm = {}
    selected_by_norm = {}
    selected_by_degree = {}
    for cid, members in comp_members.items():
        if len(members) <= 1 or all(m in isolate_set for m in members):
            central[cid] = (None, None)
            selected_by_unnorm[cid] = None
            selected_by_norm[cid] = None
            selected_by_degree[cid] = None
            continue
        sub = G.subgraph(members)
        sub_unnorm = nx.betweenness_centrality(
            sub, normalized=False, weight=None)
        sub_norm = nx.betweenness_centrality(
            sub, normalized=True, weight=None)
        for node in members:
            btw_unnorm[node] = float(sub_unnorm.get(node, 0.0))
            btw_norm[node] = float(sub_norm.get(node, 0.0))
        best_unnorm = max(members, key=lambda m: (btw_unnorm.get(m, 0.0), m))
        best_norm = max(members, key=lambda m: (btw_norm.get(m, 0.0), m))
        best_degree = max(members, key=lambda m: (G.degree(m), m))
        selected_by_unnorm[cid] = best_unnorm
        selected_by_norm[cid] = best_norm
        selected_by_degree[cid] = best_degree
        central[cid] = (
            best_unnorm,
            round(float(btw_unnorm.get(best_unnorm, 0.0)), 6),
        )
    return (central, btw_unnorm, btw_norm, selected_by_unnorm,
            selected_by_norm, selected_by_degree)


def draw_resolved_panel(ax, cohort, G, pos, node_to_comp, isolates,
                        comp_colour_map, central_nodes, anchor_nodes,
                        edges_df, node_stats_df, panel_letter, font_size=8):
    isolate_set = set(isolates)
    anchor_set = set(anchor_nodes) if anchor_nodes else set()
    central_set = {v for (v, _) in central_nodes.values() if v}
    edge_support = _edge_support_map(edges_df)

    def _colour(n):
        return comp_colour_map.get(node_to_comp.get(n, -1), ISOLATE_COLOUR)

    def _size(n):
        return max(300, 400 + 120 * G.degree(n))

    multi_edges, single_edges, fallback_edges = [], [], []
    for u, v in G.edges():
        key = tuple(sorted((u, v)))
        sd = edge_support.get(key)
        if sd is not None:
            (multi_edges if sd >= 2 else single_edges).append((u, v, sd))
        else:
            fallback_edges.append((u, v, 1.0))

    for edge_list, style in [(multi_edges, "solid"),
                              (single_edges, "dashed"),
                              (fallback_edges, "solid")]:
        for u, v, sd in edge_list:
            w = min(6.0, 1.5 + 0.5 * (sd - 1))
            alpha = min(0.85, 0.3 + 0.1 * sd)
            nx.draw_networkx_edges(
                G, pos, edgelist=[(u, v)],
                width=w, alpha=alpha, style=style,
                edge_color="#555555", ax=ax,
            )

    regular = [n for n in G.nodes() if n not in isolate_set]
    central_reg = [n for n in regular if n in central_set]
    other_reg = [n for n in regular if n not in central_set]

    if other_reg:
        nx.draw_networkx_nodes(
            G, pos, nodelist=other_reg,
            node_color=[_colour(n) for n in other_reg],
            node_size=[_size(n) for n in other_reg],
            linewidths=1.0, edgecolors="#333333", ax=ax,
        )
    if central_reg:
        nx.draw_networkx_nodes(
            G, pos, nodelist=central_reg,
            node_color=[_colour(n) for n in central_reg],
            node_size=[_size(n) * 1.05 for n in central_reg],
            linewidths=4.5, edgecolors=CENTRAL_BORDER_COLOUR, ax=ax,
        )

    iso_in_pos = [n for n in isolates if n in pos]
    if iso_in_pos:
        ax.scatter(
            [pos[n][0] for n in iso_in_pos],
            [pos[n][1] for n in iso_in_pos],
            s=[_size(n) * 0.55 for n in iso_in_pos],
            c=[_colour(n) for n in iso_in_pos],
            marker="D", edgecolors="#333333", linewidths=1.0, zorder=4,
        )

    disp_map = {}
    if (node_stats_df is not None and "cell_line" in node_stats_df.columns
            and "cell_line_display" in node_stats_df.columns):
        for _, row in node_stats_df.iterrows():
            disp_map[str(row["cell_line"]).strip()] = str(row["cell_line_display"])
    labels = {n: (disp_map.get(n) or make_display_label(n))
              for n in G.nodes() if n in pos}
    nx.draw_networkx_labels(
        G, pos, labels=labels,
        font_size=max(5, font_size - 2), font_color="#1a1a1a", ax=ax,
    )

    ax.set_title(
        f"Panel {panel_letter} — {cohort} resolved DSMZ neighbourhoods",
        fontsize=font_size + 2, pad=8,
    )
    ax.axis("off")

    handles = [
        mpatches.Patch(facecolor="#AAAAAA", edgecolor="#333333",
                       linewidth=1.0, label="Cell line node"),
        mpatches.Patch(facecolor="#AAAAAA", edgecolor=CENTRAL_BORDER_COLOUR,
                       linewidth=4.5,
                       label="Component central node "
                             "(highest unnormalised betweenness)"),
        mlines.Line2D([], [], marker="D", color="none",
                      markerfacecolor=ISOLATE_COLOUR,
                      markeredgecolor="#333333", markersize=7,
                      label="Isolate (degree 0)"),
        mlines.Line2D([], [], color="#555555", linewidth=2.0,
                      linestyle="solid", label="Multi supported edge"),
        mlines.Line2D([], [], color="#555555", linewidth=1.5,
                      linestyle="dashed", label="Single supported edge"),
    ]
    if anchor_set:
        handles.insert(2, mlines.Line2D(
            [], [], marker="*", color="#333333",
            markerfacecolor="none", markersize=10,
            linestyle="none", label="Anchor cell line (when provided)"))
    ax.legend(handles=handles, loc="lower left",
              fontsize=max(5, font_size - 2),
              framealpha=0.85, ncol=1, handlelength=1.8)


# =============================================================================
# Panel D table renderers
# =============================================================================

def _render_table(ax, columns, rows, row_colours, title, citation):
    ax.axis("off")
    ax.set_title(title, fontsize=12, pad=6)

    table_ax = ax.inset_axes([0.04, 0.40, 0.92, 0.55])
    table_ax.axis("off")
    table = table_ax.table(
        cellText=rows, colLabels=columns,
        cellColours=row_colours,
        colColours=["#DDDDDD"] * len(columns),
        loc="center", cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.6)

    note_ax = ax.inset_axes([0.04, 0.05, 0.92, 0.32])
    note_ax.axis("off")
    note_ax.text(
        0.0, 0.95, citation,
        ha="left", va="top",
        fontsize=8, color="#222222", style="italic",
        transform=note_ax.transAxes, wrap=True,
    )


def draw_similarity_summary_panel(ax, summary_rows):
    columns_keys = [
        ("cohort", "Cohort"),
        ("direction", "Winning direction"),
        ("n_nodes", "Cell lines"),
        ("n_edges", "Edges"),
        ("n_isolates", "Isolates"),
        ("n_leiden", "Leiden\ncommunities"),
        ("n_louvain", "Louvain\ncommunities"),
        ("median_sim", "Median\nsimilarity"),
        ("min_sim", "Min\nsimilarity"),
    ]
    columns = [label for _, label in columns_keys]
    cell_text = []
    row_colours = []
    for r in summary_rows:
        cell_text.append([str(r.get(k, "n/a")) for k, _ in columns_keys])
        bg = COHORT_PANEL_COLOURS.get(r["cohort"], "#FFFFFF") + "33"
        row_colours.append([bg] * len(columns))
    _render_table(
        ax, columns, cell_text, row_colours,
        title="Panel D — Similarity network construction summary",
        citation=(
            "Per cohort: winning feature-method × distance-metric direction "
            "from the unsupervised clinical patient sample neighbourhood "
            "consensus. Edges are Pearson similarities between cell lines "
            "based on tumour-neighbourhood profiles in that direction. "
            "Node colour in panels A–C is the Leiden community of that "
            "direction; Louvain is shown only as a comparison count. "
            "When n_Leiden equals n_cell_lines (e.g. BRCA, NBL, RBL), "
            "Leiden returned singleton-dominated communities at the "
            "pipeline resolution, so every cell line is its own Leiden "
            "community. "
            + CITATION_FIG1
        ),
    )


def draw_resolved_summary_panel(ax, summary_rows):
    columns_keys = [
        ("cohort", "Cohort"),
        ("n_nodes", "Nodes"),
        ("n_edges", "Edges"),
        ("n_components", "Components"),
        ("n_multi", "Multi-node\ncomponents"),
        ("n_isolates", "Isolates"),
        ("largest_component", "Largest\ncomponent"),
        ("median_degree", "Median\ndegree"),
        ("median_edge_support", "Median\nedge support"),
        ("central_nodes", "Central nodes\n(comp:cell line)"),
    ]
    columns = [label for _, label in columns_keys]
    cell_text = []
    row_colours = []
    for r in summary_rows:
        cell_text.append([str(r.get(k, "n/a")) for k, _ in columns_keys])
        bg = COHORT_PANEL_COLOURS.get(r["cohort"], "#FFFFFF") + "33"
        row_colours.append([bg] * len(columns))
    _render_table(
        ax, columns, cell_text, row_colours,
        title="Panel D — Resolved component summary",
        citation=(
            "Per cohort: graph-based consensus resolved DSMZ cell line "
            "neighbourhoods. Central node = highest unnormalised "
            "betweenness centrality node within each non-isolate component; "
            "normalised betweenness is reported only for audit/provenance. "
            + CITATION_FIG2
        ),
    )


# =============================================================================
# Per-cohort summary metrics
# =============================================================================

def compute_similarity_summary(cohort_data):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("similarity_available"):
            rows.append({
                "cohort": cohort,
                "direction": "n/a",
                "n_nodes": "n/a", "n_edges": "n/a", "n_isolates": "n/a",
                "n_leiden": "n/a", "n_louvain": "n/a",
                "median_sim": "n/a", "min_sim": "n/a",
            })
            continue
        G = d["similarity_graph"]
        sims = [data.get("similarity") for _, _, data in G.edges(data=True)
                if data.get("similarity") is not None
                and not (isinstance(data.get("similarity"), float)
                         and math.isnan(data["similarity"]))]
        leid = {_leiden_id(G.nodes[n]) for n in G.nodes()
                if _leiden_id(G.nodes[n]) != "NA"}
        louv = {_louvain_id(G.nodes[n]) for n in G.nodes()
                if _louvain_id(G.nodes[n]) != "NA"}
        isolates = [n for n in G.nodes() if G.degree(n) == 0]
        rows.append({
            "cohort": cohort,
            "direction": d.get("direction", "n/a"),
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "n_isolates": len(isolates),
            "n_leiden": len(leid),
            "n_louvain": len(louv),
            "median_sim": (round(float(np.median(sims)), 3) if sims else "n/a"),
            "min_sim":    (round(float(np.min(sims)), 3) if sims else "n/a"),
        })
    return rows


def compute_resolved_summary(cohort_data):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("resolved_available"):
            rows.append({
                "cohort": cohort,
                "n_nodes": "n/a", "n_edges": "n/a",
                "n_components": "n/a", "n_multi": "n/a",
                "n_isolates": "n/a",
                "largest_component": "n/a",
                "median_degree": "n/a",
                "median_edge_support": "n/a",
                "central_nodes": "n/a",
            })
            continue
        G = d["resolved_graph"]
        node_to_comp = d["node_to_comp"]
        isolates = d["resolved_isolates"]
        ed = d.get("resolved_edges_df")
        comp_sizes = defaultdict(int)
        for n in G.nodes():
            comp_sizes[node_to_comp[n]] += 1
        n_multi = sum(1 for s in comp_sizes.values() if s > 1)

        sds = []
        if ed is not None and "support_directions" in ed.columns:
            sds = pd.to_numeric(
                ed["support_directions"], errors="coerce").dropna().tolist()

        degrees = [G.degree(n) for n in G.nodes()]

        central = d.get("central_nodes", {})
        central_strs = []
        for cid, (cn, _bv) in sorted(central.items()):
            if cn is None:
                continue
            label = make_display_label(cn)
            central_strs.append(f"{cid}:{label}")
        central_blob = ("\n".join(central_strs)
                        if central_strs else "no multi-node component")

        rows.append({
            "cohort": cohort,
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "n_components": len(comp_sizes),
            "n_multi": n_multi,
            "n_isolates": len(isolates),
            "largest_component": (max(comp_sizes.values())
                                   if comp_sizes else 0),
            "median_degree": (round(float(np.median(degrees)), 2)
                              if degrees else 0),
            "median_edge_support": (round(float(np.median(sds)), 2)
                                    if sds else "n/a"),
            "central_nodes": central_blob,
        })
    return rows


# =============================================================================
# TSV writers
# =============================================================================

def write_validation_report(val_rows, outdir):
    path = Path(outdir) / "network_input_validation.tsv"
    pd.DataFrame(val_rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_similarity_edges(cohort_data, outdir):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("similarity_available"):
            continue
        G = d["similarity_graph"]
        direction = d.get("direction", "")
        for u, v, data in G.edges(data=True):
            sim = data.get("similarity")
            rows.append({
                "cohort": cohort,
                "direction": direction,
                "cell_line_a": u,
                "cell_line_b": v,
                "cell_line_a_display": (str(G.nodes[u].get("cell_line_display",
                                                             "")).strip()
                                         or make_display_label(u)),
                "cell_line_b_display": (str(G.nodes[v].get("cell_line_display",
                                                             "")).strip()
                                         or make_display_label(v)),
                "pearson_similarity": (round(float(sim), 6)
                                        if (isinstance(sim, (int, float))
                                            and not math.isnan(float(sim)))
                                        else ""),
                "leiden_community_a": _leiden_id(G.nodes[u]),
                "leiden_community_b": _leiden_id(G.nodes[v]),
                "louvain_community_a": _louvain_id(G.nodes[u]),
                "louvain_community_b": _louvain_id(G.nodes[v]),
            })
    path = Path(outdir) / "cell_line_similarity_network_edges.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_similarity_nodes(cohort_data, outdir):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("similarity_available"):
            continue
        G = d["similarity_graph"]
        direction = d.get("direction", "")
        for n in sorted(G.nodes()):
            attrs = G.nodes[n]
            rows.append({
                "cohort": cohort,
                "direction": direction,
                "cell_line": n,
                "cell_line_display": (str(attrs.get("cell_line_display",
                                                     "")).strip()
                                       or make_display_label(n)),
                "leiden_community_id": _leiden_id(attrs),
                "louvain_community_id": _louvain_id(attrs),
                "degree": G.degree(n),
                "is_isolate": str(G.degree(n) == 0).upper(),
                "mean_edge_similarity": str(attrs.get("mean_edge_sim", "")).strip(),
                "max_edge_similarity": str(attrs.get("max_edge_sim", "")).strip(),
                "is_outlier": str(attrs.get("is_outlier", "")).strip(),
            })
    path = Path(outdir) / "cell_line_similarity_network_nodes.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_similarity_construction_summary(summary_rows, outdir):
    path = Path(outdir) / "similarity_network_construction_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_resolved_node_summary(cohort_data, outdir):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("resolved_available"):
            continue
        G = d["resolved_graph"]
        node_to_comp = d["node_to_comp"]
        isolates = set(d["resolved_isolates"])
        anchor_set = set(d.get("anchor_nodes", set()))
        anchor_provided = d.get("anchor_path_provided", False)
        ns_df = d.get("resolved_node_stats_df")
        es_map = _edge_support_map(d.get("resolved_edges_df"))
        central_pairs = d.get("central_nodes", {})
        central_set = {v for (v, _) in central_pairs.values() if v}
        btw_unnorm = d.get("resolved_betweenness_unnormalised", {})
        btw_norm = d.get("resolved_betweenness_normalised", {})
        selected_unnorm = d.get("selected_by_unnormalised", {})
        selected_norm = d.get("selected_by_normalised", {})
        selected_degree = d.get("selected_by_degree", {})

        comp_sizes = defaultdict(int)
        comp_members = defaultdict(list)
        for n in G.nodes():
            comp_sizes[node_to_comp[n]] += 1
            comp_members[node_to_comp[n]].append(n)

        def _comp_label(cid, members):
            n = len(members)
            all_iso = all(m in isolates for m in members)
            if all_iso or n == 1:
                return "Isolate"
            if n <= 3:
                return "Small component"
            return f"Component {cid}"

        comp_labels = {cid: _comp_label(cid, members)
                       for cid, members in comp_members.items()}

        disp_map = {}
        if ns_df is not None and "cell_line" in ns_df.columns:
            for _, row in ns_df.iterrows():
                cl = str(row["cell_line"]).strip()
                if "cell_line_display" in ns_df.columns:
                    disp_map[cl] = str(row.get("cell_line_display", "")).strip()

        for node in sorted(G.nodes()):
            cid = node_to_comp[node]
            inc_sd = [
                es_map[tuple(sorted((node, nb)))]
                for nb in G.neighbors(node)
                if tuple(sorted((node, nb))) in es_map
            ]
            display = disp_map.get(node) or make_display_label(node)
            is_anchor_val = ("NA" if not anchor_provided
                             else str(node in anchor_set).upper())
            rows.append({
                "cohort": cohort,
                "cell_line_raw_id": node,
                "cell_line_display": display,
                "connected_component_id": cid,
                "connected_component_size": comp_sizes[cid],
                "resolved_component_label": comp_labels[cid],
                "is_central_node": str(node in central_set).upper(),
                "betweenness_centrality":
                    round(float(btw_unnorm.get(node, 0.0)), 6),
                "betweenness_unnormalised":
                    round(float(btw_unnorm.get(node, 0.0)), 6),
                "betweenness_normalised":
                    round(float(btw_norm.get(node, 0.0)), 6),
                "selected_by_unnormalised":
                    str(selected_unnorm.get(cid) == node).upper(),
                "selected_by_normalised":
                    str(selected_norm.get(cid) == node).upper(),
                "selected_by_degree":
                    str(selected_degree.get(cid) == node).upper(),
                "canonical_selected": str(node in central_set).upper(),
                "canonical_bridge_selected": str(node in central_set).upper(),
                "most_connected_selected":
                    str(selected_degree.get(cid) == node).upper(),
                "centrality_metric": "degree_and_betweenness",
                "centrality_normalised": "FALSE",
                "centrality_weighted": "FALSE",
                "centrality_scope": "within_component",
                "centrality_tie_break": "highest_value_then_node_id",
                "degree": G.degree(node),
                "is_isolate": str(node in isolates).upper(),
                "is_anchor": is_anchor_val,
                "median_incident_edge_support":
                    (round(float(np.median(inc_sd)), 3) if inc_sd else ""),
            })

    path = Path(outdir) / "resolved_DSMZ_neighbourhood_node_summary.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_resolved_component_annotations(cohort_data, outdir):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("resolved_available"):
            continue
        G = d["resolved_graph"]
        node_to_comp = d["node_to_comp"]
        isolates = set(d["resolved_isolates"])
        anchor_set = set(d.get("anchor_nodes", set()))
        es_map = _edge_support_map(d.get("resolved_edges_df"))
        btw_unnorm = d.get("resolved_betweenness_unnormalised", {})
        btw_norm = d.get("resolved_betweenness_normalised", {})
        selected_unnorm = d.get("selected_by_unnormalised", {})
        selected_norm = d.get("selected_by_normalised", {})
        selected_degree = d.get("selected_by_degree", {})
        central_pairs = d.get("central_nodes", {})

        comp_members = defaultdict(list)
        for n in G.nodes():
            comp_members[node_to_comp[n]].append(n)

        for cid, members in sorted(comp_members.items()):
            n = len(members)
            all_iso = all(m in isolates for m in members)
            if all_iso or n == 1:
                label = "Isolate"
            elif n <= 3:
                label = "Small component"
            else:
                label = f"Component {cid}"

            sub = G.subgraph(members)
            inc_sds = [
                es_map[tuple(sorted((m, nb)))]
                for m in members
                for nb in G.neighbors(m)
                if tuple(sorted((m, nb))) in es_map
            ]
            degrees = [G.degree(m) for m in members]
            btw_unnorm_vals = [btw_unnorm.get(m, 0) for m in members]
            btw_norm_vals = [btw_norm.get(m, 0) for m in members]
            cn, cn_btw = central_pairs.get(cid, (None, None))
            cn_norm = btw_norm.get(cn, "") if cn else ""
            most_connected = selected_degree.get(cid)
            most_connected_degree = (
                G.degree(most_connected) if most_connected else ""
            )
            changed_metric = (
                selected_unnorm.get(cid) != selected_norm.get(cid)
                if selected_unnorm.get(cid) is not None else False
            )

            rows.append({
                "cohort": cohort,
                "resolved_component_id": cid,
                "component_label": label,
                "n_cell_lines": n,
                "member_cell_lines": ";".join(sorted(members)),
                "central_node": cn or "",
                "bridge_node": cn or "",
                "most_connected_node": most_connected or "",
                "most_connected_degree": most_connected_degree,
                "central_node_betweenness":
                    (cn_btw if cn_btw is not None else ""),
                "central_node_betweenness_unnormalised":
                    (cn_btw if cn_btw is not None else ""),
                "central_node_betweenness_normalised":
                    (round(float(cn_norm), 6) if cn_norm != "" else ""),
                "central_node_selection_metric":
                    ("highest_unnormalised_betweenness"
                     if cn else "single_member_or_isolate"),
                "most_connected_selection_metric":
                    ("highest_degree" if most_connected
                     else "single_member_or_isolate"),
                "centrality_metric": "degree_and_betweenness",
                "centrality_normalised": "FALSE",
                "centrality_weighted": "FALSE",
                "centrality_scope": "within_component",
                "centrality_tie_break": "highest_value_then_node_id",
                "selected_by_normalised": selected_norm.get(cid) or "",
                "selected_by_unnormalised": selected_unnorm.get(cid) or "",
                "selected_by_degree": selected_degree.get(cid) or "",
                "degree_vs_betweenness_node_differs":
                    str(most_connected != cn).upper(),
                "normalised_vs_unnormalised_anchor_differs":
                    str(changed_metric).upper(),
                "anchor_cell_lines":
                    ";".join(sorted(set(members) & anchor_set)) or "",
                "isolate_members":
                    ";".join(m for m in members if m in isolates),
                "n_edges": sub.number_of_edges(),
                "density": round(nx.density(sub), 4),
                "median_degree": round(float(np.median(degrees)), 2),
                "max_degree": int(max(degrees)),
                "median_betweenness":
                    round(float(np.median(btw_unnorm_vals)), 6),
                "max_betweenness": round(float(max(btw_unnorm_vals)), 6),
                "median_betweenness_unnormalised":
                    round(float(np.median(btw_unnorm_vals)), 6),
                "max_betweenness_unnormalised":
                    round(float(max(btw_unnorm_vals)), 6),
                "median_betweenness_normalised":
                    round(float(np.median(btw_norm_vals)), 6),
                "max_betweenness_normalised":
                    round(float(max(btw_norm_vals)), 6),
                "median_edge_support_directions":
                    (round(float(np.median(inc_sds)), 2) if inc_sds else ""),
            })

    ann_df = pd.DataFrame(rows)
    ann_path = Path(outdir) / "resolved_DSMZ_component_annotations.tsv"
    ann_df.to_csv(ann_path, sep="\t", index=False)
    log.info("Wrote: %s", ann_path)

    supp_cols = [
        "cohort", "resolved_component_id", "component_label",
        "n_cell_lines", "member_cell_lines", "central_node",
        "bridge_node", "most_connected_node", "most_connected_degree",
        "central_node_betweenness", "central_node_selection_metric",
        "most_connected_selection_metric",
        "centrality_normalised", "centrality_weighted",
        "centrality_scope", "normalised_vs_unnormalised_anchor_differs",
        "degree_vs_betweenness_node_differs",
        "anchor_cell_lines", "n_edges", "median_degree",
        "median_edge_support_directions",
    ]
    supp_df = ann_df[[c for c in supp_cols if c in ann_df.columns]]
    supp_path = Path(outdir) / "supp_resolved_DSMZ_component_annotations.tsv"
    supp_df.to_csv(supp_path, sep="\t", index=False)
    log.info("Wrote: %s", supp_path)
    return ann_path, supp_path


def write_central_nodes(cohort_data, outdir):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("resolved_available"):
            continue
        G = d["resolved_graph"]
        btw_norm = d.get("resolved_betweenness_normalised", {})
        selected_norm = d.get("selected_by_normalised", {})
        selected_unnorm = d.get("selected_by_unnormalised", {})
        selected_degree = d.get("selected_by_degree", {})
        for cid, (cn, btw) in sorted(d.get("central_nodes", {}).items()):
            if cn is None:
                continue
            degree_node = selected_degree.get(cid)
            rows.append({
                "cohort": cohort,
                "resolved_component_id": cid,
                "central_node": cn,
                "bridge_node": cn,
                "most_connected_node": degree_node or "",
                "most_connected_degree":
                    (G.degree(degree_node) if degree_node else ""),
                "central_node_display": make_display_label(cn),
                "betweenness_centrality": btw,
                "betweenness_unnormalised": btw,
                "betweenness_normalised":
                    round(float(btw_norm.get(cn, 0.0)), 6),
                "selection_metric": "highest_unnormalised_betweenness",
                "centrality_metric": "betweenness",
                "centrality_normalised": "FALSE",
                "centrality_weighted": "FALSE",
                "centrality_scope": "within_component",
                "centrality_tie_break": "highest_value_then_node_id",
                "selected_by_normalised": selected_norm.get(cid) or "",
                "selected_by_unnormalised": selected_unnorm.get(cid) or "",
                "selected_by_degree": degree_node or "",
                "degree_vs_betweenness_node_differs":
                    str(degree_node != cn).upper(),
                "normalised_vs_unnormalised_anchor_differs":
                    str(selected_norm.get(cid) != selected_unnorm.get(cid)).upper(),
            })
    path = Path(outdir) / "resolved_component_central_nodes.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_anchor_centrality_audit(cohort_data, outdir):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("resolved_available"):
            continue
        G = d["resolved_graph"]
        node_to_comp = d["node_to_comp"]
        btw_norm = d.get("resolved_betweenness_normalised", {})
        btw_unnorm = d.get("resolved_betweenness_unnormalised", {})
        selected_norm = d.get("selected_by_normalised", {})
        selected_unnorm = d.get("selected_by_unnormalised", {})
        selected_degree = d.get("selected_by_degree", {})
        for node in sorted(G.nodes()):
            cid = node_to_comp[node]
            sel_norm = selected_norm.get(cid)
            sel_unnorm = selected_unnorm.get(cid)
            sel_degree = selected_degree.get(cid)
            rows.append({
                "cohort": cohort,
                "component_id": cid,
                "node_id": node,
                "display_label": make_display_label(node),
                "degree": G.degree(node),
                "betweenness_normalised":
                    round(float(btw_norm.get(node, 0.0)), 6),
                "betweenness_unnormalised":
                    round(float(btw_unnorm.get(node, 0.0)), 6),
                "selected_by_normalised": str(sel_norm == node).upper(),
                "selected_by_unnormalised": str(sel_unnorm == node).upper(),
                "selected_by_degree": str(sel_degree == node).upper(),
                "canonical_selected": str(sel_unnorm == node).upper(),
                "canonical_bridge_selected": str(sel_unnorm == node).upper(),
                "most_connected_selected": str(sel_degree == node).upper(),
                # Metric-explicit alias columns (schema clarity patch)
                "degree_anchor_selected": str(sel_degree == node).upper(),
                "bridge_betweenness_selected": str(sel_unnorm == node).upper(),
                "anchor_selected": str(
                    (sel_degree == node) or (sel_unnorm == node)
                ).upper(),
                "anchor_selection_reason": (
                    "degree;betweenness_unnormalised"
                    if ((sel_degree == node) and (sel_unnorm == node))
                    else "degree"
                    if (sel_degree == node)
                    else "betweenness_unnormalised"
                    if (sel_unnorm == node)
                    else ""
                ),
                "changed_relative_to_legacy":
                    str(sel_norm != sel_unnorm).upper(),
                "degree_vs_betweenness_node_differs":
                    str(sel_degree != sel_unnorm).upper(),
                "canonical_metric": "degree_and_unnormalised_betweenness",
                "alternative_metric": "normalised_betweenness",
                "centrality_scope": "within_component",
                "centrality_weighted": "FALSE",
                "centrality_tie_break": "highest_value_then_node_id",
            })
    path = Path(outdir) / "anchor_centrality_audit.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_layout_coordinates(cohort_data, outdir, layout_seed):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        for kind in ("similarity", "resolved"):
            G_key = f"{kind}_graph"
            pos_key = f"{kind}_layout"
            avail_key = f"{kind}_available"
            if not d.get(avail_key) or G_key not in d:
                continue
            G = d[G_key]
            pos = d[pos_key]
            for node in sorted(G.nodes()):
                if node not in pos:
                    continue
                attrs = G.nodes[node]
                display = (str(attrs.get("cell_line_display", "")).strip()
                            or make_display_label(node))
                rows.append({
                    "cohort": cohort,
                    "graph_type": kind,
                    "direction": d.get("direction", "") if kind == "similarity" else "",
                    "node_id": node,
                    "display_label": display,
                    "layout_method": "spring_packed_grid",
                    "layout_seed": layout_seed,
                    "x": round(float(pos[node][0]), 6),
                    "y": round(float(pos[node][1]), 6),
                })
    path = Path(outdir) / "cell_line_network_layout_coordinates.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_results_summary(cohort_data, sim_summary, res_summary, outdir):
    lines = [
        "DSMZ Cell Line Networks — Main Figures Results Summary",
        "=" * 78,
        f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "Main Figure 1 — clinical patient sample referenced cell line "
        "similarity networks",
        "-" * 78,
    ]
    sim_by_cohort = {r["cohort"]: r for r in sim_summary}
    for cohort in COHORTS:
        r = sim_by_cohort.get(cohort, {})
        if r.get("n_nodes", "n/a") == "n/a":
            lines.append(f"  {cohort}: similarity network not available.")
            continue
        lines.append(
            f"  {cohort} (winning direction = {r['direction']}): "
            f"{r['n_nodes']} cell lines, {r['n_edges']} edges, "
            f"{r['n_isolates']} isolates, "
            f"{r['n_leiden']} Leiden communities, "
            f"{r['n_louvain']} Louvain communities, "
            f"median Pearson similarity {r['median_sim']}."
        )

    lines += [
        "",
        "Main Figure 2 — graph based consensus resolved DSMZ cell line "
        "neighbourhoods",
        "-" * 78,
    ]
    res_by_cohort = {r["cohort"]: r for r in res_summary}
    for cohort in COHORTS:
        r = res_by_cohort.get(cohort, {})
        if r.get("n_nodes", "n/a") == "n/a":
            lines.append(f"  {cohort}: resolved network not available.")
            continue
        lines.append(
            f"  {cohort}: {r['n_nodes']} cell lines, {r['n_edges']} edges, "
            f"{r['n_components']} components ({r['n_multi']} multi-node, "
            f"{r['n_isolates']} isolates), "
            f"largest component = {r['largest_component']}, "
            f"median degree = {r['median_degree']}."
        )
        cn = r.get("central_nodes", "")
        if cn and cn != "n/a":
            cn_inline = "; ".join(cn.split("\n"))
            lines.append(f"    Central nodes: {cn_inline}.")

    lines += [
        "",
        "Note. Communities in Main Figure 1 are method-specific (per direction). "
        "They are reported as a comparison of Leiden vs Louvain counts only and "
        "must not be compared as raw IDs across cohorts. The final resolved "
        "components in Main Figure 2 are the cohort-level neighbourhood "
        "consensus, independent of any single direction. The separate "
        "final_consensus_all/community_stability/ directory is supplementary "
        "graph stability validation: cross-direction Leiden/Louvain "
        "co-assignment compared with the final resolved components.",
    ]

    path = Path(outdir) / "network_results_summary.txt"
    path.write_text("\n".join(lines))
    log.info("Wrote: %s", path)
    return path


def write_provenance(args, cohort_data, fig_paths, git_commit, outdir):
    input_files = []
    for cohort in ("brca", "nbl", "rbl"):
        for attr in (f"{cohort}_tn_dir",
                     f"{cohort}_resolved",
                     f"{cohort}_anchors"):
            v = getattr(args, attr, None)
            if v:
                input_files.append(str(v))

    enc_fig1 = (
        "node_colour=leiden_community_id_per_cohort_palette; "
        "node_shape=is_isolate(circle/diamond); "
        "node_size=400+120*degree; "
        "edge_width_scaled_by_pearson_similarity; "
        "edge_alpha_scaled_by_pearson_similarity"
    )
    enc_fig2 = (
        "node_colour=connected_component_id"
        "(Wong_CB_safe_palette; isolates=#AAAAAA); "
        "node_shape=is_isolate(circle/diamond); "
        "node_size=400+120*degree; "
        "central_node_border=highest_unnormalised_betweenness_within_component"
        "(thick red border, 4.5_pt); "
        "most_connected_node=highest_degree_within_component_reported_in_tsv; "
        "edge_style=support_type(solid_sd>=2/dashed_sd==1); "
        "edge_width=min(6.0,1.5+0.5*(support_directions-1))"
    )
    row = {
        "script": str(Path(__file__).resolve()),
        "git_commit": git_commit,
        "date_time": datetime.datetime.now().isoformat(),
        "input_files": "; ".join(input_files),
        "output_files": "; ".join(fig_paths),
        "layout_seed": args.layout_seed,
        "layout_method": "spring_packed_grid",
        "fig1_visual_encodings": enc_fig1,
        "fig2_visual_encodings": enc_fig2,
        "fig1_node_colour_column": "leiden_community_id",
        "fig2_node_colour_column": "connected_component_id",
        "fig2_central_node_metric":
            "highest_unnormalised_betweenness_in_component",
        "fig2_most_connected_node_metric": "highest_degree_in_component",
        "centrality_metric": "degree_and_betweenness",
        "centrality_normalised": "FALSE",
        "centrality_normalised_reported": "TRUE",
        "centrality_unnormalised_reported": "TRUE",
        "centrality_weighted": "FALSE",
        "centrality_scope": "within_component",
        "centrality_tie_break": "highest_value_then_node_id",
        "anchor_selection_script": str(Path(__file__).resolve()),
        "anchor_selection_table": "resolved_component_central_nodes.tsv",
        "anchor_audit_table": "anchor_centrality_audit.tsv",
        "canonical_anchor_metric": "degree_and_unnormalised_betweenness",
        "alternative_anchor_metric": "normalised_betweenness",
        "colour_palette": "Wong_CB_safe_8colour",
        "fig1_panel_d_source":
            "similarity_network_construction_summary.tsv;"
            "cell_line_similarity_network_nodes.tsv;"
            "cell_line_similarity_network_edges.tsv",
        "fig2_panel_d_source":
            "resolved_DSMZ_component_annotations.tsv;"
            "resolved_DSMZ_neighbourhood_node_summary.tsv;"
            "resolved_component_central_nodes.tsv;"
            "anchor_centrality_audit.tsv",
    }
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        row[f"direction_{cohort.lower()}"] = d.get("direction", "n/a")
        row[f"direction_source_{cohort.lower()}"] = (
            d.get("direction_source", "n/a"))
        if d.get("similarity_available"):
            G = d["similarity_graph"]
            row[f"sim_n_nodes_{cohort.lower()}"] = G.number_of_nodes()
            row[f"sim_n_edges_{cohort.lower()}"] = G.number_of_edges()
        if d.get("resolved_available"):
            G = d["resolved_graph"]
            row[f"res_n_nodes_{cohort.lower()}"] = G.number_of_nodes()
            row[f"res_n_edges_{cohort.lower()}"] = G.number_of_edges()

    path = Path(outdir) / "figure_provenance.tsv"
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


# =============================================================================
# Argument parser
# =============================================================================

def parse_args():
    ap = argparse.ArgumentParser(
        prog="plot_publication_cell_line_similarity_and_resolved_networks.py",
        description=(
            "Generate the two main thesis figures: Fig 1 = clinical patient "
            "sample referenced cell line similarity networks (panels A, B, C "
            "for BRCA, NBL, RBL with a network construction summary in panel "
            "D); Fig 2 = graph based consensus resolved DSMZ cell line "
            "neighbourhoods (panels A, B, C with a resolved component "
            "summary in panel D). Community stability and co-assignment "
            "outputs are handled separately as supplementary graph stability "
            "validation and are not used as main figure inputs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Strict mode (default) — all three cohorts required:
              python %(prog)s \\
                --brca-tn-dir PATH --brca-resolved PATH \\
                --nbl-tn-dir  PATH --nbl-resolved  PATH \\
                --rbl-tn-dir  PATH --rbl-resolved  PATH \\
                --rbl-anchors PATH --outdir PATH

            Diagnostic — skip absent cohorts:
              python %(prog)s ... --allow-missing-cohorts ...

            The winning direction per cohort is detected from
            <tn-dir>/final_consensus_all/winning_direction.txt; for cohorts
            without that file, the top entry of p_consensus_direction_summary
            that covers all resolved cell lines is used. Override per cohort
            with --<cohort>-direction NAME.
        """),
    )
    for cohort in ("brca", "nbl", "rbl"):
        ap.add_argument(f"--{cohort}-tn-dir", metavar="PATH",
                        help=f"{cohort.upper()} tumour_neighbourhoods directory")
        ap.add_argument(f"--{cohort}-resolved", metavar="PATH",
                        help=f"{cohort.upper()} resolved DSMZ neighbours TSV")
        ap.add_argument(f"--{cohort}-anchors", metavar="PATH",
                        help=f"{cohort.upper()} anchor components TSV (optional)")
        ap.add_argument(f"--{cohort}-direction", metavar="DIRECTION",
                        help=f"Override winning direction for {cohort.upper()}")
    ap.add_argument("--outdir", required=True, metavar="PATH",
                    help="Output directory for all figures and tables")
    ap.add_argument("--allow-missing-cohorts", action="store_true",
                    help="Continue if any cohort cannot be loaded")
    ap.add_argument("--layout-seed", type=int, default=42,
                    help="Random seed for reproducible layouts (default 42)")
    ap.add_argument("--fig-width", type=float, default=26.0,
                    help="Figure width in inches (default 26)")
    ap.add_argument("--fig-height", type=float, default=18.0,
                    help="Figure height in inches (default 18)")
    ap.add_argument("--dpi", type=int, default=300,
                    help="PNG output resolution (default 300)")
    ap.add_argument("--pad", type=float, default=8.0,
                    help="Layout padding between components (default 8.0)")
    ap.add_argument("--log-file", metavar="PATH", default=None,
                    help="Path to log file")
    return ap.parse_args()


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = (
        Path(args.log_file) if args.log_file
        else out_dir /
        "plot_publication_cell_line_similarity_and_resolved_networks.log"
    )
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(log_path, mode="w")],
        force=True,
    )
    log.info("plot_publication_cell_line_similarity_and_resolved_networks.py — "
             "start")
    log.info("Output directory: %s", out_dir.resolve())

    git_commit = get_git_commit()
    log.info("Git commit: %s", git_commit)

    cohort_args = {
        "BRCA": (args.brca_tn_dir, args.brca_resolved,
                 args.brca_anchors, args.brca_direction, None),
        "NBL":  (args.nbl_tn_dir,  args.nbl_resolved,
                 args.nbl_anchors,  args.nbl_direction,  None),
        "RBL":  (args.rbl_tn_dir,  args.rbl_resolved,
                 args.rbl_anchors,  args.rbl_direction,  DEFAULT_RBL_DIRECTION),
    }

    if not args.allow_missing_cohorts:
        missing = []
        for c, (tn, res, *_rest) in cohort_args.items():
            if tn is None or res is None:
                missing.append(c)
        if missing:
            sys.exit(
                "[ERROR] Strict mode requires both --<cohort>-tn-dir and "
                f"--<cohort>-resolved for: {', '.join(missing)}.\n"
                "  Pass --allow-missing-cohorts only for diagnostic runs."
            )
        nonexistent = []
        for c, (tn, res, anc, _direction, _fb) in cohort_args.items():
            for label, p in (("tn-dir", tn), ("resolved", res),
                              ("anchors", anc)):
                if p and not Path(p).exists():
                    nonexistent.append(f"  [{c}] --{c.lower()}-{label}: {p}")
        if nonexistent:
            sys.exit("[ERROR] Input paths not found:\n" +
                     "\n".join(nonexistent))

    cohort_data = {}
    val_rows = []

    # ------------------------------------------------------------------
    # 1. Load per-cohort
    # ------------------------------------------------------------------
    for cohort, (tn_dir, res_path, anc_path,
                 direction_override, default_dir) in cohort_args.items():
        d = {
            "cohort": cohort,
            "tn_dir": tn_dir,
            "similarity_available": False,
            "resolved_available": False,
            "messages": [],
        }

        if res_path is None or not Path(res_path).exists():
            d["messages"].append("Resolved neighbours file missing")
            cohort_data[cohort] = d
            val_rows.append({
                "cohort": cohort, "input_file": str(res_path) or "N/A",
                "scope": "resolved", "n_rows": 0,
                "validation_status": "FAIL",
                "validation_message": "Resolved file missing"})
            continue

        resolved_df, err = load_tsv(res_path)
        if resolved_df is None:
            d["messages"].append(f"Cannot load resolved: {err}")
            cohort_data[cohort] = d
            val_rows.append({
                "cohort": cohort, "input_file": str(res_path),
                "scope": "resolved", "n_rows": 0,
                "validation_status": "FAIL",
                "validation_message": err})
            continue

        # Resolved graph (Main Fig 2)
        rs_resolved, rs_node_stats, rs_edges, rs_anchors = (None, None, None, None)
        if tn_dir and Path(tn_dir).exists():
            rs_resolved, rs_node_stats, rs_edges, rs_anchors = (
                discover_resolved_files(tn_dir))
        node_stats_df, _ = (load_tsv(rs_node_stats)
                             if rs_node_stats else (None, None))
        edges_df, _      = (load_tsv(rs_edges, dtype=None)
                             if rs_edges else (None, None))

        anchor_nodes = set()
        anchor_provided = False
        anchor_source = anc_path or (str(rs_anchors) if rs_anchors else None)
        if anchor_source and Path(anchor_source).exists():
            anchor_provided = True
            anc_df, _ = load_tsv(anchor_source)
            if anc_df is not None:
                for col in ("anchor_cell_line", "anchor", "cell_line"):
                    if col in anc_df.columns:
                        anchor_nodes = (
                            set(anc_df[col].dropna().astype(str).str.strip())
                            - {"nan", ""}
                        )
                        break

        G_res, node_to_comp, isolates_res, _ = build_resolved_graph(
            resolved_df, node_stats_df)
        if G_res is None:
            d["messages"].append("Could not build resolved graph "
                                  "(check column names)")
        else:
            res_pos, _ = pack_components(
                G_res, pad=args.pad, seed=args.layout_seed,
                k_min=1.2, k_max=3.0, t_min=2.6, t_max=6.0)
            res_pos = _repel(res_pos, G_res,
                              min_d=0.18 * 2.6, iters=60, step=0.04,
                              seed=args.layout_seed)
            comp_colours = assign_component_colours(
                G_res, node_to_comp, isolates_res)
            (central_nodes, btw_unnorm, btw_norm, selected_unnorm,
             selected_norm, selected_degree) = central_node_per_component(
                G_res, node_to_comp, set(isolates_res))
            d.update({
                "resolved_available": True,
                "resolved_graph": G_res,
                "node_to_comp": node_to_comp,
                "resolved_isolates": isolates_res,
                "resolved_layout": res_pos,
                "resolved_comp_colours": comp_colours,
                "resolved_edges_df": edges_df,
                "resolved_node_stats_df": node_stats_df,
                "anchor_nodes": anchor_nodes,
                "anchor_path_provided": anchor_provided,
                "anchor_path": anchor_source,
                "central_nodes": central_nodes,
                "resolved_betweenness": btw_unnorm,
                "resolved_betweenness_unnormalised": btw_unnorm,
                "resolved_betweenness_normalised": btw_norm,
                "selected_by_unnormalised": selected_unnorm,
                "selected_by_normalised": selected_norm,
                "selected_by_degree": selected_degree,
            })
            log.info(
                "%s resolved: %d nodes, %d edges, "
                "%d components, %d isolates, %d central nodes",
                cohort, G_res.number_of_nodes(), G_res.number_of_edges(),
                len(set(node_to_comp.values())),
                len(isolates_res),
                sum(1 for v, _ in central_nodes.values() if v),
            )
            val_rows.append({
                "cohort": cohort,
                "input_file": str(res_path),
                "scope": "resolved",
                "n_rows": len(resolved_df),
                "validation_status": "PASS",
                "validation_message":
                    f"{G_res.number_of_nodes()} nodes, "
                    f"{G_res.number_of_edges()} edges",
            })

        # Similarity graph (Main Fig 1) — needs tn_dir
        if not (tn_dir and Path(tn_dir).exists()):
            d["messages"].append(
                "tn_dir not provided; similarity network skipped")
            cohort_data[cohort] = d
            val_rows.append({
                "cohort": cohort, "input_file": str(tn_dir) or "N/A",
                "scope": "similarity", "n_rows": 0,
                "validation_status": "WARN",
                "validation_message":
                    "tn_dir missing; similarity network skipped"})
            continue

        required_cells = set(
            resolved_df["cell_line"].dropna().astype(str).str.strip()
        ) - {"nan", ""}

        direction, dir_source = detect_winning_direction(
            tn_dir, override=direction_override,
            fallback=default_dir, required_cells=required_cells,
        )
        if direction is None:
            d["messages"].append("No winning direction could be determined")
            cohort_data[cohort] = d
            val_rows.append({
                "cohort": cohort, "input_file": str(tn_dir),
                "scope": "similarity", "n_rows": 0,
                "validation_status": "WARN",
                "validation_message":
                    "No winning direction available; similarity panel skipped"})
            continue

        log.info("%s direction: %s (%s)", cohort, direction, dir_source)

        sim_edges_path, sim_nodes_path, sim_iso_path = discover_similarity_files(
            tn_dir, direction)
        if sim_edges_path is None or sim_nodes_path is None:
            d["messages"].append(
                f"Similarity files missing for direction {direction}")
            cohort_data[cohort] = d
            val_rows.append({
                "cohort": cohort, "input_file": str(tn_dir),
                "scope": "similarity", "n_rows": 0,
                "validation_status": "WARN",
                "validation_message":
                    f"Similarity files missing for {direction}"})
            continue

        sim_edges_df, _ = load_tsv(sim_edges_path, dtype=None)
        sim_nodes_df, _ = load_tsv(sim_nodes_path)
        sim_iso_df, _   = (load_tsv(sim_iso_path)
                            if sim_iso_path else (None, None))

        G_sim = build_similarity_graph(sim_edges_df, sim_nodes_df, sim_iso_df)
        if G_sim is None or G_sim.number_of_nodes() == 0:
            d["messages"].append(
                f"Similarity graph empty for direction {direction}")
            cohort_data[cohort] = d
            val_rows.append({
                "cohort": cohort, "input_file": str(sim_nodes_path),
                "scope": "similarity", "n_rows": 0,
                "validation_status": "WARN",
                "validation_message": "Similarity graph empty"})
            continue

        sim_pos, _ = pack_components(
            G_sim, pad=args.pad, seed=args.layout_seed,
            k_min=1.2, k_max=3.0, t_min=2.6, t_max=6.0)
        sim_pos = _repel(sim_pos, G_sim,
                          min_d=0.18 * 2.6, iters=60, step=0.04,
                          seed=args.layout_seed)

        leid_ids = [_leiden_id(G_sim.nodes[n]) for n in G_sim.nodes()]
        sim_isolates = {n for n in G_sim.nodes() if G_sim.degree(n) == 0}
        leid_for_palette = list(dict.fromkeys(leid_ids))
        leiden_colour_map = assign_palette_colours(
            leid_for_palette, isolate_ids={"NA"}, palette=LEIDEN_PALETTE)
        node_colours = {}
        for n in G_sim.nodes():
            lid = _leiden_id(G_sim.nodes[n])
            if n in sim_isolates and lid == "NA":
                node_colours[n] = ISOLATE_COLOUR
            else:
                node_colours[n] = leiden_colour_map.get(lid, ISOLATE_COLOUR)

        d.update({
            "similarity_available": True,
            "direction": direction,
            "direction_source": dir_source,
            "similarity_graph": G_sim,
            "similarity_layout": sim_pos,
            "similarity_isolates": sim_isolates,
            "similarity_node_colours": node_colours,
            "similarity_edges_path": str(sim_edges_path),
            "similarity_nodes_path": str(sim_nodes_path),
        })
        log.info(
            "%s similarity (%s): %d nodes, %d edges, %d isolates",
            cohort, direction,
            G_sim.number_of_nodes(), G_sim.number_of_edges(),
            len(sim_isolates),
        )
        val_rows.append({
            "cohort": cohort,
            "input_file": str(sim_nodes_path),
            "scope": "similarity",
            "n_rows": G_sim.number_of_nodes(),
            "validation_status": "PASS",
            "validation_message":
                f"direction={direction}, "
                f"edges={G_sim.number_of_edges()}, "
                f"isolates={len(sim_isolates)}",
        })

        # Sanity checks
        if required_cells and G_sim.number_of_nodes() > 0:
            missing = required_cells - set(G_sim.nodes())
            if missing:
                val_rows.append({
                    "cohort": cohort,
                    "input_file": str(sim_nodes_path),
                    "scope": "similarity",
                    "n_rows": len(missing),
                    "validation_status": "WARN",
                    "validation_message":
                        f"{len(missing)} resolved cell lines absent from "
                        f"similarity graph for direction {direction}",
                })
                log.warning(
                    "%s: %d resolved cell lines absent from similarity graph "
                    "(direction %s)",
                    cohort, len(missing), direction)

        cohort_data[cohort] = d

    # ------------------------------------------------------------------
    # 2. Validate at least Fig 2 is fully covered for strict mode
    # ------------------------------------------------------------------
    if not args.allow_missing_cohorts:
        missing_res = [c for c in COHORTS
                       if not cohort_data.get(c, {}).get("resolved_available")]
        if missing_res:
            write_validation_report(val_rows, out_dir)
            sys.exit(
                "[ERROR] Resolved network could not be built for: "
                f"{', '.join(missing_res)}. Strict mode aborted."
            )

    write_validation_report(val_rows, out_dir)

    # ------------------------------------------------------------------
    # 3. Write supporting tables
    # ------------------------------------------------------------------
    sim_summary_rows = compute_similarity_summary(cohort_data)
    res_summary_rows = compute_resolved_summary(cohort_data)

    write_similarity_edges(cohort_data, out_dir)
    write_similarity_nodes(cohort_data, out_dir)
    write_similarity_construction_summary(sim_summary_rows, out_dir)
    write_resolved_node_summary(cohort_data, out_dir)
    write_resolved_component_annotations(cohort_data, out_dir)
    write_central_nodes(cohort_data, out_dir)
    write_anchor_centrality_audit(cohort_data, out_dir)
    write_layout_coordinates(cohort_data, out_dir, args.layout_seed)
    write_results_summary(cohort_data, sim_summary_rows, res_summary_rows,
                          out_dir)

    save_kw = dict(bbox_inches="tight", pad_inches=0.3,
                   facecolor="white", transparent=False)

    # ------------------------------------------------------------------
    # 4. Main Figure 1 — similarity networks
    # ------------------------------------------------------------------
    panel_letters = {"BRCA": "A", "NBL": "B", "RBL": "C"}
    fig1 = plt.figure(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)
    fig1.patch.set_facecolor("white")
    gs1 = gridspec.GridSpec(
        2, 2, figure=fig1, hspace=0.30, wspace=0.20,
        left=0.04, right=0.97, top=0.95, bottom=0.04,
    )
    panel_pos = [("BRCA", 0, 0), ("NBL", 0, 1), ("RBL", 1, 0)]
    for cohort, row, col in panel_pos:
        ax = fig1.add_subplot(gs1[row, col])
        d = cohort_data.get(cohort, {})
        if not d.get("similarity_available"):
            ax.set_facecolor("#f5f5f5")
            ax.axis("off")
            ax.text(0.5, 0.5,
                    f"{cohort} similarity network not available",
                    ha="center", va="center", fontsize=12, color="#888888",
                    transform=ax.transAxes)
            ax.set_title(f"Panel {panel_letters[cohort]} — {cohort}",
                         fontsize=11, pad=10)
            continue
        draw_similarity_panel(
            ax, cohort, d["direction"], d["similarity_graph"],
            d["similarity_layout"], d["similarity_node_colours"],
            d["similarity_isolates"], panel_letters[cohort], font_size=8,
        )
    ax_d = fig1.add_subplot(gs1[1, 1])
    draw_similarity_summary_panel(ax_d, sim_summary_rows)

    stem1 = out_dir / (
        "Fig_cell_line_similarity_networks_clinical_patient_referenced_combined")
    fig1.savefig(str(stem1) + ".pdf", **save_kw)
    fig1.savefig(str(stem1) + ".svg", **save_kw)
    fig1.savefig(str(stem1) + ".png", dpi=args.dpi, **save_kw)
    plt.close(fig1)
    log.info("Wrote Main Figure 1: %s.{pdf,svg,png}", stem1)

    # ------------------------------------------------------------------
    # 5. Main Figure 2 — resolved networks
    # ------------------------------------------------------------------
    fig2 = plt.figure(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)
    fig2.patch.set_facecolor("white")
    gs2 = gridspec.GridSpec(
        2, 2, figure=fig2, hspace=0.30, wspace=0.20,
        left=0.04, right=0.97, top=0.95, bottom=0.04,
    )
    for cohort, row, col in panel_pos:
        ax = fig2.add_subplot(gs2[row, col])
        d = cohort_data.get(cohort, {})
        if not d.get("resolved_available"):
            ax.set_facecolor("#f5f5f5")
            ax.axis("off")
            ax.text(0.5, 0.5,
                    f"{cohort} resolved network not available",
                    ha="center", va="center", fontsize=12, color="#888888",
                    transform=ax.transAxes)
            ax.set_title(f"Panel {panel_letters[cohort]} — {cohort}",
                         fontsize=11, pad=10)
            continue
        draw_resolved_panel(
            ax, cohort, d["resolved_graph"], d["resolved_layout"],
            d["node_to_comp"], d["resolved_isolates"],
            d["resolved_comp_colours"], d["central_nodes"],
            d.get("anchor_nodes", set()),
            d.get("resolved_edges_df"),
            d.get("resolved_node_stats_df"),
            panel_letters[cohort], font_size=8,
        )
    ax_d2 = fig2.add_subplot(gs2[1, 1])
    draw_resolved_summary_panel(ax_d2, res_summary_rows)

    stem2 = out_dir / "Fig_consensus_resolved_DSMZ_neighbourhoods_combined"
    fig2.savefig(str(stem2) + ".pdf", **save_kw)
    fig2.savefig(str(stem2) + ".svg", **save_kw)
    fig2.savefig(str(stem2) + ".png", dpi=args.dpi, **save_kw)
    plt.close(fig2)
    log.info("Wrote Main Figure 2: %s.{pdf,svg,png}", stem2)

    # ------------------------------------------------------------------
    # 6. Provenance and final stdout
    # ------------------------------------------------------------------
    fig_paths = []
    for stem in (stem1, stem2):
        for ext in (".pdf", ".svg", ".png"):
            fig_paths.append(str(stem) + ext)
    write_provenance(args, cohort_data, fig_paths, git_commit, out_dir)

    print()
    print("=" * 78)
    print("Cell line similarity + resolved network main figures complete.")
    print("=" * 78)
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if d.get("similarity_available") and d.get("resolved_available"):
            G_sim = d["similarity_graph"]
            G_res = d["resolved_graph"]
            print(f"  {cohort}: "
                  f"similarity {G_sim.number_of_nodes()} nodes / "
                  f"{G_sim.number_of_edges()} edges "
                  f"(direction {d['direction']}); "
                  f"resolved {G_res.number_of_nodes()} nodes / "
                  f"{G_res.number_of_edges()} edges, "
                  f"central nodes per component recorded.")
        else:
            avail = []
            if d.get("similarity_available"):
                avail.append("similarity")
            if d.get("resolved_available"):
                avail.append("resolved")
            print(f"  {cohort}: partial ({', '.join(avail) if avail else 'none'}).")

    print()
    print(f"Output directory: {out_dir.resolve()}")
    print(f"Main Figure 1   : {stem1}.{{pdf,svg,png}}")
    print(f"Main Figure 2   : {stem2}.{{pdf,svg,png}}")
    print(f"Log             : {log_path}")
    print(f"Git commit      : {git_commit}")
    print()
    print("Note: community_stability/ outputs are supplementary graph stability validation, not main figure inputs.")
    print()


if __name__ == "__main__":
    main()
