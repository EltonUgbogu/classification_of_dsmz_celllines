#!/usr/bin/env python3
"""
plot_publication_cell_line_component_networks.py

Standalone publication-quality figure and summary tables for DSMZ cell line
similarity component and community analysis across BRCA, NBL, and RBL cohorts.

Usage (complete-cohort — all three cohorts required):
    python plot_publication_cell_line_component_networks.py \
        --brca-resolved PATH --brca-node-stats PATH --brca-edges PATH \
        --nbl-resolved  PATH --nbl-node-stats  PATH --nbl-edges  PATH \
        --rbl-resolved  PATH --rbl-node-stats  PATH --rbl-edges  PATH \
        --outdir results/publication/cell_line_component_networks

Usage (relaxed — RBL not yet available):
    python plot_publication_cell_line_component_networks.py \
        --brca-resolved PATH --brca-node-stats PATH --brca-edges PATH \
        --nbl-resolved  PATH --nbl-node-stats  PATH --nbl-edges  PATH \
        --allow-missing-cohorts \
        --outdir results/publication/cell_line_component_networks
"""

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

import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import numpy as np
import pandas as pd

mpl.rcParams["text.usetex"] = False
mpl.rcParams["mathtext.default"] = "regular"

try:
    import community as community_louvain
    HAS_LOUVAIN = True
except ImportError:
    HAS_LOUVAIN = False

try:
    import leidenalg
    import igraph as ig
    HAS_LEIDEN = True
except ImportError:
    HAS_LEIDEN = False

CB_PALETTE = [
    "#4477AA", "#EE6677", "#228833", "#CCBB44",
    "#66CCEE", "#AA3377", "#BB5566", "#DDAA33",
]
ISOLATE_COLOUR = "#AAAAAA"
COHORTS = ["BRCA", "NBL", "RBL"]

log = logging.getLogger(__name__)


# =============================================================================
# Layout helpers
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
    return max(k_min, min(k_max, k_min + (math.log(n) / math.log(50)) * (k_max - k_min)))


def _dynamic_box(n, t_min=2.6, t_max=6.0):
    if n <= 1:
        return t_min
    return max(t_min, min(t_max, t_min + (math.log(n) / math.log(50)) * (t_max - t_min)))


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


def pack_components(G, pad=8.0, seed=42, k_min=1.2, k_max=3.0, t_min=2.6, t_max=6.0):
    """Spring layout per component, tiled on a grid. Returns (pos, comps_list)."""
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
            lp = _repel(lp, sub, min_d=0.55 * target, iters=140, step=0.10, seed=seed)
        row_, col_ = i // cols, i % cols
        pad_i = pad * (1.0 + 0.15 * math.log(max(n, 2)))
        offset = np.array([col_ * pad_i, -row_ * pad_i], float)
        for nd, p in lp.items():
            pos[nd] = np.array(p, float) + offset
    return pos, comps


# =============================================================================
# Utilities
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


def check_ambiguous_rbl_ids(ids):
    """Return RBL IDs missing the required replicate suffix (bare RBL_15 or RBL_20)."""
    pat = re.compile(r"^RBL_(15|20)$", re.IGNORECASE)
    return [x for x in ids if pat.match(str(x).strip())]


def cross_cohort_overlap(cohort_id_sets):
    """Return list of (cohort_a, cohort_b, n_overlap, examples) for any ID overlaps."""
    problems = []
    cohorts = sorted(cohort_id_sets)
    for i in range(len(cohorts)):
        for j in range(i + 1, len(cohorts)):
            ca, cb = cohorts[i], cohorts[j]
            shared = cohort_id_sets[ca] & cohort_id_sets[cb]
            if shared:
                problems.append((ca, cb, len(shared), sorted(shared)[:10]))
    return problems


def load_tsv(path, dtype=str):
    """Return (DataFrame, None) on success or (None, error_string) on failure."""
    if path is None:
        return None, "no path provided"
    p = Path(path)
    if not p.exists():
        return None, f"file not found: {p}"
    try:
        return pd.read_csv(p, sep="\t", dtype=dtype), None
    except Exception as e:
        return None, f"read error: {e}"


def infer_id_style(ids):
    pat = re.compile(r"^NG-\d+_")
    n_long = sum(1 for x in ids if pat.match(str(x)))
    if n_long == len(ids):
        return "long"
    if n_long == 0:
        return "short"
    return "mixed"


def make_display_label(cell_line_id):
    """Shorten long NG-XXXXX_CELLLINE_libNNNN_PLATE_REP IDs to CELLLINE_rREP."""
    m = re.match(r"^NG-\d+_(.+)_lib\d+_\d+_(\d+)$", str(cell_line_id))
    if m:
        return f"{m.group(1)}_r{m.group(2)}"
    return cell_line_id


# =============================================================================
# Per-cohort validation
# =============================================================================

def validate_cohort(cohort, resolved_df, node_stats_df, edges_df, paths):
    """
    Return list of validation row dicts for network_input_validation.tsv.
    Columns: cohort, input_file, n_rows, n_unique_nodes, id_style,
             validation_status (PASS/WARN/FAIL), validation_message.
    """
    res_path, ns_path, ed_path = paths
    rows = []

    def _row(fpath, n_rows, n_uniq, id_style, status, msg):
        rows.append(dict(
            cohort=cohort, input_file=str(fpath) if fpath else "N/A",
            n_rows=n_rows, n_unique_nodes=n_uniq, id_style=id_style,
            validation_status=status, validation_message=msg,
        ))

    # Resolved neighbours
    if resolved_df is None:
        _row(res_path, 0, 0, "unknown", "FAIL", "File not loaded")
        return rows

    has_nb = ("final_neighbors" in resolved_df.columns or
              "final_neighbours" in resolved_df.columns)
    if "cell_line" not in resolved_df.columns or not has_nb:
        _row(res_path, len(resolved_df), 0, "unknown", "FAIL",
             f"Missing required columns. Found: {list(resolved_df.columns)}")
        return rows

    ids = resolved_df["cell_line"].dropna().astype(str).str.strip().tolist()
    n_uniq = len(set(ids))
    _row(res_path, len(resolved_df), n_uniq, infer_id_style(ids), "PASS",
         "Required columns present")

    # Ambiguous RBL replicate IDs
    if cohort == "RBL":
        ambig = check_ambiguous_rbl_ids(ids)
        if ambig:
            _row(res_path, len(resolved_df), n_uniq, infer_id_style(ids), "FAIL",
                 f"Ambiguous bare RBL replicate IDs (missing suffix _1/_2/_3/_4): "
                 f"{ambig[:5]}{'...' if len(ambig) > 5 else ''}")

    # Node stats consistency
    if node_stats_df is not None and ns_path:
        if "cell_line" not in node_stats_df.columns:
            _row(ns_path, len(node_stats_df), 0, "unknown", "FAIL",
                 f"Missing 'cell_line' column. Found: {list(node_stats_df.columns)}")
        else:
            ns_ids = set(node_stats_df["cell_line"].astype(str).str.strip()) - {"nan", ""}
            res_ids = set(ids)
            only_ns = ns_ids - res_ids
            only_res = res_ids - ns_ids
            msgs = []
            if only_ns:
                msgs.append(f"{len(only_ns)} IDs only in node_stats: "
                            f"{sorted(only_ns)[:3]}")
            if only_res:
                msgs.append(f"{len(only_res)} IDs only in resolved: "
                            f"{sorted(only_res)[:3]}")
            status = "PASS" if not (only_ns or only_res) else "WARN"
            _row(ns_path, len(node_stats_df), len(ns_ids),
                 infer_id_style(list(ns_ids)), status,
                 "; ".join(msgs) if msgs else "OK")
            if "cell_line_display" in node_stats_df.columns:
                disp = node_stats_df["cell_line_display"].dropna().astype(str)
                dups = disp[disp.duplicated()].tolist()
                if dups and rows:
                    rows[-1]["validation_message"] += (
                        f"; Duplicated display labels: {dups[:5]}")
                    rows[-1]["validation_status"] = "WARN"

    # Edge file
    if edges_df is not None and ed_path:
        if not {"node1", "node2"}.issubset(edges_df.columns):
            _row(ed_path, len(edges_df), 0, "unknown", "FAIL",
                 f"Missing node1/node2 columns. Found: {list(edges_df.columns)}")
        else:
            edge_ids = (
                set(edges_df["node1"].astype(str).str.strip()) |
                set(edges_df["node2"].astype(str).str.strip())
            ) - {"nan", ""}
            res_ids = set(ids)
            extra = edge_ids - res_ids
            msgs = (
                [f"{len(extra)} edge node IDs absent from resolved"]
                if extra else []
            )
            _row(ed_path, len(edges_df), len(edge_ids),
                 infer_id_style(list(edge_ids)),
                 "WARN" if extra else "PASS",
                 "; ".join(msgs) if msgs else "OK")

    return rows


# =============================================================================
# Graph construction
# =============================================================================

def build_graph(resolved_df, node_stats_df=None):
    """
    Build resolved networkx Graph from resolved neighbours TSV.
    Returns (G, node_to_comp, isolates, nb_col).
    - node_to_comp: {node: connected_component_id}
    - isolates: list of nodes with degree 0 in the resolved graph
    """
    nb_col = next(
        (c for c in ("final_neighbours", "final_neighbors")
         if c in resolved_df.columns), None
    )
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


# =============================================================================
# Community detection
# =============================================================================

def compute_communities(G, node_stats_df=None):
    """
    Return {node: {"louvain_community_id": str, "leiden_community_id": str}}.
    Reads from node_stats columns community_louv / community_leid when present.
    Falls back to computing Louvain or Leiden if packages are available.
    """
    _NA = "not_available"
    result = {n: {"louvain_community_id": _NA, "leiden_community_id": _NA}
              for n in G.nodes()}

    has_louv_col = (node_stats_df is not None and
                    "community_louv" in node_stats_df.columns)
    has_leid_col = (node_stats_df is not None and
                    "community_leid" in node_stats_df.columns)

    if has_louv_col or has_leid_col:
        ns = node_stats_df.copy()
        ns["cell_line"] = ns["cell_line"].astype(str).str.strip()
        for _, row in ns.iterrows():
            n = str(row["cell_line"]).strip()
            if n in result:
                if has_louv_col:
                    val = str(row.get("community_louv", "")).strip()
                    result[n]["louvain_community_id"] = val if val and val != "nan" else _NA
                if has_leid_col:
                    val = str(row.get("community_leid", "")).strip()
                    result[n]["leiden_community_id"] = val if val and val != "nan" else _NA
        return result

    if HAS_LOUVAIN and G.number_of_edges() > 0:
        partition = community_louvain.best_partition(G)
        for n, cid in partition.items():
            if n in result:
                result[n]["louvain_community_id"] = str(cid)

    if HAS_LEIDEN and G.number_of_edges() > 0:
        try:
            nx_nodes = list(G.nodes())
            ig_g = ig.Graph.TupleList(
                [(str(u), str(v)) for u, v in G.edges()], directed=False
            )
            for n in nx_nodes:
                try:
                    ig_g.vs.find(name=str(n))
                except ValueError:
                    ig_g.add_vertex(name=str(n))
            part = leidenalg.find_partition(
                ig_g, leidenalg.ModularityVertexPartition)
            for cid, members in enumerate(part):
                for m in members:
                    n = ig_g.vs[m]["name"]
                    if n in result:
                        result[n]["leiden_community_id"] = str(cid)
        except Exception as exc:
            log.warning("Leiden community detection failed: %s", exc)

    return result


# =============================================================================
# Colour assignment
# =============================================================================

def assign_component_colours(G, node_to_comp, isolates):
    """Return {comp_id: hex_colour}. Single-member or all-isolate comps → ISOLATE_COLOUR."""
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


# =============================================================================
# Drawing — cohort panel
# =============================================================================

def draw_cohort_panel(ax, cohort_label, G, pos, node_to_comp, isolates,
                      comp_colour_map, anchor_nodes, edges_df, node_stats_df,
                      panel_letter, font_size=8):
    isolate_set = set(isolates)
    anchor_set = set(anchor_nodes) if anchor_nodes else set()

    # Edge support lookup
    edge_support = {}
    if (edges_df is not None and
            {"node1", "node2", "support_directions"}.issubset(edges_df.columns)):
        for _, er in edges_df.iterrows():
            n1 = str(er["node1"]).strip()
            n2 = str(er["node2"]).strip()
            try:
                edge_support[tuple(sorted((n1, n2)))] = float(
                    er["support_directions"])
            except (ValueError, TypeError):
                pass

    def _node_colour(n):
        return comp_colour_map.get(node_to_comp.get(n, -1), ISOLATE_COLOUR)

    def _node_size(n):
        return max(300, 400 + 120 * G.degree(n))

    # Classify edges
    multi_edges, single_edges, fallback_edges = [], [], []
    for u, v in G.edges():
        key = tuple(sorted((u, v)))
        if key in edge_support:
            sd = edge_support[key]
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

    # Non-isolate nodes
    regular = [n for n in G.nodes() if n not in isolate_set]
    anchor_reg = [n for n in regular if n in anchor_set]
    non_anchor_reg = [n for n in regular if n not in anchor_set]

    if non_anchor_reg:
        nx.draw_networkx_nodes(
            G, pos, nodelist=non_anchor_reg,
            node_color=[_node_colour(n) for n in non_anchor_reg],
            node_size=[_node_size(n) for n in non_anchor_reg],
            linewidths=1.0, edgecolors="#333333", ax=ax,
        )
    if anchor_reg:
        nx.draw_networkx_nodes(
            G, pos, nodelist=anchor_reg,
            node_color=[_node_colour(n) for n in anchor_reg],
            node_size=[_node_size(n) for n in anchor_reg],
            linewidths=6.0, edgecolors="#333333", ax=ax,
        )

    # Isolate nodes — diamond marker via scatter
    iso_in_pos = [n for n in isolates if n in pos]
    if iso_in_pos:
        iso_x = [pos[n][0] for n in iso_in_pos]
        iso_y = [pos[n][1] for n in iso_in_pos]
        iso_c = [_node_colour(n) for n in iso_in_pos]
        iso_s = [_node_size(n) * 0.55 for n in iso_in_pos]
        ax.scatter(iso_x, iso_y, s=iso_s, c=iso_c,
                   marker="D", edgecolors="#333333",
                   linewidths=1.0, zorder=4)

    # Labels
    disp_map = {}
    if (node_stats_df is not None and
            "cell_line" in node_stats_df.columns and
            "cell_line_display" in node_stats_df.columns):
        ns = node_stats_df.copy()
        ns["cell_line"] = ns["cell_line"].astype(str).str.strip()
        for _, row in ns.iterrows():
            disp_map[str(row["cell_line"]).strip()] = str(
                row["cell_line_display"])

    labels = {n: (disp_map.get(n) or make_display_label(n))
              for n in G.nodes() if n in pos}
    nx.draw_networkx_labels(
        G, pos, labels=labels,
        font_size=max(5, font_size - 2), font_color="#1a1a1a", ax=ax,
    )

    ax.set_title(f"Panel {panel_letter} — {cohort_label}",
                 fontsize=font_size + 2, pad=8)
    ax.axis("off")

    # Per-panel mini legend
    handles = [
        mpatches.Patch(facecolor="#AAAAAA", edgecolor="#333333",
                       linewidth=1.0, label="Non-anchor cell line (thin border)"),
        mpatches.Patch(facecolor="#AAAAAA", edgecolor="#333333",
                       linewidth=6.0, label="Anchor cell line (thick border)"),
        mlines.Line2D([], [], marker="D", color="none",
                      markerfacecolor=ISOLATE_COLOUR, markeredgecolor="#333333",
                      markersize=7, label="Isolate (diamond, degree 0)"),
        mlines.Line2D([], [], color="#555555", linewidth=2.0,
                      linestyle="solid", label="Multi supported (sd>=2, solid)"),
        mlines.Line2D([], [], color="#555555", linewidth=1.5,
                      linestyle="dashed", label="Single supported (sd=1, dashed)"),
    ]
    ax.legend(handles=handles, loc="lower left",
              fontsize=max(5, font_size - 2),
              framealpha=0.8, ncol=1, handlelength=1.8)


# =============================================================================
# Drawing — legend / summary panel
# =============================================================================

def draw_legend_panel(ax, all_comp_colours, cohort_comp_data):
    ax.axis("off")
    legend_ax = ax.inset_axes([0.0, 0.42, 1.0, 0.56])
    bar_ax = ax.inset_axes([0.05, 0.0, 0.9, 0.38])
    legend_ax.axis("off")

    handles = []

    # Component colour patches (deduplicated)
    seen = set()
    for (cohort, cid), colour in sorted(all_comp_colours.items()):
        label = (f"{cohort} isolate" if colour == ISOLATE_COLOUR
                 else f"{cohort} component {cid}")
        if label not in seen:
            handles.append(mpatches.Patch(
                facecolor=colour, edgecolor="#333333",
                linewidth=0.8, label=label))
            seen.add(label)

    # Node shape
    handles.append(mlines.Line2D(
        [], [], marker="o", color="none",
        markerfacecolor="#AAAAAA", markeredgecolor="#333333",
        markersize=8, label="Non-isolate node (circle)"))
    handles.append(mlines.Line2D(
        [], [], marker="D", color="none",
        markerfacecolor=ISOLATE_COLOUR, markeredgecolor="#333333",
        markersize=8, label="Isolate node (diamond, degree 0)"))

    # Node size
    handles.append(mlines.Line2D(
        [], [], marker="o", color="none",
        markerfacecolor="#AAAAAA", markeredgecolor="#333333",
        markersize=5, label="Low degree (small)"))
    handles.append(mlines.Line2D(
        [], [], marker="o", color="none",
        markerfacecolor="#AAAAAA", markeredgecolor="#333333",
        markersize=11, label="High degree (large)"))

    # Node border
    handles.append(mpatches.Patch(
        facecolor="#AAAAAA", edgecolor="#333333",
        linewidth=1.0, label="Non-anchor cell line (thin border)"))
    handles.append(mpatches.Patch(
        facecolor="#AAAAAA", edgecolor="#333333",
        linewidth=6.0, label="Anchor cell line (thick border)"))

    # Edge line type
    handles.append(mlines.Line2D(
        [], [], color="#555555", linewidth=2.5,
        linestyle="solid", label="Multi supported edge (sd>=2, solid)"))
    handles.append(mlines.Line2D(
        [], [], color="#555555", linewidth=1.5,
        linestyle="dashed", label="Single supported edge (sd=1, dashed)"))

    # Edge width
    handles.append(mlines.Line2D(
        [], [], color="#555555", linewidth=1.5,
        label="Narrow edge (low support)"))
    handles.append(mlines.Line2D(
        [], [], color="#555555", linewidth=4.5,
        label="Wide edge (high support)"))

    legend_ax.legend(
        handles=handles, loc="upper left",
        fontsize=6.5, ncol=2,
        title="Visual encoding guide", title_fontsize=7.5,
        framealpha=0.9, handlelength=2.0,
    )
    legend_ax.set_title("Panel D — Visual encoding and cohort summary",
                        fontsize=9, pad=6)

    # Stacked bar chart
    cohorts_present = [c for c in COHORTS if c in cohort_comp_data]
    if cohorts_present:
        for yi, cohort in enumerate(cohorts_present):
            left = 0
            for cid, n in sorted(cohort_comp_data[cohort].items()):
                colour = all_comp_colours.get((cohort, cid), ISOLATE_COLOUR)
                bar_ax.barh(yi, n, left=left, color=colour,
                            edgecolor="white", linewidth=0.5, height=0.6)
                if n >= 2:
                    bar_ax.text(left + n / 2, yi, str(cid),
                                ha="center", va="center",
                                fontsize=6, fontweight="bold", color="black")
                left += n
        bar_ax.set_yticks(list(range(len(cohorts_present))))
        bar_ax.set_yticklabels(cohorts_present, fontsize=8)
        bar_ax.set_xlabel("Number of cell lines", fontsize=8)
        bar_ax.set_title("Component composition per cohort", fontsize=8, pad=4)
        for spine in ["top", "right"]:
            bar_ax.spines[spine].set_visible(False)


# =============================================================================
# Summary table writers
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


def write_validation_report(val_rows, outdir):
    path = Path(outdir) / "network_input_validation.tsv"
    pd.DataFrame(val_rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_edge_table(cohort_data, graphs, outdir):
    """
    cell_line_network_edges.tsv — one row per edge with all plotted visual properties.
    Columns: cohort, node1, node2, support_directions, support_type,
             edge_style (solid/dashed), edge_width_plotted.
    This makes edge width and edge line type fully traceable to machine readable values.
    """
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("available") or graphs.get(cohort) is None:
            continue
        G = graphs[cohort]["G"]
        es_map = _edge_support_map(d.get("edges_df"))

        for u, v in G.edges():
            key = tuple(sorted((u, v)))
            sd = es_map.get(key)
            if sd is not None:
                support_type = "multi supported" if sd >= 2 else "single supported"
                edge_style = "solid" if sd >= 2 else "dashed"
                edge_width = round(min(6.0, 1.5 + 0.5 * (sd - 1)), 3)
            else:
                support_type = "not_available"
                edge_style = "solid"
                edge_width = 1.5
            rows.append({
                "cohort": cohort,
                "node1": u,
                "node2": v,
                "support_directions": sd if sd is not None else "not_available",
                "support_type": support_type,
                "edge_style": edge_style,
                "edge_width_plotted": edge_width,
            })

    path = Path(outdir) / "cell_line_network_edges.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_component_community_summary(cohort_data, graphs, community_maps, outdir):
    """
    cell_line_component_community_summary.tsv
    Columns clearly distinguish connected_component_id, resolved_component_id,
    louvain_community_id, and leiden_community_id.
    """
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("available") or graphs.get(cohort) is None:
            continue
        gd = graphs[cohort]
        G = gd["G"]
        node_to_comp = gd["node_to_comp"]
        isolate_set = set(gd["isolates"])
        anchor_set = set(d.get("anchor_nodes", set()))
        anchor_path_provided = d.get("anchor_path_provided", False)
        ns_df = d.get("node_stats_df")
        es_map = _edge_support_map(d.get("edges_df"))
        comm_map = community_maps.get(cohort, {})

        comp_sizes = defaultdict(int)
        for n in G.nodes():
            comp_sizes[node_to_comp[n]] += 1

        disp_map = {}
        if ns_df is not None and "cell_line_display" in ns_df.columns:
            for _, row in ns_df.iterrows():
                disp_map[str(row["cell_line"]).strip()] = str(
                    row["cell_line_display"])

        btw = nx.betweenness_centrality(G) if G.number_of_nodes() > 0 else {}

        for node in sorted(G.nodes()):
            cid = node_to_comp[node]
            inc_sd = [
                es_map[tuple(sorted((node, nb)))]
                for nb in G.neighbors(node)
                if tuple(sorted((node, nb))) in es_map
            ]
            comm = comm_map.get(node, {})
            # is_anchor: NA when no anchor file was provided; TRUE/FALSE otherwise
            if not anchor_path_provided:
                is_anchor_val = "NA"
            else:
                is_anchor_val = str(node in anchor_set).upper()
            display = disp_map.get(node) or make_display_label(node)
            rows.append({
                "cohort": cohort,
                "cell_line_raw_id": node,
                "cell_line_display": display,
                # node colour in figure maps to connected_component_id:
                "connected_component_id": cid,
                "connected_component_size": comp_sizes[cid],
                # community columns (not_available when detection not run):
                "louvain_community_id": comm.get("louvain_community_id", "not_available"),
                "leiden_community_id": comm.get("leiden_community_id", "not_available"),
                # node properties used as visual encodings:
                "degree": G.degree(node),
                "is_isolate": str(node in isolate_set).upper(),
                "is_anchor": is_anchor_val,
                "betweenness_centrality": round(btw.get(node, 0.0), 6),
                "median_incident_edge_support": (
                    round(float(np.median(inc_sd)), 3) if inc_sd else ""),
            })

    path = Path(outdir) / "cell_line_component_community_summary.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_component_annotations(cohort_data, graphs, outdir):
    """
    cell_line_component_annotations.tsv and supp_cell_line_component_annotations.tsv.
    Neutral labels only: Component N, Small component, Isolate.
    """
    ann_rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("available") or graphs.get(cohort) is None:
            continue
        gd = graphs[cohort]
        G = gd["G"]
        node_to_comp = gd["node_to_comp"]
        isolate_set = set(gd["isolates"])
        anchor_set = set(d.get("anchor_nodes", set()))
        es_map = _edge_support_map(d.get("edges_df"))
        btw = nx.betweenness_centrality(G) if G.number_of_nodes() > 0 else {}

        comp_members = defaultdict(list)
        for n in G.nodes():
            comp_members[node_to_comp[n]].append(n)

        for cid, members in sorted(comp_members.items()):
            n = len(members)
            all_iso = all(m in isolate_set for m in members)

            if all_iso or n == 1:
                label = "Isolate"
            elif n <= 3:
                label = "Small component"
            else:
                label = f"Component {cid}"

            anchor, anchor_metric = "", ""
            if n > 1:
                named = [m for m in members if m in anchor_set]
                if named:
                    anchor, anchor_metric = named[0], "explicit_anchor"
                else:
                    best = max(members, key=lambda m: btw.get(m, 0))
                    anchor = best
                    anchor_metric = f"betweenness={btw.get(best, 0):.4f}"
            elif members:
                anchor, anchor_metric = members[0], "single_member"

            sub = G.subgraph(members)
            inc_sds = [
                es_map[tuple(sorted((m, nb)))]
                for m in members
                for nb in G.neighbors(m)
                if tuple(sorted((m, nb))) in es_map
            ]
            degrees = [G.degree(m) for m in members]
            btw_vals = [btw.get(m, 0) for m in members]

            ann_rows.append({
                "cohort": cohort,
                "resolved_component_id": cid,
                "component_label": label,
                "n_cell_lines": n,
                "member_cell_lines": ";".join(sorted(members)),
                "anchor_cell_line": anchor,
                "anchor_selection_metric": anchor_metric,
                "isolate_members": ";".join(
                    m for m in members if m in isolate_set),
                "n_edges": sub.number_of_edges(),
                "density": round(nx.density(sub), 4),
                "median_degree": round(float(np.median(degrees)), 2),
                "max_degree": int(max(degrees)),
                "median_betweenness": round(float(np.median(btw_vals)), 6),
                "max_betweenness": round(float(max(btw_vals)), 6),
                "median_edge_support_directions": (
                    round(float(np.median(inc_sds)), 2) if inc_sds else ""),
                "annotation_status": "neutral_label",
            })

    ann_df = pd.DataFrame(ann_rows)
    ann_path = Path(outdir) / "cell_line_component_annotations.tsv"
    ann_df.to_csv(ann_path, sep="\t", index=False)
    log.info("Wrote: %s", ann_path)

    supp_cols = [
        "cohort", "component_label", "n_cell_lines", "member_cell_lines",
        "anchor_cell_line", "isolate_members",
        "median_degree", "median_edge_support_directions",
    ]
    supp_df = ann_df[[c for c in supp_cols if c in ann_df.columns]]
    supp_path = Path(outdir) / "supp_cell_line_component_annotations.tsv"
    supp_df.to_csv(supp_path, sep="\t", index=False)
    log.info("Wrote: %s", supp_path)

    return ann_path, supp_path


def write_layout_coordinates(cohort_data, graphs, layouts, args, outdir):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("available") or graphs.get(cohort) is None:
            continue
        G = graphs[cohort]["G"]
        pos = layouts[cohort]
        ns_df = d.get("node_stats_df")
        disp_map = {}
        if ns_df is not None and "cell_line_display" in ns_df.columns:
            for _, row in ns_df.iterrows():
                disp_map[str(row["cell_line"]).strip()] = str(
                    row["cell_line_display"])
        for node in sorted(G.nodes()):
            if node not in pos:
                continue
            rows.append({
                "cohort": cohort,
                "graph_type": "resolved",
                "node_id": node,
                "display_label": disp_map.get(node, node),
                "layout_method": "spring_packed_grid",
                "layout_seed": args.layout_seed,
                "x": round(float(pos[node][0]), 6),
                "y": round(float(pos[node][1]), 6),
            })
    path = Path(outdir) / "cell_line_network_layout_coordinates.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_comparison_table(cohort_data, graphs, outdir):
    rows = []
    for cohort in COHORTS:
        d = cohort_data.get(cohort, {})
        if not d.get("available") or graphs.get(cohort) is None:
            continue
        G = graphs[cohort]["G"]
        ed_df = d.get("edges_df")
        comps = list(nx.connected_components(G))
        iso_comps = [c for c in comps if len(c) == 1]
        sds = (
            ed_df["support_directions"].dropna().tolist()
            if ed_df is not None and "support_directions" in ed_df.columns
            else []
        )
        rows.append({
            "cohort": cohort, "graph_type": "resolved",
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "n_connected_components": len(comps),
            "n_isolate_nodes": len(iso_comps),
            "n_multi_node_components": len(comps) - len(iso_comps),
            "density": round(nx.density(G), 6),
            "median_edge_support_directions": (
                round(float(np.median(sds)), 2) if sds else ""),
        })
        # Optional consensus edges in same directory
        if d.get("edges_path"):
            ep = Path(d["edges_path"])
            for name in [
                "dsmz_cellline_similarity_consensus_edges.tsv",
                "consensus_edges.tsv",
            ]:
                cp = ep.parent / name
                if cp.exists():
                    try:
                        cons = pd.read_csv(cp, sep="\t")
                        if {"node1", "node2"}.issubset(cons.columns):
                            Gc = nx.Graph()
                            for _, er in cons.iterrows():
                                Gc.add_edge(str(er["node1"]).strip(),
                                            str(er["node2"]).strip())
                            cc = list(nx.connected_components(Gc))
                            ic = [c for c in cc if len(c) == 1]
                            csds = (
                                cons["support_directions"].dropna().tolist()
                                if "support_directions" in cons.columns else []
                            )
                            rows.append({
                                "cohort": cohort, "graph_type": "consensus",
                                "n_nodes": Gc.number_of_nodes(),
                                "n_edges": Gc.number_of_edges(),
                                "n_connected_components": len(cc),
                                "n_isolate_nodes": len(ic),
                                "n_multi_node_components": len(cc) - len(ic),
                                "density": round(nx.density(Gc), 6),
                                "median_edge_support_directions": (
                                    round(float(np.median(csds)), 2)
                                    if csds else ""),
                            })
                    except Exception:
                        pass
                    break

    path = Path(outdir) / "consensus_resolved_comparison.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


def write_provenance(args, cohort_data, graphs, fig_paths, git_commit, outdir):
    input_files = []
    for attr in [
        "brca_resolved", "brca_node_stats", "brca_edges", "brca_anchors",
        "nbl_resolved", "nbl_node_stats", "nbl_edges", "nbl_anchors",
        "rbl_resolved", "rbl_node_stats", "rbl_edges", "rbl_anchors",
    ]:
        val = getattr(args, attr, None)
        if val:
            input_files.append(str(val))

    per_cohort = {}
    for cohort in COHORTS:
        if graphs.get(cohort) is not None:
            G = graphs[cohort]["G"]
            per_cohort[cohort] = {
                "n_nodes": G.number_of_nodes(),
                "n_edges": G.number_of_edges(),
            }

    enc = (
        "node_colour=connected_component_id"
        "(Wong_CB_safe_palette; isolates=#AAAAAA); "
        "node_shape=is_isolate(circle/diamond); "
        "node_size=400+120*degree; "
        "node_border_width=is_anchor"
        "(non_anchor=1.0_pt/anchor=6.0_pt); "
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
        "n_cohorts_processed": sum(1 for c in COHORTS
                                   if graphs.get(c) is not None),
        "n_nodes_total": sum(v["n_nodes"] for v in per_cohort.values()),
        "n_edges_total": sum(v["n_edges"] for v in per_cohort.values()),
        "n_nodes_brca": per_cohort.get("BRCA", {}).get("n_nodes", 0),
        "n_edges_brca": per_cohort.get("BRCA", {}).get("n_edges", 0),
        "n_nodes_nbl": per_cohort.get("NBL", {}).get("n_nodes", 0),
        "n_edges_nbl": per_cohort.get("NBL", {}).get("n_edges", 0),
        "n_nodes_rbl": per_cohort.get("RBL", {}).get("n_nodes", 0),
        "n_edges_rbl": per_cohort.get("RBL", {}).get("n_edges", 0),
        "output_width_in": args.fig_width,
        "output_height_in": args.fig_height,
        "output_dpi": args.dpi,
        "isolate_colour": ISOLATE_COLOUR,
        "anchor_border_linewidth": 6.0,
        "non_anchor_border_linewidth": 1.0,
        "node_colour_column": "connected_component_id",
        "node_size_formula": "400+120*degree",
        "edge_width_formula": "min(6.0,1.5+0.5*(support_directions-1))",
        "edge_line_type": "solid_if_sd>=2; dashed_if_sd==1",
        "visual_encodings": enc,
        "colour_palette": "Wong_CB_safe_8colour",
    }
    path = Path(outdir) / "figure_provenance.tsv"
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
    log.info("Wrote: %s", path)
    return path


# =============================================================================
# Argument parser
# =============================================================================

def parse_args():
    ap = argparse.ArgumentParser(
        prog="plot_publication_cell_line_component_networks.py",
        description=(
            "Standalone publication-quality network figure and summary tables "
            "for DSMZ cell line similarity component and community analysis."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Complete-cohort mode (default) — all three cohorts required:
              python %(prog)s \\
                --brca-resolved PATH --brca-node-stats PATH --brca-edges PATH \\
                --nbl-resolved  PATH --nbl-node-stats  PATH --nbl-edges  PATH \\
                --rbl-resolved  PATH --rbl-node-stats  PATH --rbl-edges  PATH \\
                --outdir results/publication/cell_line_component_networks

            Relaxed mode — skip absent cohorts:
              python %(prog)s \\
                --brca-resolved PATH --brca-node-stats PATH --brca-edges PATH \\
                --allow-missing-cohorts \\
                --outdir results/publication/cell_line_component_networks
        """),
    )

    for cohort in ("brca", "nbl", "rbl"):
        ap.add_argument(f"--{cohort}-resolved", metavar="PATH",
                        help=f"{cohort.upper()} resolved DSMZ neighbours TSV")
        ap.add_argument(f"--{cohort}-node-stats", metavar="PATH",
                        help=f"{cohort.upper()} node statistics TSV (optional)")
        ap.add_argument(f"--{cohort}-edges", metavar="PATH",
                        help=f"{cohort.upper()} edge statistics TSV (optional)")
        ap.add_argument(f"--{cohort}-anchors", metavar="PATH",
                        help=f"{cohort.upper()} anchor components TSV (optional)")

    ap.add_argument("--outdir", required=True, metavar="PATH",
                    help="Output directory for all figures and tables")
    ap.add_argument(
        "--allow-missing-cohorts", action="store_true",
        help=(
            "Continue if NBL or RBL resolved files are not provided. "
            "Default: stop with an informative error."
        ),
    )
    ap.add_argument("--layout-seed", type=int, default=42,
                    help="Random seed for reproducible layout (default: 42)")
    ap.add_argument("--fig-width", type=float, default=26.0,
                    help="Figure width in inches (default: 26)")
    ap.add_argument("--fig-height", type=float, default=18.0,
                    help="Figure height in inches (default: 18)")
    ap.add_argument("--dpi", type=int, default=300,
                    help="PNG output resolution (default: 300)")
    ap.add_argument("--pad", type=float, default=8.0,
                    help="Layout padding between components (default: 8.0)")
    ap.add_argument("--log-file", metavar="PATH", default=None,
                    help="Path to log file (default: outdir/plot_publication_*.log)")
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
        else out_dir / "plot_publication_cell_line_component_networks.log"
    )
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_path),
        ],
    )
    log.info("plot_publication_cell_line_component_networks.py — start")
    log.info("Output directory: %s", out_dir.resolve())

    git_commit = get_git_commit()
    log.info("Git commit: %s", git_commit)

    # ------------------------------------------------------------------
    # 1. Complete-cohort / relaxed mode check
    # ------------------------------------------------------------------
    cohort_paths = {
        "BRCA": (args.brca_resolved, args.brca_node_stats,
                 args.brca_edges, args.brca_anchors),
        "NBL":  (args.nbl_resolved,  args.nbl_node_stats,
                 args.nbl_edges,  args.nbl_anchors),
        "RBL":  (args.rbl_resolved,  args.rbl_node_stats,
                 args.rbl_edges,  args.rbl_anchors),
    }

    if not args.allow_missing_cohorts:
        missing = [c for c, (res, *_) in cohort_paths.items() if res is None]
        if missing:
            sys.exit(
                f"[ERROR] Resolved neighbour files not provided for: "
                f"{', '.join(missing)}.\n"
                f"  Pass --{missing[0].lower()}-resolved PATH for each "
                f"missing cohort, or use --allow-missing-cohorts to skip "
                f"absent cohorts."
            )
        nonexistent = []
        for cohort, (res, ns, ed, anc) in cohort_paths.items():
            for label, path in [
                ("resolved", res), ("node-stats", ns), ("edges", ed)
            ]:
                if path and not Path(path).exists():
                    nonexistent.append(
                        f"  [{cohort}] --{cohort.lower()}-{label}: {path}")
        if nonexistent:
            sys.exit(
                "[ERROR] The following input files were not found:\n" +
                "\n".join(nonexistent)
            )

    # ------------------------------------------------------------------
    # 2. Load data
    # ------------------------------------------------------------------
    cohort_data = {}
    for cohort, (res_path, ns_path, ed_path, anc_path) in cohort_paths.items():
        if res_path is None:
            cohort_data[cohort] = {
                "cohort": cohort, "available": False,
                "resolved_df": None, "node_stats_df": None,
                "edges_df": None, "anchor_nodes": set(),
                "anchor_path_provided": anc_path is not None,
                "messages": ["No resolved file provided"],
                "resolved_path": None,
                "node_stats_path": ns_path,
                "edges_path": ed_path,
            }
            continue

        resolved_df, err = load_tsv(res_path)
        if resolved_df is None:
            if not args.allow_missing_cohorts:
                sys.exit(f"[ERROR] Cannot load {cohort} resolved file: {err}")
            cohort_data[cohort] = {
                "cohort": cohort, "available": False,
                "resolved_df": None, "node_stats_df": None,
                "edges_df": None, "anchor_nodes": set(),
                "anchor_path_provided": anc_path is not None,
                "messages": [err],
                "resolved_path": res_path,
                "node_stats_path": ns_path,
                "edges_path": ed_path,
            }
            continue

        node_stats_df, _ = load_tsv(ns_path) if ns_path else (None, None)
        edges_df, _      = load_tsv(ed_path, dtype=None) if ed_path else (None, None)

        anchor_nodes = set()
        anchor_path_provided = anc_path is not None
        if anc_path and Path(anc_path).exists():
            anc_df, _ = load_tsv(anc_path)
            if anc_df is not None:
                # Accept column named 'anchor', 'anchor_cell_line', or 'cell_line'
                for col in ("anchor_cell_line", "anchor", "cell_line"):
                    if col in anc_df.columns:
                        anchor_nodes = set(
                            anc_df[col].dropna().astype(str).str.strip()
                        ) - {"nan", ""}
                        break

        cohort_data[cohort] = {
            "cohort": cohort, "available": True,
            "resolved_df": resolved_df,
            "node_stats_df": node_stats_df,
            "edges_df": edges_df,
            "anchor_nodes": anchor_nodes,
            "anchor_path_provided": anchor_path_provided,
            "messages": [],
            "resolved_path": res_path,
            "node_stats_path": ns_path,
            "edges_path": ed_path,
        }

    available = [c for c in COHORTS if cohort_data[c]["available"]]
    if not available:
        sys.exit("[ERROR] No cohort data available. Cannot produce figure.")
    log.info("Available cohorts: %s", ", ".join(available))

    # ------------------------------------------------------------------
    # 3. Cross-cohort ID overlap (mixed cohort file detection)
    # ------------------------------------------------------------------
    cohort_id_sets = {}
    for cohort in available:
        df = cohort_data[cohort]["resolved_df"]
        if df is not None and "cell_line" in df.columns:
            cohort_id_sets[cohort] = (
                set(df["cell_line"].dropna().astype(str).str.strip())
                - {"nan", ""}
            )
        else:
            cohort_id_sets[cohort] = set()

    overlaps = cross_cohort_overlap(cohort_id_sets)
    if overlaps:
        lines = []
        for ca, cb, n, ex in overlaps:
            lines.append(
                f"  {ca} and {cb}: {n} shared ID(s) — "
                f"{', '.join(ex[:5])}{'...' if n > 5 else ''}"
            )
        sys.exit(
            "[ERROR] Mixed cohort files detected. The same cell line IDs "
            "appear in more than one cohort slot, which typically means a "
            "wrong TSV was passed to a cohort argument.\n" +
            "\n".join(lines) +
            "\n  Check that each --{cohort}-resolved path points to the "
            "correct cohort file."
        )

    # ------------------------------------------------------------------
    # 4. Per-cohort validation
    # ------------------------------------------------------------------
    val_rows = []
    for cohort in available:
        d = cohort_data[cohort]
        val_rows.extend(validate_cohort(
            cohort,
            d["resolved_df"], d["node_stats_df"], d["edges_df"],
            (d["resolved_path"], d["node_stats_path"], d["edges_path"]),
        ))

    for cohort in COHORTS:
        if not cohort_data[cohort]["available"]:
            val_rows.append({
                "cohort": cohort, "input_file": "N/A",
                "n_rows": 0, "n_unique_nodes": 0, "id_style": "N/A",
                "validation_status": "SKIP",
                "validation_message": "Cohort not provided",
            })

    write_validation_report(val_rows, out_dir)

    fail_rows = [r for r in val_rows if r["validation_status"] == "FAIL"]
    if fail_rows:
        lines = [
            f"  [{r['cohort']}] {r['input_file']}: {r['validation_message']}"
            for r in fail_rows
        ]
        sys.exit(
            "[ERROR] Input validation failed:\n" + "\n".join(lines) +
            "\n  See network_input_validation.tsv for the full report."
        )

    warn_rows = [r for r in val_rows if r["validation_status"] == "WARN"]
    for r in warn_rows:
        log.warning("[%s] %s: %s", r["cohort"], r["input_file"],
                    r["validation_message"])

    # ------------------------------------------------------------------
    # 5. Build graphs, layouts, communities
    # ------------------------------------------------------------------
    graphs = {}
    layouts = {}
    comp_colour_maps = {}
    community_maps = {}
    all_comp_colours = {}
    cohort_comp_data = {}

    for cohort in COHORTS:
        d = cohort_data[cohort]
        if not d["available"]:
            continue
        G, node_to_comp, isolates, nb_col = build_graph(
            d["resolved_df"], d.get("node_stats_df")
        )
        if G is None:
            log.warning("%s: could not build graph (check column names)", cohort)
            continue

        pos, _ = pack_components(
            G, pad=args.pad, seed=args.layout_seed,
            k_min=1.2, k_max=3.0, t_min=2.6, t_max=6.0,
        )
        pos = _repel(pos, G, min_d=0.18 * 2.6, iters=60, step=0.04,
                     seed=args.layout_seed)

        comp_colour_map = assign_component_colours(G, node_to_comp, isolates)
        comm_map = compute_communities(G, d.get("node_stats_df"))

        graphs[cohort] = {
            "G": G, "node_to_comp": node_to_comp,
            "isolates": isolates, "nb_col": nb_col,
        }
        layouts[cohort] = pos
        comp_colour_maps[cohort] = comp_colour_map
        community_maps[cohort] = comm_map

        comp_n = defaultdict(int)
        for n in G.nodes():
            comp_n[node_to_comp[n]] += 1
        cohort_comp_data[cohort] = dict(comp_n)

        for cid, colour in comp_colour_map.items():
            all_comp_colours[(cohort, cid)] = colour

        log.info("%s: %d nodes, %d edges, %d components, %d isolates",
                 cohort, G.number_of_nodes(), G.number_of_edges(),
                 len(comp_n), len(isolates))

    # ------------------------------------------------------------------
    # 5b. Append community detection status to validation rows
    # ------------------------------------------------------------------
    sample_comm = next(
        (community_maps[c] for c in COHORTS if community_maps.get(c)), {}
    )
    sample_node = next(iter(sample_comm), None)
    if sample_node is not None:
        louv_status = sample_comm[sample_node].get("louvain_community_id", "")
        leid_status = sample_comm[sample_node].get("leiden_community_id", "")
    else:
        louv_status = leid_status = "not_available"

    val_rows.append({
        "cohort": "ALL", "input_file": "community_detection",
        "n_rows": 0, "n_unique_nodes": 0, "id_style": "N/A",
        "validation_status": (
            "PASS" if louv_status != "not_available" else "WARN"),
        "validation_message": (
            f"louvain_community_id={louv_status}; "
            f"leiden_community_id={leid_status}. "
            f"Community columns absent from node_stats and packages "
            f"python-louvain/leidenalg not installed; "
            f"columns set to not_available in summary TSV."
            if louv_status == "not_available"
            else f"Community IDs populated from node_stats columns."
        ),
    })
    write_validation_report(val_rows, out_dir)   # overwrite with updated rows

    # ------------------------------------------------------------------
    # 6. Summary tables
    # ------------------------------------------------------------------
    write_component_community_summary(
        cohort_data, graphs, community_maps, out_dir)
    write_component_annotations(cohort_data, graphs, out_dir)
    write_layout_coordinates(cohort_data, graphs, layouts, args, out_dir)
    write_comparison_table(cohort_data, graphs, out_dir)
    write_edge_table(cohort_data, graphs, out_dir)

    # ------------------------------------------------------------------
    # 7. Combined 4-panel figure
    # ------------------------------------------------------------------
    panel_letters = {"BRCA": "A", "NBL": "B", "RBL": "C"}

    fig = plt.figure(figsize=(args.fig_width, args.fig_height), dpi=args.dpi)
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        hspace=0.35, wspace=0.25,
        left=0.04, right=0.97, top=0.95, bottom=0.04,
    )

    panel_positions = [("BRCA", 0, 0), ("NBL", 0, 1), ("RBL", 1, 0)]
    for cohort, row, col in panel_positions:
        ax = fig.add_subplot(gs[row, col])
        letter = panel_letters[cohort]
        if not cohort_data[cohort]["available"] or graphs.get(cohort) is None:
            ax.set_facecolor("#f5f5f5")
            ax.axis("off")
            ax.text(0.5, 0.5, f"{cohort} data not available",
                    ha="center", va="center", fontsize=12, color="#888888",
                    transform=ax.transAxes)
            ax.set_title(f"Panel {letter} — {cohort}", fontsize=11, pad=10)
            continue
        d = cohort_data[cohort]
        gd = graphs[cohort]
        draw_cohort_panel(
            ax=ax, cohort_label=cohort,
            G=gd["G"], pos=layouts[cohort],
            node_to_comp=gd["node_to_comp"],
            isolates=gd["isolates"],
            comp_colour_map=comp_colour_maps[cohort],
            anchor_nodes=d["anchor_nodes"],
            edges_df=d.get("edges_df"),
            node_stats_df=d.get("node_stats_df"),
            panel_letter=letter, font_size=8,
        )

    ax_d = fig.add_subplot(gs[1, 1])
    draw_legend_panel(ax_d, all_comp_colours, cohort_comp_data)

    stem = out_dir / "Fig_cell_line_similarity_components_combined"
    save_kw = dict(bbox_inches="tight", pad_inches=0.3,
                   facecolor="white", transparent=False)
    fig.savefig(str(stem) + ".pdf", **save_kw)
    fig.savefig(str(stem) + ".svg", **save_kw)
    fig.savefig(str(stem) + ".png", dpi=args.dpi, **save_kw)
    plt.close(fig)
    log.info("Wrote figure: %s.{pdf,svg,png}", stem)

    fig_paths = (
        [str(stem) + ext for ext in (".pdf", ".svg", ".png")] +
        [str(out_dir / f) for f in [
            "cell_line_network_edges.tsv",
            "cell_line_component_community_summary.tsv",
            "cell_line_component_annotations.tsv",
            "cell_line_network_layout_coordinates.tsv",
            "figure_provenance.tsv",
        ]]
    )

    # ------------------------------------------------------------------
    # 8. Provenance
    # ------------------------------------------------------------------
    write_provenance(args, cohort_data, graphs, fig_paths, git_commit, out_dir)

    # ------------------------------------------------------------------
    # Done
    # ------------------------------------------------------------------
    processed = [c for c in COHORTS if graphs.get(c) is not None]
    print(
        f"\nStandalone publication network figure generation complete.\n"
        f"  Cohorts processed : {len(processed)} / {len(COHORTS)}"
        f"  ({', '.join(processed)})\n"
        f"  Output directory  : {out_dir.resolve()}\n"
        f"  Figure            : {stem}.{{pdf,svg,png}}\n"
        f"  Log               : {log_path}\n"
        f"  Git commit        : {git_commit}\n"
    )


if __name__ == "__main__":
    main()
