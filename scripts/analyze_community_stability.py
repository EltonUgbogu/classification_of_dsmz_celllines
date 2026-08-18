#!/usr/bin/env python3
"""Cross-direction cell-line graph community stability assessment.

This script intentionally treats Leiden/Louvain assignments as auxiliary
assessments on per-direction cell-line similarity graphs. It does not infer
tumour communities.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ALG_COLS = {"Leiden": "community_leid", "Louvain": "community_louv"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute cross-direction Leiden/Louvain co-assignment stability for cell-line graphs."
    )
    p.add_argument("--cohort", required=True)
    p.add_argument("--tumour-nh-root", required=True)
    p.add_argument("--directions", required=True, help="Comma-separated direction list")
    p.add_argument("--resolved-tsv", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--shortnames-tsv", default="")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def split_field(value: str) -> list[str]:
    if value is None:
        return []
    s = str(value).strip()
    if not s or s.upper() in {"NA", "NAN", "."}:
        return []
    return [x.strip() for x in s.split(";") if x.strip()]


def load_short_map(path: Path | None) -> dict[str, str]:
    if not path or not path.exists():
        return {}
    df = read_tsv(path)
    if {"long_id", "short_id"}.issubset(df.columns):
        return dict(zip(df["long_id"], df["short_id"]))
    return {}


def canonical_id(x: str, short_map: dict[str, str]) -> str:
    x = str(x).strip()
    return short_map.get(x, x)


def load_resolved_components(resolved_tsv: Path, short_map: dict[str, str]) -> pd.DataFrame:
    if not resolved_tsv.exists():
        raise FileNotFoundError(f"final resolved graph TSV not found: {resolved_tsv}")

    df = read_tsv(resolved_tsv)
    if "cell_line" not in df.columns:
        raise ValueError(f"{resolved_tsv} must contain column 'cell_line'")
    nb_col = "final_neighbors" if "final_neighbors" in df.columns else "final_neighbours"
    if nb_col not in df.columns:
        raise ValueError(f"{resolved_tsv} must contain final_neighbors or final_neighbours")

    nodes = [canonical_id(x, short_map) for x in df["cell_line"] if str(x).strip()]
    nodes = sorted(set(nodes))
    adj = {n: set() for n in nodes}

    for _, row in df.iterrows():
        src = canonical_id(row["cell_line"], short_map)
        if src not in adj:
            adj[src] = set()
        for raw_tgt in split_field(row.get(nb_col, "")):
            tgt = canonical_id(raw_tgt, short_map)
            if tgt not in adj:
                adj[tgt] = set()
            if tgt != src:
                adj[src].add(tgt)
                adj[tgt].add(src)

    seen = set()
    comps = []
    for node in sorted(adj):
        if node in seen:
            continue
        stack = [node]
        seen.add(node)
        members = []
        while stack:
            cur = stack.pop()
            members.append(cur)
            for nxt in sorted(adj[cur]):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        comps.append(sorted(members))

    non_iso = [c for c in comps if len(c) > 1]
    non_iso = sorted(non_iso, key=lambda c: (-len(c), c[0]))
    comp_id = {}
    for i, comp in enumerate(non_iso, start=1):
        for node in comp:
            comp_id[node] = f"component_{i}"
    for comp in comps:
        if len(comp) == 1:
            comp_id[comp[0]] = f"isolate__{comp[0]}"

    rows = []
    for node in sorted(adj):
        rows.append(
            {
                "cell_line_id": node,
                "component_id": comp_id[node],
                "is_isolate": len(adj[node]) == 0,
                "degree_final_resolved": len(adj[node]),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError(f"no final resolved graph nodes found in {resolved_tsv}")
    return out


def expand_assignment_row(row: pd.Series, final_nodes: set[str], short_map: dict[str, str]) -> list[str]:
    candidates = []
    for col in ("cell_line", "sample_id", "cell_tech_id", "cell_line_display"):
        if col in row.index:
            candidates.extend(split_field(row.get(col, "")))

    expanded = []
    for x in candidates:
        cx = canonical_id(x, short_map)
        if cx in final_nodes:
            expanded.append(cx)

    if not expanded and "cell_line" in row.index:
        cx = canonical_id(row["cell_line"], short_map)
        if cx in final_nodes:
            expanded.append(cx)

    return sorted(set(expanded))


def normalize_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def load_assignments(
    cohort: str,
    tumour_root: Path,
    directions: list[str],
    final_nodes: set[str],
    short_map: dict[str, str],
) -> tuple[pd.DataFrame, list[dict[str, str]], list[dict[str, object]]]:
    rows = []
    missing = []
    per_dir = []

    for direction in directions:
        path = tumour_root / direction / "final_consensus" / f"cell_line_similarity_graph_node_annotations_{direction}.tsv"
        if not path.exists():
            missing.append({"direction": direction, "path": str(path)})
            continue

        df = read_tsv(path)
        needed = {"cell_line", "component", "community_louv", "community_leid", "is_outlier"}
        absent = sorted(needed - set(df.columns))
        if absent:
            raise ValueError(f"{path} missing required columns: {', '.join(absent)}")

        n_leiden = 0
        n_louvain = 0
        n_nodes = 0
        for _, row in df.iterrows():
            cell_ids = expand_assignment_row(row, final_nodes, short_map)
            if not cell_ids:
                continue
            comp = str(row.get("component", "")).strip()
            is_isolate = normalize_bool(row.get("is_outlier", "FALSE"))

            for cell_id in cell_ids:
                n_nodes += 1
                for alg, col in ALG_COLS.items():
                    raw_comm = str(row.get(col, "")).strip()
                    if not raw_comm or raw_comm.upper() in {"NA", "NAN", "."}:
                        community_id = ""
                    elif is_isolate:
                        community_id = ""
                    else:
                        prefix = f"component_{comp}" if comp else "component_unknown"
                        community_id = f"{prefix}__{alg.lower()}_{raw_comm}"
                    rows.append(
                        {
                            "cohort": cohort,
                            "direction": direction,
                            "algorithm": alg,
                            "cell_line_id": cell_id,
                            "connected_component": f"component_{comp}" if comp else "",
                            "community_id": community_id,
                            "is_isolate": bool(is_isolate),
                        }
                    )
                    if alg == "Leiden" and community_id:
                        n_leiden += 1
                    if alg == "Louvain" and community_id:
                        n_louvain += 1

        per_dir.append(
            {
                "direction": direction,
                "path": str(path),
                "n_cell_lines_loaded": n_nodes,
                "n_leiden_assignments_loaded": n_leiden,
                "n_louvain_assignments_loaded": n_louvain,
            }
        )

    if not rows:
        raise ValueError("no node-level community assignment tables were loaded")

    return pd.DataFrame(rows), missing, per_dir


def matrix_for_algorithm(assignments: pd.DataFrame, algorithm: str, ordered_nodes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    sub = assignments[assignments["algorithm"] == algorithm].copy()
    valid = sub[sub["community_id"].astype(str).str.len() > 0].copy()
    directions_used = sorted(valid["direction"].unique())
    if len(directions_used) < 2:
        raise ValueError(f"fewer than two directions contain valid {algorithm} assignments")

    n = len(ordered_nodes)
    idx = {node: i for i, node in enumerate(ordered_nodes)}
    same = np.zeros((n, n), dtype=float)
    denom = np.zeros((n, n), dtype=float)

    for direction in directions_used:
        dsub = valid[valid["direction"] == direction]
        comm_by_node = dict(zip(dsub["cell_line_id"], dsub["community_id"]))
        present = [node for node in ordered_nodes if node in comm_by_node]
        for a_i, a in enumerate(present):
            ia = idx[a]
            for b in present[a_i:]:
                ib = idx[b]
                denom[ia, ib] += 1
                denom[ib, ia] += 1
                if comm_by_node[a] == comm_by_node[b]:
                    same[ia, ib] += 1
                    same[ib, ia] += 1

    with np.errstate(divide="ignore", invalid="ignore"):
        freq = np.divide(same, denom, out=np.full_like(same, np.nan), where=denom > 0)

    return (
        pd.DataFrame(freq, index=ordered_nodes, columns=ordered_nodes),
        pd.DataFrame(denom.astype(int), index=ordered_nodes, columns=ordered_nodes),
        len(directions_used),
    )


def pair_records(freq: pd.DataFrame, nvalid: pd.DataFrame, components: pd.DataFrame) -> pd.DataFrame:
    comp = components.set_index("cell_line_id")
    nodes = list(freq.index)
    rows = []
    for i, a in enumerate(nodes):
        for b in nodes[i + 1 :]:
            val = freq.loc[a, b]
            nv = int(nvalid.loc[a, b])
            if pd.isna(val) or nv == 0:
                continue
            a_iso = bool(comp.loc[a, "is_isolate"])
            b_iso = bool(comp.loc[b, "is_isolate"])
            same_component = comp.loc[a, "component_id"] == comp.loc[b, "component_id"]
            rows.append(
                {
                    "cell_line_1": a,
                    "cell_line_2": b,
                    "coassignment": float(val),
                    "n_valid": nv,
                    "component_1": comp.loc[a, "component_id"],
                    "component_2": comp.loc[b, "component_id"],
                    "same_component": bool(same_component and not a_iso and not b_iso),
                    "involves_isolate": bool(a_iso or b_iso),
                }
            )
    return pd.DataFrame(rows)


def safe_mean(s: pd.Series) -> float:
    return float(s.mean()) if len(s) else math.nan


def safe_median(s: pd.Series) -> float:
    return float(s.median()) if len(s) else math.nan


def component_summary(cohort: str, algorithm: str, pairs: pd.DataFrame, components: pd.DataFrame) -> pd.DataFrame:
    rows = []
    comp_groups = components.groupby("component_id", sort=False)
    non_iso_nodes = set(components.loc[~components["is_isolate"], "cell_line_id"])

    for cid, cdf in comp_groups:
        members = set(cdf["cell_line_id"])
        is_isolate_comp = bool(cdf["is_isolate"].all())
        within = pairs[(pairs["component_1"] == cid) & (pairs["component_2"] == cid) & (~pairs["involves_isolate"])]
        if is_isolate_comp:
            between = pairs[
                ((pairs["cell_line_1"].isin(members)) | (pairs["cell_line_2"].isin(members)))
                & (~pairs["cell_line_1"].eq(pairs["cell_line_2"]))
            ]
        else:
            between = pairs[
                (
                    (pairs["cell_line_1"].isin(members) & pairs["cell_line_2"].isin(non_iso_nodes - members))
                    | (pairs["cell_line_2"].isin(members) & pairs["cell_line_1"].isin(non_iso_nodes - members))
                )
                & (~pairs["involves_isolate"])
            ]
        rows.append(
            {
                "cohort": cohort,
                "algorithm": algorithm,
                "component_id": cid,
                "n_cell_lines": len(members),
                "n_within_pairs": len(within),
                "mean_within_coassignment": safe_mean(within["coassignment"]),
                "median_within_coassignment": safe_median(within["coassignment"]),
                "mean_between_coassignment": safe_mean(between["coassignment"]),
                "median_between_coassignment": safe_median(between["coassignment"]),
                "within_minus_between": safe_mean(within["coassignment"]) - safe_mean(between["coassignment"]),
                "mean_valid_directions_within": safe_mean(within["n_valid"]),
                "mean_valid_directions_between": safe_mean(between["n_valid"]),
            }
        )
    return pd.DataFrame(rows)


def global_summary(
    cohort: str,
    algorithm: str,
    n_dirs: int,
    pairs: pd.DataFrame,
    components: pd.DataFrame,
) -> dict[str, object]:
    within = pairs[pairs["same_component"] & (~pairs["involves_isolate"])]
    between = pairs[(~pairs["same_component"]) & (~pairs["involves_isolate"])]
    n_components = int((~components["is_isolate"]).groupby(components["component_id"]).any().sum())
    return {
        "cohort": cohort,
        "algorithm": algorithm,
        "n_directions_used": n_dirs,
        "n_cell_lines": len(components),
        "n_components": n_components,
        "n_isolates": int(components["is_isolate"].sum()),
        "mean_within_coassignment": safe_mean(within["coassignment"]),
        "median_within_coassignment": safe_median(within["coassignment"]),
        "mean_between_coassignment": safe_mean(between["coassignment"]),
        "median_between_coassignment": safe_median(between["coassignment"]),
        "within_minus_between": safe_mean(within["coassignment"]) - safe_mean(between["coassignment"]),
    }


def ordered_component_table(components: pd.DataFrame) -> pd.DataFrame:
    out = components.copy()
    out["_size"] = out.groupby("component_id")["cell_line_id"].transform("size")
    out["_iso"] = out["is_isolate"].astype(int)
    out = out.sort_values(["_iso", "component_id", "cell_line_id"]).drop(columns=["_size", "_iso"])
    return out


def plot_heatmap(freq: pd.DataFrame, components: pd.DataFrame, algorithm: str, cohort: str, outpath: Path) -> None:
    ordered = ordered_component_table(components)
    nodes = ordered["cell_line_id"].tolist()
    mat = freq.loc[nodes, nodes].to_numpy(dtype=float)

    fig_w = max(7, min(14, 0.35 * len(nodes) + 3))
    fig_h = max(6, min(13, 0.35 * len(nodes) + 2))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(mat, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(nodes)))
    ax.set_yticks(range(len(nodes)))
    ax.set_xticklabels(nodes, rotation=90, fontsize=7)
    ax.set_yticklabels(nodes, fontsize=7)
    ax.set_title(
        f"{cohort.upper()} {algorithm} cross-direction cell-line community co-assignment frequency",
        fontsize=11,
    )

    last = None
    for i, cid in enumerate(ordered["component_id"].tolist()):
        if last is not None and cid != last:
            ax.axhline(i - 0.5, color="white", linewidth=1.5)
            ax.axvline(i - 0.5, color="white", linewidth=1.5)
        last = cid

    iso_positions = [i for i, v in enumerate(ordered["is_isolate"].tolist()) if bool(v)]
    for i in iso_positions:
        ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="red", linewidth=1.2))

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Cross-direction co-assignment frequency")
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_within_between(pairs_by_alg: dict[str, pd.DataFrame], cohort: str, outpath_by_alg: dict[str, Path]) -> None:
    for alg, pairs in pairs_by_alg.items():
        within = pairs[pairs["same_component"] & (~pairs["involves_isolate"])]["coassignment"].dropna()
        between = pairs[(~pairs["same_component"]) & (~pairs["involves_isolate"])]["coassignment"].dropna()
        fig, ax = plt.subplots(figsize=(5.5, 5))
        data = [within.to_numpy(), between.to_numpy()]
        ax.boxplot(data, labels=["Within final\ncomponents", "Between final\ncomponents"], showmeans=True)
        for x, vals in enumerate(data, start=1):
            if len(vals):
                jitter = np.linspace(-0.08, 0.08, min(len(vals), 100))
                shown = vals[: len(jitter)]
                ax.scatter(np.full(len(shown), x) + jitter, shown, s=12, alpha=0.45)
        ax.set_ylim(-0.05, 1.05)
        ax.set_ylabel("Cross-direction co-assignment frequency")
        ax.set_title(f"{cohort.upper()} {alg}: cell-line community stability vs final components", fontsize=11)
        fig.tight_layout()
        fig.savefig(outpath_by_alg[alg])
        plt.close(fig)


def write_validation_report(
    outpath: Path,
    cohort: str,
    directions: list[str],
    assignments: pd.DataFrame,
    missing: list[dict[str, str]],
    per_dir: list[dict[str, object]],
    components: pd.DataFrame,
    outputs: list[Path],
) -> None:
    lines = []
    lines.append("Community stability validation report")
    lines.append(f"cohort\t{cohort}")
    lines.append(f"scope\tcell-line graph stability assessment")
    lines.append(f"directions_discovered\t{len(directions)}")
    lines.append(f"directions_loaded\t{len(per_dir)}")
    lines.append(f"missing_assignment_files\t{len(missing)}")
    lines.append(f"final_resolved_components\t{components.loc[~components['is_isolate'], 'component_id'].nunique()}")
    lines.append(f"isolates\t{int(components['is_isolate'].sum())}")
    lines.append(f"leiden_assignments_loaded\t{int((assignments['algorithm'].eq('Leiden') & assignments['community_id'].astype(str).str.len().gt(0)).sum())}")
    lines.append(f"louvain_assignments_loaded\t{int((assignments['algorithm'].eq('Louvain') & assignments['community_id'].astype(str).str.len().gt(0)).sum())}")
    lines.append("")
    lines.append("Per-direction assignment loading")
    for row in per_dir:
        lines.append(
            f"{row['direction']}\tcell_lines={row['n_cell_lines_loaded']}\t"
            f"Leiden={row['n_leiden_assignments_loaded']}\tLouvain={row['n_louvain_assignments_loaded']}\t{row['path']}"
        )
    if missing:
        lines.append("")
        lines.append("Missing assignment files")
        for row in missing:
            lines.append(f"{row['direction']}\t{row['path']}")
    lines.append("")
    lines.append("Generated outputs")
    for path in outputs:
        lines.append(str(path))
    outpath.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    cohort = args.cohort
    tumour_root = Path(args.tumour_nh_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    directions = [d.strip() for d in args.directions.split(",") if d.strip()]
    if len(directions) < 2:
        raise ValueError("at least two directions are required")

    short_map = load_short_map(Path(args.shortnames_tsv) if args.shortnames_tsv else None)
    components = load_resolved_components(Path(args.resolved_tsv), short_map)
    components_ordered = ordered_component_table(components)
    final_nodes = set(components["cell_line_id"])

    assignments, missing, per_dir = load_assignments(cohort, tumour_root, directions, final_nodes, short_map)
    if len(per_dir) < 2:
        raise ValueError(f"fewer than two directions contain valid assignment tables for {cohort}")
    if missing:
        sys.stderr.write(f"[WARN] {cohort}: {len(missing)} direction assignment files are missing\n")

    outputs = []
    assign_path = outdir / "community_assignments_long.tsv"
    assignments = assignments.sort_values(["cohort", "direction", "algorithm", "cell_line_id"])
    assignments.to_csv(assign_path, sep="\t", index=False)
    outputs.append(assign_path)

    comp_path = outdir / "final_resolved_component_assignments.tsv"
    components_ordered.to_csv(comp_path, sep="\t", index=False)
    outputs.append(comp_path)

    ordered_nodes = components_ordered["cell_line_id"].tolist()
    comp_summaries = []
    global_summaries = []
    pairs_by_alg = {}

    for alg in ("Leiden", "Louvain"):
        freq, nvalid, n_dirs = matrix_for_algorithm(assignments, alg, ordered_nodes)
        freq_path = outdir / f"community_coassignment_matrix_{alg}.tsv"
        nvalid_path = outdir / f"community_coassignment_n_valid_{alg}.tsv"
        freq.to_csv(freq_path, sep="\t", na_rep="NA")
        nvalid.to_csv(nvalid_path, sep="\t")
        outputs.extend([freq_path, nvalid_path])

        pairs = pair_records(freq, nvalid, components_ordered)
        pairs_by_alg[alg] = pairs
        pairs_path = outdir / f"community_coassignment_pairs_{alg}.tsv"
        pairs.to_csv(pairs_path, sep="\t", index=False)
        outputs.append(pairs_path)

        comp_summaries.append(component_summary(cohort, alg, pairs, components_ordered))
        global_summaries.append(global_summary(cohort, alg, n_dirs, pairs, components_ordered))

        heatmap_path = outdir / f"Fig_community_coassignment_heatmap_{alg}.pdf"
        plot_heatmap(freq, components_ordered, alg, cohort, heatmap_path)
        outputs.append(heatmap_path)

    wb_paths = {
        "Leiden": outdir / "Fig_community_stability_within_between_Leiden.pdf",
        "Louvain": outdir / "Fig_community_stability_within_between_Louvain.pdf",
    }
    plot_within_between(pairs_by_alg, cohort, wb_paths)
    outputs.extend(wb_paths.values())

    comp_summary_path = outdir / "community_stability_component_summary.tsv"
    pd.concat(comp_summaries, ignore_index=True).to_csv(comp_summary_path, sep="\t", index=False, na_rep="NA")
    outputs.append(comp_summary_path)

    global_summary_path = outdir / "community_stability_global_summary.tsv"
    pd.DataFrame(global_summaries).to_csv(global_summary_path, sep="\t", index=False, na_rep="NA")
    outputs.append(global_summary_path)

    validation_path = outdir / "community_stability_validation_report.txt"
    write_validation_report(validation_path, cohort, directions, assignments, missing, per_dir, components_ordered, outputs)
    outputs.append(validation_path)

    print("Community stability analysis complete")
    print(f"cohort: {cohort}")
    print(f"directions loaded: {len(per_dir)} / {len(directions)}")
    print(f"final resolved components: {components_ordered.loc[~components_ordered['is_isolate'], 'component_id'].nunique()}")
    print(f"isolates: {int(components_ordered['is_isolate'].sum())}")
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        sys.stderr.write(f"[ERROR] {exc}\n")
        raise
