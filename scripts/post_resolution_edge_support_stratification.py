#!/usr/bin/env python3
"""Post-resolution edge-support stratification for cell-line graph audits."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


REQUIRED_RBL_ROWS = {"RBL_14", "RBL_7", "WERI_RB1", "Y_79", "RBL_15"}

OUTPUT_COLUMNS = [
    "Cohort",
    "Cell line",
    "Resolved-neighbour status",
    "Resolved neighbours",
    "Majority/support-threshold neighbours",
    "Union-only neighbours",
    "Highest lost-edge support count",
    "Highest lost-edge support fraction, if available",
    "Edge-support class",
    "Interpretation",
]

INTERPRETATIONS = {
    "Resolved edge": (
        "Stable resolved cell-line neighbour under the primary anti-overfitting rule."
    ),
    "Post-resolution recurrent edge": (
        "Recurrent edge with support across representations, but not accepted by the "
        "resolved-neighbour rule."
    ),
    "Union-only edge": (
        "Weak or representation-specific edge; useful for auditing connectivity, "
        "but not sufficient for stable resolved-neighbour interpretation."
    ),
    "No detected edge": (
        "No observed cell-line edge under the current feature-distance representation set."
    ),
}


def read_tsv(path: Path) -> List[dict]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def write_tsv(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=OUTPUT_COLUMNS,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def first_present(row: dict, candidates: Sequence[str]) -> str | None:
    for key in candidates:
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def canonical_pair(a: str, b: str) -> Tuple[str, str]:
    a = str(a).strip()
    b = str(b).strip()
    return tuple(sorted((a, b)))


def safe_float(value: object) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NONE"}:
        return None
    try:
        out = float(text)
    except ValueError:
        return None
    if math.isnan(out):
        return None
    return out


def format_number(value: float | None) -> str:
    if value is None:
        return "NA"
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.6g}"


def edge_support(row: dict) -> Tuple[float | None, float | None]:
    support = safe_float(first_present(row, ["support", "support_directions"]))
    n_configured = safe_float(first_present(row, ["n_configured_directions"]))
    fraction = None
    if support is not None and n_configured is not None and n_configured > 0:
        fraction = support / n_configured
    elif support is not None:
        fraction = safe_float(
            first_present(row, ["support_fraction", "support_frac", "support_weight_mean"])
        )
    return support, fraction


def read_edges(path: Path) -> Dict[Tuple[str, str], dict]:
    rows = read_tsv(path)
    out: Dict[Tuple[str, str], dict] = {}
    for row in rows:
        node1 = first_present(row, ["node1", "node_a", "source", "from"])
        node2 = first_present(row, ["node2", "node_b", "target", "to"])
        if not node1 or not node2 or node1 == node2:
            continue
        pair = canonical_pair(node1, node2)
        support, fraction = edge_support(row)
        current = out.get(pair)
        if current is None or (
            support is not None and (current["support"] is None or support > current["support"])
        ):
            out[pair] = {
                "support": support,
                "fraction": fraction,
            }
    return out


def read_nodes(path: Path, edge_maps: Iterable[Dict[Tuple[str, str], dict]]) -> List[str]:
    rows = read_tsv(path)
    nodes = set()
    for row in rows:
        node = first_present(row, ["node", "cell_line", "short_id", "cell_line_short"])
        if node:
            nodes.add(node)
    for edge_map in edge_maps:
        for a, b in edge_map:
            nodes.add(a)
            nodes.add(b)
    return sorted(nodes)


def neighbours_for(pairs: Iterable[Tuple[str, str]]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = defaultdict(list)
    for a, b in pairs:
        out[a].append(b)
        out[b].append(a)
    return {node: sorted(set(neighbours)) for node, neighbours in out.items()}


def joined(values: Sequence[str]) -> str:
    return ";".join(values) if values else "None"


def classify_cohort(args: argparse.Namespace) -> None:
    cohort = args.cohort.upper()
    union_edges = read_edges(Path(args.union_edges))
    consensus_edges = read_edges(Path(args.consensus_edges))
    resolved_edges = read_edges(Path(args.resolved_edges))
    nodes = read_nodes(Path(args.nodes), [union_edges, consensus_edges, resolved_edges])

    if cohort == "RBL":
        missing = sorted(REQUIRED_RBL_ROWS - set(nodes))
        if missing:
            raise SystemExit(
                "RBL post-resolution stratification is missing required row(s): "
                + ", ".join(missing)
            )

    union_pairs = set(union_edges)
    consensus_pairs = set(consensus_edges)
    resolved_pairs = set(resolved_edges)
    consensus_lost_pairs = consensus_pairs - resolved_pairs
    union_only_lost_pairs = union_pairs - consensus_pairs - resolved_pairs
    all_lost_pairs = (union_pairs | consensus_pairs) - resolved_pairs

    resolved_neighbours = neighbours_for(resolved_pairs)
    consensus_lost_neighbours = neighbours_for(consensus_lost_pairs)
    union_only_lost_neighbours = neighbours_for(union_only_lost_pairs)

    rows = []
    for node in nodes:
        resolved = resolved_neighbours.get(node, [])
        consensus_lost = consensus_lost_neighbours.get(node, [])
        union_only_lost = union_only_lost_neighbours.get(node, [])

        lost_supports = []
        for pair in all_lost_pairs:
            if node not in pair:
                continue
            info = union_edges.get(pair) or consensus_edges.get(pair)
            if info and info["support"] is not None:
                lost_supports.append(info)

        if lost_supports:
            best_lost = max(
                lost_supports,
                key=lambda item: (
                    item["support"] if item["support"] is not None else -1,
                    item["fraction"] if item["fraction"] is not None else -1,
                ),
            )
            highest_lost_support = format_number(best_lost["support"])
            highest_lost_fraction = format_number(best_lost["fraction"])
        else:
            highest_lost_support = "NA"
            highest_lost_fraction = "NA"

        if resolved:
            edge_class = "Resolved edge"
            status = "Resolved neighbour present"
        elif consensus_lost:
            edge_class = "Post-resolution recurrent edge"
            status = "Resolved isolate"
        elif union_only_lost:
            edge_class = "Union-only edge"
            status = "Resolved isolate"
        else:
            edge_class = "No detected edge"
            status = "Resolved isolate"

        rows.append(
            {
                "Cohort": cohort,
                "Cell line": node,
                "Resolved-neighbour status": status,
                "Resolved neighbours": joined(resolved),
                "Majority/support-threshold neighbours": joined(consensus_lost),
                "Union-only neighbours": joined(union_only_lost),
                "Highest lost-edge support count": highest_lost_support,
                "Highest lost-edge support fraction, if available": highest_lost_fraction,
                "Edge-support class": edge_class,
                "Interpretation": INTERPRETATIONS[edge_class],
            }
        )

    write_tsv(Path(args.out), rows)
    print(f"[OK] wrote {args.out} with {len(rows)} rows for {cohort}")


def combine_tables(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = False
    with out.open("w", newline="") as out_handle:
        writer = None
        for table in args.tables:
            path = Path(table)
            with path.open(newline="") as in_handle:
                reader = csv.DictReader(in_handle, delimiter="\t")
                if reader.fieldnames != OUTPUT_COLUMNS:
                    raise SystemExit(
                        f"{path} has an unexpected header: {reader.fieldnames}"
                    )
                if not wrote_header:
                    writer = csv.DictWriter(
                        out_handle,
                        delimiter="\t",
                        fieldnames=OUTPUT_COLUMNS,
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    wrote_header = True
                assert writer is not None
                for row in reader:
                    writer.writerow(row)
    print(f"[OK] wrote {args.out} from {len(args.tables)} tables")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Classify cell lines by post-resolution edge support."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    classify = subparsers.add_parser("classify")
    classify.add_argument("--cohort", required=True)
    classify.add_argument("--union-edges", required=True)
    classify.add_argument("--consensus-edges", required=True)
    classify.add_argument("--resolved-edges", required=True)
    classify.add_argument("--nodes", required=True)
    classify.add_argument("--out", required=True)
    classify.set_defaults(func=classify_cohort)

    combine = subparsers.add_parser("combine")
    combine.add_argument("--out", required=True)
    combine.add_argument("tables", nargs="+")
    combine.set_defaults(func=combine_tables)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
