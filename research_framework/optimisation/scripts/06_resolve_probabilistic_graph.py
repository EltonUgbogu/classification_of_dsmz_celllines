import argparse
import csv
from collections import defaultdict, deque
from pathlib import Path

import yaml


def read_tsv(path):
    with path.open("r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def append_report(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write("\n".join(lines) + "\n")


def disjoint_set_find(parent, node):
    while parent[node] != node:
        parent[node] = parent[parent[node]]
        node = parent[node]
    return node


def disjoint_set_union(parent, first_node, second_node):
    first_root = disjoint_set_find(parent, first_node)
    second_root = disjoint_set_find(parent, second_node)
    if first_root != second_root:
        parent[second_root] = first_root


def connected_components(nodes, edges):
    parent = {node: node for node in nodes}
    for first_node, second_node, _ in edges:
        disjoint_set_union(parent, first_node, second_node)
    grouped_nodes = defaultdict(list)
    for node in nodes:
        grouped_nodes[disjoint_set_find(parent, node)].append(node)
    return list(grouped_nodes.values())


def graph_betweenness(nodes, adjacency):
    betweenness = {node: 0.0 for node in nodes}
    for source in nodes:
        stack = []
        predecessors = {node: [] for node in nodes}
        sigma = dict.fromkeys(nodes, 0.0)
        distance = dict.fromkeys(nodes, -1)
        sigma[source] = 1.0
        distance[source] = 0
        queue = deque([source])
        while queue:
            vertex = queue.popleft()
            stack.append(vertex)
            for neighbour in adjacency[vertex]:
                if distance[neighbour] < 0:
                    queue.append(neighbour)
                    distance[neighbour] = distance[vertex] + 1
                if distance[neighbour] == distance[vertex] + 1:
                    sigma[neighbour] += sigma[vertex]
                    predecessors[neighbour].append(vertex)
        dependency = dict.fromkeys(nodes, 0.0)
        while stack:
            vertex = stack.pop()
            for predecessor in predecessors[vertex]:
                if sigma[vertex] > 0:
                    dependency[predecessor] += (sigma[predecessor] / sigma[vertex]) * (1 + dependency[vertex])
            if vertex != source:
                betweenness[vertex] += dependency[vertex]
    return {node: value / 2.0 for node, value in betweenness.items()}


def load_deterministic_isolates(results_root):
    isolate_status = {}
    for path in results_root.glob("**/isolate_cell_lines.tsv"):
        with path.open("r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                cell_line = row.get("cell_line") or row.get("cell_line_id")
                if cell_line:
                    isolate_status[cell_line] = "deterministic_isolate"
    return isolate_status


def main():
    parser = argparse.ArgumentParser(description="Resolve probabilistic cell-line graph edges.")
    parser.add_argument("--config", default="research_framework/optimisation/config/optimisation.yaml")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    config = yaml.safe_load(config_path.read_text())

    output_root = Path(config["outputs"]["optimisation_results_root"])
    docs_root = Path(config["outputs"]["optimisation_docs_root"])
    thesis_results_root = Path(config["inputs"]["thesis_results_root"])
    graph_config = config["probabilistic_graph"]
    threshold_grid = [float(value) for value in graph_config["edge_probability_threshold_grid"]]
    maximum_expected_false_edge_rate = float(graph_config["max_expected_false_edge_rate"])
    allow_no_valid_threshold = bool(graph_config.get("allow_no_valid_threshold", True))

    edge_path = output_root / "probabilistic_graphs" / "probabilistic_cellline_edges.tsv"
    if not edge_path.exists():
        raise FileNotFoundError(f"Probabilistic edge table not found: {edge_path}")

    edge_rows = read_tsv(edge_path)
    if not edge_rows:
        raise ValueError("Probabilistic edge table is empty.")

    for row in edge_rows:
        row["posterior_edge_probability_float"] = float(row["posterior_edge_probability"])

    threshold_rows = []
    valid_thresholds = []
    for threshold in threshold_grid:
        retained = [row for row in edge_rows if row["posterior_edge_probability_float"] >= threshold]
        if retained:
            expected_false_edge_rate = sum(1 - row["posterior_edge_probability_float"] for row in retained) / len(retained)
        else:
            expected_false_edge_rate = 0.0
        threshold_rows.append({
            "candidate_threshold": threshold,
            "n_retained_edges": len(retained),
            "expected_false_edge_rate": expected_false_edge_rate,
            "passes_expected_false_edge_rate": expected_false_edge_rate <= maximum_expected_false_edge_rate and len(retained) > 0,
        })
        if expected_false_edge_rate <= maximum_expected_false_edge_rate and retained:
            valid_thresholds.append(threshold)

    if valid_thresholds:
        selected_threshold = min(valid_thresholds)
        threshold_note = "configured expected false edge-rate criterion satisfied"
    elif allow_no_valid_threshold:
        selected_threshold = max(threshold_grid)
        threshold_note = "no threshold satisfied the criterion; most conservative configured threshold used"
    else:
        raise ValueError("No candidate threshold satisfied the expected false edge-rate criterion.")

    retained_rows = [row for row in edge_rows if row["posterior_edge_probability_float"] >= selected_threshold]
    expected_false_edge_rate = (
        sum(1 - row["posterior_edge_probability_float"] for row in retained_rows) / len(retained_rows)
        if retained_rows else 0.0
    )

    all_nodes = sorted(set([row["cell_line_1"] for row in edge_rows] + [row["cell_line_2"] for row in edge_rows]))
    retained_edges = [(row["cell_line_1"], row["cell_line_2"], row["posterior_edge_probability_float"]) for row in retained_rows]
    components = connected_components(all_nodes, retained_edges)

    adjacency = {node: set() for node in all_nodes}
    incident_probabilities = defaultdict(list)
    for first_node, second_node, probability in retained_edges:
        adjacency[first_node].add(second_node)
        adjacency[second_node].add(first_node)
    for row in edge_rows:
        incident_probabilities[row["cell_line_1"]].append(row["posterior_edge_probability_float"])
        incident_probabilities[row["cell_line_2"]].append(row["posterior_edge_probability_float"])

    betweenness = graph_betweenness(all_nodes, adjacency)
    deterministic_isolates = load_deterministic_isolates(thesis_results_root)

    resolved_edge_rows = []
    for row in edge_rows:
        probability = row["posterior_edge_probability_float"]
        retained = probability >= selected_threshold
        resolved_edge_rows.append({
            "cell_line_1": row["cell_line_1"],
            "cell_line_2": row["cell_line_2"],
            "posterior_edge_probability": row["posterior_edge_probability"],
            "selected_threshold": selected_threshold,
            "expected_false_edge_rate": expected_false_edge_rate,
            "retained_in_probabilistic_graph": retained,
        })

    component_rows = []
    component_id_by_node = {}
    for component_index, component_nodes in enumerate(components, start=1):
        component_edges = [edge for edge in retained_edges if edge[0] in component_nodes and edge[1] in component_nodes]
        edge_probabilities = [edge[2] for edge in component_edges]
        possible_edges = len(component_nodes) * (len(component_nodes) - 1) / 2
        component_density = len(component_edges) / possible_edges if possible_edges else 0.0
        composition = "unknown"
        for node in component_nodes:
            component_id_by_node[node] = component_index
            component_rows.append({
                "cell_line_id": node,
                "probabilistic_component_id": component_index,
                "component_size": len(component_nodes),
                "component_mean_edge_probability": sum(edge_probabilities) / len(edge_probabilities) if edge_probabilities else 0.0,
                "component_min_edge_probability": min(edge_probabilities) if edge_probabilities else 0.0,
                "component_density": component_density,
                "component_cancer_type_composition": composition,
            })

    isolate_rows = []
    for node in all_nodes:
        probabilities = incident_probabilities.get(node, [])
        posterior_isolate_probability = 1.0
        for probability in probabilities:
            posterior_isolate_probability *= 1 - probability
        probabilistic_isolate_status = len(adjacency[node]) == 0
        isolate_rows.append({
            "cell_line_id": node,
            "posterior_isolate_probability": posterior_isolate_probability,
            "deterministic_isolate_status": deterministic_isolates.get(node, "not_deterministic_isolate"),
            "probabilistic_isolate_status": probabilistic_isolate_status,
            "notes": threshold_note,
        })

    degree_by_node = {node: len(adjacency[node]) for node in all_nodes}
    nonzero_betweenness = [value for value in betweenness.values() if value > 0]
    bridge_cutoff = sorted(nonzero_betweenness)[int(0.75 * (len(nonzero_betweenness) - 1))] if nonzero_betweenness else float("inf")
    bridge_rows = []
    for node in all_nodes:
        component_connection_score = degree_by_node[node]
        posterior_bridge_score = betweenness[node] * max(incident_probabilities.get(node, [0.0]))
        bridge_rows.append({
            "cell_line_id": node,
            "posterior_bridge_score": posterior_bridge_score,
            "betweenness_centrality": betweenness[node],
            "component_connection_score": component_connection_score,
            "probabilistic_bridge_anchor_status": betweenness[node] >= bridge_cutoff and degree_by_node[node] > 1,
            "notes": "bridge-like status uses unweighted retained probabilistic graph topology",
        })

    graph_root = output_root / "probabilistic_graphs"
    write_tsv(graph_root / "probabilistic_resolved_edges.tsv", resolved_edge_rows, [
        "cell_line_1", "cell_line_2", "posterior_edge_probability", "selected_threshold",
        "expected_false_edge_rate", "retained_in_probabilistic_graph"
    ])
    write_tsv(graph_root / "probabilistic_components.tsv", component_rows, [
        "cell_line_id", "probabilistic_component_id", "component_size", "component_mean_edge_probability",
        "component_min_edge_probability", "component_density", "component_cancer_type_composition"
    ])
    write_tsv(graph_root / "probabilistic_isolates.tsv", isolate_rows, [
        "cell_line_id", "posterior_isolate_probability", "deterministic_isolate_status",
        "probabilistic_isolate_status", "notes"
    ])
    write_tsv(graph_root / "probabilistic_bridge_anchors.tsv", bridge_rows, [
        "cell_line_id", "posterior_bridge_score", "betweenness_centrality",
        "component_connection_score", "probabilistic_bridge_anchor_status", "notes"
    ])
    write_tsv(graph_root / "probabilistic_threshold_sweep.tsv", threshold_rows, [
        "candidate_threshold", "n_retained_edges", "expected_false_edge_rate",
        "passes_expected_false_edge_rate"
    ])

    threshold_report_path = docs_root / "03_threshold_validation_report.md"
    threshold_report_path.parent.mkdir(parents=True, exist_ok=True)
    threshold_report_path.write_text("\n".join([
        "# Threshold validation report",
        "",
        "## Scope",
        "",
        "The threshold validation step evaluated posterior edge-probability thresholds for the probabilistic graph.",
        "",
        "## Threshold selection rule",
        "",
        "For each candidate threshold, the expected false edge rate was computed as the sum of one minus posterior edge probability across retained edges, divided by the number of retained edges.",
        "",
        f"Configured maximum expected false edge rate: {maximum_expected_false_edge_rate}.",
        f"Selected threshold: {selected_threshold}.",
        f"Expected false edge rate at selected threshold: {expected_false_edge_rate}.",
        f"Selection note: {threshold_note}.",
        "",
        "## Output",
        "",
        "- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_threshold_sweep.tsv`",
        "",
        "## Manual inspection required",
        "",
        "- Inspect candidate thresholds with very few retained edges.",
        "- Inspect thresholds near the configured expected false edge-rate limit.",
        "- Inspect whether the selected threshold removes deterministic thesis edges that remain scientifically plausible.",
    ]) + "\n")

    append_report(docs_root / "05_probabilistic_graph_report.md", [
        "",
        "## Probabilistic graph resolution",
        "",
        f"Selected threshold: {selected_threshold}.",
        f"Expected false edge rate at selected threshold: {expected_false_edge_rate}.",
        f"Threshold note: {threshold_note}.",
        "",
        "## Resolution outputs",
        "",
        "- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_resolved_edges.tsv`",
        "- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_components.tsv`",
        "- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_isolates.tsv`",
        "- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_bridge_anchors.tsv`",
        "- `research_framework/optimisation/results/probabilistic_graphs/probabilistic_threshold_sweep.tsv`",
        "",
        "## Manual inspection required",
        "",
        "- Inspect thresholds where the expected false edge rate is close to the configured maximum.",
        "- Inspect probabilistic isolates that differ from deterministic isolate status.",
        "- Inspect bridge-like anchors with high posterior bridge scores before interpreting them biologically.",
    ])


if __name__ == "__main__":
    main()
