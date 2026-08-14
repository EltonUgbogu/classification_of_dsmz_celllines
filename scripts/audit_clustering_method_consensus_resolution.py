#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path
import yaml


def read_profiled_config(config_path: Path, profile: str):
    cfg = yaml.safe_load(config_path.read_text())
    defaults = cfg.get("defaults", {})
    profiles = cfg.get("profiles", {})
    if profile not in profiles:
        raise SystemExit(f"Profile not found in config: {profile}")
    merged = {}
    merged.update(defaults)
    merged.update(profiles[profile])
    return cfg, merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--profile", required=True)
    ap.add_argument("--directions", required=True)
    ap.add_argument("--out_tsv", required=True)
    args = ap.parse_args()

    config_path = Path(args.config)
    cfg_full = yaml.safe_load(config_path.read_text())
    profiles = cfg_full.get("profiles", {})
    prof = profiles.get(args.profile)
    if prof is None:
        raise SystemExit(f"Profile not found: {args.profile}")

    unsup_root = Path("/work/ugbogu/pipeline") / prof["paths"]["unsup_root"]
    defaults = cfg_full.get("defaults", {})
    graph_cfg = defaults.get("patient_referenced_graph")
    if not isinstance(graph_cfg, dict):
        raise SystemExit("Missing required config section: defaults.patient_referenced_graph")
    if "p_consensus_threshold" not in graph_cfg:
        raise SystemExit("Missing required config key: defaults.patient_referenced_graph.p_consensus_threshold")
    if "clustering_methods_by_distance" not in graph_cfg:
        raise SystemExit("Missing required config key: defaults.patient_referenced_graph.clustering_methods_by_distance")
    threshold = float(graph_cfg["p_consensus_threshold"])
    if not (0 < threshold <= 1):
        raise SystemExit("defaults.patient_referenced_graph.p_consensus_threshold must satisfy 0 < threshold <= 1")
    methods_by_distance = graph_cfg["clustering_methods_by_distance"]
    directions = [d.strip() for d in args.directions.split(",") if d.strip()]

    rows = []
    for direction in directions:
        distance = direction.rsplit("_", 1)[-1]
        tsv_path = unsup_root / "tumour_neighbourhoods" / direction / "final_consensus" / f"Final_consensus_tumour_neighbourhoods_{direction}.tsv"
        if not tsv_path.exists():
            raise SystemExit(f"Missing final consensus TSV: {tsv_path}")
        declared_methods = methods_by_distance.get(distance)
        if not declared_methods:
            raise SystemExit(f"No declared clustering methods for distance '{distance}'")
        neighbourhood_root = unsup_root / "tumour_neighbourhoods" / direction
        discovered_method_dirs = sorted(
            p.name for p in neighbourhood_root.iterdir()
            if p.is_dir() and p.name != "final_consensus"
        )
        declared_method_set = sorted(str(x) for x in declared_methods)
        discovered_method_set = sorted(
            method_id for method_id in discovered_method_dirs
            if list((neighbourhood_root / method_id).glob("Top_m_long_*.csv"))
        )
        missing_declared_methods = sorted(set(declared_method_set) - set(discovered_method_set))
        unexpected_discovered_methods = sorted(set(discovered_method_set) - set(declared_method_set))
        if missing_declared_methods or unexpected_discovered_methods:
            raise SystemExit(
                "Declared clustering methods do not match discovered neighbourhood files for "
                f"{direction}. missing={missing_declared_methods} unexpected={unexpected_discovered_methods}"
            )
        observed_n_methods_values = set()
        with tsv_path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for idx, row in enumerate(reader):
                if row.get("n_methods"):
                    try:
                      observed_n_methods_values.add(int(float(row["n_methods"])))
                    except ValueError:
                      pass
                if idx > 5000 and observed_n_methods_values:
                    break
        observed_n_methods = sorted(observed_n_methods_values)
        a_r = len(declared_method_set)
        possible_values = (
            ";".join(
                str(round(k / a_r, 10)).rstrip("0").rstrip(".")
                for k in range(a_r + 1)
            )
            if a_r is not None else ""
        )
        rows.append({
            "profile": args.profile,
            "representation": direction,
            "distance": distance,
            "n_clustering_methods": a_r,
            "declared_clustering_methods": ";".join(declared_method_set),
            "discovered_clustering_methods": ";".join(discovered_method_set),
            "method_file_validation": "exact_match",
            "possible_p_consensus_values": possible_values,
            "p_consensus_threshold": threshold,
            "consensus_requirement_scope": "profile_level",
            "minimum_method_count": None if a_r is None else math.ceil(threshold * a_r),
            "effective_requirement": None if a_r is None else f"{math.ceil(threshold * a_r)}/{a_r}",
            "observed_n_methods_values": ";".join(str(x) for x in observed_n_methods),
        })

    with Path(args.out_tsv).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
