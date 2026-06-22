#!/usr/bin/env python3
"""Build enrichment queries for the active ranked marker-source panel."""
import argparse
import csv
import hashlib
import re
from collections import defaultdict
from pathlib import Path

METHOD = "ranked_marker_source_pan_cancer_panel"
PREFIX = "ranked_marker_source_panel"
MIN_QUERY = 3
PROFILES = ("brca", "nbl", "rbl")
EVIDENCE = (
    "anchor_source_recurrent",
    "isolate_source_recurrent",
    "anchor_singleton_ranked",
    "ranked_nonrecurrent_marker",
)

def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path, columns, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in data:
            writer.writerow({key: row.get(key, "") for key in columns})

def clean_id(value):
    return re.sub(r"\.\d+$", "", str(value or "").strip())

def unique(values):
    seen, out = set(), []
    for value in values:
        value = clean_id(value)
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out

def split_field(value):
    return {x for x in str(value or "").split(";") if x}

def safe(value):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_").lower()

def numeric(value):
    try:
        x = float(value)
        return None if x != x else x
    except (TypeError, ValueError):
        return None

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def read_gene_file(path):
    genes = []
    for i, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines()):
        if not line.strip():
            continue
        value = line.split("\t", 1)[0].strip()
        if i == 0 and value.lower() in {"gene", "gene_id", "ensembl_gene_id"}:
            continue
        genes.append(value)
    return unique(genes)

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pipeline-root", default=Path(__file__).resolve().parents[1])
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-query-size", type=int, default=MIN_QUERY)
    p.add_argument("--background-min-normalised-count", type=float, default=1.0)
    return p.parse_args()

def main():
    args = parse_args()
    root = Path(args.pipeline_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    feature_path = root / "results/unsupervised/pan_cancer/feature_space/pan_cancer_features.tsv"
    clean_path = root / "results/unsupervised/pan_cancer/feature_space/pan_cancer_features_clean.txt"
    ranking_path = root / "results/unsupervised/pan_cancer/feature_space/ranked_marker_source_panel_ranking_components.tsv"
    feature_rows = read_tsv(feature_path)
    ranking_rows = read_tsv(ranking_path)
    clean = read_gene_file(clean_path)
    assert len(feature_rows) == len(clean) == 379
    assert [clean_id(x["clean_gene_id"]) for x in feature_rows] == clean
    assert {x["method"] for x in feature_rows} == {METHOD}

    panel_universe = unique(x["gene_id"] for x in ranking_rows)
    if not panel_universe:
        raise SystemExit("eligible marker-selection universe is empty")
    panel_set = set(clean)
    if not panel_set.issubset(set(panel_universe)):
        missing = sorted(panel_set - set(panel_universe))
        raise SystemExit(f"selected genes missing from eligible universe: {missing[:10]}")

    gene_root = outdir / "query_genes"
    bg_root = outdir / "backgrounds"
    gene_root.mkdir(parents=True, exist_ok=True)
    bg_root.mkdir(parents=True, exist_ok=True)
    manifest = []
    skipped = []
    background_manifest = []
    background_seen = {}

    def write_genes(path, genes):
        write_tsv(path, ["gene_id"], [{"gene_id": x} for x in unique(genes)])

    panel_bg_name = f"{PREFIX}__eligible_marker_selection_universe"
    panel_bg_path = bg_root / f"{panel_bg_name}.tsv"
    write_genes(panel_bg_path, panel_universe)
    background_seen[panel_bg_name] = str(panel_bg_path)
    background_manifest.append({
        "background_name": panel_bg_name,
        "background_type": "panel_eligible_marker_selection_universe",
        "cohort": "ALL", "contrast": "", "gene_count": len(panel_universe),
        "background_path": str(panel_bg_path), "sha256": sha(panel_bg_path),
        "strategy": "unique clean genes in the active ranked-marker-source ranking-component table",
    })

    def add_query(label, category, genes, *, cohort="ALL", family="", evidence="", direction="ALL", contrast="", background=None, background_name=panel_bg_name, background_strategy="active eligible marker-selection universe", source_marker_table="", rank_rows=None):
        query_id = f"{PREFIX}__{safe(label)}"
        genes = unique(genes)
        background = unique(background if background is not None else panel_universe)
        qpath = gene_root / f"{query_id}.tsv"
        bpath = Path(background_seen.get(background_name, "")) if background_name in background_seen else bg_root / f"{safe(background_name)}.tsv"
        if background_name not in background_seen:
            write_genes(bpath, background)
            background_seen[background_name] = str(bpath)
        write_genes(qpath, genes)
        rpath = gene_root / f"{query_id}_ranked.tsv"
        ranks = rank_rows or []
        if ranks:
            write_tsv(rpath, ["gene_id", "rank_stat", "rank_source"], ranks)
        else:
            write_tsv(rpath, ["gene_id", "rank_stat", "rank_source"], [])
        reason = "" if len(genes) >= args.min_query_size else f"query_size_below_{args.min_query_size}"
        row = {
            "method": METHOD, "query_id": query_id, "query_name": query_id,
            "category": category, "query_family": category, "cohort": cohort,
            "marker_source_class": family, "evidence_class": evidence,
            "contrast": contrast, "direction": direction,
            "gene_count": len(genes), "gene_list_path": str(qpath),
            "background_name": background_name, "background_count": len(background),
            "background_path": str(bpath), "background_strategy": background_strategy,
            "source_marker_table": source_marker_table,
            "ranked_genes_path": str(rpath), "ordered": "FALSE",
            "custom_background_available": "TRUE" if background else "FALSE",
            "skip": "TRUE" if reason else "FALSE", "skipped_reason": reason,
            "minimum_query_size": args.min_query_size,
        }
        manifest.append(row)
        if reason:
            skipped.append(row.copy())

    def panel_rows(cohort=None, family=None, evidence=None, direction=None):
        out = []
        for row in feature_rows:
            # Cohort-stratified panel queries use the mutually exclusive owner
            # assignment, not the non-exclusive retained-cohort membership field.
            if cohort and row["owner_profile"].lower() != cohort.lower(): continue
            if family and family not in split_field(row["marker_source_class"]): continue
            if evidence and evidence not in split_field(row["evidence_classes"]): continue
            if direction and row["direction"].upper() != direction.upper(): continue
            out.append(row)
        return out

    def ids(rows): return [x["clean_gene_id"] for x in rows]

    add_query("final_panel_all", "final_panel", clean)
    for cohort in PROFILES:
        add_query(f"cohort_owner_{cohort}", "cohort_owner", ids(panel_rows(cohort=cohort)), cohort=cohort.upper())
    for family in ("anchor", "isolate"):
        add_query(f"marker_source_class_{family}", "marker_source_class", ids(panel_rows(family=family)), family=family)
        for cohort in PROFILES:
            add_query(f"cohort_{cohort}_marker_source_class_{family}", "cohort_marker_source_class", ids(panel_rows(cohort=cohort, family=family)), cohort=cohort.upper(), family=family)
    for evidence in EVIDENCE:
        add_query(f"evidence_{evidence}", "evidence_class", ids(panel_rows(evidence=evidence)), evidence=evidence)
        for cohort in PROFILES:
            add_query(f"cohort_{cohort}_evidence_{evidence}", "cohort_evidence_class", ids(panel_rows(cohort=cohort, evidence=evidence)), cohort=cohort.upper(), evidence=evidence)
    for direction in ("UP", "DOWN", "MIXED"):
        add_query(f"direction_{direction}", "direction", ids(panel_rows(direction=direction)), direction=direction)
        for cohort in PROFILES:
            add_query(f"cohort_{cohort}_direction_{direction}", "cohort_direction", ids(panel_rows(cohort=cohort, direction=direction)), cohort=cohort.upper(), direction=direction)
        for family in ("anchor", "isolate"):
            add_query(f"marker_source_class_{family}_direction_{direction}", "marker_source_class_direction", ids(panel_rows(family=family, direction=direction)), family=family, direction=direction)
            for cohort in PROFILES:
                add_query(f"cohort_{cohort}_marker_source_class_{family}_direction_{direction}", "cohort_marker_source_class_direction", ids(panel_rows(cohort=cohort, family=family, direction=direction)), cohort=cohort.upper(), family=family, direction=direction)

    # Per-contrast sets and exact DESeq2-tested backgrounds.
    rank_by_key = {(x["cohort"], x["marker_source_class"], clean_id(x["gene_id"])) for x in ranking_rows}
    per_contrast_validation = []
    for profile in PROFILES:
        base = root / f"results/unsupervised/{profile}/deseq2_markers"
        manifest_path = base / "markers/marker_sets_manifest.tsv"
        for item in read_tsv(manifest_path):
            contrast = item["contrast"]
            family = "anchor" if contrast.startswith("anchor_") else "isolate"
            marker_path = base / item["marker_file"]
            table_path = base / item["table_file"]
            marker_genes = read_gene_file(marker_path)
            table_rows = read_tsv(table_path)
            gene_col = "gene_id" if "gene_id" in table_rows[0] else next(iter(table_rows[0]))
            norm_col = next((x for x in ("normalised_count_in_test_sample", "normalized_count_in_test_sample", "test_sample_normalised_count", "test_sample_normalized_count") if x in table_rows[0]), "")
            if not norm_col:
                raise SystemExit(f"no test-sample normalised-count column: {table_path}")
            info = {}
            background = []
            for row in table_rows:
                gene = clean_id(row.get(gene_col))
                padj = numeric(row.get("padj"))
                norm = numeric(row.get(norm_col))
                lfc = numeric(row.get("log2FoldChange"))
                if gene: info[gene] = lfc
                if gene and padj is not None and norm is not None and norm >= args.background_min_normalised_count:
                    background.append(gene)
            background = unique(background)
            mapping_missing = [g for g in marker_genes if (profile, family, g) not in rank_by_key]
            if mapping_missing:
                raise SystemExit(f"per-contrast marker mapping failed for {profile}:{contrast}: {mapping_missing[:10]}")
            bg_name = f"{PREFIX}__{profile}__{safe(contrast)}__eligible_deseq2_background"
            bg_path = bg_root / f"{bg_name}.tsv"
            write_genes(bg_path, background)
            background_seen[bg_name] = str(bg_path)
            background_manifest.append({
                "background_name": bg_name, "background_type": "per_contrast_eligible_deseq2_background",
                "cohort": profile.upper(), "contrast": contrast, "gene_count": len(background),
                "background_path": str(bg_path), "sha256": sha(bg_path),
                "strategy": f"DESeq2-tested genes with non-missing adjusted p-value and {norm_col} >= {args.background_min_normalised_count}",
            })
            for direction, query in (
                ("ALL", marker_genes),
                ("UP", [g for g in marker_genes if (info.get(g) or 0) > 0]),
                ("DOWN", [g for g in marker_genes if (info.get(g) or 0) < 0]),
            ):
                ranks = [{"gene_id": g, "rank_stat": abs(info[g]), "rank_source": "abs_log2FoldChange"} for g in query if info.get(g) is not None]
                ranks.sort(key=lambda x: (-x["rank_stat"], x["gene_id"]))
                add_query(f"contrast_{profile}_{contrast}_{direction}", "per_contrast", query,
                          cohort=profile.upper(), family=family, direction=direction, contrast=contrast,
                          background=background, background_name=bg_name,
                          background_strategy="contrast-specific eligible DESeq2-tested background",
                          source_marker_table=str(marker_path), rank_rows=ranks)
            per_contrast_validation.append({"cohort": profile.upper(), "contrast": contrast, "marker_genes": len(marker_genes), "mapped_marker_genes": len(marker_genes)-len(mapping_missing), "background_genes": len(background), "status": "PASS"})

    manifest_columns = ["method", "query_id", "query_name", "category", "query_family", "cohort", "marker_source_class", "evidence_class", "contrast", "direction", "gene_count", "gene_list_path", "background_name", "background_count", "background_path", "background_strategy", "source_marker_table", "ranked_genes_path", "ordered", "custom_background_available", "skip", "skipped_reason", "minimum_query_size"]
    write_tsv(outdir / f"{PREFIX}_enrichment_query_manifest.tsv", manifest_columns, manifest)
    write_tsv(outdir / f"{PREFIX}_enrichment_background_manifest.tsv", ["background_name", "background_type", "cohort", "contrast", "gene_count", "background_path", "sha256", "strategy"], background_manifest)
    write_tsv(outdir / f"{PREFIX}_enrichment_skipped_queries.tsv", manifest_columns, skipped)
    write_tsv(outdir / f"{PREFIX}_per_contrast_background_validation.tsv", ["cohort", "contrast", "marker_genes", "mapped_marker_genes", "background_genes", "status"], per_contrast_validation)
    write_tsv(outdir / f"{PREFIX}_enrichment_preparation_metrics.tsv", ["metric", "value"], [
        {"metric": "method", "value": METHOD},
        {"metric": "selected_panel_genes", "value": len(clean)},
        {"metric": "eligible_panel_universe_genes", "value": len(panel_universe)},
        {"metric": "minimum_query_size", "value": args.min_query_size},
        {"metric": "total_queries", "value": len(manifest)},
        {"metric": "runnable_queries", "value": sum(x["skip"] == "FALSE" for x in manifest)},
        {"metric": "skipped_queries", "value": len(skipped)},
        {"metric": "per_contrast_backgrounds", "value": len(per_contrast_validation)},
    ])
    print(f"queries={len(manifest)} runnable={sum(x['skip']=='FALSE' for x in manifest)} skipped={len(skipped)}")
    print(f"panel_background={len(panel_universe)} per_contrast={len(per_contrast_validation)}")

if __name__ == "__main__":
    main()
