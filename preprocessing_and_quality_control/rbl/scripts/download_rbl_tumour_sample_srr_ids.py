#!/usr/bin/env python3
# =============================================================================
# download_rbl_tumour_sample_srr_ids.py
#
# SRR accession lists for retinoblastoma (RBL) RNA-seq cohorts (NCBI SRA / GEO).
#
# Supported targets:
#   GSE111168  — primary tumour samples identified from GEO titles (tumor1–3)
#   GSE196420  — bulk primary retinoblastoma RNA-seq (title filter)
#   GSE268136  — primary tumour RB_* sample titles from GEO
#
# Outputs per dataset include all_primary_tumour_srr.txt, all_primary_tumour_pe.txt,
# metadata/SraRunInfo*.csv, and GEO gsm/primary GSM tables where applicable.
#
# Usage:
#   python download_rbl_tumour_sample_srr_ids.py GSE111168|GSE196420|GSE268136|all
# =============================================================================

import csv
import gzip
import logging
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[3]
RBL_DATA_ROOT = Path(os.environ.get("RBL_DATA_ROOT", REPO_ROOT / "data" / "rbl")).expanduser()


@dataclass(frozen=True)
class DatasetConfig:
    query: str
    outdir: Path
    mode: str  # "geo_title_filter" or "runinfo_filter"
    title_pattern: Optional[str] = None
    runinfo_gsm_fields: tuple[str, ...] = ("LibraryName", "SampleName")
    require_strategy: str = "RNA-Seq"
    require_layout: str = "PAIRED"
    require_species: str = "Homo sapiens"


DATASETS: dict[str, DatasetConfig] = {
    "GSE111168": DatasetConfig(
        query="PRJNA436090",
        outdir=RBL_DATA_ROOT / "GSE111168_primary_tumours",
        mode="geo_title_filter",
        # GEO uses American spelling (tumor1, tumor2, tumor3), not tumour*.
        title_pattern=r"^tumor[123]$",
        runinfo_gsm_fields=("LibraryName", "SampleName"),
    ),
    "GSE196420": DatasetConfig(
        query="PRJNA804875",
        outdir=RBL_DATA_ROOT / "GSE196420_primary_tumours",
        mode="geo_title_filter",
        title_pattern=r"Primary Retinoblastoma - RNA-Seq",
        runinfo_gsm_fields=("LibraryName", "SampleName"),
    ),
    "GSE268136": DatasetConfig(
        query="PRJNA1114744",
        outdir=RBL_DATA_ROOT / "GSE268136_primary_tumours",
        mode="geo_title_filter",
        title_pattern=r"^RB_\d+$",
        runinfo_gsm_fields=("LibraryName", "SampleName"),
    ),
}

REQUIRED_RUNINFO_COLUMNS = {
    "Run",
    "LibraryLayout",
    "LibraryStrategy",
    "ScientificName",
}


def run_shell(cmd: str) -> str:
    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        text=True,
        capture_output=True,
        executable="/bin/bash",
    )
    return result.stdout


def is_valid_runinfo_csv(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with path.open(newline="") as handle:
            header = set(next(csv.reader(handle)))
        return REQUIRED_RUNINFO_COLUMNS.issubset(header)
    except Exception:
        return False


def download_runinfo(query: str, out_csv: Path) -> None:
    stdout = run_shell(f'esearch -db sra -query "{query}" | efetch -format runinfo')
    out_csv.write_text(stdout)
    if not is_valid_runinfo_csv(out_csv):
        raise RuntimeError(f"Invalid or empty RunInfo CSV: {out_csv}")


def geo_family_soft_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/soft/{gse}_family.soft.gz"


def ensure_family_soft(gse: str, outdir: Path) -> Path:
    family_soft = outdir / f"{gse}_family.soft.gz"
    if family_soft.exists() and family_soft.stat().st_size > 0:
        log.info("Using cached SOFT file: %s", family_soft)
        return family_soft

    url = geo_family_soft_url(gse)
    log.info("Downloading SOFT file: %s", url)
    subprocess.run(["wget", "-O", str(family_soft), url], check=True)
    return family_soft


def extract_gsm_titles(family_soft: Path, gsm_titles_tsv: Path) -> dict[str, str]:
    gsm_titles: dict[str, str] = {}
    current_gsm: str | None = None

    with gzip.open(family_soft, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                current_gsm = line.split(" = ", 1)[1]
            elif current_gsm and line.startswith("!Sample_title = "):
                gsm_titles[current_gsm] = line.split(" = ", 1)[1]

    with gsm_titles_tsv.open("w") as out:
        for gsm, title in gsm_titles.items():
            out.write(f"{gsm}\t{title}\n")

    log.info("Extracted %d GSM titles from %s", len(gsm_titles), family_soft.name)
    return gsm_titles


def select_primary_gsms(
    gsm_titles: dict[str, str],
    title_pattern: str,
    primary_gsms_tsv: Path,
) -> set[str]:
    pattern = re.compile(title_pattern, flags=re.IGNORECASE)
    selected = {gsm: title for gsm, title in gsm_titles.items() if pattern.search(title)}

    with primary_gsms_tsv.open("w") as out:
        for gsm, title in selected.items():
            out.write(f"{gsm}\t{title}\n")

    log.info(
        "Selected %d primary GSMs matching pattern '%s'",
        len(selected),
        title_pattern,
    )
    return set(selected.keys())


def resolve_gsm(row: dict[str, str], candidate_fields: tuple[str, ...]) -> str:
    for field_name in candidate_fields:
        value = row.get(field_name, "").strip()
        if value:
            return value
    return ""


def write_runinfo_filter_outputs(
    cfg: DatasetConfig,
    runinfo_csv: Path,
    out_all_txt: Path,
    out_all_csv: Path,
    out_pe_txt: Path,
    out_pe_csv: Path,
) -> tuple[int, int]:
    all_hits: list[str] = []
    pe_hits: list[str] = []

    with runinfo_csv.open(newline="") as f, \
         out_all_csv.open("w", newline="") as all_handle, \
         out_pe_csv.open("w", newline="") as pe_handle:

        reader = csv.DictReader(f)
        writer_all = csv.DictWriter(all_handle, fieldnames=reader.fieldnames)
        writer_pe = csv.DictWriter(pe_handle, fieldnames=reader.fieldnames)
        writer_all.writeheader()
        writer_pe.writeheader()

        for row in reader:
            if (
                row.get("LibraryStrategy", "").strip() == cfg.require_strategy
                and row.get("ScientificName", "").strip() == cfg.require_species
            ):
                run = row.get("Run", "").strip()
                if not run:
                    continue

                all_hits.append(run)
                writer_all.writerow(row)

                if row.get("LibraryLayout", "").strip() == cfg.require_layout:
                    pe_hits.append(run)
                    writer_pe.writerow(row)

    with out_all_txt.open("w") as out:
        for run in all_hits:
            out.write(run + "\n")

    with out_pe_txt.open("w") as out:
        for run in pe_hits:
            out.write(run + "\n")

    return len(all_hits), len(pe_hits)


def write_geo_title_filter_outputs(
    cfg: DatasetConfig,
    runinfo_csv: Path,
    primary_gsms: set[str],
    out_all_txt: Path,
    out_all_csv: Path,
    out_pe_txt: Path,
    out_pe_csv: Path,
) -> tuple[int, int]:
    all_hits: list[str] = []
    pe_hits: list[str] = []

    with runinfo_csv.open(newline="") as f, \
         out_all_csv.open("w", newline="") as all_handle, \
         out_pe_csv.open("w", newline="") as pe_handle:

        reader = csv.DictReader(f)
        writer_all = csv.DictWriter(all_handle, fieldnames=reader.fieldnames)
        writer_pe = csv.DictWriter(pe_handle, fieldnames=reader.fieldnames)
        writer_all.writeheader()
        writer_pe.writeheader()

        for row in reader:
            gsm = resolve_gsm(row, cfg.runinfo_gsm_fields)

            if (
                row.get("LibraryStrategy", "").strip() == cfg.require_strategy
                and row.get("ScientificName", "").strip() == cfg.require_species
                and gsm in primary_gsms
            ):
                run = row.get("Run", "").strip()
                if not run:
                    continue

                all_hits.append(run)
                writer_all.writerow(row)

                if row.get("LibraryLayout", "").strip() == cfg.require_layout:
                    pe_hits.append(run)
                    writer_pe.writerow(row)

    with out_all_txt.open("w") as out:
        for run in all_hits:
            out.write(run + "\n")

    with out_pe_txt.open("w") as out:
        for run in pe_hits:
            out.write(run + "\n")

    return len(all_hits), len(pe_hits)


def process_dataset(name: str) -> int:
    if name not in DATASETS:
        log.error("Unsupported dataset '%s'. Supported: %s", name, ", ".join(DATASETS))
        return 1

    cfg = DATASETS[name]
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    metadata_dir = cfg.outdir / "metadata"
    metadata_dir.mkdir(exist_ok=True)

    runinfo_csv = metadata_dir / "SraRunInfo.csv"
    out_all_txt = cfg.outdir / "all_primary_tumour_srr.txt"
    out_all_csv = metadata_dir / "SraRunInfo_primary.csv"
    out_pe_txt = cfg.outdir / "all_primary_tumour_pe.txt"
    out_pe_csv = metadata_dir / "SraRunInfo_primary_PE.csv"

    if is_valid_runinfo_csv(runinfo_csv):
        log.info("Using cached RunInfo: %s", runinfo_csv)
    else:
        log.info("Downloading RunInfo for %s (%s)", name, cfg.query)
        download_runinfo(cfg.query, runinfo_csv)

    if cfg.mode == "runinfo_filter":
        total_rna, total_pe = write_runinfo_filter_outputs(
            cfg=cfg,
            runinfo_csv=runinfo_csv,
            out_all_txt=out_all_txt,
            out_all_csv=out_all_csv,
            out_pe_txt=out_pe_txt,
            out_pe_csv=out_pe_csv,
        )
        log.info("%s", name)
        log.info("Tumour RNA-seq runs      : %d", total_rna)
        log.info("Paired-end tumour RNA-seq: %d", total_pe)
        log.info("Wrote: %s", out_all_txt)
        log.info("Wrote: %s", out_pe_txt)
        return 0

    if cfg.mode == "geo_title_filter":
        if cfg.title_pattern is None:
            raise RuntimeError(f"title_pattern is required for {name}")

        family_soft = ensure_family_soft(name, cfg.outdir)
        gsm_titles_tsv = metadata_dir / "gsm_titles.tsv"
        primary_gsms_tsv = metadata_dir / "primary_gsms.tsv"

        gsm_titles = extract_gsm_titles(family_soft, gsm_titles_tsv)
        primary_gsms = select_primary_gsms(
            gsm_titles=gsm_titles,
            title_pattern=cfg.title_pattern,
            primary_gsms_tsv=primary_gsms_tsv,
        )

        if not primary_gsms:
            raise RuntimeError(f"No primary GSMs matched title_pattern for {name}")

        total_rna, total_pe = write_geo_title_filter_outputs(
            cfg=cfg,
            runinfo_csv=runinfo_csv,
            primary_gsms=primary_gsms,
            out_all_txt=out_all_txt,
            out_all_csv=out_all_csv,
            out_pe_txt=out_pe_txt,
            out_pe_csv=out_pe_csv,
        )

        log.info("%s", name)
        log.info("Primary GSMs matched     : %d", len(primary_gsms))
        log.info("Tumour RNA-seq runs      : %d", total_rna)
        log.info("Paired-end tumour RNA-seq: %d", total_pe)
        log.info("Wrote: %s", out_all_txt)
        log.info("Wrote: %s", out_pe_txt)
        return 0

    raise RuntimeError(f"Unknown mode: {cfg.mode}")


def main() -> int:
    if len(sys.argv) != 2:
        print(
            f"Usage: {Path(sys.argv[0]).name} "
            f"{' | '.join(DATASETS)} | all"
        )
        return 1

    target = sys.argv[1]

    if target == "all":
        rc = 0
        for name in DATASETS:
            log.info("=== Processing %s ===", name)
            try:
                ret = process_dataset(name)
                if ret != 0:
                    rc = ret
            except Exception as exc:
                log.error("%s: %s", name, exc)
                rc = 1
        return rc

    try:
        return process_dataset(target)
    except Exception as exc:
        log.error("%s: %s", target, exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
