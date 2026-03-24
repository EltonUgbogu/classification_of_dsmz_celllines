#!/bin/bash
#SBATCH --job-name=target_nbl_counts
#SBATCH --chdir=/work/ugbogu/pipeline/preprocessing_and_quality_control/scripts/nbl
#SBATCH --output=/work/ugbogu/pipeline/data/nbl/target_nbl/logs/target_nbl_counts_%j.out
#SBATCH --error=/work/ugbogu/pipeline/data/nbl/target_nbl/logs/target_nbl_counts_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00

set -euo pipefail

echo "== $(date) :: START =="
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"

BASE=/work/ugbogu/pipeline/data/nbl/target_nbl
LOG_DIR=/work/ugbogu/pipeline/data/nbl/target_nbl/logs
PROJECT_ID="TARGET-NBL"
DATA_TYPE="Gene Expression Quantification"
WORKFLOW_TYPE="STAR - Counts"

mkdir -p "${BASE}" "${LOG_DIR}"
cd "${BASE}"

echo "Working directory: ${BASE}"

source /home/ugbogu/miniforge3/etc/profile.d/conda.sh
conda activate gdc39

echo "Conda env: ${CONDA_DEFAULT_ENV:-none}"
echo "python: $(command -v python)"
python --version
echo "gdc-client: $(command -v gdc-client)"
gdc-client --version || true

for cmd in curl python gdc-client; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: Missing required command: ${cmd}" >&2
    exit 1
  fi
done

cat > payload_counts.json <<JSON
{
  "filters": {
    "op": "and",
    "content": [
      {
        "op": "=",
        "content": {
          "field": "cases.project.project_id",
          "value": "${PROJECT_ID}"
        }
      },
      {
        "op": "=",
        "content": {
          "field": "files.data_type",
          "value": "${DATA_TYPE}"
        }
      },
      {
        "op": "=",
        "content": {
          "field": "files.analysis.workflow_type",
          "value": "${WORKFLOW_TYPE}"
        }
      }
    ]
  },
  "format": "TSV",
  "fields": "analysis.workflow_type,cases.project.project_id,cases.samples.portions.analytes.aliquots.submitter_id,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.tissue_type,cases.submitter_id,data_category,data_type,experimental_strategy,file_id,file_name,id",
  "size": 2000
}
JSON

echo "== $(date) :: Querying GDC metadata table =="
curl --fail --silent --show-error \
  --request POST \
  --header "Content-Type: application/json" \
  --data @payload_counts.json \
  "https://api.gdc.cancer.gov/files" \
  > target_nbl_star_counts_mapping.tsv

if [[ ! -s target_nbl_star_counts_mapping.tsv ]]; then
  echo "ERROR: target_nbl_star_counts_mapping.tsv is missing or empty." >&2
  exit 1
fi

cat > payload_counts_manifest.json <<JSON
{
  "filters": {
    "op": "and",
    "content": [
      {
        "op": "=",
        "content": {
          "field": "cases.project.project_id",
          "value": "${PROJECT_ID}"
        }
      },
      {
        "op": "=",
        "content": {
          "field": "files.data_type",
          "value": "${DATA_TYPE}"
        }
      },
      {
        "op": "=",
        "content": {
          "field": "files.analysis.workflow_type",
          "value": "${WORKFLOW_TYPE}"
        }
      }
    ]
  },
  "format": "TSV",
  "fields": "file_id,file_name"
}
JSON

echo "== $(date) :: Querying GDC manifest =="
curl --fail --silent --show-error \
  --request POST \
  --header "Content-Type: application/json" \
  --data @payload_counts_manifest.json \
  "https://api.gdc.cancer.gov/files?return_type=manifest" \
  > gdc_manifest_target_nbl_star_counts.txt

if [[ ! -s gdc_manifest_target_nbl_star_counts.txt ]]; then
  echo "ERROR: gdc_manifest_target_nbl_star_counts.txt is missing or empty." >&2
  exit 1
fi

echo "== $(date) :: Building file-to-aliquot mapping =="
python - <<'PY'
import csv
import sys

inp = "target_nbl_star_counts_mapping.tsv"
out = "target_nbl_file_to_aliquot.tsv"
aliquot_list = "target_nbl_aliquot_ids.txt"
aliquot_unique = "target_nbl_aliquot_ids_unique.txt"

with open(inp, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames or []

    required = [
        "file_id",
        "file_name",
        "cases.0.submitter_id",
        "cases.0.samples.0.submitter_id",
        "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id",
    ]
    missing = [x for x in required if x not in fieldnames]
    if missing:
        sys.stderr.write("ERROR: Missing expected columns:\n")
        sys.stderr.write("\n".join(missing) + "\n")
        sys.stderr.write("Observed columns:\n")
        sys.stderr.write("\n".join(fieldnames) + "\n")
        sys.exit(1)

    rows = list(reader)

with open(out, "w", newline="") as g:
    writer = csv.writer(g, delimiter="\t")
    writer.writerow([
        "file_id",
        "file_name",
        "case_submitter_id",
        "sample_submitter_id",
        "aliquot_submitter_id"
    ])

    aliquots = []
    for row in rows:
        aliquot = row["cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"]
        writer.writerow([
            row["file_id"],
            row["file_name"],
            row["cases.0.submitter_id"],
            row["cases.0.samples.0.submitter_id"],
            aliquot
        ])
        if aliquot:
            aliquots.append(aliquot)

with open(aliquot_list, "w") as h:
    for a in aliquots:
        h.write(a + "\n")

with open(aliquot_unique, "w") as h:
    for a in sorted(set(aliquots)):
        h.write(a + "\n")

print(f"rows={len(rows)}")
print(f"non_empty_aliquot_ids={len(aliquots)}")
print(f"unique_aliquot_ids={len(set(aliquots))}")
PY

echo "== $(date) :: Summary =="
echo "Mapping rows      : $(tail -n +2 target_nbl_file_to_aliquot.tsv | wc -l)"
echo "Unique aliquots   : $(wc -l < target_nbl_aliquot_ids_unique.txt)"
echo "Manifest entries  : $(tail -n +2 gdc_manifest_target_nbl_star_counts.txt | wc -l)"
head -n 5 target_nbl_file_to_aliquot.tsv || true

echo "== $(date) :: Downloading TARGET-NBL STAR counts =="
gdc-client download \
  -n 1 \
  --log-file download.log \
  -m gdc_manifest_target_nbl_star_counts.txt

echo "== $(date) :: DONE =="
echo "Outputs in: ${BASE}"