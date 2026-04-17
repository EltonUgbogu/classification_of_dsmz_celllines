#!/usr/bin/env bash
set -euo pipefail

# Portable runner for the pipeline
# Works for: local execution OR HPC via Snakemake profile (e.g., slurm)

usage() {
  cat <<'EOF'
Usage:
  ./run_pipeline.sh --profile <name> --executor <local|slurm> [--cores N] [--jobs N] [target...]

Examples:
  ./run_pipeline.sh --profile brca --executor local --cores 8
  ./run_pipeline.sh --profile brca --executor slurm --jobs 200
  ./run_pipeline.sh --profile brca --executor local --cores 2 smoke_test

Notes:
  - Creates run/<profile>/ as the Snakemake working directory
  - Uses profiles/<executor>/config.yaml
  - Passes --config pipeline_profile=<profile> into Snakemake (mandatory)
EOF
}

PROFILE=""
EXECUTOR=""
CORES=""
JOBS=""
TARGETS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile) PROFILE="${2:-}"; shift 2;;
    --executor) EXECUTOR="${2:-}"; shift 2;;
    --cores) CORES="${2:-}"; shift 2;;
    --jobs) JOBS="${2:-}"; shift 2;;
    -h|--help) usage; exit 0;;
    --) shift; TARGETS+=("$@"); break;;
    *) TARGETS+=("$1"); shift;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "${PROFILE}" || -z "${EXECUTOR}" ]]; then
  echo "[ERROR] --profile and --executor are required." >&2
  echo "" >&2
  echo "Available profiles (from config/config.yaml):" >&2
  if [[ -f "${REPO_ROOT}/config/config.yaml" ]]; then
    # Extract profile keys (lines like "  brca:" under "profiles:" section)
    awk '/^profiles:/{flag=1; next} /^[a-z]/ && flag==1 {flag=0} flag==1 && /^  [a-z]+:/ {gsub(/:/, "", $1); print "  -", $1}' "${REPO_ROOT}/config/config.yaml" 2>/dev/null | head -10 || echo "  (could not parse config)" >&2
  else
    echo "  (config/config.yaml not found)" >&2
  fi
  echo "" >&2
  usage
  exit 2
fi

if [[ "${EXECUTOR}" != "local" && "${EXECUTOR}" != "slurm" ]]; then
  echo "[ERROR] --executor must be 'local' or 'slurm'."
  exit 2
fi
RUN_DIR="${REPO_ROOT}/run/${PROFILE}"
SNAKEFILE="${REPO_ROOT}/Snakefile"
CONFIG_YAML="${REPO_ROOT}/config/config.yaml"
PROFILE_DIR="${REPO_ROOT}/profiles/${EXECUTOR}"

# Validate profile exists in config before proceeding
if [[ -f "${CONFIG_YAML}" ]]; then
  # Check if profile exists under profiles: section
  if ! awk '/^profiles:/{flag=1; next} /^[a-z]/ && flag==1 {flag=0} flag==1 && /^  '"${PROFILE}"':/ {found=1} END {exit !found}' "${CONFIG_YAML}" 2>/dev/null; then
    echo "[ERROR] Profile '${PROFILE}' not found in config/config.yaml" >&2
    echo "" >&2
    echo "Available profiles:" >&2
    awk '/^profiles:/{flag=1; next} /^[a-z]/ && flag==1 {flag=0} flag==1 && /^  [a-z]+:/ {gsub(/:/, "", $1); print "  -", $1}' "${CONFIG_YAML}" 2>/dev/null | head -10
    exit 2
  fi
fi

mkdir -p "${RUN_DIR}"

# Setup conda (portable - tries multiple common locations)
CONDA_SETUP=false
if command -v conda &>/dev/null; then
  CONDA_SETUP=true
elif [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  CONDA_SETUP=true
elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
  source "${HOME}/anaconda3/etc/profile.d/conda.sh"
  CONDA_SETUP=true
elif [[ -n "${CONDA_PREFIX:-}" ]]; then
  # Already in a conda environment
  CONDA_SETUP=true
fi

# Try to activate snakemake environment if conda is available
if [[ "${CONDA_SETUP}" == true ]]; then
  if conda env list 2>/dev/null | grep -qE "^smk\s"; then
    conda activate smk 2>/dev/null || true
  fi
fi

# Auto defaults
if [[ -z "${CORES}" && "${EXECUTOR}" == "local" ]]; then
  if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    CORES="${SLURM_CPUS_PER_TASK}"
  elif command -v nproc >/dev/null 2>&1; then
    CORES="$(nproc)"
  else
    CORES="4"
  fi
fi

if [[ -z "${JOBS}" && "${EXECUTOR}" == "slurm" ]]; then
  JOBS="200"
fi

# Build snakemake args
SMK_ARGS=(--snakefile "${SNAKEFILE}"
          --directory "${RUN_DIR}"
          --configfile "${CONFIG_YAML}"
          --config "pipeline_profile=${PROFILE}")

# Add executor profile if it exists
if [[ -d "${PROFILE_DIR}" ]]; then
  SMK_ARGS+=(--profile "${PROFILE_DIR}")
fi

# Add conda support if available
if [[ "${CONDA_SETUP}" == true ]]; then
  SMK_ARGS+=(--use-conda --conda-frontend conda)
fi

# Provide cores/jobs only if relevant / set
if [[ "${EXECUTOR}" == "local" ]]; then
  SMK_ARGS+=(--cores "${CORES}")
else
  SMK_ARGS+=(--jobs "${JOBS}")
fi

# Default target if none provided
if [[ ${#TARGETS[@]} -eq 0 ]]; then
  TARGETS=(all)
fi

# Export for R scripts
export PIPELINE_PROFILE="${PROFILE}"
export SNAKEMAKE_PROFILE="${PROFILE}"

echo "[INFO] Repo root: ${REPO_ROOT}"
echo "[INFO] Run dir:   ${RUN_DIR}"
echo "[INFO] Executor:  ${EXECUTOR}"
echo "[INFO] Profile:   ${PROFILE}"
if [[ "${EXECUTOR}" == "local" ]]; then
  echo "[INFO] Cores:     ${CORES}"
else
  echo "[INFO] Jobs:      ${JOBS}"
fi
echo "[INFO] Targets:   ${TARGETS[*]}"
echo ""

exec snakemake "${SMK_ARGS[@]}" "${TARGETS[@]}"
