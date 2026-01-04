#!/bin/bash

# run_snakemake_exprtf.sh

#SBATCH --job-name=smk_exprtf
#SBATCH --partition=cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/work/ugbogu/expr_transformer/logs/%x_%j.out
#SBATCH --error=/work/ugbogu/expr_transformer/logs/%x_%j.err

set -euo pipefail

# Purge any Spack modules that might interfere with conda
module purge || true

# -------------------------------------------------------------------
# Snakemake controller (CPU) launching GPU rules via Slurm profile
# -------------------------------------------------------------------

REPO_ROOT="${REPO_ROOT:-/work/${USER}/expr_transformer}"
SNAKEFILE="${SNAKEFILE:-${REPO_ROOT}/src/Snakefile}"
WF_PROFILE="${WF_PROFILE:-${REPO_ROOT}/profiles/slurm}"

LOGDIR="${LOGDIR:-${REPO_ROOT}/logs}"
mkdir -p "$LOGDIR"

# Conda activation for the controller job
CONDA_SH="${CONDA_SH:-/work/${USER}/software/miniforge3/etc/profile.d/conda.sh}"
CONDA_ENV="${SMK_CONDA_ENV:-smk}"   # controller env must have snakemake

# Snakemake knobs
CORES="${CORES:-4}"
JOBS="${JOBS:-50}"
LATENCY_WAIT="${LATENCY_WAIT:-300}"  # Increased for HPC filesystem latency
RESTART_TIMES="${RESTART_TIMES:-1}"

# Optional per-rule conda environments (default ON to match profile config)
USE_CONDA="${USE_CONDA:-1}"
# Use dedicated conda prefix outside any conda env tree to avoid nesting errors
CONDA_PREFIX="${CONDA_PREFIX:-/work/${USER}/snakemake_conda_envs/expr_transformer}"
mkdir -p "$CONDA_PREFIX"

# Force shared TMPDIR so cluster jobs can see Snakemake temp files
export TMPDIR="${REPO_ROOT}/.snakemake/tmp"
mkdir -p "$TMPDIR"

TARGETS="${TARGETS:-}"

die(){ echo "ERROR: $*" >&2; exit 1; }

[[ -d "$REPO_ROOT" ]] || die "REPO_ROOT not found: $REPO_ROOT"
[[ -f "$SNAKEFILE" ]] || die "Snakefile not found: $SNAKEFILE"
[[ -d "$WF_PROFILE" ]] || die "Workflow profile not found: $WF_PROFILE"
[[ -f "$CONDA_SH" ]] || die "conda.sh not found: $CONDA_SH"

# Temporarily disable nounset for conda activation (MKL hooks may reference unset vars)
set +u
# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate "$CONDA_ENV"
set -u

# --- HARD PIN conda to /work miniforge (avoid Spack conda) ---
MINIFORGE_ROOT="/work/${USER}/software/miniforge3"

export PATH="${MINIFORGE_ROOT}/condabin:${MINIFORGE_ROOT}/bin:${PATH}"
hash -r

# conda.sh creates a shell function that uses $CONDA_EXE internally
export CONDA_EXE="${MINIFORGE_ROOT}/condabin/conda"

# (optional, but helpful if present)
if [[ -x "${MINIFORGE_ROOT}/bin/mamba" ]]; then
  export MAMBA_EXE="${MINIFORGE_ROOT}/bin/mamba"
fi

# Diagnostic prints to verify environment
echo "[controller] type conda:"
type conda || true
echo "[controller] command -v conda: $(command -v conda || true)"
echo "[controller] CONDA_EXE=${CONDA_EXE}"
"${CONDA_EXE}" --version || true
echo ""
echo "[controller] CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-NA}"
echo "[controller] Conda envs available:"
conda info --envs | head -n 30 || true
echo "[controller] Snakemake location:"
which snakemake || true
echo "[controller] Snakemake version:"
snakemake --version || true

command -v snakemake >/dev/null 2>&1 || die "snakemake not found in env '$CONDA_ENV'"

cd "$REPO_ROOT"

SMK_LOG="${LOGDIR}/snakemake_${SLURM_JOB_ID}.log"

echo "=== Snakemake controller ==="
echo "JobID:      $SLURM_JOB_ID"
echo "Host:       $(hostname)"
echo "Partition:  ${SLURM_JOB_PARTITION:-NA}"
echo "QOS:        ${SLURM_JOB_QOS:-NA}"
echo "Walltime:   ${SLURM_TIMELIMIT:-NA}"
echo "Repo:       $REPO_ROOT"
echo "Snakefile:  $SNAKEFILE"
echo "Profile:    $WF_PROFILE"
echo "Env:        $CONDA_ENV"
echo "Cores:      $CORES"
echo "Jobs:       $JOBS"
echo "USE_CONDA:  $USE_CONDA"
echo "Targets:    ${TARGETS:-<default>}"
echo "Log:        $SMK_LOG"
echo "============================"

SMK_ARGS=(
  --snakefile "$SNAKEFILE"
  --profile "$WF_PROFILE"
  --cores "$CORES"
  --jobs "$JOBS"
  --latency-wait "$LATENCY_WAIT"
  --restart-times "$RESTART_TIMES"
  --rerun-incomplete
  --keep-going
  --printshellcmds
  --reason
)

if [[ "$USE_CONDA" == "1" ]]; then
  SMK_ARGS+=( --use-conda --conda-prefix "$CONDA_PREFIX" )
fi

if [[ -n "$TARGETS" ]]; then
  # shellcheck disable=SC2206
  SMK_ARGS+=( $TARGETS )
fi

# Unlock working directory if locked (safe to run even if no lock exists)
echo "[controller] Checking for Snakemake lock..."
snakemake --snakefile "$SNAKEFILE" --profile "$WF_PROFILE" --unlock 2>&1 | grep -v "Nothing to unlock" || true
echo ""

# Run Snakemake
snakemake "${SMK_ARGS[@]}" 2>&1 | tee "$SMK_LOG"

