#!/usr/bin/env bash
# =============================================================================
# jarvis_hpc_mac_commit.sh
#
# This script automates the synchronisation of pipeline source files from a
# remote HPC cluster (jarvis) to a local Git repository on macOS, followed
# by an optional commit and push to the remote origin.
#
# It performs three sequential operations:
#   1. Rsync  — transfers pipeline source files from the HPC working directory
#               to the local repository, applying exclude rules to suppress
#               large, generated, or environment-specific files.
#   2. Commit — stages any changes detected under pipeline/ and records them
#               as a Git commit with a configurable message.
#   3. Push   — propagates the commit to the remote Git origin (typically
#               GitHub), keeping the repository in sync with HPC development.
#
# All three stages are configurable via environment variables, allowing the
# script to be run in dry-run or no-push mode without modifying source code.
# =============================================================================

set -euo pipefail
# set -e        : exits immediately if any command returns a non-zero status.
# set -u        : treats references to undefined variables as errors.
# set -o pipefail: a pipeline fails if any component command fails, not just
#                  the last one — critical for chains involving rsync and git.


# =============================================================================
# Section 1 — Configuration
#
# All paths and behavioural flags are resolved from environment variables,
# with sensible defaults applied via the ${VAR:-default} pattern. This allows
# the script to be invoked without modification across different user
# environments or CI contexts.
# =============================================================================
REMOTE_HOST="${REMOTE_HOST:-jarvis}"                                         # SSH alias for the HPC login node.
REMOTE_DIR="${REMOTE_DIR:-/work/ugbogu/pipeline/}"                           # Absolute path to the pipeline directory on the HPC.
LOCAL_REPO="${LOCAL_REPO:-$HOME/classification_of_dsmz_celllines}"           # Root of the local Git repository.
LOCAL_PIPELINE_DIR="${LOCAL_PIPELINE_DIR:-$LOCAL_REPO/pipeline}"             # Subdirectory within the repo that mirrors the HPC pipeline.
DEFAULT_MSG="${DEFAULT_MSG:-Sync pipeline changes from HPC}"                 # Fallback commit message when none is supplied as a positional argument.

# Behavioural flags; both default to disabled (0).
DRY_RUN="${DRY_RUN:-0}"   # When set to 1, rsync previews changes without writing any files to disk.
NO_PUSH="${NO_PUSH:-0}"   # When set to 1, the script commits locally but skips git push.


# =============================================================================
# Section 2 — Safety Checks
#
# Before any file transfer or Git operation is attempted, the script confirms
# that the local repository is a valid Git working tree. This prevents
# accidental rsync or commits into an unversioned directory.
# =============================================================================
if [[ ! -d "$LOCAL_REPO/.git" ]]; then
    echo "ERROR: Not a git repo: $LOCAL_REPO"
    exit 1
fi

# The local pipeline directory is created if it does not yet exist, ensuring
# rsync has a valid destination on the first run of the script.
mkdir -p "$LOCAL_PIPELINE_DIR"


# =============================================================================
# Section 3 — Rsync Exclude Rules
#
# A persistent exclude file governs which files and directories are suppressed
# during transfer. If the file does not exist it is written with a default
# rule set covering the following categories of non-source artefacts:
#
#   • .snakemake/      — Snakemake's runtime cache; environment-specific and
#                        potentially very large.
#   • results/, logs/, count_data/ — Generated outputs, runtime logs, and
#                        count matrices that should not be mirrored into Git.
#   • *_family.soft.gz  — GEO family soft archives.
#   • data/reference/gencode_v44/ and pipeline/data/reference/gencode_v44/
#     plus GRCh38.primary_assembly.genome.fa and generic FASTA assets
#     (*.fa, *.fasta, *.fna).
#   • data/reference/star_*/ and pipeline/data/reference/star_*/
#                        STAR index/reference folders.
#   • data/nbl/target_nbl/manifests/*.json,
#     data/nbl/target_nbl/metadata/payload*.json and their
#     pipeline/data/... equivalents
#                        target-nbl payload and manifest JSON files.
#   • envs/.conda/     — Generated Conda payloads; YAML environment specs and
#                        setup scripts elsewhere under envs/ remain tracked.
#   • *.fastq.gz, *.fq.gz, *.bam, *.bai, *.sra — Large sequencing/alignment
#                        artefacts that should not be synced/committed.
#   • *.rds / *.Rds    — Serialised R objects; binary files that change
#                        frequently and are reproducible from source.
#   • *.gtf / *.gtf.gz — Reference genome annotation files; large static
#                        resources that belong in dedicated data stores.
#   • Reports / artifacts — PDFs, HTML, Parquet, and common image outputs.
#   • Python / R junk  — Bytecode caches (__pycache__/, *.pyc), R session
#                        files (.Rhistory, .RData), and backup files (*.bak).
#   • !*.txt / !*.out / !*.err — Explicitly kept for mirrored text outputs and
#                        scheduler logs.
#
# A file-based exclude list (--exclude-from) is preferred over inline rsync
# flags because it is auditable, versionable alongside the repository, and
# straightforward to extend without modifying this script.
# =============================================================================
RSYNC_EXCLUDE_FILE="$LOCAL_REPO/.rsyncignore"
if [[ ! -f "$RSYNC_EXCLUDE_FILE" ]]; then
    cat >"$RSYNC_EXCLUDE_FILE" <<'EOF'
# Snakemake runtime (huge)
.snakemake/
**/.snakemake/

# Big / generated directories inside pipeline
results/
logs/
count_data/
*_family.soft.gz
data/reference/gencode_v44/
data/reference/**/GRCh38.primary_assembly.genome.fa
data/reference/star_*/
pipeline/data/reference/gencode_v44/
pipeline/data/reference/**/GRCh38.primary_assembly.genome.fa
pipeline/data/reference/star_*/
*.fa
*.fasta
*.fna
data/nbl/target_nbl/manifests/*.json
data/nbl/target_nbl/metadata/payload*.json
pipeline/data/nbl/target_nbl/manifests/*.json
pipeline/data/nbl/target_nbl/metadata/payload*.json


# Generated conda cache/env internals — exclude these, but keep env spec files
envs/.conda/
**/envs/.conda/

# Large sequencing/alignment artifacts
*.fastq.gz
*.fq.gz
*.bam
*.bai
*.sra

# Python / R junk
__pycache__/
*.pyc
.Rhistory
.RData
*.rds
*.Rds
*.bak

# Reports / artifacts
*.pdf
*.html
*.parquet
*.png
*.jpg
*.jpeg

# Annotations (too big for GitHub)
*.gtf
*.gtf.gz

# Keep text and scheduler logs in sync
!*.txt
!*.out
!*.err
!*.tsv
!*.csv
EOF
    echo "Created default rsync exclude file at $RSYNC_EXCLUDE_FILE"
fi


# =============================================================================
# Section 4 — Rsync Transfer
#
# Pipeline source files are synchronised from the HPC to the local repository
# using rsync over SSH. The flags applied are:
#
#   -a (archive)         : preserves permissions, timestamps, and symlinks,
#                          and recursively transfers all subdirectories.
#   -v (verbose)         : logs each transferred file, supporting auditability.
#   --delete             : removes local files that no longer exist on the HPC,
#                          keeping the local copy an accurate mirror of source.
#   --delete-excluded    : also removes locally cached files that match the
#                          exclude rules, preventing stale artefact accumulation.
#   --prune-empty-dirs   : suppresses creation of empty directories that arise
#                          after exclusions are applied.
#   --exclude-from       : loads the exclusion pattern file defined in Section 3.
# =============================================================================
echo "==> Rsync from ${REMOTE_HOST}:${REMOTE_DIR} -> ${LOCAL_PIPELINE_DIR}/"
RSYNC_ARGS=(
    -av
    --delete
    --delete-excluded
    --prune-empty-dirs
    --exclude-from="$RSYNC_EXCLUDE_FILE"
)

if [[ "$DRY_RUN" == "1" ]]; then
    # In dry-run mode --dry-run is appended so rsync reports what would change
    # without modifying any files on disk.
    RSYNC_ARGS+=(--dry-run)
    echo "==> DRY_RUN enabled (no files will be changed)"
fi

rsync "${RSYNC_ARGS[@]}" \
    "${REMOTE_HOST}:${REMOTE_DIR}" \
    "${LOCAL_PIPELINE_DIR}/"

# Hard post-transfer safety check: if any Snakemake runtime directory was
# transferred despite the exclude rules (e.g. due to a misconfigured pattern),
# the script aborts before any Git operation. Committing nested .snakemake/
# directories would introduce environment-specific cache files into the repo.
if [[ -n "$(find "$LOCAL_PIPELINE_DIR" -type d -name .snakemake -print -quit)" ]]; then
    echo "ERROR: One or more .snakemake directories exist after rsync. Refusing to commit."
    echo "Fix your excludes, then remove them manually."
    exit 2
fi


# =============================================================================
# Section 5 — Git Commit
#
# This section stages all changes under the pipeline/ subdirectory and records
# them as a Git commit. Two idempotency guards prevent unnecessary commits:
#
#   Guard 1 — git status --porcelain returns an empty string when the working
#             tree is clean; the script exits without error in this case.
#   Guard 2 — git diff --cached --quiet returns true when staging produced no
#             changes (e.g. all changed files were suppressed by .gitignore);
#             the script again exits cleanly.
#
# The commit message is taken from the first positional argument ($1) when
# supplied, falling back to DEFAULT_MSG for unattended or scheduled invocations.
# =============================================================================
cd "$LOCAL_REPO"

echo "==> Git status (porcelain):"
git status --porcelain

# Guard 1: exit cleanly if the working tree contains no changes.
if [[ -z "$(git status --porcelain)" ]]; then
    echo "==> No changes to commit."
    exit 0
fi

echo "==> Staging changes under pipeline/"
git add -A pipeline

# Guard 2: exit cleanly if .gitignore suppressed all changed files after staging.
if git diff --cached --quiet; then
    echo "==> Nothing staged after git add (maybe everything ignored)."
    exit 0
fi

MSG="${1:-$DEFAULT_MSG}"
echo "==> Committing: $MSG"
git commit -m "$MSG"


# =============================================================================
# Section 6 — Git Push
#
# The local commit is pushed to the configured remote (origin/main). This
# step is skipped when NO_PUSH=1, which is useful for testing the rsync and
# commit stages in isolation without affecting the remote repository.
# =============================================================================
if [[ "$NO_PUSH" == "1" ]]; then
    echo "==> NO_PUSH enabled; skipping git push."
    exit 0
fi

echo "==> Pushing to origin/main"
git push
echo "==> Done."
