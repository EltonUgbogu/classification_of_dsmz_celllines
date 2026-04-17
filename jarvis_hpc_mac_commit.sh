#!/usr/bin/env bash
set -euo pipefail

REMOTE_HOST="jarvis"
REMOTE_DIR="/work/ugbogu/pipeline/"
LOCAL_DIR="$HOME/classification_of_dsmz_celllines"

cd "$LOCAL_DIR"

echo "==> Rsync from ${REMOTE_HOST}:${REMOTE_DIR} -> ${LOCAL_DIR}/"

RSYNC_ARGS=(
  -avh
  --progress
  --exclude ".git/"
)

if [[ -f ".rsyncignore" ]]; then
  RSYNC_ARGS+=(--exclude-from ".rsyncignore")
fi

rsync "${RSYNC_ARGS[@]}" "${REMOTE_HOST}:${REMOTE_DIR}" "${LOCAL_DIR}/"

echo
echo "==> Git status"
git status

if ! git diff --quiet || ! git diff --cached --quiet || [[ -n "$(git ls-files --others --exclude-standard)" ]]; then
  COMMIT_MSG="${1:-Sync pipeline changes from HPC}"
  echo
  echo "==> Staging changes"
  git add -A

  echo
  echo "==> Committing"
  git commit -m "$COMMIT_MSG" || echo "Nothing to commit."

  echo
  echo "==> Pushing to origin/main"
  git push origin main
else
  echo
  echo "No changes detected after rsync."
fi
