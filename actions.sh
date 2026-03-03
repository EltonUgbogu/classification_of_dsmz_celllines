#!/usr/bin/env bash


cat > ~/.gh_env <<'EOF'
export GH_OWNER="EltonUgbogu"
export GH_REPO="classification_of_dsmz_celllines"
export GH_TOKEN="$(< ~/.secrets/github_token)"
EOF

# when needed
source ~/.gh_env
