#!/usr/bin/env bash
# 02_apply_to_repo.sh
#
# Patch the tcrtrifold-experiments repo with the new entries discovered by
# step 01. Run with --dry_run first to see what will change before committing.

set -euo pipefail
cd "$(dirname "$0")/.."

: "${REPO_ROOT:?Please set REPO_ROOT to the tcrtrifold-experiments-main checkout}"

echo "==> [02] Patching repo at $REPO_ROOT"
python scripts/02_apply_to_repo.py \
    --repo_root "$REPO_ROOT" \
    "$@"

echo "==> [02] done."
