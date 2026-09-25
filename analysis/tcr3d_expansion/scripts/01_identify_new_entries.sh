#!/usr/bin/env bash
# 01_identify_new_entries.sh
#
# Identify the TCR3d entries that need to be added to the Figure 1 pipeline.
# This must run on a machine with network access to RCSB (data.rcsb.org and
# www.rcsb.org), but does NOT need to be on Gemini — running locally is fine
# and saves intermediate calls in data/rcsb_cache/.

set -euo pipefail
cd "$(dirname "$0")/.."

: "${REPO_ROOT:?Please set REPO_ROOT to the tcrtrifold-experiments-main checkout (e.g. ./tcrtrifold-experiments-main)}"
: "${TCR3D_CLASSI:=data/tcr3d_classI_complexes.csv}"
: "${TCR3D_CLASSII:=data/tcr3d_classII_complexes.csv}"

echo "==> [01] Identifying new TCR3d entries vs $REPO_ROOT"
python scripts/01_identify_new_entries.py \
    --tcr3d_classI  "$TCR3D_CLASSI" \
    --tcr3d_classII "$TCR3D_CLASSII" \
    --repo_root     "$REPO_ROOT" \
    "$@"

echo "==> [01] done."
