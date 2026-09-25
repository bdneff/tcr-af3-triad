#!/usr/bin/env bash
# 03_analyze_results.sh
#
# Step 3 of the cross-reactivity controls pipeline.
#
# Run this AFTER AF3 inference + feature extraction have finished on Gemini.
# It joins the AF3 PTI-PAE results with the class labels and BLOSUM62
# distances from step 02, writes a clean CSV, and produces the scatter plot.
#
# Inputs:   $1 = path to the AF3 feature-extraction output (.parquet or .csv)
#                e.g. data/cross_reactivity_controls.conf_af3.parquet
#           data/cross_reactivity_controls_pipeline.csv  (from step 02)
# Outputs:  data/cross_reactivity_controls_with_PTIPAE.csv
#           data/cross_reactivity_controls_plot.png
#
# Usage:
#   bash scripts/03_analyze_results.sh path/to/conf_af3.parquet

set -euo pipefail
cd "$(dirname "$0")/.."

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <results.parquet|results.csv>" >&2
    exit 1
fi
RESULTS="$1"

echo "==> [03] Joining results with labels and plotting"
python scripts/03_analyze_results.py \
    --results "${RESULTS}" \
    --raw     data/cross_reactivity_controls_pipeline.csv \
    --out_csv data/cross_reactivity_controls_with_PTIPAE.csv \
    --out_png data/cross_reactivity_controls_plot.png

echo "==> [03] done"
