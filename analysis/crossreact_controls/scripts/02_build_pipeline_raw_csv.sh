#!/usr/bin/env bash
# 02_build_pipeline_raw_csv.sh
#
# Step 2 of the cross-reactivity controls pipeline.
#
# Converts the BLOSUM62-controlled CSV from step 01 into the column format
# the tcrtrifold-experiments pipeline expects (adds `cognate` and MD5 `name`
# columns, reorders to match cross_reactivity_raw.csv).
#
# Inputs:   data/cross_reactivity_controls_raw.csv
# Outputs:  data/cross_reactivity_controls_pipeline.csv
#
# After this completes, place the pipeline CSV at:
#   tcrtrifold-experiments/data/cross_reactivity_controls/raw/cross_reactivity_controls_raw.csv
# Then run the clean -> MSA -> inference -> feature-extract workflow exactly
# as for the prior 750-triad cross_reactivity job.
#
# Usage:
#   bash scripts/02_build_pipeline_raw_csv.sh

set -euo pipefail
cd "$(dirname "$0")/.."

echo "==> [02] Building pipeline-ready raw CSV"
python scripts/02_build_pipeline_raw_csv.py \
    --in  data/cross_reactivity_controls_raw.csv \
    --out data/cross_reactivity_controls_pipeline.csv

echo "==> [02] done"
echo
echo "Next: copy data/cross_reactivity_controls_pipeline.csv into"
echo "      tcrtrifold-experiments/data/cross_reactivity_controls/raw/ and run"
echo "      the standard 01_clean_and_inf -> 02_extract pipeline on Gemini."
