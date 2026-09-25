#!/usr/bin/env bash
# 04_replot_figure1.sh
#
# Regenerate Figure 1 (panels b, c, d) using the expanded RMSD parquets.
# Run this AFTER the pipeline (step 03) completes on Gemini and the staged
# parquets have been copied back to a local checkout, OR run on Gemini
# directly pointing at the production paths.

set -euo pipefail
cd "$(dirname "$0")/.."

: "${REPO_ROOT:?Please set REPO_ROOT to the tcrtrifold-experiments-main checkout}"

AF3_RMSD="${AF3_RMSD:-$REPO_ROOT/data/pdb/triad/staged/pdb_triad.af3_rmsd.parquet}"
BOLTZ_RMSD="${BOLTZ_RMSD:-$REPO_ROOT/data/pdb/triad/staged/pdb_triad.boltz_rmsd.parquet}"
OUT_DIR="${OUT_DIR:-results/fig1_v2}"

echo "==> [04] Regenerating Figure 1 panels b/c/d"
echo "    AF3 RMSD:   $AF3_RMSD"
echo "    Boltz RMSD: $BOLTZ_RMSD"
echo "    Out dir:    $OUT_DIR"

python scripts/04_replot_figure1.py \
    --af3_rmsd   "$AF3_RMSD" \
    --boltz_rmsd "$BOLTZ_RMSD" \
    --out_dir    "$OUT_DIR"

echo "==> [04] done. Open $OUT_DIR/ for results."
