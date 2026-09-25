#!/usr/bin/env bash
# 01_generate_random_peptides.sh
#
# Step 1 of the cross-reactivity controls pipeline.
#
# Generates randomized "pink" peptide controls spanning BLOSUM62 distance
# from the original peptide, paired with each class's fixed TCR + MHC.
#
# Inputs:   data/I.csv, data/II.csv
# Outputs:  data/cross_reactivity_controls_raw.csv
#
# Configurable via environment variables (defaults shown):
#   N_PER_CLASS=400   number of randomized peptides per class
#   SEED=42           random seed for reproducibility
#
# Usage:
#   bash scripts/01_generate_random_peptides.sh
#   N_PER_CLASS=300 SEED=7 bash scripts/01_generate_random_peptides.sh

set -euo pipefail
cd "$(dirname "$0")/.."

N_PER_CLASS="${N_PER_CLASS:-5000}"
SEED="${SEED:-42}"

echo "==> [01] Generating randomized controls (n=${N_PER_CLASS} per class, seed=${SEED})"
python scripts/01_generate_random_peptides.py \
    --i_csv  data/I.csv  \
    --ii_csv data/II.csv \
    --out    data/cross_reactivity_controls_raw.csv \
    --n_per_class "${N_PER_CLASS}" \
    --seed "${SEED}"

echo "==> [01] done"
