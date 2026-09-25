#!/bin/bash
# Rerun the class I comparison with v2 analyzer to ALSO get TCR-centric AUCs.
# (Antigen-centric numbers will match what you already saw; TCR-centric is new.)
# No re-running of TCRdock needed - reuses user_output_w_pae.tsv from class I run.
set -euo pipefail
WORKDIR="${WORKDIR:-/scratch/bneff/tcrdock_validation}"
cd "$WORKDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$WORKDIR/env"

python analyze_tcrdock_vs_af3_v2.py \
    --tcrdock_tsv  user_output_w_pae.tsv \
    --meta_xlsx    tcrdock_input_class_I_full.xlsx \
    --af3_val_xlsx supplementaryTable3_updated.xlsx \
    --mhc_class    I \
    --prefix       classI
