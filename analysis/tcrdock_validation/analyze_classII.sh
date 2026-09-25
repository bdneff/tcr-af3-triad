#!/bin/bash
# Class II AUC comparison: TCRdock vs AF3 PTI-PAE, both antigen-centric and TCR-centric.
# Run after finalize_classII.sh produces user_output_classII_w_pae.tsv.
set -euo pipefail
WORKDIR="${WORKDIR:-/scratch/bneff/tcrdock_validation}"
cd "$WORKDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$WORKDIR/env"

python analyze_tcrdock_vs_af3_v2.py \
    --tcrdock_tsv  user_output_classII_w_pae.tsv \
    --meta_xlsx    tcrdock_input_class_II_full.xlsx \
    --af3_val_xlsx supplementaryTable3_updated.xlsx \
    --mhc_class    II \
    --prefix       classII
