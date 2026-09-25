#!/bin/bash
#SBATCH --job-name=tcrdock_arr
#SBATCH --output=logs/tcrdock_%A_%a.log
#SBATCH --error=logs/tcrdock_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --array=0-9         # 10 chunks of ~20 triads each (197 total)
# --- EDIT THESE for Gemini ---
#SBATCH --partition=gpu-a100
##SBATCH --account=your_account
##SBATCH --constraint=a100
# -----------------------------
#
# Parallel TCRdock runner via SLURM job array.
#
# How to use:
#   1) Run setup_env.sh once (as before) to build env + grab AF params
#   2) Tune the --array=0-9 line above to your taste:
#        --array=0-9         -> 10 chunks, ~20 triads each, max 10 GPUs
#        --array=0-19        -> 20 chunks, ~10 triads each, max 20 GPUs
#        --array=0-9%4       -> 10 chunks but at most 4 running concurrently
#                              (use if you only have a small GPU quota)
#   3) Submit:
#        sbatch submit_tcrdock_array.sh
#   4) After all array tasks finish, run:
#        bash finalize.sh
#      to concatenate per-chunk outputs into one user_output_w_pae.tsv
#
# Re-running a subset (e.g., if chunks 3 and 7 failed):
#   sbatch --array=3,7 submit_tcrdock_array.sh

set -euo pipefail

WORKDIR="${WORKDIR:-$PWD}"
INPUT_TSV="${INPUT_TSV:-tcrdock_input_class_I.tsv}"
NCHUNKS="${NCHUNKS:-10}"
CHUNK_ID="${SLURM_ARRAY_TASK_ID:-0}"

cd "$WORKDIR"
mkdir -p logs chunks runs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$WORKDIR/env"

# --- 1. Split input TSV into chunks (idempotent; every task does this, but
#        it's a fraction of a second and ensures the chunk files exist)
python3 - <<EOF
import pandas as pd, math, os
df = pd.read_csv("$INPUT_TSV", sep="\t")
n_per = math.ceil(len(df) / $NCHUNKS)
for i in range($NCHUNKS):
    chunk = df.iloc[i*n_per:(i+1)*n_per]
    if len(chunk):
        chunk.to_csv(f"chunks/chunk_{i:03d}.tsv", sep="\t", index=False)
print(f"Split {len(df)} rows into {$NCHUNKS} chunks of up to {n_per} rows")
EOF

CHUNK_TSV="chunks/chunk_$(printf '%03d' $CHUNK_ID).tsv"
if [ ! -f "$CHUNK_TSV" ]; then
    echo "Chunk $CHUNK_ID has no rows (probably $NCHUNKS > #triads); exiting cleanly"
    exit 0
fi

echo "Chunk $CHUNK_ID: $(wc -l < $CHUNK_TSV) rows (incl. header) from $CHUNK_TSV"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"

# --- 2. Activate env
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$WORKDIR/env"

# --- 3. Per-chunk working dir to avoid any output collisions between tasks
CHUNK_DIR="$WORKDIR/runs/chunk_$(printf '%03d' $CHUNK_ID)"
mkdir -p "$CHUNK_DIR"
cp "$WORKDIR/$CHUNK_TSV" "$CHUNK_DIR/input.tsv"

cd "$WORKDIR/TCRdock"

# --- 4. Setup, predict, extract PAE — all scoped to this chunk's dir
python setup_for_alphafold.py \
    --targets_tsvfile "$CHUNK_DIR/input.tsv" \
    --output_dir      "$CHUNK_DIR/setup" \
    --new_docking \
    --benchmark

python run_prediction.py --verbose \
    --targets        "$CHUNK_DIR/setup/targets.tsv" \
    --outfile_prefix "$CHUNK_DIR/output" \
    --model_names    model_2_ptm_ft4 \
    --data_dir       "$WORKDIR/alphafold_params/" \
    --model_params_files "$WORKDIR/alphafold_params/params/tcrpmhc_run4_af_mhc_params_891.pkl"

python add_pmhc_tcr_pae_to_tsvfile.py \
    --infile  "$CHUNK_DIR/output_final.tsv" \
    --outfile "$CHUNK_DIR/output_w_pae.tsv"

echo "=== Chunk $CHUNK_ID done (results in $CHUNK_DIR/output_w_pae.tsv) ==="
