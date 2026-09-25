#!/bin/bash
#SBATCH --job-name=tcrdock_cII
#SBATCH --output=logs_classII/tcrdock_%A_%a.log
#SBATCH --error=logs_classII/tcrdock_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu-a100
#SBATCH --array=0-19         # 20 chunks of ~72 triads each
#
# Class II TCRdock job array. Run with:
#   sbatch submit_classII_array.sh

set -euo pipefail

WORKDIR="${WORKDIR:-/scratch/bneff/tcrdock_validation}"
INPUT_TSV="${INPUT_TSV:-tcrdock_input_class_II.tsv}"
NCHUNKS="${NCHUNKS:-20}"
CHUNK_ID="${SLURM_ARRAY_TASK_ID:-0}"

cd "$WORKDIR"
mkdir -p logs_classII chunks_classII runs_classII

# Activate env BEFORE any python (needed for pandas in the chunker)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$WORKDIR/env"

# Split input TSV into chunks (idempotent across array tasks)
python3 - <<EOF
import pandas as pd, math
df = pd.read_csv("$INPUT_TSV", sep="\t")
n_per = math.ceil(len(df) / $NCHUNKS)
for i in range($NCHUNKS):
    chunk = df.iloc[i*n_per:(i+1)*n_per]
    if len(chunk):
        chunk.to_csv(f"chunks_classII/chunk_{i:03d}.tsv", sep="\t", index=False)
print(f"Split {len(df)} rows into {$NCHUNKS} chunks of up to {n_per} rows")
EOF

CHUNK_TSV="chunks_classII/chunk_$(printf '%03d' $CHUNK_ID).tsv"
[ ! -f "$CHUNK_TSV" ] && { echo "No chunk $CHUNK_ID; exiting"; exit 0; }

echo "Chunk $CHUNK_ID: $(wc -l < $CHUNK_TSV) lines from $CHUNK_TSV"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo unknown)"

# Per-chunk working dir to keep outputs from colliding
CHUNK_DIR="$WORKDIR/runs_classII/chunk_$(printf '%03d' $CHUNK_ID)"
mkdir -p "$CHUNK_DIR"
cp "$WORKDIR/$CHUNK_TSV" "$CHUNK_DIR/input.tsv"

cd "$WORKDIR/TCRdock"

# 1. Setup (with --benchmark to exclude pre-AF3-cutoff templates, --new_docking for 3x speedup)
python setup_for_alphafold.py \
    --targets_tsvfile "$CHUNK_DIR/input.tsv" \
    --output_dir      "$CHUNK_DIR/setup" \
    --new_docking \
    --benchmark

# 2. Predict
python run_prediction.py --verbose \
    --targets        "$CHUNK_DIR/setup/targets.tsv" \
    --outfile_prefix "$CHUNK_DIR/output" \
    --model_names    model_2_ptm_ft4 \
    --data_dir       "$WORKDIR/alphafold_params/" \
    --model_params_files "$WORKDIR/alphafold_params/params/tcrpmhc_run4_af_mhc_params_891.pkl"

# 3. Extract pmhc_tcr_pae
python add_pmhc_tcr_pae_to_tsvfile.py \
    --infile  "$CHUNK_DIR/output_final.tsv" \
    --outfile "$CHUNK_DIR/output_w_pae.tsv"

echo "=== Class II chunk $CHUNK_ID done ==="
