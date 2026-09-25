#!/bin/bash
# One-shot env setup for TCRdock on an HPC node.
# Run this ONCE on a Gemini login node (no GPU needed for setup itself, but
# the env it creates is GPU-capable).
#
# Usage:
#   bash setup_env.sh /path/to/scratch/tcrdock_work
#
# After it finishes, the env lives at $WORKDIR/env and AF params at
# $WORKDIR/alphafold_params/. Both get reused by submit_tcrdock.sh.

set -euo pipefail

WORKDIR="${1:-$PWD/tcrdock_work}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"
echo "Working in: $WORKDIR"

# --- 1. Conda env (mirrors what the TCRdock Colab does, minus the Colab cruft)
# We try mamba first (much faster); falls back to conda if not present.
SOLVER=$(command -v mamba >/dev/null 2>&1 && echo mamba || echo conda)
echo "Using $SOLVER to build env"

# If conda is not on PATH at all, load the system module (adjust to your site).
# Common module names: anaconda3, miniconda3, conda. Edit if your cluster differs.
if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found - trying 'module load anaconda3'"
  module load anaconda3 || { echo "ERROR: install/load conda first"; exit 1; }
fi

# Accept Anaconda ToS preemptively (silences the July 2025 ToS prompt issue)
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r    2>/dev/null || true

if [ ! -d env ]; then
  $SOLVER create -y -p ./env -c conda-forge -c bioconda \
    python=3.10 openmm=8.0.0 pdbfixer cudatoolkit=11.8.0 cudnn=8.9 \
    biopython numpy pandas matplotlib scipy py3dmol \
    hmmer blast
else
  echo "env/ already exists - skipping conda env creation"
fi

# Activate and install pip-only bits
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ./env

# --- 2. TCRdock (its modified AlphaFold lives inside this repo)
if [ ! -d TCRdock ]; then
  git clone --depth 1 https://github.com/phbradley/TCRdock.git
fi
pip install -r TCRdock/requirements_colab_af232.txt
# (jax with CUDA wheels are listed inside that requirements file)

# --- 3. Set up BLAST databases that TCRdock needs for gene parsing
(cd TCRdock && python download_blast.py)

# --- 4. AlphaFold parameter files (~3 GB)
PARAMS_DIR="$WORKDIR/alphafold_params/params"
mkdir -p "$PARAMS_DIR"
for url in \
  https://www.dropbox.com/s/e3uz9mwxkmmv35z/params_model_2_ptm.npz \
  https://www.dropbox.com/s/jph8v1mfni1q4y8/tcrpmhc_run4_af_mhc_params_891.pkl
do
  fname=$(basename "$url" | sed 's/?.*//')
  if [ ! -f "$PARAMS_DIR/$fname" ]; then
    wget -O "$PARAMS_DIR/$fname" "$url?dl=1"
  fi
done

echo
echo "=== Setup done ==="
echo "Env:    $WORKDIR/env"
echo "Repo:   $WORKDIR/TCRdock"
echo "Params: $WORKDIR/alphafold_params/params/"
echo
echo "Next: copy your input TSV here and submit:"
echo "  cp tcrdock_input_class_I.tsv $WORKDIR/"
echo "  sbatch submit_tcrdock.sh"
