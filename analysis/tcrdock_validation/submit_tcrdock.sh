#!/bin/bash
#SBATCH --job-name=tcrdock_classI
#SBATCH --output=tcrdock_%j.log
#SBATCH --error=tcrdock_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
# --- EDIT THESE for Gemini ---
##SBATCH --partition=gpu          # uncomment and set to your GPU partition
##SBATCH --account=your_account   # uncomment if your site requires it
##SBATCH --constraint=a100         # if your cluster lets you request A100 specifically
# -----------------------------

set -euo pipefail

# Resolve the working directory we set up earlier.
# If you put this script next to setup_env.sh's WORKDIR, this default works.
WORKDIR="${WORKDIR:-$PWD}"
INPUT_TSV="${INPUT_TSV:-tcrdock_input_class_I.tsv}"

cd "$WORKDIR"
echo "Working in: $WORKDIR"
echo "Input: $INPUT_TSV"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"

# --- Activate the env we built in setup_env.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ./env

cd TCRdock
# Make sure the input TSV is visible in TCRdock's working directory
ln -sf "../$INPUT_TSV" .

# --- 1. Generate per-target AlphaFold inputs.
#    --benchmark excludes pre-AF3-cutoff PDB templates (critical: prevents
#                TCRdock from using our test-set structures as templates)
#    --new_docking uses the fine-tuned/faster docking-geometry approach
python setup_for_alphafold.py \
    --targets_tsvfile "$INPUT_TSV" \
    --output_dir user_output \
    --new_docking \
    --benchmark

# --- 2. Run prediction on all 197 targets using the fine-tuned model
python run_prediction.py --verbose \
    --targets user_output/targets.tsv \
    --outfile_prefix user_output \
    --model_names model_2_ptm_ft4 \
    --data_dir   "$WORKDIR/alphafold_params/" \
    --model_params_files "$WORKDIR/alphafold_params/params/tcrpmhc_run4_af_mhc_params_891.pkl"

# --- 3. Extract the pmhc_tcr_pae metric we need
python add_pmhc_tcr_pae_to_tsvfile.py \
    --infile  user_output_final.tsv \
    --outfile user_output_w_pae.tsv

# --- 4. Copy the file we care about into the workdir root for easy retrieval
cp user_output_w_pae.tsv "$WORKDIR/"
echo
echo "=== Done. Results at: $WORKDIR/user_output_w_pae.tsv ==="
echo "=== Per-triad PDBs:    $WORKDIR/TCRdock/user_output_T*_model_*.pdb ==="
