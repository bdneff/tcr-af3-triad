#!/usr/bin/env bash
# =============================================================================
# af3_array.sbatch — fold the post-AF3-cutoff class II triads with AlphaFold 3.
#
# Submit FROM THE REPO ROOT, one array task per JSON:
#     cd /scratch/$USER/tcr-af3-triad
#     sbatch --array=1-9 analysis/classII_supp3c/af3_array.sbatch
#
# Purpose: reviewer round 2, minor point 1 — extend Supp Fig 3c's PTI-PAE vs CDR-RMSD
# benchmark with class II structures released after the AF3 training cutoff.
#
# RUN CONFIGURATION IS NOT A CHOICE HERE — it is fixed by the manuscript so these
# predictions sit alongside the published ones. Methods, "AlphaFold 3 configuration":
# MSAs and templates ON, no custom templates or alignments, seeds 1-5, ONE diffusion
# sample, default 10 recycles.
#   * seeds live in the JSON (staged by scripts/01_stage_af3.py)
#   * --num_diffusion_samples 1  (NOT 5 — the paper used 1)
#   * --norun_data_pipeline is ABSENT on purpose. The JSONs carry no unpairedMsa/
#     pairedMsa/templates, so AF3 builds real MSAs and searches templates. The JSON and
#     the flag are a matched pair; adding the flag would fold from a bare sequence and
#     the numbers would not be comparable to the published triads.
#
# Resume: each prediction stamps _DONE only once a *_model.cif exists, so re-submitting
# skips finished work and a job that exits 0 without a structure is caught rather than
# silently counted. Retry specific ones with --array=3,7
# =============================================================================
#SBATCH -J af3_classII
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o analysis/classII_supp3c/logs/slurm_%A_%a.out
#SBATCH -e analysis/classII_supp3c/logs/slurm_%A_%a.err

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
SUB="analysis/classII_supp3c"
if   [[ -d "$ROOT/$SUB/af3_inputs" ]]; then BASE="$ROOT/$SUB"
elif [[ -d "$ROOT/af3_inputs"      ]]; then BASE="$ROOT"
else
    echo "FATAL: cannot find $SUB/af3_inputs from $ROOT"
    echo "  submit from the repo root:  sbatch --array=1-9 $SUB/af3_array.sbatch"
    exit 2
fi

IN_DIR="$BASE/af3_inputs"
OUT_ROOT="$BASE/af3_outputs"
mkdir -p "$OUT_ROOT" "$BASE/logs"

SIF=/tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif
MODELS=/ref_genomes/alphafold/alphafold3/models
DBS=/ref_genomes/alphafold/alphafold3/

N_JSON=$(ls -1 "$IN_DIR"/*.json 2>/dev/null | wc -l | tr -d ' ')
TASK="${SLURM_ARRAY_TASK_ID:?submit with --array=1-$N_JSON}"
JSON="$(ls -1 "$IN_DIR"/*.json | sed -n "${TASK}p")"
[[ -n "${JSON:-}" ]] || { echo "No JSON at index $TASK in $IN_DIR (have $N_JSON)"; exit 2; }

PRED_ID="$(basename "$JSON" .json)"
OUT_DIR="$OUT_ROOT/$(echo "$PRED_ID" | tr '[:upper:]' '[:lower:]')"
DONE="$OUT_DIR/_DONE"

if [[ -f "$DONE" ]]; then
    echo "Already done: $PRED_ID — nothing to do."
    exit 0
fi

echo "=== $PRED_ID | task $TASK/$N_JSON | $(date) | $(hostname) ==="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo unknown)"
echo "JSON: $JSON"

module load singularity

singularity exec --nv -B /home,/scratch,/tgen_labs,/ref_genomes "$SIF" \
    python /app/alphafold/run_alphafold.py \
        --json_path "$JSON" \
        --model_dir "$MODELS" \
        --db_dir    "$DBS" \
        --output_dir "$OUT_ROOT" \
        --num_diffusion_samples 1

n_cif=$(find "$OUT_DIR" -name '*_model.cif' 2>/dev/null | wc -l | tr -d ' ')
if [[ "$n_cif" -eq 0 ]]; then
    echo "FATAL: $PRED_ID produced no *_model.cif — not marking done."
    exit 1
fi

touch "$DONE"
echo "=== $PRED_ID done: $n_cif model(s) | $(date) ==="
