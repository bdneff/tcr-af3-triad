#!/bin/bash
#
# 03_run_pipeline.sh
# ------------------
# Run the existing pdb.nf workflow (AF3 + Boltz-2 + RMSD) over the expanded
# set of triads. The pipeline will:
#   - Use the patched table_S1 (130 + 178 new = 308 pre-cutoff entries)
#   - Use the patched post_af3_cutoff_valid_pdbs (22 + 49 new = 71 entries)
#   - Plus the outlier_pdb entries for STCRDab-missing PDBs
#
# Run from the tcrtrifold-experiments repo root. The pipeline is -resume-able,
# so existing AF3 and Boltz-2 predictions are NOT re-run. Only the new triads
# will go through inference.
#
# Note: this assumes the patch has been applied via 02_apply_to_repo.py.
#       If you want to verify before submitting, do a dry run there first.

#SBATCH --job-name=pdb_tcr3d_expansion
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USER@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=10-00:00:00
#SBATCH --output=logs/pdb_tcr3d_expansion/slurm.%j.log

set -euo pipefail

# env vars
export NXF_LOG_FILE=logs/pdb_tcr3d_expansion/.nextflow.log
export NXF_CACHE_DIR=logs/pdb_tcr3d_expansion/.nextflow

mkdir -p logs/pdb_tcr3d_expansion

# --- AF3 + RMSD (uses existing pdb.nf) -----------------------------------
# This is the same call as scripts/pdb/01_clean_and_inf.sh, unchanged.
# With -resume, existing predictions are preserved; only new triads are run.
echo "==> [AF3 + RMSD] Running pdb.nf with expanded inputs..."
conda run -n nf-core --live-stream  nextflow run \
    ./workflows/pdb.nf \
    --input_replication "data/pdb/raw/table_S1_structure_benchmark_complexes.csv" \
    --input_stcr "data/pdb/raw/db_summary.dat" \
    -output-dir data/pdb \
    -profile gemini \
    -resume

# --- Boltz-2 on the new triads --------------------------------------------
# pdb.nf already runs Boltz-2 inline (via BOLTZ_FROM_TRIAD_PARQUET), so
# this is also handled by the -resume above. The 03_boltz_with_template.sh
# script in scripts/pdb/ is for the *template-conditioned* Boltz variant
# used in a different supplementary figure — NOT needed for Fig 1b/c/d.

echo "==> Done. Outputs are in data/pdb/triad/staged/"
echo "    Inspect:"
echo "      pdb_triad.af3_rmsd.parquet     (AF3 RMSDs, new + old entries)"
echo "      pdb_triad.boltz_rmsd.parquet   (Boltz-2 RMSDs)"
echo "      pdb_triad.af3_tcrdock.parquet  (docking geometries)"
echo
echo "Next: from a local checkout, run scripts/04_replot_figure1.py"
echo "      with these parquets to regenerate Figure 1b/c/d."
