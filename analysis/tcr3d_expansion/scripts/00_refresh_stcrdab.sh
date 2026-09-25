#!/usr/bin/env bash
# 00_refresh_stcrdab.sh
#
# Fetch fresh per-PDB summary files from STCRDab and concatenate into a
# replacement for db_summary.dat. This is recommended before running step 01
# because the on-disk db_summary.dat is a snapshot that misses many recent
# entries (~140 PDBs as of mid-2026).
#
# Reduces the "needs classifier" set in step 01 from ~46 to a handful or zero.

set -euo pipefail
cd "$(dirname "$0")/.."

: "${REPO_ROOT:?Please set REPO_ROOT to the tcrtrifold-experiments-main checkout}"

OLD_SUMMARY="$REPO_ROOT/data/pdb/raw/db_summary.dat"

echo "==> [00] Fetching fresh STCRDab per-PDB summaries"
python scripts/00_refresh_stcrdab.py \
    --tcr3d_classI  data/tcr3d_classI_complexes.csv \
    --tcr3d_classII data/tcr3d_classII_complexes.csv \
    --also_keep_old "$OLD_SUMMARY" \
    --out           data/db_summary_fresh.dat \
    "$@"

echo "==> [00] done."
echo
echo "Next: pass the fresh summary explicitly to step 01:"
echo "    python scripts/01_identify_new_entries.py \\"
echo "        --tcr3d_classI  data/tcr3d_classI_complexes.csv \\"
echo "        --tcr3d_classII data/tcr3d_classII_complexes.csv \\"
echo "        --repo_root     \"\$REPO_ROOT\" \\"
echo "        --stcrdab_summary data/db_summary_fresh.dat"
