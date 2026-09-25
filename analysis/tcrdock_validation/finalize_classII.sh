#!/bin/bash
# Concatenate class II chunk outputs after submit_classII_array.sh finishes.
# Run from /scratch/bneff/tcrdock_validation/ on the login node.
set -euo pipefail
WORKDIR="${WORKDIR:-/scratch/bneff/tcrdock_validation}"
INPUT_TSV="${INPUT_TSV:-tcrdock_input_class_II.tsv}"
cd "$WORKDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$WORKDIR/env"

python3 - <<EOF
import pandas as pd, glob, os, sys
expected = len(pd.read_csv("$INPUT_TSV", sep="\t"))
dirs = sorted(glob.glob("runs_classII/chunk_*"))
done   = [d for d in dirs if os.path.exists(f"{d}/output_w_pae.tsv")]
missing= [d for d in dirs if d not in done]

print(f"Expected: {expected} rows from $INPUT_TSV")
print(f"Found:    {len(done)}/{len(dirs)} chunk dirs have output_w_pae.tsv")

if missing:
    ids = ",".join(str(int(d.rsplit('_',1)[1])) for d in missing)
    print(f"\nMissing chunks. To retry:")
    print(f"  sbatch --array={ids} submit_classII_array.sh")

if done:
    dfs = [pd.read_csv(f"{d}/output_w_pae.tsv", sep="\t") for d in done]
    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv("user_output_classII_w_pae.tsv", sep="\t", index=False)
    print(f"\nWrote user_output_classII_w_pae.tsv ({len(combined)} rows)")
    if len(combined) != expected:
        print(f"  WARN: {expected-len(combined)} rows still missing")
else:
    sys.exit("Nothing to concatenate.")
EOF
