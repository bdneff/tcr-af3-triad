#!/bin/bash
# Run this after `sbatch submit_tcrdock_array.sh` completes (login node is fine,
# no GPU needed). It:
#   1) concatenates per-chunk output_w_pae.tsv files into one combined TSV
#   2) reports any chunks that didn't finish (so you can re-submit them)
#   3) confirms the row count matches the input

set -euo pipefail
WORKDIR="${WORKDIR:-$PWD}"
INPUT_TSV="${INPUT_TSV:-tcrdock_input_class_I.tsv}"
cd "$WORKDIR"

python3 - <<EOF
import pandas as pd, glob, os, sys

expected = len(pd.read_csv("$INPUT_TSV", sep="\t"))

# Find which chunk dirs exist vs which have finished outputs
chunk_dirs = sorted(glob.glob("runs/chunk_*"))
done   = [d for d in chunk_dirs if os.path.exists(f"{d}/output_w_pae.tsv")]
missing= [d for d in chunk_dirs if d not in done]

print(f"Expected: {expected} rows from $INPUT_TSV")
print(f"Found:    {len(done)}/{len(chunk_dirs)} chunk dirs have output_w_pae.tsv")

if missing:
    print("\nChunks missing output_w_pae.tsv (re-submit these IDs):")
    for d in missing:
        cid = int(d.rsplit("_",1)[1])
        print(f"  chunk {cid}: {d}/")
    print("\nTo retry only these chunks:")
    ids = ",".join(str(int(d.rsplit('_',1)[1])) for d in missing)
    print(f"  sbatch --array={ids} submit_tcrdock_array.sh")

if not done:
    sys.exit("\nNothing to concatenate.")

dfs = [pd.read_csv(f"{d}/output_w_pae.tsv", sep="\t") for d in done]
combined = pd.concat(dfs, ignore_index=True)
combined.to_csv("user_output_w_pae.tsv", sep="\t", index=False)
print(f"\nConcatenated {len(done)} chunks -> user_output_w_pae.tsv ({len(combined)} rows)")

if len(combined) != expected:
    print(f"  WARN: {expected - len(combined)} rows are missing from the combined output")
EOF
