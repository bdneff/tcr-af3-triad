# Class II TCRdock comparison — quick run instructions

You should already have everything set up from the class I run (the conda env
at `/scratch/bneff/tcrdock_validation/env`, TCRdock cloned, AF params downloaded).
This adds class II on top of that without touching any class I files.

All commands are run from `/scratch/bneff/tcrdock_validation/`.

## Files in this bundle
- `build_tcrdock_input_classII.py` — builds the input TSV from supplementaryTable3
- `submit_classII_array.sh`        — SLURM array job (20 chunks)
- `finalize_classII.sh`            — combines chunk outputs
- `analyze_tcrdock_vs_af3_v2.py`   — analyzer with antigen + TCR-centric AUCs
- `analyze_classII.sh`             — wrapper that runs the analyzer with class II args
- `analyze_classI_v2.sh`           — same but for class I (gives you TCR-centric AUCs without rerunning TCRdock)

## Steps

```bash
cd /scratch/bneff/tcrdock_validation
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/bneff/tcrdock_validation/env

# --- 1. Install ANARCI + hmmer one time (needed by the input builder)
pip install anarci
conda install -y -c bioconda hmmer

# --- 2. Build class II input TSV from supplementaryTable3_updated.xlsx
python build_tcrdock_input_classII.py
#   Expected: ~1433 rows, 8 antigens, "0 issues" on gene validation.
#   If gene validation reports unrecognized genes, STOP and report back.

# --- 3. Submit the 20-way job array
sbatch submit_classII_array.sh
squeue -u bneff

# --- 4. Wait. ~2-3 hours wall clock once chunks start running (~72 triads
#        per chunk, ~1.5-2 min per triad on A100, 20-way parallelism).
#        Walk away.

# --- 5. Once squeue is empty, combine chunks then analyze
bash finalize_classII.sh
bash analyze_classII.sh
```

Step 5 produces:
- `classII_per_antigen_AUC.csv` — 8 antigens × {AF3, TCRdock}
- `classII_per_TCR_AUC.csv`     — 205 cognate TCRs × {AF3, TCRdock}
- Paired Wilcoxon p-values for both views, printed to stdout

## Bonus: TCR-centric AUCs for class I (no new compute)

The class I TCRdock run is already done. To get the TCR-centric numbers too:

```bash
bash analyze_classI_v2.sh
```

Produces `classI_per_antigen_AUC.csv` and `classI_per_TCR_AUC.csv` so you
have all 4 panels (class × {antigen-centric, TCR-centric}) for Figure 4.

## If a chunk fails
- `finalize_classII.sh` will tell you which chunks are missing
- It also prints the exact `sbatch --array=X,Y` retry command
- Rerun that, wait, run `finalize_classII.sh` again

## If you hit GPU memory errors
Open `submit_classII_array.sh`, add this line right after the `conda activate` line:
```
export XLA_PYTHON_CLIENT_MEM_FRACTION=0.6
```
Class II peptides can be longer (up to 15-mer) and MHCs have 2 polymorphic
chains, so memory pressure is higher than class I. Only add this if needed.
