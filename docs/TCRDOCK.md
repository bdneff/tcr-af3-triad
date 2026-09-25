# Running AF-TCRdock — the working recipe

AF-TCRdock (Bradley 2023) is the **baseline** the manuscript compares AF3 against
(Figure 4c, `pmhc_tcr_pae` vs AF3's PTI-PAE). Getting it to run on Gemini was the
fiddliest part of the paper, almost entirely because it is pinned to an
AlphaFold-2.3.2-era software stack that no longer matches a modern GPU node.

This file records **what actually worked**, recovered from the scripts that ran it.

> **Read this before `README_gemini.md`.** That file is the *generic pre-flight*
> version written before the runs happened — it has placeholder partitions
> (`##SBATCH --account=your_account`), tells you to go find your GPU partition
> with `sinfo`, and quotes **197 triads**. The real runs used `gpu-a100` and
> **145 class I / 1433 class II** triads. Where the two disagree, this file wins.

---

## 1. The thing that makes it hard: the CUDA pin

TCRdock runs a modified AlphaFold 2.3.2 on JAX. That JAX is built against
**CUDA 11.8**, and the environment therefore pins:

```
cudatoolkit=11.8.0
cudnn=8.9
python=3.10
openmm=8.0.0
```

This is the crux. A newer `cudatoolkit` gets you `Unable to initialize backend
'gpu'` or an XLA/cuDNN version error, and JAX silently falls back to CPU — at
which point a prediction that should take 1–2 minutes takes hours. The failure
is *slow*, not loud, so watch the per-target timings in the log rather than
waiting for a crash.

The constraint is one-directional: **the driver must be at least as new as the
toolkit**. `nvidia-smi` on a GPU node reports the maximum CUDA the driver
supports; pick a `cudatoolkit` less than or equal to that. On Gemini's A100
nodes, 11.8 was satisfied and did not need changing.

Two more setup snags worth knowing:

- **Anaconda ToS prompt.** Since July 2025 `conda create` can block on a
  terms-of-service acceptance. `setup_env.sh` pre-accepts for the `main` and
  `r` channels; without it, setup hangs waiting on stdin inside a batch job.
- **`conda: command not found`.** Setup falls back to `module load anaconda3`.
  Run `module avail conda` if your site names it differently.

---

## 2. One-time setup (~15–20 min, login node)

```bash
mkdir -p /scratch/$USER/tcrdock_work
cd       /scratch/$USER/tcrdock_work
bash setup_env.sh /scratch/$USER/tcrdock_work
```

`setup_env.sh` is idempotent — it skips `env/` and `TCRdock/` if they already
exist — and does four things:

1. **Conda env** at `./env` (mamba if available, else conda) with the pinned
   stack above plus biopython/numpy/pandas/scipy/matplotlib, and `hmmer` +
   `blast` (TCRdock needs them for V/J gene parsing, not just for MSAs).
2. **Clones TCRdock** (`github.com/phbradley/TCRdock`, shallow) and installs
   `TCRdock/requirements_colab_af232.txt` — the JAX CUDA wheels are listed
   *inside that requirements file*, which is why the conda `cudatoolkit` pin and
   the pip JAX wheel have to agree.
3. **BLAST databases** via `python download_blast.py`, run from inside the
   TCRdock checkout.
4. **AlphaFold parameters (~3 GB)** into `alphafold_params/params/`:
   - `params_model_2_ptm.npz` — stock AF2 weights
   - `tcrpmhc_run4_af_mhc_params_891.pkl` — **TCRdock's fine-tuned TCR:pMHC
     weights**, the ones that make it a TCR docking model rather than plain AF2

Both are fetched from Dropbox. If those links rot, this is the step that breaks,
and the `.pkl` is the irreplaceable one — keep a copy on `/tgen_labs`.

---

## 3. Running

Sequential is simplest; the array is what we actually used.

```bash
sbatch submit_tcrdock_array.sh        # class I
sbatch submit_classII_array.sh        # class II
```

**Class I** (`submit_tcrdock_array.sh`): `gpu-a100`, 1 GPU, 4 CPUs, 32 GB,
4 h, `--array=0-9` → 10 chunks.
**Class II** (`submit_classII_array.sh`): same, but **48 GB and 8 h**,
`--array=0-19` → 20 chunks. Class II peptides are longer and the groove is open
at both ends, so both memory and time go up — this is the one parameter change
between the arms.

Each array task:

1. Splits the input TSV into `chunks/chunk_NNN.tsv` (idempotent — every task
   redoes it in a fraction of a second, which guarantees the file exists no
   matter which task lands first).
2. Works in its **own** `runs/chunk_NNN/` directory. This matters: TCRdock
   writes fixed-name outputs, so without per-chunk dirs concurrent tasks
   overwrite each other.
3. Runs the three-stage pipeline from inside the TCRdock checkout:

```bash
python setup_for_alphafold.py \
    --targets_tsvfile "$CHUNK_DIR/input.tsv" \
    --output_dir      "$CHUNK_DIR/setup" \
    --new_docking --benchmark

python run_prediction.py --verbose \
    --targets        "$CHUNK_DIR/setup/targets.tsv" \
    --outfile_prefix "$CHUNK_DIR/output" \
    --model_names    model_2_ptm_ft4 \
    --data_dir       "$WORKDIR/alphafold_params/" \
    --model_params_files "$WORKDIR/alphafold_params/params/tcrpmhc_run4_af_mhc_params_891.pkl"

python add_pmhc_tcr_pae_to_tsvfile.py \
    --infile  "$CHUNK_DIR/output_final.tsv" \
    --outfile "$CHUNK_DIR/output_w_pae.tsv"
```

Three flags carry real meaning:

- `--new_docking` — the TCR-docking-geometry mode; the whole point of TCRdock.
- `--benchmark` — templates from the deposited structure are **excluded**, so
  the model is not handed its own answer. Without it the comparison against AF3
  is meaningless.
- `--model_names model_2_ptm_ft4` — must match the fine-tuned `.pkl`; the name
  and the weights file are a matched pair.

The third stage is not optional: `pmhc_tcr_pae`, the column the entire Figure 4c
comparison rests on, is added by `add_pmhc_tcr_pae_to_tsvfile.py`, not by
`run_prediction.py`.

### Chunking

Each chunk pays a **one-time ~5–10 min JAX compile**, then ~1–2 min per target.
So finer chunking is not monotonically faster — past a point you are just paying
more compiles. 5–10 chunks is the sweet spot. Useful variants:

| `--array=` | effect |
|---|---|
| `0-9` | 10 chunks (class I default) |
| `0-19` | 20 chunks (class II default) |
| `0-9%4` | 10 chunks, at most 4 concurrent — for a tight GPU quota |
| `3,7` | **retry only chunks 3 and 7** |

---

## 4. Finalizing — and why it is not just `cat`

```bash
bash finalize.sh                                  # class I
INPUT_TSV=tcrdock_input_class_II.tsv bash finalize.sh   # class II
```

`finalize.sh` concatenates `runs/chunk_*/output_w_pae.tsv` into
`user_output_w_pae.tsv`, and — importantly — **counts**:

- reports how many chunk dirs exist vs how many produced output,
- prints the exact `sbatch --array=3,7 ...` line to retry the missing ones,
- warns if the combined row count differs from the input.

That last check is the one that matters. A chunk that dies leaves a *smaller*
combined TSV that is otherwise perfectly well-formed, and every downstream AUC
would be computed on a silently truncated set. Always confirm the final row
count equals the input count before analysing.

**For the record, both arms came out whole:**

| arm | input triads | output rows |
|---|---|---|
| class I | 145 | 145 |
| class II | 1433 | 1433 |

(The 197 in `README_gemini.md` is from an earlier, larger input set.)

---

## 5. Analysis (local, no GPU)

```bash
python analyze_tcrdock_vs_af3.py \
    --tcrdock_tsv user_output_w_pae.tsv \
    --meta_xlsx   tcrdock_input_class_I_full.xlsx
```

Emits the per-antigen AUC table (TCRdock `pmhc_tcr_pae` vs AF3 PTI-PAE) and the
paired Wilcoxon p-value. Results land in `classI_per_antigen_AUC.csv`,
`classI_per_TCR_AUC.csv`, `classI_per_triad_long.csv` (and `classII_*`), which
feed `figures/Figure_4_validation_and_leakage/Fig4_c_AF3_vs_AFTCRdock.ipynb`.

`analyze_classI_v2.sh` / `analyze_classII.sh` are the thin wrappers that were
actually invoked; `analyze_tcrdock_vs_af3_v2.py` is the revision-pass version.

---

## 6. Troubleshooting, in the order things actually go wrong

| symptom | cause | fix |
|---|---|---|
| Predictions take hours each | JAX fell back to CPU | CUDA/cuDNN mismatch — check the log for a backend warning; re-pin `cudatoolkit` ≤ driver CUDA |
| `Unable to initialize backend 'gpu'` | same | as above |
| Setup hangs with no output | Anaconda ToS prompt on stdin | the `conda tos accept` lines in `setup_env.sh` |
| `conda: command not found` | no conda module | edit the `module load anaconda3` fallback |
| Job hits the walltime partway | too many targets per chunk | just resubmit — `run_prediction.py` **skips targets whose output already exists**, so it resumes |
| Combined TSV short of the input | a chunk OOM'd or timed out | `finalize.sh` names the chunks and prints the retry command |
| OOM on one target | long peptide | raise `--mem` (class II already uses 48 GB); with 32 GB + A100 40 GB this did not occur for class I peptides ≤ 11 residues |

---

## 7. What is tracked here vs what is not

Tracked in this repo: every driver script (`setup_env.sh`, the two submit
scripts, `finalize.sh`, the `build_tcrdock_input*.py` and `analyze_*` scripts),
the input TSVs/XLSXs, and the result CSVs.

**Not tracked** (see `.gitignore` and `docs/PROVENANCE.md`): the conda `env/`
(7.7 GB), the `TCRdock/` checkout (1.5 GB, upstream), `alphafold_params/`
(0.75 GB, redistributable weights), and `runs/` + `runs_classII/` (4.8 GB of
per-target PDBs and intermediates). Rebuild the first three with `setup_env.sh`;
the run outputs are reproduced by resubmitting the arrays.
