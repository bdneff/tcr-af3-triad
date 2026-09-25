# Running TCRdock on Gemini — short version

You need three things in a scratch directory: the env setup script, the SLURM
submission script, and the input TSV. After running the setup once, every
future job is just `sbatch submit_tcrdock.sh`.

## Step 1 — copy files to scratch

```bash
mkdir -p /scratch/$USER/tcrdock_work
cd       /scratch/$USER/tcrdock_work
# copy these three files in from wherever you downloaded them:
#   setup_env.sh
#   submit_tcrdock.sh
#   tcrdock_input_class_I.tsv
```

## Step 2 — set up the environment (once, ~15-20 min)

On a Gemini login node:

```bash
bash setup_env.sh /scratch/$USER/tcrdock_work
```

This does everything in one shot: builds a conda env in `./env`, clones
TCRdock, installs all dependencies, downloads BLAST databases, and pulls
the two AlphaFold parameter files (~3 GB). If your site doesn't have
conda available by default, the script tries `module load anaconda3` —
edit line 24 of `setup_env.sh` if your cluster uses a different module name.

## Step 3 — open the submission script and check the SLURM headers

Open `submit_tcrdock.sh`. Lines 11-13 are commented out (`##SBATCH`) — these
are partition/account/GPU constraints that vary by site. Uncomment them and
fill in whatever Gemini requires. If you're not sure:

```bash
sinfo                  # lists partitions
scontrol show partition <gpu_partition_name>   # shows what GPUs it has
```

A safe starting point if Gemini doesn't enforce accounts:

```
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
```

## Step 4 — submit

You have two options. Pick one.

### Option A: Sequential (simplest, one job, ~6-7h on A100)

```bash
sbatch submit_tcrdock.sh
```

That's it. Single SLURM job, single GPU, processes all 197 triads in series.

### Option B: Parallel job array (fastest, ~30-50 min wall clock)

```bash
sbatch submit_tcrdock_array.sh
```

This submits 10 independent jobs that SLURM will schedule across whatever
GPUs Gemini gives you. Each handles ~20 triads. The total wall time is just
the slowest chunk: ~30-50 min if all 10 run truly concurrently, longer if
your queue is busy.

**After all array tasks finish**, concatenate the chunk outputs into one file:

```bash
bash finalize.sh
```

`finalize.sh` also reports any chunks that didn't finish, with the exact
`sbatch` command to retry just those:

```
Chunks missing output_w_pae.tsv (re-submit these IDs):
  chunk 3: runs/chunk_003/
To retry only these chunks:
  sbatch --array=3 submit_tcrdock_array.sh
```

### Tuning the parallelism

Edit `#SBATCH --array=0-9` near the top of `submit_tcrdock_array.sh`:

- `--array=0-9` (default): 10 chunks of ~20 triads each
- `--array=0-19`: 20 chunks of ~10 triads. Faster if you have ≥20 GPUs available, but each chunk pays the ~5 min JAX compile cost so finer ≠ always faster.
- `--array=0-9%4`: 10 chunks, max 4 running at once. Use if your GPU quota is tight.
- `--array=0-1`: 2 chunks of ~100 each. Closer to the "batches of 100" pattern, lighter on the scheduler.

The sweet spot is usually 5-10 chunks: enough parallelism to matter, not so fine that compilation overhead dominates.

## Timing

**Sequential** (`submit_tcrdock.sh`): on an A100 with `--new_docking`:
- First target: ~5-10 min (one-time JAX compilation)
- Each subsequent target: ~1-2 min
- Total for 197 triads: ~6-7 hours

**Parallel array** (`submit_tcrdock_array.sh` with 10 chunks): each chunk
pays its own ~5-10 min compile cost, then does ~20 triads at ~1-2 min each
= ~30-50 min per chunk. With 10-way concurrency you finish in roughly the
time of one chunk, ~40 min wall clock. Caveat: if Gemini's queue is busy
and you get scheduled in waves of 2-3, expect closer to 1.5-2h.

Either way it fits in a 24h SLURM allocation.

## What you get back

Two things end up in `/scratch/$USER/tcrdock_work/`:

- `user_output_w_pae.tsv` — the one file you need. Has 197 rows with all the
  TCRdock metrics including `pmhc_tcr_pae` (the column we use for the AUC
  comparison). Hand this to `analyze_tcrdock_vs_af3.py`.
- `TCRdock/user_output_T*.pdb` — the predicted structures, one per triad.

## Final analysis (local, no GPU needed)

```bash
python analyze_tcrdock_vs_af3.py \
    --tcrdock_tsv user_output_w_pae.tsv \
    --meta_xlsx   tcrdock_input_class_I_full.xlsx
```

Prints the per-antigen AUC table (TCRdock vs AF3 PTI-PAE side by side) and
the paired Wilcoxon p-value. Drop the resulting CSVs into whatever
plotting code you'd use for Figure 4a.

## If something breaks

- **`conda: command not found`** when running setup — edit line 24 of
  `setup_env.sh` to match your cluster's conda module (try `module avail conda`
  on a login node to see what's available).
- **CUDA/JAX errors in the SLURM log** — the env was built for CUDA 11.8.
  If Gemini's GPU drivers are newer, you may need to swap `cudatoolkit=11.8.0`
  in `setup_env.sh` for a matching version (`nvidia-smi` on a GPU node shows
  what CUDA the drivers support; pick a `cudatoolkit` ≤ that).
- **Job times out partway through** — `run_prediction.py` skips targets whose
  output already exists, so just `sbatch submit_tcrdock.sh` again and it'll
  pick up where it left off.
- **Out of memory on a particular target** — TCRdock will fail that one and
  move on. Check `user_output_w_pae.tsv` row count at the end; if it's < 197,
  the missing rows are the OOMs. With 64 GB system RAM and an A100 (40 GB
  HBM) this should not happen for class I peptides ≤ 11 residues.
