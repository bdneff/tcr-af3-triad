# Environment — the Gemini cluster

Everything compute-heavy runs on TGen's Gemini cluster. This is the reference for
partitions, modules, containers, and where things live.

## Filesystems

| Path | Nature | Use for |
|---|---|---|
| `/tgen_labs/altin/alphafold3/workspace` (`$WS`) | durable lab workspace | anything long-lived: pipeline checkouts, weights, staged data |
| `/scratch/$USER` | **ephemeral, swept ~30 days** | repo checkouts, run outputs in flight |
| `/ref_genomes/alphafold/alphafold3/` | read-only reference | AF3 model weights + genetic databases |
| `/home/$USER` | small, backed up | not for run data |

A `/scratch` checkout is disposable by design — clone it again whenever. Anything you
would be sorry to lose must be copied to `$WS`.

### The tcrtrifold pipeline's data layout (verified 2026-09-09)

Everything the paper's predictions and metrics come from lives under
`$WS/tcrtrifold-experiments/data/pdb/triad/`. Written down because it took several wrong
guesses to find, and `find` over that tree is slow enough to be unusable interactively —
bound the depth or go straight to the path.

| Path (under `data/pdb/triad/`) | Contents |
|---|---|
| `cleaned_pdb/` | **288** canonicalised crystal structures, `<pdb>.pdb`, one per triad — this is `rmsd.py`'s `clean_pdb_dir` |
| `staged/` | all the metric parquets (`pdb_triad.af3_rmsd.parquet`, `pdb_triad.conf_af3.parquet`, `pdb_validation_triad.*`, …) |
| `predictions/` | inference outputs, directories named by a **content hash**, not by PDB ID. Contents are `<hash>_model_0..4.cif`, `confidence_<hash>_model_N.json`, `pae_*.npz`, `pde_*.npz`. The confidence schema (`ligand_iptm`, `complex_pde`, `complex_ipde`, `pair_chains_iptm`) and the presence of `pde_*` are **Boltz**, not AF3 — AF3 emits no PDE. So this directory is probably the Boltz arm; locate the AF3 outputs before assuming. |
| `predictions_with_templates/` | a second, separate inference run — check which one a number came from before quoting template settings |
| `inference/`, `raw/` | inputs and unprocessed downloads |

**No input JSON is retained** in a prediction directory, so chain order cannot be recovered from
there. **PAE is not in the confidence JSON** either — it is in the sibling `pae_*.npz`, which means
`mean_p_tcr_interface_pae` is computed by pipeline code from the full PAE matrix rather than read
from a field. Find that code before recomputing it for new predictions.

The cleaner itself is at
`$WS/tcrtrifold-experiments/workflows/subworkflows/local/cleaning/resources/usr/bin/clean_pdb.py`
(copied into `docs/provenance/` for reference).

`/tgen_labs` is **not mounted on the login node** — `ls` there fails with
`No such file or directory` and `mkdir` with `Permission denied`. Grab an allocation first
(`srun -p compute -c 1 --mem=2G -t 00:10:00 --pty bash -i`) for anything touching it, including
copies out of `/scratch`.

### Canonical chain convention

`clean_pdb.py` emits five chains in a fixed order, the same for class I and class II:

| chain | class I | class II |
|---|---|---|
| **A** | **peptide** | **peptide** |
| **B** | MHC heavy | MHC α |
| **C** | β₂-microglobulin | MHC β |
| **D** | TCR α | TCR α |
| **E** | TCR β | TCR β |

Verified by CA counts: `8enh` (class I) → A=9, B=275, C=99, D=205, E=245; `8vd0` (class II)
→ A=16, B=183, C=178, D=202, E=240.

**The peptide is chain A, not chain C.** This is easy to get backwards — some tools and inputs
use A=MHC, B=β₂m, C=peptide, the opposite way round — and `rmsd.py` selects `segid D`/`segid E`
for CDR RMSD, so a swapped A/B/C survives *that* metric while corrupting anything touching the
MHC frame or the peptide. Match this convention in any new AF3 input JSON.

**Tethered class II peptides are already handled.** `8vd0`, `8vq8` and `9aud` fuse the peptide
to the MHC β chain through a linker, so the deposit has four polymer entities where a prediction
has five. `clean_pdb.py` splits them correctly: `8vd0`'s cleaned A is exactly the 16-residue
peptide and its C is the β chain without it. No special handling is needed downstream.

### AF3 template search is date-capped by default

`max_template_date` defaults to **`2021-09-30`** (`run_alphafold.py` ~line 194, *"By default,
use the date from the AlphaFold 3 paper"*), read into the search config at ~line 718. Nothing in
the pipeline or in our sbatch overrides it. So "templates on" cannot reach anything deposited
after Sept 2021 — including every structure in the post-training-cutoff validation sets. Verify
with:

```bash
singularity exec /tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif \
    sed -n "188,202p" /app/alphafold/run_alphafold.py
```

## Partitions and scheduling

- **`gpu-a100`** — the GPU partition used throughout (A100, 40 GB HBM). 4-day walltime cap.
- Typical GPU request: `--gres=gpu:1 --cpus-per-task=8 --mem=32G` (raise memory for class II
  and for AF3 with MSAs, which we run at 64 GB).
- Submit from the directory the script expects; scripts anchor to `$SLURM_SUBMIT_DIR`.
- Job arrays: `--array=0-9`, throttle with `%4`, retry a subset with `--array=3,7`.
  Chain stages with `--dependency=afterok:$JOBID` (capture the ID with `sbatch --parsable`).
- Never paste a `<placeholder>` into a shell — `<` is read as a redirect. Substitute the
  real job IDs.

## AlphaFold 3 (lab container)

```bash
module load singularity
singularity exec --nv -B /home,/scratch,/tgen_labs,/ref_genomes \
    /tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif \
    python /app/alphafold/run_alphafold.py \
        --json_path "$JSON" \
        --model_dir /ref_genomes/alphafold/alphafold3/models \
        --db_dir    /ref_genomes/alphafold/alphafold3/ \
        --output_dir "$OUT"
```

Paths are in `configs/paths.py`.

**`--norun_data_pipeline` is the difference that matters.** With it, AF3 folds from whatever
MSA the JSON carries (fast, no databases, ~10 min). Without it, AF3 searches the genetic
databases itself (slow, hours, and the only mode comparable to the manuscript's triads).
The JSON and the flag are a **matched pair**: JSONs staged with `af3_stage.py --mode single`
embed a self-MSA and need the flag; `--mode pipeline` omits the MSA fields and must not have
it. Mismatching them does not error — it quietly folds from a bare sequence.

Note that every earlier AF3 job in these repos used `--norun_data_pipeline`, so `--db_dir`
was passed but never read. A job that runs the data pipeline is the first to actually
exercise the genetic databases; if they are absent, that is where it will show.

Output: `<out>/<job_name_lowercased>/*_model.cif` — AF3 lowercases the job name.

## Python on the cluster

**`conda activate tcrtrifold-experiments`** — the working env for analysis on Gemini
It is the same env as the
upstream prediction pipeline, so analysis runs where the predictions were made.

Do **not** run from `base` — it is the shared Mamba install at `/packages/Mamba/`, not a
project environment.

Note this is distinct from **`$WORKDIR/env`** on scratch, which is the CUDA-11.8-pinned
AlphaFold-2.3.2 environment built by `setup_env.sh` **solely** for AF-TCRdock. Do not use it
for general analysis and do not install into it — its pins are load-bearing (see
[`TCRDOCK.md`](TCRDOCK.md)).

## AF-TCRdock

A self-built conda environment rather than a container, pinned to CUDA 11.8 / cuDNN 8.9.
Full recipe and failure modes: **[`TCRDOCK.md`](TCRDOCK.md)**.

## Getting the repo there

```bash
cd /scratch/$USER
git clone https://github.com/bdneff/tcr-af3-triad.git
cd tcr-af3-triad
```

Afterwards, `git pull` to pick up changes.
