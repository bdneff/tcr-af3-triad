# TCR3d expansion for Figure 1

Add the TCR3d-curated TCR-pMHC triads that aren't already in the Figure 1
pipeline, then re-run AF3 + Boltz-2 + RMSD on the union, and replot Figure 1.

## What's added vs the current Figure 1

| | Current Fig 1 | After expansion |
|---|---|---|
| **Pre-cutoff** (`replication=True` in pipeline) | 130 | **130 + 178 = 308** |
| **Post-cutoff** (`post_af3_cutoff_valid_pdbs`) | 22 | **22 + 49 = 71** |
| Outlier (6l9l) | 1 | 1 |
| **Total** | 153 | **380** |

The 178 + 49 new entries come from cross-checking the TCR3d class I + II
peptide-complex CSVs against the current pipeline's inputs.

## Why this is needed (Reviewer #1 comment #4)

The current Figure 1 set was built around the **Bradley 2023 *eLife*** input
(`table_S1_structure_benchmark_complexes.csv`), which has 130 pre-cutoff
entries. TCR3d's pre-AF3-cutoff peptide-complex set has 308 entries — 178 of
which are absent from the current Fig 1, and 84% of those are from
2020–2022 (the years the Bradley set stopped collecting). The 49 new
post-cutoff additions extend the held-out benchmark to a broader sample of
recently-deposited structures.

## Pipeline architecture (no new workflow needed)

The repo's existing `pdb.nf` does AF3 + Boltz-2 + all four flavors of RMSD
for whatever is in the input PDB list. **We don't build a new workflow.**
We just patch the inputs:

1. **Append** new pre-cutoff rows to
   `data/pdb/raw/table_S1_structure_benchmark_complexes.csv`.
2. **Append** new PDB IDs to the `post_af3_cutoff_valid_pdbs` list in
   `clean_pdb.py`.
3. For the small subset of TCR3d entries that haven't yet been picked up by
   STCRDab (`db_summary.dat`), **append** rows to `clean_pdb.py`'s
   `outlier_pdb` polars DataFrame, with chain segids and species filled in
   by our sequence-based classifier (validated against STCRDab on the 152
   entries that overlap).

Then `nextflow run ./workflows/pdb.nf -resume` picks up where it left off:
existing predictions are kept, only the new triads run through AF3 +
Boltz-2 + RMSD.

## Run order

### Prerequisites

- Local checkout of `tcrtrifold-experiments-main` (the public repo)
- This package directory (`tcr3d_expansion/`)
- Both TCR3d peptide-complex CSVs in `data/`:
    - `data/tcr3d_classI_complexes.csv`
    - `data/tcr3d_classII_complexes.csv`
- Network access to RCSB (steps 01 & validation only) — Gemini works,
  local also works
- Gemini access for step 03 (the actual pipeline run)

### Step 00 (recommended) — Refresh STCRDab

The on-disk `db_summary.dat` in the repo is a snapshot. STCRDab itself has
been updated since (634 PDBs vs ~496 in the snapshot). The bulk
`db_summary.dat` download URL is no longer published, but per-PDB summaries
still work, so this step assembles a fresh `db_summary.dat` by querying the
public per-PDB endpoint for every TCR3d PDB.

```bash
export REPO_ROOT=/path/to/tcrtrifold-experiments-main
bash scripts/00_refresh_stcrdab.sh
```

This produces `data/db_summary_fresh.dat`. Pass it to step 01 with
`--stcrdab_summary` (see below). With the fresh STCRDab, ~140 additional
PDBs become covered, dropping the "needs classifier" count from ~46 to
near zero.

Per-PDB fetches are cached in `data/stcrdab_cache/` so reruns are fast.

### Step 01 — Identify new entries

Run from anywhere with network access. Hits RCSB to fetch chain sequences
for STCRDab-missing entries; cached on disk (`data/rcsb_cache/`) so reruns
are fast.

```bash
export REPO_ROOT=/path/to/tcrtrifold-experiments-main
bash scripts/01_identify_new_entries.sh
```

Outputs (all in `data/`):
- `tcr3d_new_pre_cutoff.csv`        — rows to append to table_S1
- `tcr3d_new_post_cutoff_pdbs.txt`  — IDs to add to the post-cutoff list
- `tcr3d_new_outlier_rows.csv`      — segid annotations for STCRDab-missing PDBs
- `classifier_validation.csv`       — sanity check vs STCRDab (152 entries)

**Look at `classifier_validation.csv` before going further.** The classifier
should reproduce STCRDab on ≥95% of the 152 cross-checked entries. If
accuracy is lower, the motif library in `src/chain_classifier.py` may need
tightening before trusting the predictions on STCRDab-missing PDBs.

If you want to refresh `db_summary.dat` first (recommended — STCRDab updates
weekly), download it from `http://opig.stats.ox.ac.uk/webapps/stcrdab/` and
pass it explicitly:

```bash
python scripts/01_identify_new_entries.py \
    --tcr3d_classI  data/tcr3d_classI_complexes.csv \
    --tcr3d_classII data/tcr3d_classII_complexes.csv \
    --repo_root     "$REPO_ROOT" \
    --stcrdab_summary /path/to/fresh_db_summary.dat
```

A fresh STCRDab summary will likely cover almost all of the 45 currently-
missing post-cutoff entries (since they're mostly 2024-2025 PDBs that
weekly-updated STCRDab will have picked up by now). That means the classifier
becomes a fallback rather than the main path.

### Step 02 — Patch the repo

**Dry run first:**

```bash
export REPO_ROOT=/path/to/tcrtrifold-experiments-main
bash scripts/02_apply_to_repo.sh --dry_run
```

This shows what will change without modifying anything. If it looks right:

```bash
bash scripts/02_apply_to_repo.sh
```

Three files are modified, each with a timestamped backup written next to it
(`*.tcr3d_backup.YYYYMMDD-HHMMSS`):

- `data/pdb/raw/table_S1_structure_benchmark_complexes.csv`
- `workflows/subworkflows/local/cleaning/resources/usr/bin/clean_pdb.py`
- `workflows/subworkflows/local/cleaning/resources/usr/bin/clean_pdb.py` (twice — `post_af3_cutoff_valid_pdbs` list and `outlier_pdb` DataFrame)

The patch is idempotent: re-running skips PDBs already present. Use
`--min_confidence high` to only inject classifier predictions tagged
high-confidence (default is medium).

### Step 03 — Run the pipeline on Gemini

```bash
# On Gemini, inside the patched repo:
cd /path/to/tcrtrifold-experiments-main
mkdir -p logs/pdb_tcr3d_expansion

# (update USER@tgen.org in 03_run_pipeline.sh first)
sbatch /path/to/tcr3d_expansion/scripts/03_run_pipeline.sh
```

This kicks off the existing `pdb.nf` pipeline with `-resume`. Only the new
triads run through AF3 + Boltz-2 + RMSD; existing 153 entries reuse cached
predictions.

Approximate cost (new triads × 2 models, MSA-deduped by unique chains):
- AF3: a few hours (similar density to the cross_reactivity job that took
  ~2 h for 750 triads)
- Boltz-2: significantly heavier per-triad with 5 recycle × 5 diffusion;
  budget 1–2 days
- RMSD: minutes (CPU work)

### Step 04 — Replot Figure 1

Once Gemini has produced the updated staged parquets
(`data/pdb/triad/staged/pdb_triad.af3_rmsd.parquet` and `..._boltz_rmsd...`),
either copy them locally or run step 04 directly on Gemini:

```bash
export REPO_ROOT=/path/to/tcrtrifold-experiments-main
bash scripts/04_replot_figure1.sh
```

Outputs (in `results/fig1_v2/`):
- `fig1b_cdr_rmsd_compare.{png,pdf}` — AF-TCRdock / AF3 / Boltz-2 boxplots,
  one row per MHC class, blue=class I and red=class II palette
- `fig1c_rmsd_by_chain.{png,pdf}`    — per-chain breakdown (peptide / MHC / TCR)
- `fig1d_confidence_correlations.{png,pdf}` — ipTM and docking-Z vs CDR RMSD

## File map

```
tcr3d_expansion/
├── README.md                              (you are here)
├── data/
│   ├── tcr3d_classI_complexes.csv         (TCR3d class I peptide complexes)
│   ├── tcr3d_classII_complexes.csv        (TCR3d class II peptide complexes)
│   ├── rcsb_cache/                        (auto-created; cached RCSB API responses)
│   ├── tcr3d_new_pre_cutoff.csv           (step 01 output)
│   ├── tcr3d_new_post_cutoff_pdbs.txt     (step 01 output)
│   ├── tcr3d_new_outlier_rows.csv         (step 01 output)
│   └── classifier_validation.csv          (step 01 output, sanity check)
├── src/
│   ├── chain_classifier.py                (motif-based chain ID; works on raw RCSB FASTA)
│   └── rcsb_fetch.py                      (RCSB FASTA + entity-metadata wrappers, with caching)
└── scripts/
    ├── 00_refresh_stcrdab.py / .sh         (build a fresh db_summary.dat from STCRDab per-PDB endpoints)
    ├── 01_identify_new_entries.py / .sh   (cross-check vs current Fig 1, classify missing chains)
    ├── 02_apply_to_repo.py / .sh          (patch table_S1.csv + clean_pdb.py; idempotent + backed up)
    ├── 03_run_pipeline.sh                 (SLURM submission; runs pdb.nf -resume on Gemini)
    └── 04_replot_figure1.py / .sh         (replot Fig 1b/c/d with the new data; blue/red palette)
```

## Notes and gotchas

- **The pipeline run uses `-resume`.** This is intentional and means the
  pipeline reuses every existing AF3 and Boltz-2 prediction. If you've made
  other changes upstream that should invalidate the cache, drop `-resume`
  from `03_run_pipeline.sh`. Otherwise leave it: rerunning ~153 triads of
  AF3+Boltz wastefully would add days.

- **The 22 hardcoded post-cutoff PDBs in `clean_pdb.py`** can drift over
  time. If Lawson updates that list before this script runs, the constant
  `EXISTING_POST_AF3` in `01_identify_new_entries.py` should be re-synced.
  Otherwise step 01 will treat newly-added IDs as "new" and step 02 will
  silently skip them (idempotency saves us).

- **The classifier's motif library** in `src/chain_classifier.py` is
  designed for full RCSB FASTA sequences (V+J+C TCR chains, mature MHC
  chains). It will NOT work correctly on the trimmed/cleaned sequences that
  live in the pipeline's `*.cleaned.parquet`. This is fine because the
  classifier only ever sees raw RCSB FASTA.

- **Confidence thresholds.** Step 02 defaults to `--min_confidence medium`,
  meaning predictions tagged "low" are skipped. Examine
  `classifier_validation.csv` (step 01) to see whether "low" predictions
  are actually wrong on the validation set — if they're consistently
  correct, you can re-run with `--min_confidence low` to accept everything.

- **Plot tweaks.** `scripts/04_replot_figure1.py` keeps the manuscript's
  three-panel structure (b/c/d) but cleans up the palette to blue=class I /
  red=class II throughout, with darker shades for post-cutoff. The
  function is parameterized — adjust `CLASS_COLORS` at the top of the file
  to retune.
