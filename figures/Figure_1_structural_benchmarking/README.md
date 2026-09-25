# Figure 1 — Structural benchmarking

Self-contained reproduction of Figure 1 of *General prediction of T cell receptor
antigen specificity from sequence using AlphaFold 3* (Woods et al., *Science Advances*).

## What's here

```
.
├── data/
│   ├── pdb_triad.af3_rmsd.parquet           — AF3 + AF-TCRdock RMSDs (318 rows, 287 unique PDBs)
│   ├── pdb_triad.boltz_rmsd.parquet         — Boltz-2 RMSDs (287 unique PDBs)
│   ├── pdb_triad.conf_af3.parquet           — AF3 confidence summaries for all 287 PDBs
│   └── table_S1_structure_benchmark_complexes.csv
│                                              — Lawson's original input list of 130 PDBs
├── make_figure1.py                          — standalone script (canonical source)
├── make_figure1.ipynb                       — same code as the script, split into cells
├── figures/                                 — 6 PNGs (panel b/c/d × {full, lawson130})
└── README.md
```

The script is the canonical source; the notebook just splits it into cells with
markdown headers and re-runs the same panel functions. Both produce identical
output.

## Quickstart

```bash
pip install pandas pyarrow matplotlib scipy
python3 make_figure1.py
```

This regenerates all six PNGs in `./figures/`:

```
figure1_panel_b.png            — 3-method CDR RMSD comparison (full TCR3d set)
figure1_panel_c.png            — Per-component AF3 RMSDs (full TCR3d set)
figure1_panel_d.png            — 4-subpanel scatter: ipTM, docking Z, peptide length, species
figure1_panel_b_lawson130.png  — same panels restricted to Lawson's original ~152-PDB set
figure1_panel_c_lawson130.png
figure1_panel_d_lawson130.png
```

## What each panel shows

**Panel B — 3-method CDR RMSD comparison.** AF-TCRdock vs AF3 vs Boltz-2, CDR
RMSD of predicted vs crystal (best of 5 diffusion samples). Boxes split by
whether each PDB was deposited **before** the model's training cutoff
(pre-training, filled) or **after** (post-training, hollow). Per-method
cutoffs (verified against Lawson's notebook):

| Method     | Cutoff      |
|------------|-------------|
| AF-TCRdock | 2018-05-01  | (AF2-Multimer's PDB date cutoff)
| AF3        | 2023-01-12  | (AF3 paper)
| Boltz-2    | 2023-06-01  | (Boltz-2 paper)

Pre/post sets differ across methods because each method has its own training
cutoff — this is intended (per Lawson's manuscript).

The panel carries 4 staggered Mann-Whitney brackets per MHC class comparing
adjacent methods at matched pre/post status:

- AF-TCRdock-pre ↔ AF3-pre
- AF-TCRdock-post ↔ AF3-post
- AF3-pre ↔ Boltz-2-pre
- AF3-post ↔ Boltz-2-post

Tests are **one-tailed** (`alternative='greater'` in `scipy.stats.mannwhitneyu`,
testing H₁: the older model has higher RMSD than the newer). p-values are
always shown numerically (no `n.s.` substitution).

**Panel C — per-component AF3 RMSDs.** For AF3 only, RMSDs of individual
components (peptide, MHC, TCR) of the predicted vs crystal structure. Shows
that individual chain folding is solved (sub-1 Å) and that the docking
arrangement is the harder problem (which panel B's CDR RMSD captures).

Three pre-vs-post Mann-Whitney brackets per MHC class, one per component.
Tests are **one-tailed** (`alternative='less'`, testing H₁: pre RMSD < post
RMSD — AF3 expected to do worse on unseen data).

**Panel D — confidence/feature vs accuracy.** Four scatter sub-panels relating
AF3 prediction features to CDR RMSD:
- ipTM (AF3 internal confidence) — strong negative correlation, AF3's
  confidence really does track accuracy
- docking-geometry Z-score from the consensus PDB distribution — weak positive
- peptide length — flat (no length dependence)
- MHC species (human vs mouse) — also flat

Spearman / Mann-Whitney shown in-panel.

## Two flavors of every panel

- **No suffix** — uses the full TCR3d-expanded dataset (~269 unique PDBs). This
  is the dataset for the revised Figure 1.
- **`_lawson130` suffix** — subsets the same data to ~152 PDBs matching the
  manuscript's original Figure 1 (130 from Bradley 2023's `table_S1` + ~22
  STCRDab post-cutoff additions). Kept around for sanity-checking against the
  earlier manuscript figure.

## Statistical output (printed to stdout)

`make_figure1.py` prints Mann-Whitney one-tailed p-values to the console as a
sanity check, in addition to drawing them on the figures. Sample:

```
  Mann-Whitney U-test (one-tailed, H1: a > b), pre vs post:
    AF-TCRdock  class I:  n_pre= 67 med= 4.41 | n_post= 22 med= 3.85 | p=8.25e-01
    AF-TCRdock  class II: n_pre= 30 med= 3.02 | n_post= 10 med= 3.18 | p=3.04e-01
    AF3         class I:  n_pre=196 med= 1.49 | n_post= 18 med= 3.18 | p=3.78e-03
    AF3         class II: n_pre= 63 med= 1.78 | n_post=  7 med= 2.35 | p=4.75e-01
    Boltz-2     class I:  n_pre=203 med= 0.55 | n_post= 11 med= 4.10 | p=1.41e-07
    Boltz-2     class II: n_pre= 63 med= 0.74 | n_post=  7 med= 2.50 | p=2.83e-02
```

What this says:
- **AF-TCRdock** doesn't degrade on post-training entries for either class (its
  curated chain templates generalize fine).
- **AF3** shows a real degradation on post-training class I (p=4e-3) — a real
  but small generalization gap.
- **Boltz-2** shows the biggest degradation on post-training entries for both
  classes (p<0.03), so it leans more on training-data memorization than AF3.
- Pre-training, **Boltz-2 ≈ AF3 > AF-TCRdock**; post-training, **AF3 holds up
  best** — consistent with the manuscript's "AF3 generalizes better" point.

The on-figure brackets in panel B show the cross-method versions of these
comparisons (older model vs newer at matched pre/post).

## Notes on the dataset

The parquets here were produced by running the
[tcrtrifold-experiments](https://github.com/AltinLab/tcrtrifold-experiments)
Nextflow pipeline with the TCR3d expansion described in the manuscript
revisions. The repo checkout on gemini is at
`/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/`. The pipeline:

1. Resolves all currently-available TCR:p:MHC triads from STCRDab/TCR3d
2. Filters to αβ-TCR + classical class I/II MHC + peptide antigens
3. Runs MSA, then AF3 and Boltz-2 inference on each
4. Computes RMSDs against the deposited crystal structure
5. Extracts AF3 confidence metrics

The current parquets contain **287 unique PDBs**. The earlier "12-13 PDBs
missing due to an over-aggressive null-sequence filter in clean_pdb.py" caveat
(mouse class I H2-Db/H2-Ld plus a few class II) is **resolved** — those entries
are now in the parquets. Panel D also runs on the full set now (no more 139-row
caveat from missing confidence metrics).

## How to extend or restyle

The script is structured as three independent functions (`panel_b`, `panel_c`,
`panel_d`) plus a `subset_to_lawson_original` helper. Edit `make_figure1.py`
directly for colors, layout, or new panels — it's a single self-contained
file. The notebook re-imports the same functions, so the script is the only
thing that needs to change.

## Revision pass — what changed from the original Figure 1

- **Panel B**: added 4 staggered cross-method Mann-Whitney brackets per MHC
  class (AF-TCRdock-pre ↔ AF3-pre, AF-TCRdock-post ↔ AF3-post, AF3-pre ↔
  Boltz-2-pre, AF3-post ↔ Boltz-2-post). Tests are now **one-tailed** and the
  numeric p is always shown. Boxes widened (0.7 → 0.85), method groups
  compressed (3.0 → 2.2 spacing), figure trimmed (6.4" → 5.8" wide), ylim
  raised to clear the long class-II AF3-pre whisker.
- **Panel C**: added 3 within-component pre-vs-post Mann-Whitney brackets per
  MHC class (Peptide / MHC / TCR), also **one-tailed**. Boxes widened to 0.85;
  the per-box `n=X` annotations were removed (n values are in the manuscript
  text).
- **Panel D**: axis labels and subplot titles bumped to publication weights
  (axes 16, subplot titles 16, suptitle 15, ticks 13, legends 11, stat boxes
  12). Species sub-panel boxes widened to 0.85; figure widened to 11.5×9.0".
- **Output format**: PNGs only (no PDFs).
