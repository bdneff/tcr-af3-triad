# Figure 1 reproducibility package

Self-contained package for regenerating Figure 1 of *General prediction of T
cell receptor antigen specificity from sequence using AlphaFold 3*
(Woods et al., *Science Advances*).

## What's in this package

```
.
├── data/
│   ├── pdb_triad.af3_rmsd.parquet           — AF3 + AF-TCRdock RMSDs (one row per PDB)
│   ├── pdb_triad.boltz_rmsd.parquet         — Boltz-2 RMSDs (one row per PDB)
│   ├── pdb_triad.conf_af3.parquet           — AF3 confidence summaries (ipTM, pTM, ranking)
│   └── table_S1_structure_benchmark_complexes.csv
│                                              — Lawson's original input list of 130 PDBs
├── make_figure1.py                          — single script that builds all three panels
└── README.md                                — this file
```

## Quickstart

```bash
pip install pandas pyarrow matplotlib scipy

python3 make_figure1.py
```

Output goes to `./figures/`:

```
figures/
├── figure1_panel_b.pdf / .png            — 3-method CDR RMSD comparison (full dataset)
├── figure1_panel_c.pdf / .png            — Per-component AF3 RMSDs (full dataset)
├── figure1_panel_d.pdf / .png            — AF3 ipTM vs CDR RMSD (full dataset)
├── figure1_panel_b_lawson130.pdf / .png  — same panels, subsetted to Lawson's original ~152 PDBs
├── figure1_panel_c_lawson130.pdf / .png
└── figure1_panel_d_lawson130.pdf / .png
```

## What each panel shows

**Panel B (3-method comparison)**

For each method, CDR RMSD of predicted vs crystal structure (best of 5
diffusion samples). Boxes split by whether each PDB was deposited
**before** the model's training cutoff date (pre-training, filled) or
**after** (post-training, hollow).

Per-method cutoffs (verified against Lawson's notebook
`notebooks/fig_1/3_method_comparison_rmsd.ipynb`):

| Method     | Cutoff      |
|------------|-------------|
| AF-TCRdock | 2018-05-01  | (AF2-Multimer's PDB date cutoff)
| AF3        | 2023-01-12  | (AF3 paper)
| Boltz-2    | 2023-06-01  | (Boltz-2 paper)

Pre/post sets differ across methods because each method has its own
training cutoff. This is the intended behavior (per Lawson's manuscript).

**Panel C (per-component RMSDs)**

For AF3 only, RMSDs of individual components (peptide, MHC, TCR) of
the predicted vs crystal structure. Shows that individual chain folding
is solved (sub-1 Å), and that the docking arrangement is the harder
problem (captured by CDR RMSD in panel B).

**Panel D (ipTM vs CDR RMSD)**

Scatter of AF3's ipTM confidence metric vs AF3's CDR RMSD, colored by
MHC class. Demonstrates that AF3's internal confidence predicts
structural accuracy. Annotated with Spearman correlation.

## Two flavors of every panel

- **No suffix** — uses the full TCR3d-expanded dataset (~269 unique PDBs).
  This is the dataset used in the revised Figure 1.
- **`_lawson130` suffix** — subsets the same data to ~152 unique PDBs
  matching the manuscript's original Figure 1 (130 from Bradley 2023's
  `table_S1` + ~22 STCRDab post-cutoff additions). Provided for
  reproducibility — should closely match the original manuscript figure.

## Statistical comparisons

The script also runs **Mann-Whitney U-tests** (two-sided, pre vs post) for
each method/class combination and prints results to stdout. Sample output:

```
Mann-Whitney U-test (two-sided), pre vs post:
  AF-TCRdock  class I: n_pre= 58 med= 4.64 | n_post= 22 med= 3.85 | p=8.25e-01
  AF-TCRdock  class II: n_pre= 28 med= 2.90 | n_post= 10 med= 3.15 | p=3.04e-01
  AF3         class I: n_pre=182 med= 1.48 | n_post= 17 med= 3.18 | p=7.57e-03
  AF3         class II: n_pre= 60 med= 1.77 | n_post=  7 med= 2.35 | p=4.75e-01
  Boltz-2     class I: n_pre=188 med= 0.54 | n_post= 11 med= 4.10 | p=2.81e-07
  Boltz-2     class II: n_pre= 60 med= 0.74 | n_post=  7 med= 2.50 | p=2.83e-02
```

Key observations consistent with the manuscript:

- **AF-TCRdock** shows no significant pre/post difference for either
  class (p>0.3). This is expected — AF-TCRdock uses curated chain
  templates rather than relying heavily on training-data memorization, so
  it generalizes well to unseen structures.
- **AF3** shows a significant degradation on post-training class I
  structures (p≈0.008) — a smaller but real generalization gap.
- **Boltz-2** shows the largest degradation on post-training entries
  for both classes (p<0.03), suggesting more reliance on training data.
- **Pre-training, Boltz-2 ≈ AF3 > AF-TCRdock** in accuracy; **post-training,
  AF3 outperforms Boltz-2**, consistent with the manuscript's argument
  that AF3 generalizes somewhat better.

## Notes on the dataset

The parquets here were produced by running the
[tcrtrifold-experiments](https://github.com/AltinLab/tcrtrifold-experiments)
Nextflow pipeline with the TCR3d expansion described in the manuscript
revisions. The pipeline does:

1. Resolves all currently-available TCR:p:MHC triads from STCRDab/TCR3d
2. Filters to αβ-TCR + classical class I/II MHC + peptide antigens
3. Runs MSA, then AF3 and Boltz-2 inference on each
4. Computes RMSDs against the deposited crystal structure
5. Extracts AF3 confidence metrics

The current parquets contain **269 unique PDBs** (300 rows because some
multi-chain crystals appear multiple times; the script deduplicates
internally).

### Known caveats

- ~12-13 of Lawson's original 152 PDBs (mostly mouse class I with
  H2-Db/H2-Ld alleles, plus a few class II) are missing from the current
  parquets due to an over-aggressive null-sequence filter in
  `clean_pdb.py` from the TCR3d expansion. A patched pipeline run is
  expected to restore these. They are all pre-cutoff and don't
  qualitatively change the figure.
- Panel D currently uses 139 entries (intersection of AF3 RMSDs and
  AF3 confidence metrics). The confidence parquet contains Lawson's
  original 152; the new TCR3d additions don't yet have confidence
  metrics extracted (needs one more pipeline run of `EXTRACT_CONF_FEAT`).

## Reproducing the figure

```bash
python3 make_figure1.py
# Look in ./figures/ for the outputs
```

That's it. No cluster access needed, no Nextflow, no MSAs to compute —
just the four data files and a Python interpreter with the standard
scientific stack.

## Modifying / extending

The script is structured as three independent functions (`panel_b`,
`panel_c`, `panel_d`) plus a `subset_to_lawson_original` helper. To
change colors, layout, or add new panels, edit `make_figure1.py`
directly — it's intentionally a single self-contained file.
