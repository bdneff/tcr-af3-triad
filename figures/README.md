# Manuscript plots — reproducibility package

Plots, data, and notebooks for the Woods et al. AF3 TCR:antigen specificity
manuscript, organized by manuscript figure number.

## Axis convention (applies to every figure)

All plots that show the AF3 interface-PAE feature on an axis use the inverted
quantity **`AF3 PTI-PAE` = `31 − mean_p_tcr_interface_pae`**, so that **higher =
more confident interface / better prediction** everywhere. Per the manuscript
definition, "PTI-PAE" *is* this inverted (31 − raw) quantity, so axes are
labelled **`AF3 PTI-PAE`**. The raw feature `mean_p_tcr_interface_pae` is
low-is-better; the inversion makes it high-is-better. Where an AUC is computed
from the feature, the same inversion is applied internally so AUC > 0.5 means
the classifier ranks cognates above non-cognates.

> Label note: figure axes (Fig 3c, 4a–d, Supp Fig 2, Supp Fig 3 panels) all
> read `AF3 PTI-PAE`. **Fig 4e** (`Fig4_e_*`) still reads `31 - PTI-PAE` —
> it depends on RMSD inputs that weren't re-rendered in the current pass.
> Doesn't change any numbers; flagged for a later sweep.

## Structure

```
manuscript_plots/
├── Figure_1_structural_benchmarking/      ← Manuscript Figure 1 (b, c, d)
│   ├── make_figure1.ipynb / .py            ← 3 panels × {full TCR3d set, lawson130 subset}
│   ├── data/                               ← 3 parquets + table_S1 CSV
│   ├── figures/                            ← 6 PNGs
│   └── README.md
│
├── Figure_3_affinity_correlation/         ← Manuscript Figure 3c
│   ├── Fig3c_affinity_correlation.ipynb
│   ├── corr_spr_upstream.ipynb            ← upstream IEDB DB → fig3c_data.csv derivation
│   ├── fig3c_data.csv
│   └── README.md
│
├── Figure_4_validation_and_leakage/       ← Manuscript Figure 4 (a–e)
│   ├── Fig4_ab_validation_ROCs.ipynb       ← panels a, b (ROC curves)
│   ├── Fig4_c_AF3_vs_AFTCRdock.ipynb       ← panel c (AUC violins: antigen + TCR-centric)
│   ├── Fig4_d_decomposed_leakage.ipynb     ← panel d (peptide / MHC / TCR distance vs AUC)
│   ├── Fig4_e_PTIPAE_vs_RMSD_validation.ipynb ← panel e (PTI-PAE vs CDR RMSD)
│   ├── *.csv                               ← per-antigen / per-TCR AUCs + distances
│   └── README.md
│
├── reviewer_responses/
│   ├── cross_reactivity/                   ← Supp Fig 2 (4 PNGs: violin/scatter × BLOSUM62/Hamming)
│   ├── structural_analysis/                ← Supp Fig 3 (RMSD + mean dihedral distance vs PTI-PAE)
│   └── study_design_flowchart/             ← study-design schematic
│
├── supplementaryTable3_updated.xlsx       ← source of truth for per-triad PTI-PAE
└── README.md                              ← this file
```

## What changed during the revision pass

A consistent set of styling and content updates was applied across every
figure:

- **One-tailed statistical tests where directional.** Figure 1 panels B and C
  switched their Mann-Whitney tests from two-sided to one-tailed:
  - **Panel B**: `alternative='greater'`, testing H₁ that the older method
    (AF-TCRdock vs AF3, AF3 vs Boltz-2) has higher RMSD than the newer.
  - **Panel C**: `alternative='less'`, testing H₁ that AF3's pre-cutoff RMSD
    is lower than its post-cutoff RMSD (AF3 expected to do worse on unseen
    data).
  - **Figure 4d** Spearman: `alternative='less'`, so the p-value directly
    tests the leakage hypothesis (negative correlation between distance and
    AUC). All six p-values come back ≥ 0.58, so the leakage hypothesis isn't
    supported.

  Numeric p-values are now always shown (no `n.s.` substitution).

- **Stat brackets on boxplots.** Panels B and C of Figure 1 carry on-figure
  Mann-Whitney brackets per MHC class — 4 staggered cross-method brackets in
  B (pre and post × adjacent method pairs), 3 within-component brackets in C
  (Peptide / MHC / TCR).

- **Wider boxes, compressed spacing, bigger fonts.** Figure 1 panels B/C/D
  figure widths trimmed, box widths bumped from 0.55–0.70 to 0.85. All axis
  labels, tick labels, sub-panel titles, legends and stat-box text in panel
  D raised to publication weights (axes 16, ticks 13, subplot titles 16,
  suptitle 15, legends 11).

- **Figure 3c re-styled.** Decade boxes span the full log interval; the
  non-cognate box is red (matching its legend); decade boxes are colored
  green or blue by majority class of the cognate points inside them. Legend
  and Spearman stats now live as a compact stacked inset in the empty
  bottom-right corner of the cognate panel, so the figure width drops from
  8.6" to 7.0".

- **Figure 4 panel D simplified.** Distance axes are labelled
  `Peptide distance to nearest PDB triad`, `MHC distance to nearest PDB triad`,
  and `TCRdist to nearest TCR3d entry`. Stat lines read
  `class I  n=14  r=+0.62, p=0.99` (the test is Spearman, one-tailed
  `alternative='less'`; documented in the Figure_4 README).

- **Supp Fig 2 (cross-reactivity).** All four PNGs now use the `AF3-PTI-PAE`
  label. The MW caption and right-side Spearman inset were removed; per-bin
  MW annotation boxes are staggered vertically so adjacent bins do not
  overlap.

- **Supp Fig 3 (structural).** The CDR-RMSD-vs-PTI-PAE plot has its stats box
  moved into the figure top margin (outside the axes), with a fixed 0–31
  x-axis. The previous Total Variation / KDE metric was replaced with **mean
  per-residue dihedral distance (°)** for the second structural plot —
  defined as the average over CDR residues of √(Δφ² + Δψ²) with angle
  differences wrapped into [-180°, 180°]. Much simpler to describe in
  methods, and tracks the same local backbone fidelity.

## Per-figure summary

### Figure 1 — structural benchmarking
- **Panel B**: AF-TCRdock vs AF3 vs Boltz-2 CDR RMSD, split by MHC class and
  pre- vs post-model-cutoff PDB deposit. Includes a `_lawson130` variant
  restricting to Lawson's original 130-triad subset (sanity check against the
  pre-revision manuscript numbers).
- **Panel C**: per-component AF3 RMSDs (peptide, MHC α/β, TCR V-region) on the
  same pre/post split. The class-I post drop in TCR RMSD (p=2e-04) and
  peptide RMSD (p=0.02) is consistent with the manuscript's claim that AF3
  generalizes somewhat less for class-I TCR docking than for MHC backbone.
- **Panel D**: 2×2 grid relating CDR RMSD to (a) ipTM, (b) docking-geometry
  Z-score, (c) peptide length, (d) MHC species. Spearman / Mann-Whitney
  shown in-panel; species boxes widened to 0.85.

### Figure 3 — affinity correlation
Single panel pairing the noncognate PTI-PAE distribution (red, n=31,590)
against log10-Kd-decade boxes for the cognate triads (n=63 natural + 8
engineered sub-µM with measured Kds in IEDB). Spearman: all cognate r=−0.40
(p=4.7e-4); natural-only (>1 µM) r=−0.26 (p=0.044).

### Figure 4 — validation + leakage
- **Panel a/b**: ROC curves for the AF3 PTI-PAE classifier on (a) 14 class I
  PDB validation antigens (16 cognate / 181 noncognate) and (b) 8 class II
  CRESTA antigens (205 cognate / 1228 noncognate). Median AUCs: I = 0.92,
  II = 0.81.
- **Panel c**: AUC violin distributions comparing AF3 PTI-PAE to AF-TCRdock
  across 4 strata (MHC class × antigen-centric / TCR-centric). All four
  paired Wilcoxon comparisons: AF3 beats AF-TCRdock (p ≤ 0.008).
- **Panel d**: three scatter panels showing per-antigen / per-TCR AUC against
  distance to the nearest reference in (i) peptide identity, (ii) MHC
  identity, (iii) TCRdist. **One-tailed Spearman directly testing the
  leakage hypothesis (r<0): all six p ≥ 0.58, so the data don't support a
  leakage explanation.** Class I peptide and MHC distances show
  significantly *positive* correlations under a two-sided test — i.e., AUC
  is a bit higher for antigens further from the training set, the opposite
  of leakage.
- **Panel e**: per-triad PTI-PAE vs CDR RMSD on the 16 class I validation
  triads. AF3's confidence really does track its own structural accuracy on
  held-out data (r=−0.70 on 5-seed median RMSD).

### Reviewer responses
- **`cross_reactivity/`**: Supp Fig 2. Pink (randomized) vs blue (annotated
  cross-reactive) PTI-PAE distributions binned by BLOSUM62 or Hamming
  distance to the original peptide.
- **`structural_analysis/`**: Supp Fig 3. CDR RMSD and mean per-residue
  dihedral distance vs PTI-PAE on the 16 class I validation triads.
- **`study_design_flowchart/`**: schematic of the study design.

## Colors

```python
CI_BLUE = "#3a78b8"   # MHC class I
CII_RED = "#c83737"   # MHC class II
```

## Output formats

PNG only — no PDFs anywhere in the package.

## Data provenance

The raw data behind every CSV/figure in this package lives on the TGen
**gemini** cluster, primarily in Lawson's repo checkout:

    /tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/

with the AF3 confidence/RMSD outputs and intermediate parquets under its
`data/` tree, and the IEDB affinity source database at:

    /tgen_labs/altin/alphafold3/IEDB/2025-04-15/iedb_public.db

Some derived/scratch files (cross-reactivity controls, PDB RMSD parquets)
live under `/scratch/bneff/tcrtrifold/`. Each figure's own README has a
"Data provenance" section pointing to the specific file(s) behind its CSVs —
these are pointers for verification only; the local CSVs already in this
package are everything the notebooks need to re-run.

## How to re-run

Each subdirectory has a self-contained notebook that produces the figure(s)
in that directory. The Figure 1 directory also has `make_figure1.py` as a
runnable script that produces all three panels (and their `_lawson130`
variants) in one shot.

```
cd Figure_1_structural_benchmarking && python make_figure1.py
cd Figure_3_affinity_correlation && jupyter nbconvert --execute Fig3c_*.ipynb
cd Figure_4_validation_and_leakage && jupyter nbconvert --execute Fig4_*.ipynb
cd reviewer_responses/cross_reactivity && jupyter nbconvert --execute 02_*.ipynb
cd reviewer_responses/structural_analysis && jupyter nbconvert --execute 06_*.ipynb
cd reviewer_responses/study_design_flowchart && jupyter nbconvert --execute 03_*.ipynb
```

For Figure 1, the canonical source is `make_figure1.py`; `make_figure1.ipynb`
splits the same code across cells (panel B, C, D, helpers) and is kept in
sync with the .py. For every other figure, the notebook *is* the canonical
source.

## Status

| Figure / panel | Status | Notes |
|---|---|---|
| Fig 1 b | ✅ complete | 4 staggered MW brackets (pre+post × TCRdock↔AF3 / AF3↔Boltz); one-tailed `greater`; numeric p always |
| Fig 1 c | ✅ complete | 3 pre-vs-post MW brackets per MHC class; one-tailed `less`; numeric p always |
| Fig 1 d | ✅ complete | 2×2 grid, wider species boxes (0.85), publication-weight fonts |
| Fig 3 c | ✅ complete | Compact in-panel legend/stats; non-cognate box red; decade boxes green/blue by majority class |
| Fig 4 a, b | ✅ complete | ROC curves; class I (16/181) median AUC 0.92; class II (205/1228) median AUC 0.81 |
| Fig 4 c | ✅ complete | Paired Wilcoxon: AF3 > AF-TCRdock across all 4 strata (p ≤ 0.008) |
| Fig 4 d | ✅ complete | One-tailed Spearman (r<0); all 6 p ≥ 0.58 → no leakage. Two versions: identity-distance (manuscript) + BLOSUM62-distance (alt, methodologically uniform with TCRdist). |
| Fig 4 e | ⚠ stale label | PTI-PAE vs CDR RMSD, r=−0.70 on 5-seed median; axis still says "31 - PTI-PAE" |
| Supp Fig 2 | ✅ complete | 4 PNGs, AF3-PTI-PAE label, staggered per-bin annotation boxes |
| Supp Fig 3 (RMSD) | ✅ complete | Stats moved out of axes; 0–31 x-axis fixed |
| Supp Fig 3 (dihedral) | ✅ complete | Mean per-residue dihedral distance (°) replaces TV/KDE |
| Study-design flowchart | ✅ complete | Runs as-is |
