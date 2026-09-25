# Structural analysis — Supp Fig 3

Two scatter plots characterizing AF3's structural prediction quality on the
16 class I post-AF3-cutoff validation triads (the ones used in Fig 4a). Both
plot AF3 PTI-PAE on the x-axis (higher = more confident interface) and a
structural-accuracy metric on the y-axis.

The point of this figure is to address Reviewer 1 Major Points #1 and #2:
- **#1**: stratify by PTI-PAE and show structural accuracy in each bin. Is
  PTI-PAE a proxy for model quality?
- **#2**: report a CDR dihedral-angle metric in addition to Cartesian RMSD
  (per North & Dunbrack 2011; McMaster & Koohy 2024).

## Files in this directory

```
06_structural.ipynb                          — single notebook that produces both plots
ramachandran_all_classI.csv                  — per-residue (φ,ψ) for all 16 triads
                                                (16 triads × ~60 CDR residues = 996 rows)
ramachandran_distances_per_triad.csv         — derived per-triad metrics
                                                (mean dihedral distance, TV, JS — 16 rows)
merged_tv_rmsd_pae.csv                       — per-triad CDR RMSD + PTI-PAE, with mhc_class
                                                (16 class I + 10 class II = 26 rows; see below)
rama_7q9a_topo.png                           — example per-triad Ramachandran (7q9a topology)
plot_rmsd_vs_pae.png                         — CDR RMSD vs PTI-PAE (Supp Fig 3, RMSD panel) — ORIGINAL 16
plot_ramachandran_vs_pae.png                 — L² mean dihedral distance vs PTI-PAE (Supp Fig 3, dihedral panel)
plot_dihedral_L1_vs_pae.png                  — L¹ mean dihedral distance vs PTI-PAE (alt metric, full 0–360° y-axis)
plot_dihedral_L1_vs_rmsd.png                 — L¹ mean dihedral distance vs CDR RMSD (two structural metrics)

# class II extension (reviewer round 2, minor point 1) — see the last section
plot_supp3c.py                               — the extended panel; --check guards class I
plot_rmsd_vs_pae_classI_II.png               — Supp Fig 3c, 16 class I + 10 class II
build_classII_rows.py                        — appends the class II rows to merged_tv_rmsd_pae.csv
classII_pti_pae.csv                          — PTI-PAE (incl. 8enh/8es9 controls + excluded 8tr*)
classII_cdr_rmsd.csv                         — CDR RMSD for the five newly folded triads
classII_summary_for_review.csv               — the 10 triads with dates + provenance, for John
make_review_pdf.py                           — builds classII_summary_for_review.pdf from the CSVs
classII_summary_for_review.pdf               — one-page summary + figure, for verification
```

## How to run

```bash
jupyter nbconvert --execute 06_structural.ipynb
```

The notebook reads the local CSVs and writes both PNGs. No external data
needed — `merged_tv_rmsd_pae.csv` and `ramachandran_distances_per_triad.csv`
already contain everything.

## Plots

### Plot 1 — CDR RMSD vs AF3 PTI-PAE
- 16 dots, one per validation triad. AF3 PTI-PAE on x (0 – 31 axis, fixed),
  CDR RMSD on y.
- PDB code annotated next to each dot.
- Stats box (n, Pearson r, Spearman ρ) sits in the figure top margin —
  outside the axes — so it doesn't overlap any data point.

### Plot 2 — Mean per-residue dihedral distance vs AF3 PTI-PAE
- Same 16 triads, same x-axis.
- y-axis = the average over CDR residues of √(Δφ² + Δψ²), with Δφ and Δψ
  wrapped into [-180°, 180°] (so a 178° → -178° crystal-to-prediction jump
  reads as 4°, not 356°). Units: degrees.
- This metric was chosen because (a) it's straightforward to define in one
  sentence, (b) it's independent of rigid-body docking (so it isolates the
  local backbone accuracy), and (c) it's the dihedral-angle measure the
  reviewers asked about.



### Plot 3 — L¹ dihedral distance vs AF3 PTI-PAE (alternative metric, full range)
- Same data, same x-axis, but y-axis is the **L¹** version of the dihedral
  distance — mean |Δφ| + |Δψ| (theoretical max = 360°).
- y-axis is fixed at 0–365 so the absolute scale of error is visible. All 16
  triads sit in the bottom ~15% of the possible range.
- Pearson r = −0.20, p = 0.47; Spearman ρ = −0.09, p = 0.73. Same story as the
  L² version (no correlation between dihedral fidelity and PTI-PAE).

### Plot 4 — L¹ dihedral distance vs CDR RMSD (two structural metrics)
- Compares the two structural-accuracy metrics directly. CDR RMSDs are AF3
  model 0 (seed-0); 7q9b was computed standalone because the parquet pipeline
  failed on it.
- Pearson r = −0.35, p = 0.19; Spearman ρ = −0.21, p = 0.44. Essentially
  uncorrelated — and the four reverse-dock failures (8gom, 8gon, 7q9b, 8es9)
  sit at high RMSD (17–24 Å) but **low** dihedral distance (~10–60° out of
  360° max). That's the headline: AF3's local backbone is right even on the
  triads where it docks the loops in the wrong place.

## Headline correlations

| comparison | n | Pearson r | p | Spearman ρ | p |
|---|---|---|---|---|---|
| CDR RMSD vs AF3 PTI-PAE | 16 | −0.53 | 0.036 | −0.60 | 0.014 |
| L² mean dihedral distance vs AF3 PTI-PAE | 16 | −0.17 | 0.54 | −0.09 | 0.75 |
| L¹ mean dihedral distance vs AF3 PTI-PAE | 16 | −0.20 | 0.47 | −0.09 | 0.73 |
| L¹ mean dihedral distance vs CDR RMSD    | 16 | −0.35 | 0.19 | −0.21 | 0.44 |

**Interpretation.** CDR RMSD tracks AF3 PTI-PAE meaningfully: a more confident
interface goes with a more accurate Cartesian docking pose. The dihedral
metric doesn't — because it measures *local* backbone fidelity, which AF3
gets right even on the triads where the global docking pose is wrong (the
"reverse-dock" failures like 8gom and 8gon have low dihedral distance but
high CDR RMSD). Together the two metrics tell complementary stories:

- **CDR RMSD**: AF3's failure mode on this set is rigid-body misplacement,
  and PTI-PAE flags those cases.
- **Mean dihedral distance**: AF3 gets local CDR torsion right even when the
  global pose is off.

## Per-triad table (sorted by AF3 PTI-PAE, high = most confident)

In `ramachandran_distances_per_triad.csv` (column `mean_dihedral_deg` is the
new metric; `tv` and `js` are the older KDE-on-(φ,ψ) metrics kept for
backward-compatibility):

| pdb | peptide | PTI-PAE | mean dihedral (°) | n CDR res |
|---|---|---|---|---|
| 8enh | LPFEKSTIM  | 28.81 | 45.5 | 66 |
| 8eo8 | LPFDKATIM  | 28.79 | 36.6 | 66 |
| 8dnt | LLLDRLNQL  | 28.77 | 45.7 | 60 |
| 7q9b | EAAGIGILTV | 28.35 | 26.9 | 58 |
| 7q9a | LLLGIGILVL | 28.32 |  7.7 | 60 |
| 8en8 | LPFDKSTIM  | 28.26 | 44.4 | 66 |
| 7q99 | NLSALGIFST | 27.93 | 12.7 | 60 |
| 8gon | RLQSLQIYV  | 25.05 |  9.1 | 63 |
| 8i5d | VVGAVGVGK  | 24.88 | 27.0 | 61 |
| 8wul | VVGAVGVGK  | 24.46 | 38.0 | 63 |
| 8i5c | VVGAVGVGK  | 23.86 | 30.5 | 59 |
| 8ye4 | NYNYLYRLF  | 23.59 | 37.1 | 62 |
| 8f5a | TSTLQEQIGW | 22.93 | 46.4 | 62 |
| 8gom | RLQSLQTYV  | 22.63 | 11.6 | 63 |
| 8wte | VVGAVGVGK  | 22.33 | 32.6 | 63 |
| 8es9 | GVYDGREHTV | 16.78 | 48.8 | 64 |

(8gom and 8gon — the two reverse-dock failures — sit at PTI-PAE ≈ 23 and 25
but have **low** dihedral distance: AF3 modeled the loops right, it just put
them in the wrong place.)

## Data provenance — where this came from on gemini

Generated as part of the "Dihedral angle plots for protein structure
validation" analysis on May 27:

```
# Inputs (on gemini):
/scratch/bneff/tcrtrifold/extract_rama_all.py            → ramachandran_all_classI.csv
/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/data/pdb/triad/staged/
  pdb_validation_triad.af3_rmsd.parquet                  → cdr_rmsd_af3_0
  pdb_validation_triad.conf_af3.parquet                  → mean_p_tcr_interface_pae
```

The CSVs in this directory are mirrors of the gemini outputs; everything the
notebooks need is already local.

## Class II extension (reviewer round 2, minor point 1)

The reviewer accepted the PTI-PAE/accuracy relationship but noted it rests on only 16 post-training
structures, and asked for a larger benchmark if more are available. Rather than expand the *class I*
set — which is load-bearing across Fig 1, Fig 4a and the pre/post-cutoff splits, so growing it would
propagate through the whole paper — the answer adds **class II** structures, which touch nothing
else because class II validation elsewhere uses CRESTA.

**`plot_supp3c.py`** produces the extended panel. It is a script rather than a new notebook cell so
the figure has a named, runnable command; `06_structural.ipynb` is **kept unedited** as the record
of the original 16-point version, per this repo's rule that superseded artifacts are archived rather
than overwritten.

```bash
cd figures/reviewer_responses/structural_analysis
python plot_supp3c.py --check      # -> plot_rmsd_vs_pae_classI_II.png
```

Everything about the class I panel is deliberately unchanged — same CSV, same `31 - pti_pae`
transform, same fixed 0–31 axis, same `#3a78b8`, same margin stats box — so a reader comparing the
two sees the original sixteen points unmoved with new ones added, not a redrawn figure. Class II
uses `#A32D2D`, the class II red already used throughout Figure 1, so the two figures agree about
what red means.

`--check` is a regression guard: the class I subset must still give **n=16, Pearson −0.53,
Spearman −0.60**. If adding rows ever perturbs those, the script exits non-zero rather than quietly
publishing a moved result. It currently passes.

`merged_tv_rmsd_pae.csv` gained an **`mhc_class`** column (all original rows `I`). Three
correlations are reported — class I alone, class II alone, and combined — because the reviewer will
want to see that the original class I claim did not move.

### Where the class II numbers come from

**`build_classII_rows.py`** appends them. Both metrics are taken from the **published pipeline's own
tables** and neither is recomputed here:

| column | source |
|---|---|
| `cdr_rmsd` | `cdr_rmsd_af3_0` in `data/exports/af3_rmsd.csv` (mirror of `pdb_triad.af3_rmsd.parquet`) — or, for the five newly folded, `classII_cdr_rmsd.csv` from `hla2_rmsd.nf`, which runs the same four processes as `pdb.nf`. A PDB in both is an error. |
| `pti_pae` | `mean_p_tcr_interface_pae` from `workflows/bin/extract_triad_conf_feat.py` |

`cdr_rmsd_af3_0` is **seed index 0**, which is how the original 16 class I points were built (see
"Data provenance" above). The per-seed columns must not be averaged: `pred_sample_rank_N` ties —
`8pjg` has ranks `1,2,3,3,3` with three byte-identical RMSD rows — so a mean weights one model three
times.

```bash
python build_classII_rows.py --pae classII_pti_pae.csv
python plot_supp3c.py --check
```

Guards: every requested triad lands in exactly one bucket (added / no RMSD / no PTI-PAE) and the
buckets must sum to the number requested; the 16 class I rows are diffed and the script aborts
rather than write a file whose published half moved; an already-present triad is refused rather
than duplicated.

### Final result (2026-09-09)

| | n | Pearson r | p | Spearman ρ | p |
|---|---|---|---|---|---|
| class I (unchanged) | 16 | −0.53 | 0.036 | −0.60 | 0.014 |
| **class II** | **10** | **−0.75** | **0.012** | **−0.73** | **0.016** |
| combined | 26 | −0.50 | 0.010 | −0.59 | 0.001 |

Class II independently reproduces the class I relationship, and more strongly. Every class II triad
is post-cutoff (2023-01-12; earliest release is `8pjg` on 2024-06-19). Five were already folded in
the published `pdb.nf` run (`8pjg` `8vcx` `8vcy` `8vd0` `8vd2`) and needed only PTI-PAE. Five were
folded new (`8vq8` `9aud` `9yaf` `9yag` `27eb`). The PTI-PAE extraction was checked against
published values with `8enh` and `8es9`, which are the two class I rows in `classII_pti_pae.csv`.

**`27eb` is the point to understand.** It has the lowest confidence (PTI-PAE 17.5) and the worst
accuracy (22.2 Å) in class II, and it sits with the class I reverse-dock failures `8gom`/`8gon`.
PTI-PAE picked it out before its RMSD was computed. It is also the one triad not on John's TSV.

**Excluded: `8trl` and `8trr`** are citrullinated epitopes, and AF3 folds citrulline as arginine.
They are excluded on chemistry, not quality (both predicted well), as John decided. See
`EXCLUDED_CITRULLINATED` in `build_classII_rows.py`. Taking them out made every statistic stronger.

The full pipeline account, the run order, and the problems hit along the way are in
[`../../../analysis/classII_supp3c/README.md`](../../../analysis/classII_supp3c/README.md). The
standalone staging built first is archived under `analysis/classII_supp3c/superseded/` and is **not
the source of any row here**. Its AF3 JSONs put the peptide at chain C (the opposite convention)
instead of the pipeline's chain A.
