# Triad analysis — the main AF3 arm

The core analysis behind the paper's claim that prediction *confidence* tracks binding
*specificity*, plus the controls that defend it against the obvious objection: that the
signal is memorization of training data rather than a learned interaction.

Predictions themselves are produced upstream by
[`AltinLab/tcrtrifold-experiments`](https://github.com/AltinLab/tcrtrifold-experiments)
(see `docs/PIPELINE.md`). What lives here is everything done *to* those predictions.

---

## Leakage — is the signal just memorization?

The central worry: AF3 saw many of these complexes during training, so high confidence on a
cognate triad might reflect recall rather than modelled physics. Each version asks the same
question with a stricter similarity metric — **does per-antigen/per-TCR AUC degrade as
validation triads get further from the reference set?** If it does not, similarity is not
what is driving the signal.

| Script | Reference set | Similarity metric |
|---|---|---|
| `leakage_v3_tcr3d.py` | 308 pre-AF3-cutoff triads (226 class I + 82 class II) from TCR3d | peptide: sliding-window Hamming; TCR: **CDR3-only** TCRdist (BLOSUM62 + length penalty, weight 3 on CDR3), reimplemented locally so no V-gene call is needed |
| `leakage_v4_fulltcrdist.py` | same | **full TCRdist** via the published `tcrdist3` — CDR1+CDR2+CDR2.5+CDR3 weighted per Dash et al. 2017. Needs V-gene names, obtained through ANARCI's germline assignment (BLAST against IMGT germlines) |
| `leakage_v5_directcrs.py` | direct cross-reactivity framing | |
| `leakage_v6_blosum_triad.py` | same | **one uniform metric across all three components** — BLOSUM62 + length-gap distance with sliding-offset registration. Three panels: per-antigen AUC vs peptide distance, vs MHC distance (class II sums α+β), and per-TCR AUC vs TCR distance |

`v6` is the version to read first: applying a single metric to peptide, MHC, and TCR alike
removes the objection that each component was scored on its own favourable scale. The
earlier versions are kept because figures in the submitted manuscript cite them.

Outputs: `leakage_fig*.png`, `antigen_centric*.csv`, `tcr_centric*.csv`
(`_unstrat` = not stratified by class).

## Structural validation — good geometry, or good geometry in the wrong place?

A low RMSD can hide a badly-formed backbone, and a well-formed backbone can sit in entirely
the wrong pose. These separate the two:

- `extract_rama.py`, `extract_rama_all.py` — per-residue CDR φ/ψ for class I cognate
  validation triads → `ramachandran_all_classI.csv`
- `dihedral_analysis.py` — chain matching by *sequence* rather than by chain ID, with an
  ANARCI mouse fallback (the naive chain-ID join silently mismatches on non-human systems)
  → `dihedral_validation_classI{,_v2}.csv`
- `rmsd.py` — docking / CDR / peptide / MHC / TCR RMSDs against the deposited structure

The instructive case is **8gom / 8gon**: high interface PAE but *low* dihedral distance —
AF3 modelled the loops correctly and placed them wrongly.

## Supporting

- `dump_per_component_identities.py` → `per_component_identities.csv` — sequence identity per
  triad component, the input to the leakage plots
- `pmhc_distance_patch.py`, `apply_patch.py` — patches applied to the upstream pipeline
- `study_design_figure_v9.py` — the study-design schematic
- `feat_auc_weighted.ipynb` — feature-weighted AUC exploration

## Data here

`pdb_triad.af3_rmsd.parquet`, `pdb_triad.boltz_rmsd.parquet`,
`pdb_validation_triad.annotated.parquet` (pre/post-cutoff annotated),
`tcr3d_triad_sequences.csv`, `iedb_kd_records.csv`, `ternary_templates_v2.tsv`,
`supplementaryTable3*.xlsx`.

> **`tcr3d_triad_sequences.csv` has two known defects** — a corrupt peptide for 3D3V/3D39
> (`peptide_seq` drops the P5 anchor Phe; the table's own `epitope_listed` disagrees) and a
> faked 8-mer for 7BYD's lipopeptide. Details in `docs/PROVENANCE.md`. Take sequences from
> the crystals where it matters.

Reading the parquets needs `pyarrow`; without it pandas raises "Unable to find a usable
engine" rather than a missing-package error. CSV mirrors are in `data/exports/`.

## Versioned scripts

`_v2`…`_v6` suffixes are successive revision passes, not dead code — an earlier figure may
still cite an earlier version. Nothing here is superseded silently; see `docs/REORG.md`.
