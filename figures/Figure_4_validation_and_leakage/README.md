# Figure 4 — Validation ROCs, model comparison, and leakage analysis

Four notebooks producing the panels of the new Figure 4:

| Notebook | Panel | What it produces |
|---|---|---|
| `Fig4_ab_validation_ROCs.ipynb` | **a, b** | Per-antigen ROC curves: class I (14 antigens, 16 cognate / 181 noncognate) on left, class II (8 antigens, 205 cognate / 1228 noncognate) on right. AUCs computed from `../supplementaryTable3_updated.xlsx`. |
| `Fig4_c_AF3_vs_AFTCRdock.ipynb` | **c** | AUC violin distributions for AF3 (PTI-PAE) vs AF-TCRdock across MHC class × orientation (antigen-centric vs TCR-centric). Connecting lines between paired triads; medians + n + paired Wilcoxon p on each pair. |
| `Fig4_d_decomposed_leakage.ipynb` | **d** | Three side-by-side panels: per-antigen / per-TCR AUC vs peptide distance, MHC distance, and paired-chain TCRdist — testing whether AUC is explained by closeness to AF3's training data ("leakage"). Identity-distance version. |
| `Fig4_d_decomposed_leakage_BLOSUM62.ipynb` | **d (alt)** | Same leakage analysis but with BLOSUM62-based distances on the peptide and MHC axes. Methodologically uniform with TCRdist on the TCR axis (which is itself BLOSUM62-based). |
| `Fig4_e_PTIPAE_vs_RMSD_validation.ipynb` | **e** | AF3 PTI-PAE vs CDR RMSD for the 16 class I post-cutoff triads (5-seed median + seed-0). |

## How to run

```bash
cd Figure_4_validation_and_leakage
jupyter nbconvert --execute Fig4_ab_validation_ROCs.ipynb
# (or any of the four notebooks)
```

Run order doesn't matter — each notebook reads only its local CSVs.

## Data files

| CSV | Used by | Contents |
|---|---|---|
| `antigen_centric.csv` | 4c, 4d | Per-antigen AUC + composite pMHC distance (n=22: 14 class I + 8 class II) |
| `antigen_component_identities.csv` | **4d** | Per-antigen AUC + separate `pep_identity` and `mhc_identity` (the decomposed columns John asked for) |
| `antigen_blosum62_distances.csv` | **4d (alt)** | Per-antigen AUC + `pep_blosum62_dist` and `mhc_blosum62_dist`. Derivation cell at the bottom of `Fig4_d_decomposed_leakage_BLOSUM62.ipynb`. |
| `tcr_centric.csv` | 4c, 4d | Per-TCR AUC + paired-chain TCRdist (n=205: 13 class I + 192 class II after ANARCI filtering) |
| `per_antigen_auc_AF3.csv` | 4c | Same AUCs as `antigen_centric.csv['auc']` but recomputed directly from raw PTI-PAE scores in supp table 3 (sanity check) |
| `fig4c_tcrdock_per_antigen.csv` | 4c | Paired AF3 vs AF-TCRdock per-antigen AUC, antigen-centric (12 class I + 8 class II) |
| `fig4c_tcrdock_summary.csv` | 4c | All 4 AF3-vs-AF-TCRdock subgroups: medians, win counts, paired Wilcoxon p |
| `classII_per_TCR_AUC.csv` | 4c | 205 class II per-TCR AUCs (AF3 + AF-TCRdock) for the class II TCR-centric violin |

## Panel-by-panel notes

### Panel a, b — validation ROCs
- Class I uses a blue family (`#3a78b8`), class II red family (`#c83737`).
- Class I step curves are jittered (tiny ±0.012 random offset on both axes) so
  overlapping curves stay individually visible.
- Median AUCs: **class I = 0.92**, **class II = 0.81** (the class II number was
  0.78 in the earlier draft; recomputed against the updated supp table).

### Panel c — AF3 (PTI-PAE) vs AF-TCRdock
- Full symmetric violins for both methods, paired with connecting lines per
  triad (12 class I + 8 class II antigens, plus 205 class II per-TCR).
- Paired Wilcoxon (one-sided where appropriate): **AF3 beats AF-TCRdock in all
  four subgroups** (p ≤ 0.008).
- All four subgroups (class I/II × antigen-centric/TCR-centric) are now
  populated. The earlier "PENDING AFdock" placeholder is resolved.

### Panel d — decomposed leakage

The reviewer asked whether per-antigen / per-TCR AUC is explained by
similarity to structures AF3 could have seen in training (the "leakage"
worry). The question is: do validation antigens that are *close* to AF3's
training set get artificially-high AUCs?

There are two parallel versions of this plot, using different distance
metrics on the x-axis. Both make the same point (no leakage), and both use
**Spearman with `alternative='less'`** so the p-value directly tests the
leakage hypothesis (negative correlation between distance and AUC). A small p
(≤0.05) would support leakage; large p (≫0.05) says the data don't support
a leakage explanation.

#### Version 1 — sliding-window identity distance (`Fig4_d_decomposed_leakage.ipynb`)

This is the version that appears in the manuscript figure. Three side-by-side
panels:

- **Peptide distance** (`pep_identity` in `antigen_component_identities.csv`,
  plotted as `1 − pep_identity`): for each validation peptide, slide the
  shorter sequence ungapped along every pre-AF3-cutoff PDB peptide of the
  matching MHC class. For each register, compute identity = (matched
  residues) / (length of the shorter sequence). The maximum identity across
  all registers and references is `pep_identity`; the distance is
  `1 − pep_identity`. **0 = a near-identical peptide exists in training; 1 =
  no shared residues with anything in training.**
- **MHC distance** (`mhc_identity`, plotted as `1 − mhc_identity`): same
  recipe but on the MHC chain sequences. Because MHCs are similar in length
  across alleles, the sliding-window collapses to a near-direct comparison.
- **Paired-chain TCRdist** (`tcr_dist` in `tcr_centric.csv`): TCRdist
  (Dash et al., Nature 2017; `tcrdist3`) between the validation TCR and its
  nearest neighbour in the TCR3d reference set, summed over both chains'
  CDR1/CDR2/CDR2.5/CDR3 loops. Under the hood, TCRdist *is* BLOSUM62-based
  scoring on the CDR loops with chain-specific weighting — so this axis is
  philosophically the same kind of metric as the BLOSUM62 versions of the
  peptide/MHC axes in version 2 below.

Results:

| comparison | class | n | Spearman r | p (H₁: r < 0) |
|---|---|---|---|---|
| AUC vs peptide distance | I  | 14  | +0.62 | 0.99 |
| AUC vs peptide distance | II | 8   | +0.09 | 0.58 |
| AUC vs MHC distance     | I  | 14  | +0.50 | 0.97 |
| AUC vs MHC distance     | II | 8   | +0.12 | 0.61 |
| AUC vs TCRdist          | I  | 13  | +0.09 | 0.61 |
| AUC vs TCRdist          | II | 192 | +0.04 | 0.71 |

#### Version 2 — BLOSUM62 distance (`Fig4_d_decomposed_leakage_BLOSUM62.ipynb`)

Same analysis but with the peptide and MHC distances re-computed using
**BLOSUM62** instead of plain identity matching. The point of this version:
TCRdist (the third panel of version 1) is itself BLOSUM62-based, so doing
the peptide and MHC distances in the same currency makes the three axes
methodologically uniform.

**What is BLOSUM62?** A 20×20 substitution matrix derived from observed amino
acid replacement rates in conserved blocks of related proteins, where each
cell `B62[a, b]` is a log-odds score for substituting `a` with `b` during
evolution. Diagonal entries (`B62[a, a]`, an amino acid matched to itself)
are positive and large (A=4, W=11, etc.). Conservative substitutions are
mildly positive (I↔L=2, S↔T=1), and non-conservative substitutions are
negative (W↔E=−3). The matrix is the workhorse for sequence-similarity
scoring in protein bioinformatics, and is the substitution scoring TCRdist
uses on CDR loops.

**Why we use it here.** Plain identity matching treats every substitution as
equally distant — an A→I and an A→W look the same. BLOSUM62 weighs them by
biological similarity, so a conservative substitution looks "closer" than a
disruptive one. Using BLOSUM62 across peptide, MHC, and TCR axes makes the
three leakage panels apples-to-apples (instead of identity-vs-identity-vs-
TCRdist).

**How the distance is computed:**

For a validation query `q` and a set of reference sequences `{r₁, r₂, …}`,
the BLOSUM62 score of one register is `Σ_i B62[q_i, r_i]`. The query's
self-score is `self_score(q) = Σ_i B62[q_i, q_i]` — the score q would get
aligned to itself, no gaps, perfect match. Then:

    raw_score = max over (refs and registers) of Σ_i B62[q_i, r_i]
    distance  = 1 − raw_score / self_score(q)

This gives **distance ∈ [0, 1]**: 0 when a perfect match exists in the
reference set (raw_score = self_score), 1 when the best reference scores
zero. Normalizing by the query's self-score (rather than by length) accounts
for amino-acid composition — a tryptophan-rich sequence has higher
self-score than an alanine-rich one of the same length, and the
normalization handles that.

For peptides, "registers" means the ungapped sliding-window registers
(same as the identity version). For MHC, full-length pairwise global
alignment with affine gaps (open=−10, extend=−1) using Biopython's
`PairwiseAligner` — necessary because MHC alleles can differ in length.

**One caveat on the MHC panel.** Because MHC sequences are ~180–360 aa and
HLA polymorphism is concentrated in a few binding-groove positions, the
self-score-normalized MHC distances are very small (0–0.07) — far smaller
than the peptide distances (0–0.95). The two panels share the 0–1 x-axis,
so MHC dots all cluster at the left edge, which is meant to drive home the
point that **the MHCs are essentially all "seen"** by AF3 even though the
specific allele combinations are novel triads. The peptides have real
variation in distance to training.

Results (BLOSUM62 version):

| comparison | class | n | Spearman r | p (H₁: r < 0) |
|---|---|---|---|---|
| AUC vs peptide BLOSUM62 distance | I  | 14 | +0.51 | 0.97 |
| AUC vs peptide BLOSUM62 distance | II | 8  | −0.14 | 0.37 |
| AUC vs MHC BLOSUM62 distance     | I  | 14 | −0.23 | 0.21 |
| AUC vs MHC BLOSUM62 distance     | II | 8  | +0.13 | 0.62 |

#### What both versions say

Every p across both versions is well above 0.05, so the negative-correlation
hypothesis (leakage) isn't supported in any of the analyses. In two cases —
class I peptide identity-distance and class I peptide BLOSUM62 distance —
the observed correlation is significantly **positive**: AUC is actually a
bit *higher* for antigens further from the training set, the opposite of
leakage. The MHC and TCR axes show no signal in either direction. The
conclusion is robust to whether we use identity-distance or BLOSUM62-distance
as the leakage metric.

### Panel e — PTI-PAE vs CDR RMSD
Per-triad AF3 PTI-PAE vs CDR RMSD on the 16 class I validation triads. The
5-seed median version gives r=−0.70 — AF3's confidence really does track its
own structural accuracy on held-out data.

> ⚠️ Axis label note: panel e still says `31 - PTI-PAE` instead of the
> harmonized `AF3 PTI-PAE` label. Doesn't change any numbers; will be swept on
> a later pass when the 5-seed RMSD pipeline updates.

### Why panel c shows n=14 but panel d shows n=13 (class I, TCR-centric)
Panel c per-TCR AUCs come straight from supp table 3 — 14 unique class I
cognate TCRs. Panel d's TCR-centric x-axis is the paired-chain TCRdist, which
needs IMGT numbering of *both* chains via ANARCI; one class I TCR has a
non-standard variable domain that ANARCI couldn't number, so it drops out of
panel d (n=13) but stays in panel c (n=14). The antigen-centric counts match
(n=14 both panels) because peptide/MHC distances don't go through ANARCI.

## Data provenance — where to verify these CSVs

All CSVs come from Lawson's repo checkout on gemini:

    /tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/

- **`antigen_centric.csv`, `tcr_centric.csv`, `antigen_component_identities.csv`,
  `per_antigen_auc_AF3.csv`** — AUCs computed from the PTI-PAE scores in
  `../supplementaryTable3_updated.xlsx` (held-out validation set: PDB class I
  + CRESTA class II). Underlying AF3 confidence values trace back to the
  validation triad parquets in the repo `data/` tree.
- **`fig4c_tcrdock_per_antigen.csv`, `fig4c_tcrdock_summary.csv`,
  `classII_per_TCR_AUC.csv`** — from the AF-TCRdock run on the held-out
  antigens (the AF-TCRdock pipeline = Bradley et al. eLife 2023,
  `pmhc_tcr_pae` feature). AF-TCRdock can only model the antigens without
  cross-species KRAS chains, hence n=12 (not 14) for class I antigen-centric.
  Class I antigen- and TCR-centric counts are equal because each of those 12
  antigens has exactly one cognate TCR.

For Fig 4e, the PTI-PAE and CDR RMSD values for the 16 class I post-cutoff
triads come from supp table 3 + the PDB triad RMSD parquet
(`pdb_triad.af3_rmsd.parquet`, in the repo `data/pdb/` tree and mirrored in
`/scratch/bneff/tcrtrifold/`).
