# Cross-reactivity — Supp Fig 2

Class I cross-reactivity panel for the response-to-reviewers. One fixed TCR
recognizing the cognate peptide KAFSPEVIPMF (HLA-B\*57:03), scored against:
- the **cognate** (n=1) — original peptide, green dot
- 21 experimentally-validated **cross-reactive** peptides — blue
- 4,537 sequence-distance-matched **randomized** controls — pink

All scored by AF3, plotted as **AF3-PTI-PAE** = 31 − raw `mean_p_tcr_interface_pae`
(higher = more confident interface = more cognate-like).

## Files in this directory

```
02_cross_reactivity.ipynb                          — single notebook produces all 4 PNGs
classI_results.csv                                 — source data (4,559 rows)
plot_cross_reactivity_violin_BLOSUM62.png          — binned violin (BLOSUM62 bins of width 5)
plot_cross_reactivity_violin_hamming.png           — binned violin (Hamming bins of width 1)
plot_cross_reactivity_scatter_BLOSUM62.png         — un-binned scatter + LOWESS through randomized cloud
plot_cross_reactivity_scatter_hamming.png          — same for Hamming distance
```

## How to run

```bash
jupyter nbconvert --execute 02_cross_reactivity.ipynb
# or just open the notebook and run all
```

The notebook contains two plotting functions (`plot_violin`, `plot_scatter`),
each called twice — once for BLOSUM62 and once for Hamming.

## Which reviewer point this addresses

Reviewer 1 Major #7 (deep mutational scanning of peptide variants for a fixed
TCR) and the editor's cross-reactivity mandate. The fixed-TCR design directly
maps onto the "many pMHCs : one TCR" framing the reviewers asked about — for
one TCR, can AF3 distinguish actual cross-reactive peptides from sequence-
distance-matched random controls?

## Source data

`classI_results.csv` (4,559 rows):

| column | type | notes |
|---|---|---|
| `peptide` | str | the 11-mer peptide sequence |
| `peptide_source` | str | one of `original` (n=1), `cross-reactive` (n=21), `randomized` (n=4537) |
| `blosum62_distance_from_original` | int | BLOSUM62 distance from KAFSPEVIPMF |
| `hamming_distance_from_original` | int | Hamming distance (number of substitutions) |
| `mean_p_tcr_interface_pae` | float | raw AF3 interface PAE — **lower = more cognate-like** |

The notebook inverts the PAE column to `AF3-PTI-PAE = 31 − raw` so that higher
= better throughout.

Original peptide has raw PAE = 2.11 → AF3-PTI-PAE ≈ 28.9. Cross-reactive
median raw PAE = 2.42; randomized median ~12, with a rising trend in raw PAE
as BLOSUM62 distance grows.

## Plot details

### Violins
- BLOSUM62 panel uses bin width = 5, centered on tick labels (5, 10, 15, ..., 90),
  so bin "45" = [42.5, 47.5). Hamming panel uses bin width = 1, centered on
  integers.
- Original (distance = 0) is shown separately as a green dot at the cognate
  PTI-PAE level, with a dashed horizontal reference line.
- Pink (randomized) violins drawn where ≥5 randomized peptides fall in the bin.
- Blue (cross-reactive) violins drawn where ≥3 cross-reactive peptides fall in
  the bin (so singleton or doubleton bins don't get spurious violin shapes).
- **All 21 cross-reactive peptides are shown as dots regardless of bin
  membership** (scatter threshold is min_n=1), so the singletons at BLOSUM62
  distances 25 and 30 and the doubleton at Hamming distance 9 still appear
  even though they don't get a violin. Blue dots use widened jitter
  (0.85 × violin width in x, ±0.35 in y) so that overlapping points separate
  visually.
- **Per-bin Mann-Whitney annotations** (alt='less', testing blue raw-PAE <
  pink raw-PAE at matched distance, i.e. cross-reactive more cognate-like
  than randomized) shown only where both groups have ≥3 points. Annotation
  boxes are staggered vertically so adjacent bins stay readable. r-values
  were dropped per John's last comment — only p is shown.

### Scatter
- Randomized: small low-alpha pink with tiny x-jitter (BL62 distances are
  integer-valued, so overplotting at each integer would be massive).
- LOWESS smoother through the randomized cloud (frac=0.35) to show the rising
  background trend without binning.
- Cross-reactive: larger edged blue dots.
- Original: large green dot with edge + horizontal reference line.

## Headline numbers (per-bin Mann-Whitney, blue vs pink)

| BLOSUM62 bin | n_blue | p          |
|--------------|--------|------------|
| 35 (= [32.5, 37.5)) | 4 | 2 × 10⁻⁷ |
| 40 (= [37.5, 42.5)) | 9 | 6 × 10⁻⁶ |
| 45 (= [42.5, 47.5)) | 6 | 1 × 10⁻⁴ |

| Hamming bin | n_blue | p       |
|-------------|--------|---------|
| 6           | 3      | 0.047   |
| 7           | 7      | 2 × 10⁻⁵ |
| 8           | 8      | 1 × 10⁻⁶ |

Spearman across the full randomized cloud: ρ = +0.52, p < 10⁻¹⁰, n = 4537
(raw PAE rises with BLOSUM62 distance — AF3 reads the sequence-similarity
signal even on non-binders). Spearman across cross-reactive peptides: ρ = +0.28,
p = 0.22, n = 21 (no clear trend, as expected — they all bind).

The story: at matched BLOSUM62 distances, cross-reactive peptides cluster
*well above* the randomized cloud on the AF3-PTI-PAE axis. AF3 distinguishes
genuine cross-reactivity from sequence-similarity effects alone.

## Suggested manuscript text

> "We applied AF3 to 4,537 randomized control peptides at BLOSUM62 distances
> 5-90 from the cognate (KAFSPEVIPMF), and to 21 experimentally-validated
> cross-reactive peptides at BLOSUM62 distances 25-47 from the same cognate.
> AF3 PTI-PAE fell with peptide distance from the cognate among randomized
> controls (Spearman ρ = +0.52 in raw PAE; p < 10⁻¹⁰), confirming that AF3
> captures the sequence-similarity signal. Importantly, at matched BLOSUM62
> distances, cross-reactive peptides had significantly higher PTI-PAE than
> randomized controls (Mann-Whitney p ≤ 10⁻⁴ at BLOSUM62 distances 35, 40, and
> 45; rank-biserial r = +0.86 to +0.99), demonstrating that AF3 distinguishes
> genuine cross-reactivity from sequence-similarity effects alone."

## Data provenance

Source AF3 outputs and the scan inputs live on gemini at:

    /scratch/bneff/tcrtrifold/crossreact-controls/        (also crossreact-controls.zip)

with the broader pipeline at `/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/`.
The `mean_p_tcr_interface_pae` column in `classI_results.csv` is the raw
interface PAE; the notebook inverts it to AF3-PTI-PAE = 31 − raw.
