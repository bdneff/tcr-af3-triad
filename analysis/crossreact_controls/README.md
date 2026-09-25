# cross-reactivity controls (BLOSUM62 axis)

Generates randomized "pink" peptide controls spanning BLOSUM62 distance, to
sit alongside the green (original) and blue (published cross-reactive)
peptides in the cross-reactivity analysis. This addresses reviewer 1 point 7
and reviewer 2 point 2: showing that AF3 PTI-PAE can distinguish bona-fide
cross-reactive peptides from background random peptides at matched distance,
and ideally that cross-reactive peptides score better than randoms even
*closer* than them in sequence space.

## Datasets

The two source studies (one per HLA class), each holding the TCR + MHC fixed
and varying the peptide:

| class | HLA          | peptide           | length | source CSV |
| ----- | ------------ | ----------------- | ------ | ---------- |
| I     | B*57:03      | KAFSPEVIPMF       | 11mer  | `data/I.csv`  (eLife 2020) |
| II    | DRB3*03:01   | PPQIAANRSQLISLV   | 15mer  | `data/II.csv` (PMC7541396) |

Each input CSV contains one `original` row and 21–23 published
`cross-reactive` rows. Their AF3 PTI-PAEs are already computed (green +
blue points in the planned plot). What's missing is the `randomized` set
(pink) — that's what this package generates.

## BLOSUM62 distance (the x-axis)

For two equal-length peptides `p` and `p_ref`:

    d(p, p_ref) = sum_i BLOSUM62(p_ref_i, p_ref_i) - sum_i BLOSUM62(p_i, p_ref_i)

so `d = 0` when `p == p_ref`, and `d` grows as substitutions accumulate
(weighted by BLOSUM62 dissimilarity rather than raw position differences).
This matches the convention from `leakage_v3_tcr3d.py` for manuscript
consistency.

For reference: the published cross-reactive peptides occupy `d ≈ 19–47`
(class I and II). A fully-random equal-length peptide has expected
`d ≈ 69` (class I) and `≈ 86` (class II). The randomized set spans
`d = 5` out to roughly `1.5 × expected-random`, covering everything from
"one conservative substitution" to "deep in random space".

## Pipeline overview

Three steps, runnable independently or chained:

```
data/I.csv ──┐
data/II.csv ─┴─> 01_generate_random_peptides ─> data/cross_reactivity_controls_raw.csv
                                                    │
                                                    ▼
                          02_build_pipeline_raw_csv ─> data/cross_reactivity_controls_pipeline.csv
                                                                            │
                                                                            ▼
                                  (drop into tcrtrifold-experiments, run AF3, get back .parquet)
                                                                            │
                                                                            ▼
                                      03_analyze_results ─> CSV + scatter plot PNG
```

## Step 1: generate randomized peptides locally

Pure-Python, no external dependencies. Run anywhere:

```bash
bash scripts/01_generate_random_peptides.sh
# or, with custom settings:
N_PER_CLASS=500 SEED=7 bash scripts/01_generate_random_peptides.sh
```

Output: `data/cross_reactivity_controls_raw.csv` with ~787 triads:

- 1 original + 21 cross-reactive + ~380 randomized for class I  (~402 total)
- 1 original + 23 cross-reactive + ~360 randomized for class II (~385 total)

Acceptance-rejection sampling distributes the randomized peptides roughly
uniformly across the BLOSUM62 distance axis (default 20 bins per class).

## Step 2: build the pipeline-ready CSV

Adds `name` (MD5-derived job ID) and `cognate` columns, reorders columns to
match `cross_reactivity_raw.csv` from the prior 750-triad run.

```bash
bash scripts/02_build_pipeline_raw_csv.sh
```

Output: `data/cross_reactivity_controls_pipeline.csv`.

## Step 3: run AF3 on Gemini

Drop the pipeline CSV into tcrtrifold-experiments and run the standard
pipeline. Recommended location matches the prior `cross_reactivity` job:

```bash
# on Gemini, in the tcrtrifold-experiments repo:
mkdir -p data/cross_reactivity_controls/raw
cp /path/to/this_package/data/cross_reactivity_controls_pipeline.csv \
   data/cross_reactivity_controls/raw/cross_reactivity_controls_raw.csv

# then run the same clean -> MSA -> inference -> feature-extract scripts
# you used for the prior cross_reactivity job, pointing at the new subdir.
# (cf. scripts/cross_reactivity/01_clean_and_inf.sh in the experiments repo)
```

Expected runtime is similar to the prior 750-triad run (~6 hours on Gemini
for inference, plus feature extraction).

Output of interest: `data/cross_reactivity_controls/triad/cross_reactivity_controls_triad.conf_af3.parquet`.

## Step 4: analyze and plot

Once the `.parquet` is back, run locally:

```bash
bash scripts/03_analyze_results.sh path/to/conf_af3.parquet
```

Outputs:

- `data/cross_reactivity_controls_with_PTIPAE.csv` — one row per triad with
  AF3 PTI-PAE, class label, BLOSUM62 distance, MHC class. Ready for any
  further plotting / stats.
- `data/cross_reactivity_controls_plot.png` — two-panel scatter (class I
  left, class II right) of AF3 PTI-PAE vs. BLOSUM62 distance, color-coded
  green / blue / pink, with the original's PTI-PAE marked by a dashed line.

If you want polars-backed parquet reading: `pip install polars`. If you only
have CSV output, the script accepts CSV too.

## File map

```
.
├── README.md                                    (this file)
├── data/
│   ├── I.csv                                    (class I source, KAFSPEVIPMF + 21 cross-reactive)
│   ├── II.csv                                   (class II source, PPQIAANRSQLISLV + 23 cross-reactive)
│   └── cross_reactivity_controls_*.csv          (generated; see steps)
├── src/
│   ├── blosum62.py                              (BLOSUM62 matrix + distance function)
│   └── io_utils.py                              (CSV readers, fixed-context extractor)
└── scripts/
    ├── 01_generate_random_peptides.py / .sh     (BLOSUM62-controlled peptide generation)
    ├── 02_build_pipeline_raw_csv.py / .sh       (add name + cognate columns, reorder)
    └── 03_analyze_results.py / .sh              (join + plot, runs after AF3 returns)
```

## Reproducibility

Both `01` and `02` are deterministic given seed + input CSVs. Same seed
(default `42`), same I/II.csv => bit-identical output CSVs. Re-running the
package on a different machine will produce the same 787 triads.

## What the planned plot is supposed to show

Per John's framing (1 → 3, in increasing strength):

1. **Green + blue dominate red.** Both the original and the cross-reactive
   peptides score with substantially lower AF3 mean_p_tcr_interface_pae
   (higher PTI-PAE) than the randomized peptides at any distance bin —
   AF3 can recognize cross-reactivity, not just exact cognate identity.

2. **Cross-reactive vs matched-distance randoms.** At each BLOSUM62 distance
   in the blue range, the blue cross-reactive peptides score higher PTI-PAE
   than the pink randoms at the same distance. AF3 isn't just rewarding
   sequence similarity — it's identifying specific peptides as plausible
   binders within the broader sequence neighborhood.

3. **Cross-reactive beats *closer* randoms.** The strongest result: even
   blue peptides at d=25–47 may score higher than pink peptides at d=5–20.
   If true, AF3 is identifying biologically meaningful cross-reactivity
   independent of raw sequence distance.
