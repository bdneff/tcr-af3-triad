# Figure 3c — AF3 PTI-PAE vs TCR:antigen affinity (Kd)

`Fig3c_affinity_correlation.ipynb` builds the panel from `fig3c_data.csv` (one
row per triad: `group`, `mhc_class`, `Kd_M`, `pti_pae_raw`). The plot shows
**AF3 PTI-PAE** (high = more confident interface = predicted to bind tightly)
vs measured Kd on a log10 axis.

## Result (matches the manuscript)
- 67 class I + 4 class II cognate triads with measured Kd.
- 8 engineered sub-µM TCRs (blue), 63 natural TCRs (green), 31,590 matched
  non-cognate controls (red box on the left).
- Spearman (AF3 PTI-PAE vs Kd):
  - **all cognate**: r = −0.40, p = 4.7e-4
  - **natural (>1 µM only)**: r = −0.26, p = 0.044
  - Negative because tighter binders (lower Kd) get higher PTI-PAE — exactly
    the right direction.

## Plot layout (per John's revision feedback)
- Decade boxes span the **full** log10(Kd) interval (10^e to 10^(e+1)) instead
  of a fixed narrow width.
- The non-cognate box on the left is red-outlined to match its legend entry.
- Decade boxes are colored by the majority population in that decade —
  green where natural TCRs dominate, blue where engineered sub-µM dominate.
- Legend + Spearman stats are stacked in the empty bottom-right corner of the
  cognate panel (no wasted right margin); figure width is 7.0" (down from
  8.6" before).

## Files

- `Fig3c_affinity_correlation.ipynb` — runnable; reads `fig3c_data.csv`.
- `fig3c_data.csv` — the extracted plotting data.
- `Fig3c_affinity_correlation.png` — pre-rendered figure.
- `corr_spr_upstream.ipynb` — Lawson's original derivation pipeline (full DB
  joins; needs the repo + the gemini IEDB DB to run, not needed to reproduce
  the figure).

## How to run

```bash
jupyter nbconvert --execute Fig3c_affinity_correlation.ipynb
# or just open it
```

## Data provenance — where `fig3c_data.csv` comes from

Extracted on gemini by joining three sources in Lawson's repo checkout
`/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/`:

| piece of the CSV | source file (under repo `data/`) |
|---|---|
| `pti_pae_raw` (cognate + non-cognate) | `data/iedb_I/triad/iedb_I_triad.conf_af3.parquet`, `data/iedb_II/triad/iedb_II_triad.conf_af3.parquet` (joined to `…receptor_reference.parquet`) |
| assay-type join key | `data/iedb_meta/assay_type.parquet` |
| PDB-overlap exclusion | `data/iedb_{I,II}_full/triad/iedb_{I,II}_triad.annotated.parquet` |
| `Kd_M` | `data/iedb_meta/triad_affinity.parquet` |

`triad_affinity.parquet` is in turn derived from the IEDB database dump:

    /tgen_labs/altin/alphafold3/IEDB/2025-04-15/iedb_public.db   (Kd in nM; ÷1e9 → M)

via `data/iedb_meta/extract_affinity.py`. The join logic is reproduced in
`corr_spr_upstream.ipynb`.
