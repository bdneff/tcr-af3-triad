# Study-design flowchart

The clean version of the study-design flowchart used in the response-to-
reviewers. Three lanes:

1. **Top — HELD-OUT VALIDATION SETS** — PDB (orange/brown, class I + II) and
   CRESTA (purple, class II only). Kept untouched during classifier
   development; used only for the final test in Fig 4.
2. **Middle — TRAINING SOURCE** — IEDB (blue), shown as a 5-step pipeline
   ending in PTI-PAE feature selection.
3. **Bottom — VALIDATION** — AF3 PTI-PAE classifier (dark navy) → PDB AUC 0.92
   (Fig 4a) → CRESTA AUC 0.78 (Fig 4b).

Dashed arrows show: (a) held-out sets consulted to filter training (top →
middle), and (b) validation triads scored by the trained classifier (top →
bottom, L-shaped around the outside of IEDB).

## Files

```
03_flowchart.ipynb            — single notebook that produces the PNG
study_design_flowchart.png    — pre-rendered output
```

## How to run

```bash
jupyter nbconvert --execute 03_flowchart.ipynb
```

No external data — all counts are hardcoded from the manuscript Methods
section.

## Hardcoded counts (from Methods)

| Step | Class I | Class II | Total |
|---|---|---|---|
| Start (all eligible IEDB human triads) | 8,104 | 1,688 | 9,792 |
| After Filter 1 (drop AF3-trained look-alikes) | 7,711 | 1,625 | 9,336 |
| After Filter 2 (drop validation look-alikes) | 7,702 | 1,625 | 9,327 |
| With non-cognate controls | 7,702 / 157,025 | 1,625 / 33,628 | 9,327 / 190,653 |

Validation sets:
- **PDB (Fig 4a)**: 16 class I cognate + 181 non-cognate; 14 antigens; median AUC 0.92
- **CRESTA (Fig 4b)**: 205 class II cognate + 1,228 non-cognate; 8 antigens; median AUC 0.78

## Color palette

```python
# PDB lane (orange/brown family)
PDB_DARK = "#8a4116"   PDB_MID = "#c46a1f"

# CRESTA lane (purple family)
CR_DARK  = "#5b2a7e"   CR_MID  = "#7a3da3"

# IEDB lane (blue family)
IEDB_DARK = "#1f3a68"  IEDB_MID = "#2e5da0"

# Non-cognate step (step 4)
NEG_DARK = "#a23a3a"

# Classifier
CLF_DARK = "#1c2735"
```

These intentionally differ from the manuscript's class I/II palette
(`#3a78b8` / `#c83737`) — the flowchart is about data **provenance**, not
class membership, so it uses a separate scheme aligned to the three data
sources.

## To tweak

The notebook uses a 160 × 100 coordinate system. Major levers:
- `pdb_box_*` and `cr_box_*` control the top row geometry
- `iedb_x, iedb_y, iedb_w, iedb_h` plus `step_xs, step_w, step_h` control the
  middle row
- `clf_x, val1_x, val2_x` control the bottom row
- The four arrows are: (a) step → step inside IEDB, (b) dashed held-out →
  IEDB filter, (c) IEDB step 5 → AF3 classifier (L-shape going down then
  left), (d) dashed AF3 outputs → held-out sources (L-shape around the
  outside).
