# Superseded: do not use as a source for any figure

This is the standalone staging path for the class II triads, built before John's TSV was recognised
as already being in the tcrtrifold pipeline's input schema. The Supp Fig 3c numbers come from
`../pipeline/`, not from here. The full account is in `../README.md`, under "`superseded/`".

In short:

- `af3_inputs/*.json` use the **wrong chain order**, the opposite convention (A = MHC α, B = MHC β,
  C = peptide), where the pipeline uses A = peptide.
- `af3_array.sbatch` picks its JSON by **position**, so adding a file silently shifts every array
  index.
- `scripts/04_prepare_structures.py` is a bespoke crystal cleaner and was never adopted.
- **Paths inside these files no longer resolve.** They were written for
  `analysis/classII_supp3c/` before the move, and are kept as a record, not as something to run.

Kept rather than deleted because older commits and notes cite them, and because the ten predictions
this path produced on Gemini (job 37446788) remain an independent cross-check.
