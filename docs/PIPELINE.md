# Pipeline — end to end

How a TCR:pMHC triad becomes a number in a figure. Four arms feed the manuscript.
Each stage names the code that runs it and what it emits.

```
   STCRDab / TCR3d / IEDB / RCSB
              │
              ▼
   ┌──────────────────────┐
   │ 1. triad resolution  │   tcrtrifold-experiments (Nextflow, upstream)
   │    + filtering       │   → triad table, sequences, chain mapping
   └──────────┬───────────┘
              ▼
   ┌──────────────────────┐
   │ 2. MSA               │   MSA_WORKFLOW:RUN_MSA
   └──────────┬───────────┘
              ▼
   ┌──────────────────────┬──────────────────────┐
   │ 3a. AF3 inference    │ 3b. Boltz-2          │   ┌────────────────────┐
   └──────────┬───────────┴──────────┬───────────┘   │ 3c. AF-TCRdock     │
              │                      │               │  docs/TCRDOCK.md   │
              ▼                      ▼               └─────────┬──────────┘
   ┌──────────────────────────────────────────────────────────────────────┐
   │ 4. metrics: RMSD vs crystal · confidence (PTI-PAE) · pmhc_tcr_pae     │
   └──────────┬───────────────────────────────────────────────────────────┘
              ▼
   ┌──────────────────────┐
   │ 5. analysis + figures│   leakage · cross-reactivity · dihedral · AUC
   └──────────────────────┘
```

---

## 1–2. Triad resolution and MSAs

**Upstream: [`AltinLab/tcrtrifold-experiments`](https://github.com/AltinLab/tcrtrifold-experiments)**
(cluster checkout at `configs.paths.TRIFOLD`). Clone it; do not vendor it. **Its internals — data
layout, input schema, canonical chain order, how to run a single stage without Nextflow, and how a
new dataset gets added — are documented in [`TCRTRIFOLD.md`](TCRTRIFOLD.md).** The pipeline:

1. resolves currently-available TCR:p:MHC triads from STCRDab/TCR3d,
2. filters to αβ-TCR + classical class I/II MHC + peptide antigens,
3. **runs MSA**, then AF3 and Boltz-2 inference on each,
4. computes RMSDs against the deposited structure,
5. extracts AF3 confidence metrics.

Output parquets land in `configs.paths.TRIAD_PARQUETS` and are mirrored into `data/`.
Current set: **287 unique PDBs**.

> **MSAs are not optional for comparability.** `docs/provenance/nextflow_run.log` records
> `MSA_WORKFLOW:RUN_MSA`, `SEQ_LIST_TO_FASTA`, `FILT_FORMAT_MSA`. Any *new* prediction meant
> to sit alongside the published ones must run the data pipeline too — i.e. **not**
> `--norun_data_pipeline`. AF3's accuracy on this interface leans on the alignments.

**Expansion (revision).** `analysis/tcr3d_expansion/` answers Reviewer #1 comment #4, that
the original Figure 1 set inherited the narrower Bradley 2023 *eLife* inputs. It refreshes
STCRDab, cross-checks the TCR3d class I/II complex listings against the pipeline inputs, and
re-runs on the union: **153 → 380 triads** (pre-cutoff 130→308, post-cutoff 22→71). Steps
`00_refresh_stcrdab` → `01_identify_new_entries` → `02_apply_to_repo` → `03_run_pipeline`
→ `04_replot_figure1`.

## 3c. The AF-TCRdock baseline

Separate stack, separate recipe — **[`docs/TCRDOCK.md`](TCRDOCK.md)**. Pinned to
CUDA 11.8 / cuDNN 8.9 against an AlphaFold-2.3.2-era JAX; a newer toolkit fails *slowly*
(JAX falls back to CPU) rather than loudly. `setup_env.sh` → `submit_tcrdock_array.sh` →
`finalize.sh` → `analyze_tcrdock_vs_af3.py`. Completed 145/145 class I, 1433/1433 class II.

## 4. Metrics

| Metric | Column | Meaning |
|---|---|---|
| Docking RMSD | `docking_rmsd_af3_0` | placement of the TCR on the pMHC. **~35 Å with `tcr_rmsd` ≈ 0.5 Å means a misdock**: folded right, placed wrong |
| Component RMSDs | `cdr_`, `peptide_`, `mhc_`, `tcr_rmsd_af3_0` | per-part accuracy |
| PTI-PAE | `mean_p_tcr_interface_pae` | the paper's headline confidence feature |
| TCRdock PAE | `pmhc_tcr_pae` | the baseline's comparable score |

## 5. Analysis and figures

- **Leakage** (`analysis/triad_analysis/leakage_v*.py`) — does similarity to training data
  explain the signal? `v3` TCR3d-based, `v4` full TCRdist, `v6` BLOSUM62 triad distance.
- **Cross-reactivity controls** (`analysis/crossreact_controls/`, `01`→`03`) — scrambled and
  random peptides as negatives.
- **Structural validation** (`dihedral_analysis.py`, `extract_rama*.py`) — Ramachandran and
  dihedral checks; distinguishes *bad geometry* from *good geometry, wrong pose*.
- **Figures** — `figures/Figure_*/`, each with its own README, notebook, and data.

Most figure notebooks run **locally with no GPU**: the derived tables are tracked.
