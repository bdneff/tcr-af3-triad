# TCR:pMHC triad prediction with AlphaFold 3 — manuscript companion

Code, data, and documentation behind:

> **General prediction of T cell receptor antigen specificity from sequence using AlphaFold 3**
> Lawson J. Woods, **Brandon Neff**, Kamel Lahouel, Antoine Goursaud, Mete Mulazimoglu,
> Kameron Bates, Cristian Tomasetti, John A. Altin
> *(TGen / City of Hope / University of Sydney — accepted, Science Advances)*
>
> Preprint: [doi:10.64898/2026.06.02.729478](https://doi.org/10.64898/2026.06.02.729478) ·
> Archived release of this repository: [doi:10.5281/zenodo.22967485](https://doi.org/10.5281/zenodo.22967485)

This repository exists so a reviewer or collaborator can go from **a claim in the
paper** to **the code and data that produced it**, without cluster access. It is
organized around the manuscript's figures rather than around our directory
history.

---

## Start here

| I want to… | Go to |
|---|---|
| Reproduce a figure | [`figures/`](#figure-map) — each figure has its own README, notebook, and data |
| Understand the whole pipeline end to end | [`docs/PIPELINE.md`](docs/PIPELINE.md) |
| Re-run the AF-TCRdock baseline | [`docs/TCRDOCK.md`](docs/TCRDOCK.md) — **the hard-won recipe** |
| Understand or re-run the upstream prediction pipeline | [`docs/TCRTRIFOLD.md`](docs/TCRTRIFOLD.md) — its layout, input schema, chain convention, and how to run a stage without Nextflow |
| Know where a data file came from | [`docs/PROVENANCE.md`](docs/PROVENANCE.md) |
| Work on the cluster | [`docs/ENVIRONMENT_gemini.md`](docs/ENVIRONMENT_gemini.md) |
| See what moved during reorganization | [`docs/REORG.md`](docs/REORG.md) |

> **Code lives in git; bulk data and model weights live on the cluster.**
> Everything ignored is regenerable, third-party, or too large — and each item is
> located through `configs/paths.py`, the single file to edit per environment.
> The excluded bulk totals ~28 GB; what is tracked here is ~40 MB.

---

## What the paper does

AF3 is applied to MHC:peptide:TCR triads at scale — >9,000 TCRs against >1,000
distinct epitopes across >70 class I and class II alleles — to ask whether
*structure prediction confidence* can stand in for *binding specificity*. The
headline metric is the **PTI-PAE** (predicted TCR-interface PAE): cognate triads
are predicted more confidently than non-cognate ones, and a model built on those
features reaches median AUC 0.81–0.92 on triads unseen by AF3 and unused in
feature selection.

The claims that needed the most defending, and where they live:

1. **AF3 predicts triad structures well** — benchmarked against deposited
   crystals, and against Boltz-2 and AF-TCRdock (`Figure_1`, `Figure_4c`).
2. **The signal is not memorization** — pre/post AF3-training-cutoff splits, and
   sequence-similarity leakage analyses (`Figure_4`, `analysis/triad_analysis/`).
3. **The signal is not trivially explained by peptide chemistry** — cross-reactivity
   negative controls with scrambled/random peptides (`analysis/crossreact_controls/`).

---

## Figure map

Each directory is self-contained: notebook, script, its own `data/`, and a README
explaining what the panel shows and how it was made.

| Directory | Manuscript element | Contents |
|---|---|---|
| `figures/Figure_1_structural_benchmarking/` | Fig. 1 | AF3 vs Boltz-2 vs AF-TCRdock RMSDs; pre/post-cutoff generalization (287 unique PDBs) |
| `figures/Figure_3_affinity_correlation/` | Fig. 3c | Prediction confidence vs measured affinity (SPR/IEDB Kd) |
| `figures/Figure_4_validation_and_leakage/` | Fig. 4a–d | Validation ROCs, AF3 vs AF-TCRdock AUCs, leakage/similarity controls |
| `figures/reviewer_responses/` | revision | `cross_reactivity/`, `structural_analysis/` (Ramachandran + dihedral validation; **Supp Fig 3c class II extension**), `study_design_flowchart/` |

---

## Analysis arms

Three independent arms feed the figures. Each keeps its own scripts, inputs, and
result tables.

### `analysis/triad_analysis/` — the main AF3 triad analysis
The core tables (`pdb_triad.af3_rmsd.parquet`, `pdb_triad.conf_af3.parquet`,
`tcr3d_triad_sequences.csv`) plus the leakage series (`leakage_v3_tcr3d.py`,
`leakage_v4_fulltcrdist.py`, `leakage_v6_blosum_triad.py`), dihedral/Ramachandran
validation, and per-component identity extraction.

Predictions themselves come from the **upstream Nextflow pipeline**,
[`AltinLab/tcrtrifold-experiments`](https://github.com/AltinLab/tcrtrifold-experiments)
(cluster checkout: `/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/`),
which resolves triads from STCRDab/TCR3d, filters to αβ-TCR + classical MHC,
runs MSA → AF3 and Boltz-2 inference, computes RMSDs, and extracts confidence
metrics. It is upstream code — clone it, do not vendor it.

### `analysis/tcrdock_validation/` — the AF-TCRdock baseline
The comparison method for Figure 4c. Fiddly to build because it is pinned to an
AlphaFold-2.3.2/CUDA-11.8 stack; **[`docs/TCRDOCK.md`](docs/TCRDOCK.md) is the
working recipe**, including the CUDA pin, the fine-tuned weights, the
`--benchmark` flag that withholds the answer, and the completeness check that
catches a silently truncated result set. Runs completed whole: 145/145 class I,
1433/1433 class II.

### `analysis/crossreact_controls/` — cross-reactivity negative controls
A numbered pipeline (`01_generate_random_peptides` → `02_build_pipeline_raw_csv`
→ `03_analyze_results`) testing whether confidence separates cognate from
deliberately non-cognate peptides.

### `analysis/tcr3d_expansion/` — the Figure 1 expansion (revision)
Answers **Reviewer #1 comment #4**, that the original Figure 1 set was built
around the Bradley 2023 *eLife* inputs and was therefore narrower than the
available structural data. Refreshes STCRDab, cross-checks TCR3d's class I and II
complex listings against the pipeline's inputs, and re-runs inference on the
union: **153 → 380 triads** (pre-cutoff 130→308, post-cutoff 22→71). Numbered
`00_refresh_stcrdab` → `01_identify_new_entries` → `02_apply_to_repo` →
`03_run_pipeline` → `04_replot_figure1`.

### `analysis/classII_supp3c/` — the Supp Fig 3c class II extension (revision 2)
Answers **Reviewer round 2, minor point 1**, which asked for a benchmark larger than
the 16 post-cutoff class I structures behind Supp Fig 3c. Adds **10 class II** triads
without changing the class I ones: class II alone gives r = −0.75, ρ = −0.73. Five had
already been folded in the published run, and five were folded new through the
pipeline's own processes (`pipeline/hla2.nf`, `pipeline/hla2_rmsd.nf`). The two
citrullinated epitopes are excluded because AF3 folds citrulline as arginine. The
figure and a review PDF are in `figures/reviewer_responses/structural_analysis/`.

---

## Layout

```
.
├── README.md                  ← you are here
├── configs/paths.py           # the ONLY file with absolute cluster paths
├── docs/
│   ├── PIPELINE.md            # end to end: triads → inference → metrics → figures
│   ├── TCRDOCK.md             # ★ how the AF-TCRdock baseline was actually made to run
│   ├── PROVENANCE.md          # every data file: origin, cluster path, regeneration
│   ├── ENVIRONMENT_gemini.md  # partitions, modules, containers, weights
│   └── REORG.md               # what moved where, and what was deliberately dropped
├── figures/                   # one directory per manuscript figure, each with a README
├── analysis/
│   ├── triad_analysis/        # main AF3 analysis: leakage, dihedral, identities
│   ├── tcrdock_validation/    # AF-TCRdock baseline (Fig 4c)
│   ├── crossreact_controls/   # cross-reactivity negatives
│   └── classII_supp3c/        # Supp Fig 3c class II extension (revision 2)
├── data/                      # small tracked inputs (sequence tables, metric parquets)
└── structures/                # crystal PDBs used across arms
```

---

## Reproducing

Most figures need **no cluster and no GPU** — the derived metric tables are
tracked, so the notebooks run locally:

```bash
git clone <this repo> && cd tcr
python -m pip install -r requirements.txt
jupyter lab figures/Figure_1_structural_benchmarking/make_figure1.ipynb
```

Regenerating the *predictions* does need the cluster: see
[`docs/PIPELINE.md`](docs/PIPELINE.md) for AF3/Boltz-2 (upstream Nextflow)
and [`docs/TCRDOCK.md`](docs/TCRDOCK.md) for the baseline.

---

## Running on another system

The figure notebooks need only the tracked tables and run anywhere. Regenerating
predictions or baselines needs the original compute environment's data and model paths.

- Environment-specific locations are collected in **[`configs/paths.py`](configs/paths.py)**:
  edit that one file to point at your own copies of the model weights, databases and
  prediction outputs.
- In this release, **16 scripts still hardcode those locations** directly rather than
  reading `configs/paths.py`. Update the paths in these files as well:

- `analysis/classII_supp3c/pipeline/01_stage_triad_parquet.py`
- `analysis/classII_supp3c/superseded/af3_array.sbatch`
- `analysis/tcrdock_validation/analyze_classII.sh`
- `analysis/tcrdock_validation/analyze_classI_v2.sh`
- `analysis/tcrdock_validation/build_tcrdock_input_classII.py`
- `analysis/tcrdock_validation/finalize_classII.sh`
- `analysis/tcrdock_validation/submit_classII_array.sh`
- `analysis/triad_analysis/apply_patch.py`
- `analysis/triad_analysis/dihedral_analysis.py`
- `analysis/triad_analysis/dump_per_component_identities.py`
- `analysis/triad_analysis/extract_rama.py`
- `analysis/triad_analysis/extract_rama_all.py`
- `analysis/triad_analysis/leakage_v3_tcr3d.py`
- `analysis/triad_analysis/leakage_v4_fulltcrdist.py`
- `analysis/triad_analysis/leakage_v5_directcrs.py`
- `analysis/triad_analysis/leakage_v6_blosum_triad.py`

---

## Conventions

- **Every figure is regenerated by a named script**, and its provenance recorded
  in that figure's README. No figure is hand-edited.
- **Count before concluding.** Several analyses here fan out over chunks or
  windows; a partial result is usually well-formed and silently wrong. Inputs and
  outputs are reconciled explicitly (e.g. `finalize.sh`).
- **Superseded code is archived, not deleted** — a `_v2`/`_v3` suffix means an
  earlier version is still on disk and still cited by an older figure.
- **Absolute cluster paths live in exactly one file**, `configs/paths.py`.

---

## License

- **Code** is released under the [MIT License](LICENSE).
- **Data and figures** in this repository (`data/`, the tables and plots under `figures/`
  and `analysis/`) are released under
  [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

Third-party inputs keep their own terms and are not redistributed here: AlphaFold 3
model parameters, the full IEDB, and raw prediction outputs are listed, with how to
obtain or regenerate each, in [`docs/PROVENANCE.md`](docs/PROVENANCE.md).

## Citation

If you use this code or data, please cite the paper and the archived release
(see [`CITATION.cff`](CITATION.cff)):

> Woods LJ, Neff B, Lahouel K, Goursaud A, Mulazimoglu M, Bates K, Tomasetti C, Altin JA.
> General prediction of T cell receptor antigen specificity from sequence using
> AlphaFold 3. *Science Advances* (2026).
> Code and data: [doi:10.5281/zenodo.22967485](https://doi.org/10.5281/zenodo.22967485)
