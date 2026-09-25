# Provenance — where every input came from

For each tracked data file: what produced it, what it contains, and what depends
on it. Anything *not* tracked is listed at the bottom with the command that
regenerates it.

Cluster paths are symbolic here; the literal values live in `configs/paths.py`.

---

## Upstream sources

| Source | What it gives us |
|---|---|
| **RCSB / PDB** | crystal structures (`structures/*.pdb`) — the ground truth all RMSDs are measured against |
| **STCRDab** | curated TCR:pMHC complex annotations; the initial triad set |
| **TCR3d** | the curated triad expansion added in revision (`analysis/tcr3d_expansion/`) |
| **IEDB** (`2025-04-15`) | epitope/affinity records; `iedb_kd_records.csv`, Figure 3 affinity correlation |
| **[AltinLab/tcrtrifold-experiments](https://github.com/AltinLab/tcrtrifold-experiments)** | the Nextflow pipeline that runs MSA → AF3 + Boltz-2 → RMSD → confidence extraction. **All prediction parquets come from here.** |
| **[phbradley/TCRdock](https://github.com/phbradley/TCRdock)** | the AF-TCRdock baseline (Figure 4c) |

> **Templates were on, but capped at 2021-09-30.** AF3's `max_template_date` defaults to
> `'2021-09-30'` (the date from the AF3 paper) and nothing in the pipeline overrides it, so the
> template search cannot reach anything deposited after that — including every structure in the
> post-training-cutoff validation set. Verified 2026-09-09 by reading `run_alphafold.py` inside
> `alphafold_3.0.1.sif`. This matters because the Methods sentence *"run with MSAs and templates
> turned on"* reads, on its own, as though a post-cutoff benchmark could be handed its own answer.
> It cannot.
>
> **The manuscript's predictions used real MSAs.** `docs/provenance/nextflow_run.log`
> records `MSA_WORKFLOW:RUN_MSA`, `SEQ_LIST_TO_FASTA`, and `FILT_FORMAT_MSA` stages.
> Any new prediction intended to be comparable must run the data pipeline too —
> i.e. **not** `--norun_data_pipeline`.

---

## Prediction + metric tables (`data/`)

Produced by the tcrtrifold Nextflow pipeline; staged on the cluster at
`$TRIAD_PARQUETS` and mirrored here because the figure notebooks need them and
they are small.

| File | Contents |
|---|---|
| `pdb_triad.af3_rmsd.parquet` | AF3 RMSDs vs deposited structure — `docking_rmsd_af3_0`, `cdr_rmsd_af3_0`, `peptide/mhc/tcr_rmsd_af3_0`. 287 unique PDBs. |
| `pdb_triad.boltz_rmsd.parquet` | same, Boltz-2 |
| `pdb_triad.conf_af3_full.parquet` | AF3 confidence metrics, incl. `mean_p_tcr_interface_pae` (PTI-PAE) — the paper's headline feature |
| `pdb_validation_triad.annotated.parquet` | validation split with pre/post AF3-cutoff annotation |
| `cross_reactivity_controls_triad.conf_af3.parquet` | confidence for the scrambled/random-peptide negatives |
| `data/exports/{af3_rmsd,boltz_rmsd,conf_af3}.csv` | CSV exports of the above, for the Figure 1 notebook |

`docking_rmsd` is the column that detects a **misdock**: a value of ~35 Å with a
`tcr_rmsd` near 0.5 Å means the TCR folded correctly and was placed wrongly.

## Sequence + triad tables

| File | Contents | Caveat |
|---|---|---|
| `tcr3d_triad_sequences.csv` | 379 rows: per-PDB peptide / MHC / TCR sequences, allele, dates, chain mapping | **Two known defects — see below** |
| `ternary_templates_v2.tsv` | template triads for the pipeline | |
| `classI_complexes.csv`, `classII_complexes.csv` | TCR3d complex listings behind the expansion | |
| `hla2_9_triads.tsv` | 9 post-AF3-cutoff **class II** ternary complexes staged for prediction, to extend Supp Fig 3c (reviewer round 2, minor point 1) | **Hand-assembled by John Altin, 2026-09-08.** In the tcrtrifold pipeline's `SEQ_STRUCT` input schema, so it feeds the pipeline directly (it is *not* the schema of `tcr3d_triad_sequences.csv`). `*_chain` columns hold role labels ("alpha"/"beta"), not deposited chain IDs; see below |
| `hla2_derived_extra.tsv` | 1 row, `27eb`: the third HLA-DQ2.5 gliadin structure John described but the TSV omitted | **Derived here from RCSB** by `analysis/classII_supp3c/scripts/03_derive_27eb.py`; kept separate so John's file still says what he sent |
| `10TCRs.csv`, `75antigens.csv` | the validation subsets | |

> **Defect 1 — corrupt peptide, 3D3V and 3D39.** `peptide_seq` reads `LLFGPVYV`
> (8-mer); both crystals are the 9-mer `LLFGFPVYV`. The Phe at P5 — a primary
> anchor — is dropped. The same table's `epitope_listed` column agrees with the
> crystal, so the file is internally inconsistent and `peptide_seq` is the wrong
> field. Anything folded from `peptide_seq` for those two used the wrong antigen.
>
> **Defect 2 — 7BYD is not a peptide system.** `peptide_seq` reads `GGAIGGAI`,
> which is the crystal's 4-mer `GGAI` repeated to fake an 8-mer. It is a
> *lipopeptide* on a macaque MHC (Mamu-B\*05104) whose myristoyl anchor (`MYR`)
> has no backbone, so it is not a valid peptide:MHC triad for peptide-based analyses.
>
> Both were found by cross-checking the table against the crystal structures; take
> sequences **from the crystals** rather than from `peptide_seq`.

### `hla2_9_triads.tsv` — origin, and a discrepancy to resolve

Assembled by John Altin (Slack, 2026-09-08) to answer **reviewer round 2, minor point 1**: Supp
Fig 3c rests on 16 post-AF3-cutoff class I structures and the reviewer asked for a larger benchmark.
Expanding the *class I* set (~60 now exist) would propagate through Fig 1, Fig 4a and the
pre/post-cutoff splits, so the chosen route adds **class II** structures to Supp Fig 3c alone —
which touches nothing else, because class II validation elsewhere in the paper uses CRESTA and the
paper currently contains no PDB class II triads at all.

His stated selection: all MHC-II:peptide:TCR complexes released on/after 2023-01-12 (RCSB, cross-
checked against TCR3d class II and STCRDab) = 17 → drop `9EJH` (CLIP invariant-chain fragment, not a
specific peptide) = 16 → drop `8TRL`, `8TRQ`, `8TRR`, `9NIG` (citrullinated peptides; the AlphaFold
**Server** cannot model non-standard residues) = 12 → `8PJG` failed server-side = 11 → set aside
`9EJG`, `9EJI` (peptide-*independent* recognition of HLA-DQ2: the TCR engages the MHC α-helices
alone and sits ~12 Å from the peptide, so a peptide-referenced interface PAE and an MHC-frame CDR
RMSD do not measure the same quantity there) = **9**.

**⚠ The file does not match that prose, in two places.** Recorded rather than reconciled, per this
repo's convention:

| John's description | the actual TSV |
|---|---|
| three HLA-DQ2.5 (gliadin) | **two** — `9yaf`, `9yag` |
| `8PJG` dropped (server-side failure) | **present** as row 1 |

Both still total 9. The likely explanation is that `8pjg` was reinstated once the plan moved to
running on Gemini, where a *server*-side failure does not apply, displacing the third gliadin entry
— but that is inference, not something he stated.

**Resolution (2026-09-09).** The missing third gliadin entry was recovered from RCSB as **`27eb`**
(released 2026-08-26, five separate entities) → `hla2_derived_extra.tsv`. With `8pjg` kept, the
set is **10**. Both points are listed for John to confirm in `analysis/classII_supp3c/README.md`.

**The citrullinated entries stay excluded, but for a different reason than his.** His stated
reason (the AlphaFold *Server* cannot model non-standard residues) does not apply to a local run.
Adding them via the CCD code for citrulline was proposed and **rejected (Brandon, 2026-09-08)**:
every other prediction in the paper uses standard residues only, so those points would come from a
different modelling regime. `8trl` and `8trr` were already folded in the published run, though,
*as arginine*, which is how citrulline canonicalises. On 2026-09-09 they were briefly added to the
panel as extra finds, without noticing that John's selection above had already excluded them. He
excluded them again the same day (*"I was excluding citrullinated epitopes since they contain a
modified AA that is not represented by AF3"*). Folding arginine where the crystal has citrulline
models a different molecule, and for DR4 epitopes the modification is what creates the
neo-epitope. They are listed in `EXCLUDED_CITRULLINATED` in
`figures/reviewer_responses/structural_analysis/build_classII_rows.py`.

**Chain mapping.** The `*_chain` columns are role labels, so the deposited chain IDs that
`format_true_pdbs.py` needs are derived by `analysis/classII_supp3c/pipeline/02_derive_segids.py`,
using sequence match to RCSB entities plus a CA-contact check to pick one complex from multi-copy
deposits. `mean_p_tcr_interface_pae` is **empty in all 9 rows** of the TSV. It is a schema
placeholder, not data; the values were computed by `extract_triad_conf_feat.py`.

### Class II predictions and metrics: cluster outputs, not tracked

Everything below lives on Gemini under
`$W = /tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments` (durable storage). It is
excluded bulk, and the tracked record is the metric CSVs in
`figures/reviewer_responses/structural_analysis/`.

| Path | Contents | Regenerate with |
|---|---|---|
| `$W/data/hla2/triad/staged/hla2_triad.staged.parquet` | the 5 newly folded triads, pipeline schema, `job_name` via `generate_job_name` | `analysis/classII_supp3c/pipeline/01_stage_triad_parquet.py` |
| `$W/data/hla2/triad/staged/hla2_triad.segid.parquet` | the same, plus deposited chain IDs | `pipeline/02_derive_segids.py` |
| `$W/data/hla2/triad/inference/<job_name>/` | AF3 predictions, 5 seeds each, with MSA-embedded input JSONs (~296 MB) | `pipeline/hla2.nf` (GPU) |
| `$W/data/hla2/triad/cleaned_pdb/` | canonicalised crystals, the RMSD reference | `pipeline/hla2_rmsd.nf` (`FORMAT_TRUE_PDBS`) |

The five already folded in the published run (`8pjg` `8vcx` `8vcy` `8vd0` `8vd2`) live in
`$W/data/pdb/triad/inference/` like every other published prediction. Their PTI-PAE is regenerated
by running `workflows/bin/extract_triad_conf_feat.py` on `pdb_triad.cleaned.parquet`, and it needs
no GPU.

**Nextflow's work directory is `/scratch/$USER/work`**, which is swept every ~30 days. AF3 writes
into task directories there (`--output_dir=.`). Predictions that were not copied to `$W` by an
explicit `--triad_inf_dir` are lost after a sweep. See `analysis/classII_supp3c/README.md` for how
this happened once and how it was recovered.

## AUC / analysis outputs

Produced by `analysis/tcrdock_validation/analyze_tcrdock_vs_af3*.py` and the
leakage scripts in `analysis/triad_analysis/`:

`classI_per_antigen_AUC.csv`, `classI_per_TCR_AUC.csv`, `classI_per_triad_long.csv`,
`classII_*` equivalents, `per_antigen_auc_TCRdock_vs_AF3.csv`,
`per_component_identities.csv`, `merged_tv_rmsd_pae.csv`,
`dihedral_validation_classI{,_v2}.csv`, `antigen_centric*`, `tcr_centric*`.

`_v2`/`_v3`/`_v4` suffixes are revision-pass regenerations; the unsuffixed file is
the original submission's version. Both are kept — an older figure may still cite
the older table.

## TCRdock baseline (`analysis/tcrdock_validation/`)

| File | Contents |
|---|---|
| `tcrdock_input_class_I.tsv` / `_class_II.tsv` | inputs: **145** / **1433** triads |
| `user_output_w_pae.tsv` / `user_output_classII_w_pae.tsv` | outputs: **145** / **1433** rows — complete, nothing dropped |
| `tcrdock_input_class_I_full.xlsx` | metadata joined to the input, for the AUC analysis |
| `classI_per_*.csv`, `classII_per_*.csv` | the AUC tables feeding Figure 4c |

Full recipe: **`docs/TCRDOCK.md`**.

## Evidence (`docs/provenance/`)

`docs/provenance/` holds evidence rather than prose: `nextflow_run.log` (the
prediction pipeline's own record — this is what established that the manuscript's
predictions ran real MSAs, `MSA_WORKFLOW:RUN_MSA`) and `cleanpdb_test/`
(chain-cleaning QC).

---

## Not tracked — and how to get it back

| Item | Size | Regenerate with |
|---|---|---|
| conda env for TCRdock | 7.7 GB | `bash analysis/tcrdock_validation/setup_env.sh <workdir>` |
| `iedb_public.db` | 13.5 GB | `configs/paths.iedb()` |
| TCRdock run outputs | 4.8 GB | resubmit the arrays, `docs/TCRDOCK.md` |
| `TCRdock/` checkout | 1.5 GB | `git clone https://github.com/phbradley/TCRdock` |
| AlphaFold2 + fine-tuned weights | 0.75 GB | `setup_env.sh` step 4 |
| tcrtrifold Nextflow pipeline | 283 MB | `git clone https://github.com/AltinLab/tcrtrifold-experiments` |
