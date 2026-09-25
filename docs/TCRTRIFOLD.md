# The tcrtrifold pipeline — how it actually works, and how to run parts of it

Every prediction and every metric in this paper comes from
[`AltinLab/tcrtrifold-experiments`](https://github.com/AltinLab/tcrtrifold-experiments) (Lawson's
Nextflow pipeline). This repo clones it rather than vendoring it, which is right — but it means the
pipeline's shape is invisible from here, and working that out from scratch costs an afternoon.

This is the map. Companion to `TCRDOCK.md`, which does the same job for the AF-TCRdock baseline.

Cluster checkout: **`/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments`** (`$WS/...`,
`configs.paths.TRIFOLD`). Everything below was verified against it on 2026-09-09.

---

## The one thing to know

**Nextflow is orchestration; the logic lives in ordinary Python scripts.** Each process shells out:

```groovy
process VALIDATION_STANDARDIZE {
  script:
  """
  validation_standardize.py --pdb_parquet ${pdb_pq} \
    -ot pdb_validation_triad.cleaned.parquet -op pdb_validation_pmhc.cleaned.parquet
  """
}
```

So **you can run a stage directly without the workflow engine, and it is the same code.** That
matters because adding a new input set to the Nextflow means editing a repo that is mid-refactor —
`CLEAN_PDB_VALIDATION_WORKFLOW` and `MSA_WORKFLOW` are currently commented out in
`workflows/subworkflows/local/cleaning/main.nf`, and `conf/` holds only execution configs
(`af3.config`, `boltz.config`, `local.config`), not input-set declarations.

Calling `clean_pdb.py` / `rmsd.py` directly is **not** reinventing the pipeline. Writing your own
cleaner or your own RMSD is.

The scripts live under `workflows/subworkflows/local/<stage>/resources/usr/bin/`. The cleaning
stage holds one per dataset, which is the clearest statement of how the pipeline is organised:

```
clean_pdb.py            clean_cresta.py        clean_iedb_I.py    clean_expansion.py
clean_cross_reactivity.py  clean_iedb_II.py    format_true_pdbs.py  validation_standardize.py
```

There is **no local copy of these scripts** — read them on the cluster. (`docs/provenance/` holds
only `nextflow_run.log` and `cleanpdb_test/run.log`.)

**`docs/provenance/nextflow_run.log` is the most useful offline artifact in this repo**: a complete
27-May-2026 `pdb.nf` execution trace, so the whole DAG and every process name can be recovered
without cluster access —

```
CLEAN_PDB → FORMAT_TRUE_PDBS
  → MSA_FROM_TRIAD_PARQUET:{SEQ_LIST_TO_FASTA, MSA_WORKFLOW:{FILT_FORMAT_MSA, RUN_MSA}}
  → UNBATCHED_INFERENCE_…:{NOOP_DEP, SEQ_LIST_TO_FASTA, COMPOSE_INFERENCE_JSON,
                            INFERENCE, CLEAN_INFERENCE_DIR}      # 288 staged → 136 folded
  → DEPEND_TRIAD_ALL_ON_INFERENCE → BOLTZ_FROM_TRIAD_PARQUET → DEPEND_TRIAD_ALL_ON_BOLTZ
  → TCRDOCK_GEOM_FROM_{PDB, AF3_INFERENCE, BOLTZ_INFERENCE}
  → COMPUTE_RMSD_{AF3, BOLTZ}
```

Two things follow from it. **`COMPOSE_INFERENCE_JSON` is the pipeline's own AF3 JSON builder** — so
chain order for a new set is its job, not something to hand-roll. And the **validation set is a
second arm of the same run**, not a separate workflow:

```
VALIDATION_STANDARDIZE → EXCLUDE_AF3_TRAINING_DATA:{SPLIT_QUERY, BLAST_MHC_1/2, BLAST_TCR_1/2}
  → ANNOTATE_TRIAD_FROM_PDB → GEN_NEGATIVES
  → UNBATCHED_INFERENCE_VALIDATION (243 staged → 44 folded) → TCRDOCK_GEOM_VALIDATION
```

**Membership in the validation set is therefore not a date filter.** It runs BLAST against AF3's
training data, which is why 48 post-cutoff class I structures yield only the **16** in Supp Fig 3c.
Any class II points added to that panel must clear the same gate, or the two colours mean different
things. Note also that `pdb.nf` contains **no PTI-PAE process** — that is the standalone step below.

**They must be run inside the pipeline's environment**: `conda activate tcrtrifold-experiments`.
Outside it they fail at the first import (`ModuleNotFoundError: No module named 'tcrtrifold'`) —
`tcrtrifold` is both the package and the env name, which is easy to misread as a missing install.
Never run from `base` (it is the shared Mamba install).

### Adding a dataset — there is a precedent, do not invent one

`$WS/tcrtrifold-experiments/scripts/` holds **one directory per input set**:

```
pdb  pdb_validation  cresta  cresta_new  cresta_no_msa  expansion
iedb_I  iedb_I_full  iedb_II  iedb_II_corrected  iedb_II_full  iedb_meta
cross_reactivity  cross_reactivity_controls
```

**`pdb/` + `pdb_validation/` are the template for a structural set — not `expansion/`.**
An earlier version of this file claimed `expansion/` was the TCR3d structural expansion
(Reviewer #1 comment #4, 153 → 380 triads). It is not: `expansion/` is a clonotype /
HLA-genotyping cohort, and `clean_expansion.py` cleans *that*. Copying it would have produced a
pipeline shaped for the wrong kind of input. `pdb_validation/` is the closest analogue to a
post-cutoff structural benchmark, which is exactly what a Supp Fig 3c addition is.

The two directories are split by **stage**, not by dataset:

| | script | what it runs | conda env |
|---|---|---|---|
| clean + infer | `scripts/pdb/01_clean_and_inf.sh` | `nextflow run ./workflows/pdb.nf` | **`nf-core-tcr3d`** |
| confidence metrics | `scripts/pdb_validation/02_extract_triad_conf_feat.sh` | `workflows/bin/extract_triad_conf_feat.py` | **`tcrtrifold-experiments`** |

**Two different environments**, and confusing them is the first thing that goes wrong. Nextflow
runs in `nf-core-tcr3d`; the standalone Python steps run in `tcrtrifold-experiments`.

`workflows/bin/` holds only `extract_triad_conf_feat.py` and `extract_pmhc_conf_feat.py` — RMSD is
not there, it lives in the `extract_feat` subworkflow's `bin/` (this repo keeps a copy at
`analysis/triad_analysis/rmsd.py`).

### Entering the workflow downstream of discovery

`workflows/pdb.nf`'s entry workflow is **unnamed** (`workflow {`, line 33), so no other file can
`include` it — and it begins from *discovery* inputs:

```groovy
rep_csv  = Channel.fromPath(params.input_replication)   // table_S1_structure_benchmark_complexes.csv
stcr_csv = Channel.fromPath(params.input_stcr)          // STCRDab db_summary.dat
```

It resolves triads from those, then cleans, folds, and measures them. A set whose triads are
**already resolved** — `data/hla2_9_triads.tsv` is precisely that — does not want that front end.

The named subworkflows it composes *are* importable, and they are the downstream entry points:

```groovy
include { MSA_FROM_TRIAD_PARQUET; … }                                  from './subworkflows/local/af3_adapter'
include { CLEAN_PDB; FORMAT_TRUE_PDBS; VALIDATION_STANDARDIZE }        from './subworkflows/local/cleaning'
include { COMPUTE_RMSD; … }                                            from './subworkflows/local/extract_feat'
```

So a class II addition is shaped like this: write John's TSV to a triad parquet in `SEQ_STRUCT`
schema — a **format conversion, not a transformation**, since the columns already match — then
drive `CLEAN_PDB` → `MSA_FROM_TRIAD_PARQUET` → AF3 → `COMPUTE_RMSD` and
`extract_triad_conf_feat.py`. Check each subworkflow's `take:` before wiring it; do not assume the
channel shapes.

## Where the data is

Under **`$WS/tcrtrifold-experiments/data/pdb/triad/`**:

| Path | Contents |
|---|---|
| `cleaned_pdb/` | **288** canonicalised crystals, `<pdb>.pdb`. This is `rmsd.py`'s `clean_pdb_dir` |
| `staged/` | the metric parquets — `pdb_triad.af3_rmsd.parquet`, `pdb_triad.conf_af3.parquet`, `pdb_validation_triad.*`, … (mirrored into this repo's `data/`) |
| `inference/` | **AF3** outputs. Per job: `<hash>.h5`, `<hash>_model.cif.gz`, `seed-1_sample-0/` … `seed-5_sample-0/` |
| `predictions/` | **Boltz** outputs. Per job: `<hash>_model_0..4.cif`, `confidence_<hash>_model_N.json`, `pae_*.npz`, `pde_*.npz` |
| `predictions_with_templates/` | a second inference run — establish which one a number came from before quoting template settings |
| `raw/` | unprocessed downloads |

**Directories are named by a content hash, not by PDB ID** (`generate_job_name` in
`tcrtrifold.utils`). You cannot `ls` your way to a structure; go through the parquets.

Telling the two inference arms apart: the Boltz confidence schema has `ligand_iptm`, `complex_pde`,
`complex_ipde`, `pair_chains_iptm` and there are `pde_*.npz` files. **AF3 emits no PDE**, and AF3
output is organised as `seed-N_sample-M/` subdirectories.

## Input schema

`clean_pdb.py` declares exactly what a triad row must carry:

```python
SEQ_STRUCT = pl.Struct({
    "peptide": pl.String, "mhc_1_seq": pl.String, "mhc_2_seq": pl.String,
    "tcr_1_seq": pl.String, "tcr_2_seq": pl.String,
})
```

plus `pdb`, `mhc_class`, `job_name`, and the `mhc_*_chain` / `mhc_*_species` / `mhc_*_name`
descriptors. Class I leaves `mhc_2_seq` null (guarded at lines ~547/565); class II fills both.

**`data/hla2_9_triads.tsv` is in this format** — John built it to be fed straight in. Note it is
*not* the schema of `data/tcr3d_triad_sequences.csv`, which is the TCR3d expansion's own table;
comparing against that file and concluding the TSV was ad-hoc is a mistake that has been made here
once already.

## Canonical chain convention

`clean_pdb.py` emits five chains in a fixed order, identical for both MHC classes:

| chain | class I | class II |
|---|---|---|
| **A** | **peptide** | **peptide** |
| **B** | MHC heavy | MHC α |
| **C** | β₂-microglobulin | MHC β |
| **D** | TCR α | TCR α |
| **E** | TCR β | TCR β |

Verified by CA counts — `8enh` (I): A=9, B=275, C=99, D=205, E=245; `8vd0` (II): A=16, B=183,
C=178, D=202, E=240.

**The peptide is chain A.** Some tools and earlier inputs use the opposite order (A=MHC, B=β₂m, C=peptide),
so anything moving between the two conventions must be re-lettered. This is a quiet failure rather than a loud
one: `rmsd.py` selects `segid D`/`segid E` for CDR RMSD, so a swapped A/B/C leaves *that* metric
correct while corrupting the MHC-frame alignment and every peptide metric.

**Tethered class II peptides are handled.** Some class II deposits fuse the peptide to the MHC β
chain through a linker, giving four polymer entities where a prediction has five (`8vd0`, `8vq8`,
`9aud`). `clean_pdb.py` splits them: `8vd0`'s cleaned chain A is exactly the 16-residue peptide and
chain C is the β chain without it. No downstream special-casing is needed.

## AF3 run configuration

From the manuscript Methods, *"AlphaFold 3 configuration"*: MSAs **on**, templates **on**, no
custom templates or alignments, seeds **1–5**, **1** diffusion sample, default **10** recycles.

- MSAs are not optional for comparability — `docs/provenance/nextflow_run.log` records
  `MSA_WORKFLOW:RUN_MSA`. **But the pipeline still passes `--norun_data_pipeline`**, and an earlier
  version of this file wrongly said it must not. The two stages are separate: `RUN_MSA` computes the
  alignments first and `COMPOSE_INFERENCE_JSON` embeds them in the JSON, so the inference call has
  nothing left to search. Verified from a real task's `.command.sh`:

  ```bash
  python /app/alphafold/run_alphafold.py \
      --json_path=<hash>.json \
      --model_dir=/ref_genomes/alphafold/alphafold3/models \
      --db_dir=/ref_genomes/alphafold/alphafold3/ \
      --output_dir=. \
      --norun_data_pipeline \
      --num_diffusion_samples=1
  ```

  So the test of "did this prediction get MSAs" is **whether the JSON carries `unpairedMsa`/
  `pairedMsa`**, not whether the flag is present. A standalone run that sets the flag *without*
  precomputed MSAs in its JSON folds from bare sequence — the failure this note originally meant to
  warn about, stated the wrong way round.

- **`--output_dir=.`** — AF3 writes into the Nextflow task directory, and `workDir` is
  `/scratch/$USER/work` (`conf/local.config`). Predictions therefore land on **ephemeral** storage
  and must be copied to `/tgen_labs` before the ~30 day sweep.
- **Templates are date-capped by default.** `max_template_date` defaults to `'2021-09-30'`
  (`run_alphafold.py` ~line 194, *"the date from the AlphaFold 3 paper"*), read at ~line 718, and
  nothing overrides it. So "templates on" cannot reach anything deposited after Sept 2021 —
  including every structure in the post-training-cutoff validation sets. Check with:
  ```bash
  singularity exec /tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif \
      sed -n "188,202p" /app/alphafold/run_alphafold.py
  ```
  This is worth knowing because the Methods sentence alone reads as though a post-cutoff benchmark
  could be handed its own answer. It cannot.

## Metrics

`analysis/triad_analysis/rmsd.py` (a copy of the pipeline's) computes them, importing
`tcrtrifold.tcr` and `mdaf3` — **so it only runs in the pipeline's conda environment on Gemini**
(`conda activate tcrtrifold-experiments`), never locally.

- **CDR RMSD** (`true_pred_cdr_rmsd`) superposes both structures on the MHC frame, sequence-aligns
  the TCR chains, IMGT-numbers them via anarci, and takes CA RMSD over CDR1/2/2.5 (weight 1) and
  CDR3 including the conserved C and F (weight 3). Gaps are fine — it aligns internally, so a
  cleaned chain missing disordered residues is not a problem.
- **PTI-PAE** is `mean_p_tcr_interface_pae`, and it is **not** in `pdb_triad.conf_af3.parquet` —
  that table carries only AF3's own summary (`iptm`, `ptm`, `chain_pair_pae_min`, `chain_ptm`,
  `ranking_score`, `fraction_disordered`, `has_clash`). PTI-PAE is derived downstream by
  **`workflows/bin/extract_triad_conf_feat.py`**, and lands in the *validation* table:

  ```bash
  conda run -n tcrtrifold-experiments --live-stream python \
      ./workflows/bin/extract_triad_conf_feat.py \
      --input_parquet data/pdb/triad/staged/pdb_validation_triad.neg.parquet \
      --inference_type af3 \
      --inference_dir data/pdb/triad/inference \
      --output_path data/pdb/triad/pdb_validation_triad.conf_af3.parquet
  ```

  It reads existing inference output, so it needs **no GPU and no refold** — a triad already
  predicted can be given a PTI-PAE cheaply. This is the authoritative definition; never
  reimplement the peptide↔TCR interface mean from a `pae_*.npz`, since a plausible-but-different
  definition would put new points on a subtly wrong scale.

- **Two `conf_af3` tables exist and they are not interchangeable.** `pdb_triad.conf_af3.parquet`
  (152 rows) is the AF3 summary for the main PDB set; `pdb_validation_triad.conf_af3.parquet` is
  the post-cutoff validation set *with* PTI-PAE. Supp Fig 3c reads the latter. Note the extraction
  script writes to `data/pdb/triad/`, while `figures/reviewer_responses/structural_analysis/README.md`
  cites the copy under `staged/` — presumably the Nextflow publish target. Confirm which one a
  number came from before quoting it.

- **`cdr_rmsd` is null throughout; the live column is `cdr_rmsd_af3_<N>` for N = 0…4**, one per
  seed, each with a companion `pred_sample_rank_N`. **The figures use index 0**
  (`structural_analysis/README.md`, provenance section). Do not average across N: ranks tie
  (`8pjg` ranks `1,2,3,3,3` with three byte-identical RMSD rows), so a mean would silently
  weight one structure three times.
- Figures plot `31 − mean_p_tcr_interface_pae` so that higher reads as more confident (verified:
  `8enh` 31 − 2.19 = 28.81; `8es9` 31 − 14.22 = 16.78).

## Gotchas

- **`/tgen_labs` is not mounted on the login node.** `ls` fails with `No such file or directory`,
  `mkdir` with `Permission denied`. Grab any allocation first —
  `srun -p compute -c 1 --mem=2G -t 00:10:00 --pty bash -i` — including just to copy files out of
  `/scratch`.
- **`find` over `$WS/tcrtrifold-experiments` hangs.** `data/` holds thousands of hash directories.
  Always bound with `-maxdepth`, or grep only `workflows/`, `conf/`, `src/`, `scripts/`.
- **`polars` is not in `base`.** Activate an env that has it before reading a parquet.
- **`/scratch` is swept ~30 days** and takes `.git` with it, leaving a directory that looks like a
  repo and is not (`fatal: not a git repository` with `.git` present). Re-clone; do not `git init`.
  Copy anything expensive to `$WS` — AF3 predictions cost roughly an hour of A100 per triad.

## Repo layout, for orientation

```
main.nf  nextflow.config  conf/{af3,boltz,local}.config
workflows/subworkflows/local/<stage>/main.nf                    process definitions
workflows/subworkflows/local/<stage>/resources/usr/bin/*.py     the actual logic
src/  scripts/  envs/  notebooks/  build/  results/  logs/  tmp/
data/pdb/triad/{cleaned_pdb,staged,inference,predictions,...}
```
