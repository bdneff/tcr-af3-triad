# Class II extension of Supp Fig 3c

**Reviewer round 2, minor point 1.** Supp Fig 3c shows AF3 PTI-PAE tracking structural accuracy
(CDR RMSD) on structures released *after* the AF3 training cutoff. It rested on 16 class I entries
(Pearson −0.53, Spearman −0.60), and the reviewer asked for a larger benchmark if more post-cutoff
structures exist.

**Result (2026-09-09): 10 class II triads added.** Class II alone gives **r = −0.75 (p = 0.012),
ρ = −0.73 (p = 0.016)**, which is stronger than class I on an independent set. The combined 26 points
give ρ = −0.59 (p = 0.001). The class I points are unchanged, and `plot_supp3c.py --check` enforces
that. The figure, the per-triad table and a review PDF are in
[`../../figures/reviewer_responses/structural_analysis/`](../../figures/reviewer_responses/structural_analysis/README.md).

**Why class II rather than more class I.** About 60 post-cutoff class I triads now exist, but the class I
validation set carries weight across Fig 1, Fig 4a and the pre/post-cutoff splits, and expanding it
at final revision would ripple through the whole paper. Class II only adds: class II validation
elsewhere uses CRESTA, so no PDB class II triad appears anywhere else in the paper. It also arguably
answers the comment better, by showing the relationship holds across MHC classes.

## The ten, and where each came from

| | triads | prediction | metrics |
|---|---|---|---|
| **5** | `8pjg` `8vcx` `8vcy` `8vd0` `8vd2` | already folded in the **published `pdb.nf` run** (Figure 1's benchmark) | CDR RMSD already computed; PTI-PAE extracted from existing inference, **no GPU** |
| **5** | `8vq8` `9aud` `9yaf` `9yag` `27eb` | **newly folded** by `pipeline/hla2.nf`, 2026-09-09 | `pipeline/hla2_rmsd.nf` + `extract_triad_conf_feat.py` |

The first five were never missing. Supp Fig 3c drew from the `pdb_validation` arm, which is 100%
class I in all nine of its tables, so class II predictions that already existed were simply never
plotted. The PTI-PAE extraction was checked by running `8enh` and `8es9` alongside them, and both
reproduced their published values to the last digit.

Every triad is post-cutoff. The cutoff is **2023-01-12**, which `clean_pdb.py` takes from the AF3
paper's data-availability statement, and the earliest release here (`8pjg`, 2024-06-19) is 524 days
past it.

### Excluded: citrullinated epitopes

`8trl` (α-enolase) and `8trr` (vimentin), both on HLA-DR4, are post-cutoff and were already folded.
They are **excluded on chemistry** (John Altin, 2026-09-09: *"I was excluding citrullinated epitopes
since they contain a modified AA that is not represented by AF3"*). Citrulline has no standard
one-letter code, so RCSB canonicalises it to arginine and AF3 folds a different molecule from the
deposit. For these epitopes that matters a lot: DR4's shared-epitope P4 pocket repels arginine, so
citrullination is what creates the neo-epitope. Both predicted well (3.34 and 1.22 Å), so this is
not a quality cut. The same rule was applied earlier, when adding them through a CCD
modified-residue path was rejected (Brandon, 2026-09-08). Taking them out made every statistic
stronger. They are listed in `EXCLUDED_CITRULLINATED` in `build_classII_rows.py`, not deleted.

### For John to confirm

- **`27eb` is not on his TSV.** We derived it from RCSB (`scripts/03_derive_27eb.py` →
  `../../data/hla2_derived_extra.tsv`) as the third HLA-DQ2.5 gliadin structure he described. It is
  the **least confident (PTI-PAE 17.5) and least accurate (22.2 Å)** class II triad, flagged as
  suspicious on confidence *before* its RMSD was known. That makes it the most informative point on
  the class II side, and also the one a reviewer is most likely to question. All five of its chains
  match deposited entities exactly.
- **`8pjg` is present** although his note said it had been dropped because of a server-side failure.
  It was probably put back once the plan moved to Gemini, where that failure does not apply.

## Run order: the path that produced the figure

Everything goes through the tcrtrifold pipeline's own processes. **Stage here, run on Gemini.** Two
conda envs are involved: Nextflow runs in **`nf-core-tcr3d`** and the Python steps run in
**`tcrtrifold-experiments`**. `$W` is `/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments`.

```bash
# ── local: inputs the pipeline reads ──────────────────────────────────────────
python analysis/classII_supp3c/scripts/00_fetch_metadata.py     # RCSB dates; hard-fails if pre-cutoff
python analysis/classII_supp3c/scripts/02_xcheck_sequences.py   # staged vs deposited: 42/45 exact
python analysis/classII_supp3c/scripts/03_derive_27eb.py        # -> data/hla2_derived_extra.tsv

# ── Gemini, conda activate tcrtrifold-experiments ────────────────────────────
# 1. stage: TSVs -> triad parquet, job names from the pipeline's own generate_job_name
python analysis/classII_supp3c/pipeline/01_stage_triad_parquet.py \
    --repo . --out $W/data/hla2/triad/staged/hla2_triad.staged.parquet

# 2. fold (GPU): MSA -> COMPOSE_INFERENCE_JSON -> INFERENCE      [env: nf-core-tcr3d]
cp analysis/classII_supp3c/pipeline/hla2.nf $W/workflows/ && cd $W
nextflow run ./workflows/hla2.nf \
    --input_triad   $PWD/data/hla2/triad/staged/hla2_triad.staged.parquet \
    --triad_inf_dir $PWD/data/hla2/triad/inference -profile gemini -resume
ls data/hla2/triad/inference | wc -l                           # COUNT: want 5 dirs

# 3. PTI-PAE (no GPU)
python ./workflows/bin/extract_triad_conf_feat.py \
    --input_parquet data/hla2/triad/staged/hla2_triad.staged.parquet \
    --inference_type af3 --inference_dir data/hla2/triad/inference \
    --output_path /tmp/hla2.conf.parquet

# 4. deposited chain IDs, then CDR RMSD (no GPU)                 [env: nf-core-tcr3d for step 4b]
python analysis/classII_supp3c/pipeline/02_derive_segids.py \
    --in  $W/data/hla2/triad/staged/hla2_triad.staged.parquet \
    --out $W/data/hla2/triad/staged/hla2_triad.segid.parquet
nextflow run ./workflows/hla2_rmsd.nf \
    --input_triad   $PWD/data/hla2/triad/staged/hla2_triad.segid.parquet \
    --triad_inf_dir $PWD/data/hla2/triad/inference -profile gemini -resume

# ── local: assemble and plot ─────────────────────────────────────────────────
cd figures/reviewer_responses/structural_analysis
python build_classII_rows.py --pae classII_pti_pae.csv
python plot_supp3c.py --check
python make_review_pdf.py
```

## Things that bit, recorded so they don't happen again

**Nextflow reports success without writing the output.** The first `hla2.nf` derived its inference
directory from `$workflow.outputDir`. `-output-dir` only feeds the `output {}` publishing mechanism,
which the file did not define, so every process showed 100% and `Succeeded: 67` while
`data/hla2/triad/inference` did not exist. AF3 runs with `--output_dir=.` and had written into its
task directories under `/scratch/$USER/work` (see `conf/local.config`), which is swept every ~30 days.
The predictions were recovered through the log:

```bash
grep "WORKFLOW:INFERENCE" logs/hla2/.nextflow.log | grep -oE "workDir: [^]]+"
```

Both workflows now take `--triad_inf_dir` explicitly. **Count outputs; never trust the summary.**

**`--norun_data_pipeline` does not mean "no MSAs."** The pipeline passes it on purpose. `RUN_MSA`
runs as a separate stage and `COMPOSE_INFERENCE_JSON` embeds the alignments in the JSON. The real
check is whether the JSON contains `unpairedMsa`/`pairedMsa`, and all five do. The danger is that
flag combined with a JSON that has **no** MSA, which would fold from bare sequence.

**Crystal chain IDs have to be derived, not assumed.** `format_true_pdbs.py` selects
`segid <X>` from the deposit, and John's TSV stores role labels ("alpha"/"beta") rather than chain
letters. `02_derive_segids.py` matches each staged sequence to a deposited entity. Two conventions
come from the published table:

- **A tethered peptide shares MHC β's chain** (`8vd0`: `peptide_segid == mhc_2_segid == C`).
  `8vq8` and `9aud` are tethered and handled the same way.
- **Index order is not copy order.** The published table pairs `8trl`'s MHC `A`/`B` with TCR `I`/`J`.
  `9yaf`, `9yag` and `27eb` each hold two complexes, so each candidate copy was checked against the
  deposited CA coordinates, and every selected chain sits 4–7 Å from the MHC α chain. Picking by index
  could pair an MHC with a TCR from the other copy. That structure never existed, and its RMSD would
  look entirely reasonable.

**Canonical chain order for predictions** (set by `COMPOSE_INFERENCE_JSON`, identical for both
classes): **A = peptide, B = MHC-1, C = MHC-2, D = TCR α, E = TCR β.** An alternative convention uses
the opposite order for A/B/C, which is how the superseded staging below went wrong.

## Templates cannot reach these structures

The Methods say templates were on, which sounds alarming for a post-cutoff benchmark. But AF3's
`max_template_date` defaults to **2021-09-30** (`run_alphafold.py` ~line 194), nothing overrides it,
and that date comes before the 2023-01-12 cutoff and every deposition here. No structure could have
been given to AF3 as its own template, and the same applies to the published 16. Verified in the lab
container:

```bash
singularity exec /tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif \
    sed -n "188,202p" /app/alphafold/run_alphafold.py
```

## Layout

```
pipeline/                        THE PATH USED — every figure number comes from here
  01_stage_triad_parquet.py      TSVs -> triad parquet; job names via generate_job_name
  02_derive_segids.py            deposited chain IDs by sequence match + contact geometry
  hla2.nf                        MSA + AF3 inference (copy into $W/workflows/)
  hla2_rmsd.nf                   FORMAT_TRUE_PDBS -> TCRDOCK geometry -> COMPUTE_RMSD_AF3
scripts/                         inputs the pipeline reads
  00_fetch_metadata.py           RCSB release dates -> data/pdb_metadata.csv
  02_xcheck_sequences.py         staged sequences vs deposited entities
  03_derive_27eb.py              rebuilds the 27eb row -> data/hla2_derived_extra.tsv
data/pdb_metadata.csv            RCSB metadata (release date, resolution, entity count)
superseded/                      the standalone staging path — see below; NOT a figure source
logs/                            SLURM logs (git-ignored)
```

Cluster outputs, not tracked here (regeneration in `../../docs/PROVENANCE.md`):
`$W/data/hla2/triad/{staged,inference,cleaned_pdb}/`.

## `superseded/`: the standalone path, archived rather than deleted

Before we realised John's TSV was already in the pipeline's input schema, a direct-invocation path
was built: `scripts/01_stage_af3.py` → `af3_inputs/*.json` → `af3_array.sbatch`, plus a bespoke
crystal cleaner, `scripts/04_prepare_structures.py`, which was never adopted. It folded all ten on
Gemini (job 37446788, ~1h10m each). **None of it feeds the figure.**

- **The JSONs have the wrong chain order.** They use the opposite convention (A = MHC α, B = MHC β,
  C = peptide) instead of the pipeline's A = peptide. CDR RMSD selects `segid D`/`segid E` and would
  have survived, while the MHC-frame superposition it aligns on would not have. A metric that happens
  to be immune is not a reason to keep a wrong input.
- **`af3_array.sbatch` chooses its JSON by position** (`ls | sed -n "${TASK}p"`). When `27eb.json`
  was added it sorted *first* and shifted every index, so an `--array=1-9` submitted beforehand would
  have run nine jobs, succeeded on all of them, and silently skipped `9yag`.
- `data/manifest.csv` and `structures/_cif/` (deposited mmCIFs) were that path's outputs and inputs.

They stay because an older commit or note may cite them, and because the ten standalone predictions
are an independent cross-check if one is ever needed.
