# Reorganization — what moved, what was dropped, and why

2026-07-24. Before this, the project lived in two places: a local working
directory that had grown organically, and `/home/bneff/tcr_project_to_organize`
on gemini. This records the merge into a single manuscript companion repo, so
that anything you remember existing can still be found.

**Nothing was deleted without being checked first.** Commit `3900cc3` is a
baseline snapshot of the project exactly as found, so every move below is
reversible with `git diff 3900cc3`.

---

## The guiding principle

Structure maps to **the manuscript**, not to our directory history. A reviewer
should be able to go from a figure to the code and data that made it. Consequently
`manuscript_plots/` — which already had one directory per figure, each with its
own README — was promoted to `figures/` and became the spine of the repo.

Second principle: **code lives in git; bulk data,
weights, and run outputs live on the cluster**, reachable through
`configs/paths.py`.

---

## Moves

| From | To | Note |
|---|---|---|
| `manuscript_plots/` | `figures/` | the figure map; per-figure READMEs kept |
| `figures/` (loose PNGs) | `scratch_figs/` | exploratory plots, not manuscript figures; not included in the public release |
| `tcrtrifold/` + `scripts/` | `analysis/triad_analysis/` | merged; the local copies were a subset of gemini's |
| gemini `tcrdock_validation/` | `analysis/tcrdock_validation/` | AF-TCRdock baseline, Figure 4c |
| `crossreact-controls/` + `cross_reactivity/` | `analysis/crossreact_controls/` | merged |
| gemini `tcr3d_expansion.tar.gz` | `analysis/tcr3d_expansion/` | see *near miss* below |
| gemini `fig1_data_export/` | `data/exports/` | CSV exports behind Figure 1 |
| gemini `nextflow_run.log` | `docs/provenance/` | force-added past `*.log`; it is evidence |
| gemini `cleanpdb_test/` | `docs/provenance/cleanpdb_test/` | chain-cleaning QC record |
| `*_expanded.html` | `scratch_figs/` | interactive exploratory plots; not included in the public release |

Merges used `rsync --ignore-existing`, so where a file existed both locally and on
gemini the **local copy won**. The two were spot-checked as identical for the
overlapping parquets and CSVs.

---

## A near miss worth recording

`tcr3d_expansion.tar.gz` sat at the top level of the gemini directory alongside
four genuinely redundant nested archives (`crossreact-controls.zip`,
`tcrdock_validation.tar.gz`, and duplicate `.zip`/`.tar.gz` pairs), all of which
were re-archives of content already extracted. It looked like a fifth.

It is not. It is the `00_refresh_stcrdab` → `01_identify_new_entries` →
`02_apply_to_repo` → `03_run_pipeline` → `04_replot_figure1` pipeline that answers
**Reviewer #1 comment #4**, growing Figure 1's set from 153 to 380 triads
(pre-cutoff 130→308, post-cutoff 22→71). Nothing else in the repo reproduces it.

It was recovered by extracting it from the source tarball and reading its README
before deleting. The general lesson, and the reason this section exists: *an
archive that looks redundant because of where it sits is not redundant.* Check
the contents, not the neighbourhood.

---

## Deliberately not imported

The gemini directory was 28.2 GB uncompressed; ~36 MB of it was science. What was
left behind, and how to get it back:

| Item | Size | Why not | How to regenerate |
|---|---|---|---|
| `tcrdock_validation/env/` | 7.7 GB | a conda environment | `bash analysis/tcrdock_validation/setup_env.sh` |
| `tcrtrifold/iedb_public.db` | 13.5 GB | database | `configs/paths.IEDB_DB` on cluster; 90 MB pruned copy local |
| `tcrdock_validation/runs*/` | 4.8 GB | per-target outputs | resubmit the arrays (`docs/TCRDOCK.md`) |
| `tcrdock_validation/TCRdock/` | 1.5 GB | upstream code | `github.com/phbradley/TCRdock` |
| `alphafold_params/` | 0.75 GB | redistributable weights | `setup_env.sh` step 4 |
| `archives/tcrtrifold-experiments-main.zip` | 283 MB | upstream pipeline | `github.com/AltinLab/tcrtrifold-experiments` |
| `tcr_project_to_organize.tar.gz` | 7.4 GB | the staging tarball itself | still on gemini at `/home/bneff/` |

All are gitignored and all are declared in `configs/paths.py` or
`docs/PROVENANCE.md`.

> The staging tarball is still in the working directory and still gitignored.
> Delete it once you are satisfied nothing else is needed from it — its sha256 is
> `ef8d3404acac269d9521f2845ff7149339c311a2b352d528e2a4b65b5ac8aeb1`.

---

## Corrections made along the way

Two documents were found to disagree with the artifacts they describe. Both are
kept, with the discrepancy noted rather than silently overwritten:

- **`README_gemini.md` says 197 triads.** The actual class I run was **145**
  (and class II **1433**). It is the *pre-flight* draft, written before the runs,
  and still contains placeholder SLURM directives (`##SBATCH --account=your_account`)
  where the real jobs used `gpu-a100`. `docs/TCRDOCK.md` supersedes it.

- **`data/tcr3d_triad_sequences.csv` has a corrupt peptide for 3D3V and 3D39.**
  `peptide_seq` reads `LLFGPVYV` (8-mer) where both crystals are the 9-mer
  `LLFGFPVYV` — the Phe at P5, an anchor residue, is missing. The table's own
  `epitope_listed` column agrees with the crystal. Anything staged from
  `peptide_seq` for those two entries folded the wrong peptide; take sequences from the
  crystals instead.
