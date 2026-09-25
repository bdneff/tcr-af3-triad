#!/usr/bin/env python3
"""01_stage_triad_parquet.py — build a triad parquet for the class II triads that
the published pipeline never folded.

**Runs on gemini, inside `tcrtrifold-experiments`** — it imports the pipeline's own
`generate_job_name` so the job hashes match what the inference and feature-extraction
stages expect. It cannot run locally.

    conda activate tcrtrifold-experiments
    python 01_stage_triad_parquet.py --repo <path to tcr-af3-triad> \
        --out /tmp/hla2_triad.staged.parquet

Context. Twelve class II triads are deposited after the pipeline's AF3 training cutoff
of **2023-01-12** (`clean_pdb.py`, from the AF3 paper's data-availability statement).
Seven were already folded in the published `pdb.nf` run and need nothing. These five —
8vq8, 9aud, 9yaf, 9yag, 27eb — appear in **neither** discovery input
(`table_S1_structure_benchmark_complexes.csv`, `db_summary.dat`), so they were never
candidates and must be staged explicitly.

Why not hand these to `clean_pdb.py`? Its CLI takes the *discovery* inputs
(`--raw_csv_path`, `--raw_stcr_path`) and re-derives triads from RCSB; it has no entry
for an already-resolved triad table. The sequences here come from John's TSV, which is
already in the pipeline's `SEQ_STRUCT` shape. Cleaning the *crystals* is a separate
step and is only needed for the RMSD reference, not for inference.

**`job_name` is regenerated, not taken from the TSV.** John's TSV puts the PDB ID in
that column, but the pipeline names inference directories by a content hash. Passing
the ID through would produce directories that `extract_triad_conf_feat.py` cannot find,
and the lookup would come back empty rather than loudly wrong.
"""
import argparse
import sys
from pathlib import Path

import polars as pl

# The five class II triads absent from the published run. Kept explicit rather than
# derived, so this script cannot silently widen its own scope.
TARGETS = ["8vq8", "9aud", "9yaf", "9yag", "27eb"]

SEQ_FIELDS = ["peptide", "mhc_1_seq", "mhc_2_seq", "tcr_1_seq", "tcr_2_seq"]
AA = set("ACDEFGHIKLMNPQRSTVWY")


# generate_job_name(df, cols, name="job_name") concatenates `cols` and md5s the result.
# The pipeline's triad column set is this one, used at every call site in
# src/tcrtrifold/neg_creation.py (lines 217, 600, 713):
JOB_NAME_COLS = ["peptide", "mhc_1_seq", "mhc_2_seq", "tcr_1_seq", "tcr_2_seq"]

# NOTE — this does NOT reproduce the job_names in pdb_triad.cleaned.parquet. clean_pdb.py
# never calls generate_job_name (it does not appear in that file), so the published PDB
# table is named by a different mechanism. Recorded rather than papered over.
#
# That discrepancy does not affect us. These five triads have no pre-existing inference
# directory, so nothing needs to match a published hash; what matters is that the name
# is stable and unique, since inference writes <job_name>/ and
# extract_triad_conf_feat.py reads <job_name>/ — both from THIS parquet. Using the
# pipeline's own function with the pipeline's own triad columns keeps the convention.
#
# The gate that does matter is collision: a name equal to an existing directory would
# make inference appear "already done" and silently attach our metrics to someone
# else's structure. That is checked below against the published table.


def load_job_namer():
    """The pipeline's own hash function. A missing gate must fail, not warn."""
    try:
        from tcrtrifold.utils import generate_job_name
        return generate_job_name
    except ImportError as e:
        sys.exit(
            f"cannot import tcrtrifold.utils.generate_job_name ({e}).\n"
            "Run this inside the pipeline env on gemini:\n"
            "    conda activate tcrtrifold-experiments\n"
            "Do NOT substitute a local hash — the job name determines the inference\n"
            "directory, and a different one silently detaches the metrics from the runs."
        )


def assert_no_collision(names, reference_parquet):
    """A generated name must not equal an existing one.

    A collision would let inference see the directory as already present and skip it,
    attaching our metrics to a different structure — the silent-hole failure this
    project keeps meeting. Absence of the reference is itself an error: a check that
    cannot run is not a reason to continue.
    """
    ref = pl.read_parquet(reference_parquet)
    clash = sorted(set(names) & set(ref["job_name"].to_list()))
    if clash:
        sys.exit(f"job_name collides with published rows: {clash}")
    print(f"no collision with the {ref.height} published job_names")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", required=True, help="path to the tcr-af3-triad checkout")
    ap.add_argument("--out", required=True, help="output parquet")
    ap.add_argument("--targets", nargs="*", default=TARGETS)
    ap.add_argument("--reference",
                    default="/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments"
                            "/data/pdb/triad/staged/pdb_triad.cleaned.parquet",
                    help="published table used to calibrate the job_name hash")
    a = ap.parse_args()

    repo = Path(a.repo)
    tsvs = [repo / "data/hla2_9_triads.tsv", repo / "data/hla2_derived_extra.tsv"]
    frames = []
    for t in tsvs:
        if not t.exists():
            sys.exit(f"missing input: {t}")
        frames.append(pl.read_csv(t, separator="\t", infer_schema_length=None))
    # Align on the shared columns; the derived row carries the same SEQ_STRUCT fields.
    shared = [c for c in frames[0].columns if all(c in f.columns for f in frames)]
    df = pl.concat([f.select(shared) for f in frames], how="vertical")

    # John's TSV stores the PDB ID under 'job_name'; make that explicit before
    # regenerating the real job name.
    if "pdb" not in df.columns:
        df = df.rename({"job_name": "pdb"})
    df = df.with_columns(pl.col("pdb").str.to_lowercase())

    want = [t.lower() for t in a.targets]
    sub = df.filter(pl.col("pdb").is_in(want))

    found = sub["pdb"].to_list()
    missing = [t for t in want if t not in found]
    print(f"requested {len(want)}: found {len(found)} {sorted(found)}, missing {missing}")
    if missing:
        sys.exit(f"not in the TSVs: {missing}")
    assert len(found) == len(set(found)), f"duplicate rows: {found}"

    # Gates. Every one of these has burned this project before, so they fail rather
    # than warn.
    bad_class = sub.filter(pl.col("mhc_class") != "II")
    if bad_class.height:
        sys.exit(f"expected class II only, got: {bad_class['pdb'].to_list()}")
    for f in SEQ_FIELDS:
        empty = sub.filter(pl.col(f).is_null() | (pl.col(f).str.len_chars() == 0))
        if empty.height:
            sys.exit(f"{f} empty for {empty['pdb'].to_list()} — class II must fill both MHC chains")
        odd = [(p, s) for p, s in zip(sub["pdb"], sub[f]) if set(s.upper()) - AA]
        if odd:
            sys.exit(f"{f}: non-standard residues in {[p for p, _ in odd]}; "
                     "this pipeline does not model them, and they are omitted everywhere else")

    meta = repo / "analysis/classII_supp3c/data/pdb_metadata.csv"
    if meta.exists():
        m = (pl.read_csv(meta)
             .select(pl.col("pdb_id").str.to_lowercase().alias("pdb"), "release_date"))
        sub = sub.join(m, on="pdb", how="left")
        stale = sub.filter(pl.col("release_date").is_null()
                           | (pl.col("release_date") <= "2023-01-12"))
        if stale.height:
            sys.exit("these are not post-cutoff (2023-01-12) or have no release date: "
                     f"{stale['pdb'].to_list()}")
        print("all rows deposited after the 2023-01-12 AF3 training cutoff")

    gen = load_job_namer()
    sub = gen(sub.drop("job_name", strict=False), JOB_NAME_COLS, name="job_name")
    names = sub["job_name"].to_list()
    if len(set(names)) != len(names):
        sys.exit(f"duplicate job_name within our own rows: {names}")
    assert_no_collision(names, a.reference)

    sub.write_parquet(a.out)
    print(f"\nwrote {a.out} — {sub.height} rows")
    for p, j in zip(sub["pdb"], names):
        print(f"  {p}  {j}")


if __name__ == "__main__":
    main()
