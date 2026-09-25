#!/usr/bin/env python3
"""build_classII_rows.py — add the class II triads to `merged_tv_rmsd_pae.csv` for Supp Fig 3c.

Answers **reviewer round 2, minor point 1** (the 16-structure benchmark is small). Adds class II
post-training-cutoff triads to the Supp Fig 3c scatter, leaving the class I panel untouched.

Both metrics come from the **published pipeline's own tables**, never recomputed here:

    cdr_rmsd  <- cdr_rmsd_af3_0            in data/exports/af3_rmsd.csv   (mirror of
                                              pdb_triad.af3_rmsd.parquet)
    pti_pae   <- mean_p_tcr_interface_pae  from workflows/bin/extract_triad_conf_feat.py

`cdr_rmsd_af3_0` is **seed index 0**, matching how the original 16 class I points were built
(see this directory's README, "Data provenance"). The per-seed columns must NOT be averaged:
`pred_sample_rank_N` ties (8pjg ranks 1,2,3,3,3) repeat a structure, so a mean silently
weights one model three times.

    # 1. on gemini, produce PTI-PAE from predictions that already exist (no GPU):
    #    python workflows/bin/extract_triad_conf_feat.py \
    #        --input_parquet data/pdb/triad/staged/pdb_triad.cleaned.parquet \
    #        --inference_type af3 --inference_dir data/pdb/triad/inference \
    #        --output_path /tmp/pdb_triad.conf_af3_pti.parquet
    #    then export pdb + mean_p_tcr_interface_pae to CSV and scp it back.
    # 2. locally:
    python build_classII_rows.py --pae classII_pti_pae.csv
    python plot_supp3c.py --check

Guards, in the spirit of "count before concluding":
  * every requested triad lands in exactly one bucket (added / missing rmsd / missing pae),
    and the buckets must sum to the number requested;
  * the 16 class I rows must come through **bit-identical** — the script diffs them and
    aborts rather than writing a file whose published half moved;
  * a triad already present is refused rather than silently duplicated.
"""
import argparse
import os
import sys

import pandas as pd

MERGED = "merged_tv_rmsd_pae.csv"
RMSD = "../../../data/exports/af3_rmsd.csv"
RMSD_COL = "cdr_rmsd_af3_0"
PAE_COL = "mean_p_tcr_interface_pae"

# Class II triads deposited after the pipeline's AF3 training cutoff, **2023-01-12**
# (`clean_pdb.py`, citing the AF3 paper's data-availability statement -- not the
# 2021-09-30 `max_template_date`, which is the separate template cap).
#
#   8pjg 8vcx 8vcy 8vd0 8vd2   John's list, already folded in the published pdb.nf run
#   8vq8 9aud 9yaf 9yag 27eb   absent from the discovery inputs entirely (verified: no hit
#                              in table_S1_structure_benchmark_complexes.csv or
#                              db_summary.dat), so these needed cleaning and inference
TRIADS = ["8pjg", "8vcx", "8vcy", "8vd0", "8vd2",
          "8vq8", "9aud", "9yaf", "9yag", "27eb"]

# EXCLUDED, deliberately — citrullinated epitopes (John Altin, 2026-09-09: "I was
# excluding citrullinated epitopes since they contain a modified AA that is not
# represented by AF3"). Both are post-cutoff class II and were already folded in the
# published run, and both predict well (8trl 3.34 A, 8trr 1.22 A) — they are omitted on
# chemistry, not on quality.
#
# Citrullination converts an arginine side chain to neutral citrulline. That charge change
# is the whole point of these epitopes: DR4's shared-epitope P4 pocket repels arginine, so
# the unmodified peptide is a poor binder and the modification is what creates the
# neo-epitope (the mechanism behind RA autoimmunity). Citrulline has no standard
# one-letter code, so RCSB canonicalises it back to R and AF3 folded plain arginine —
# a different chemistry from the deposit. This matches the rule already applied elsewhere
# in the project: non-standard residues omitted everywhere else are not adopted here.
#
# Kept listed rather than deleted so the decision is visible and reversible; their metrics
# remain in classII_pti_pae.csv and data/exports/af3_rmsd.csv.
EXCLUDED_CITRULLINATED = ["8trl", "8trr"]


def load_pae(path):
    df = pd.read_csv(path)
    col = PAE_COL if PAE_COL in df.columns else None
    if col is None:
        sys.exit(f"{path}: no '{PAE_COL}' column; found {list(df.columns)}")
    if "pdb" not in df.columns:
        sys.exit(f"{path}: no 'pdb' column; found {list(df.columns)}")
    df = df[["pdb", col]].rename(columns={col: "pti_pae"})
    df["pdb"] = df["pdb"].str.lower()
    return df.dropna(subset=["pti_pae"]).drop_duplicates("pdb")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pae", required=True,
                    help="CSV with pdb + mean_p_tcr_interface_pae (from extract_triad_conf_feat.py)")
    ap.add_argument("--rmsd", default=RMSD, help=f"CDR RMSD table (default {RMSD})")
    ap.add_argument("--extra-rmsd", default="classII_cdr_rmsd.csv",
                    help="CDR RMSD for triads absent from the published run "
                         "(from hla2_rmsd.nf); default classII_cdr_rmsd.csv")
    ap.add_argument("--merged", default=MERGED)
    ap.add_argument("--out", default=None, help="output CSV (default: overwrite --merged)")
    ap.add_argument("--triads", nargs="*", default=TRIADS)
    a = ap.parse_args()
    out = a.out or a.merged

    merged = pd.read_csv(a.merged)
    before = merged[merged["mhc_class"] == "I"].copy()

    rm = pd.read_csv(a.rmsd)
    rm["pdb"] = rm["pdb"].str.lower()
    rm = rm[["pdb", "mhc_class", RMSD_COL]].drop_duplicates("pdb")

    # The five triads absent from the published run get their RMSD from our own
    # hla2_rmsd.nf, which runs the SAME four processes in the same order
    # (FORMAT_TRUE_PDBS -> TCRDOCK_GEOM_FROM_PDB + TCRDOCK_GEOM_FROM_AF3_INFERENCE ->
    # COMPUTE_RMSD_AF3), so the column means the same thing. Appended rather than
    # merged, and a PDB present in both tables is an error: two sources for one number
    # is how a figure ends up quietly mixing provenances.
    if a.extra_rmsd and os.path.exists(a.extra_rmsd):
        ex = pd.read_csv(a.extra_rmsd)
        ex["pdb"] = ex["pdb"].str.lower()
        ex = ex[["pdb", "mhc_class", RMSD_COL]].drop_duplicates("pdb")
        both = sorted(set(ex["pdb"]) & set(rm["pdb"]))
        if both:
            sys.exit(f"{a.extra_rmsd} and {a.rmsd} both define: {both}")
        rm = pd.concat([rm, ex], ignore_index=True)
        print(f"+{len(ex)} rows from {a.extra_rmsd}")
    pae = load_pae(a.pae)

    want = [t.lower() for t in a.triads]
    have = set(merged["pdb"].str.lower())
    dup = [t for t in want if t in have]
    if dup:
        sys.exit(f"already present in {a.merged}, refusing to duplicate: {dup}")

    rows, no_rmsd, no_pae = [], [], []
    for t in want:
        r = rm[rm.pdb == t]
        p = pae[pae.pdb == t]
        if r.empty or pd.isna(r.iloc[0][RMSD_COL]):
            no_rmsd.append(t)
            continue
        if p.empty:
            no_pae.append(t)
            continue
        cls = r.iloc[0]["mhc_class"]
        if cls != "II":
            sys.exit(f"{t}: expected class II, table says {cls!r}")
        rows.append({"pdb": t, "cdr_rmsd": float(r.iloc[0][RMSD_COL]),
                     "pti_pae": float(p.iloc[0]["pti_pae"]), "mhc_class": "II"})

    assert len(rows) + len(no_rmsd) + len(no_pae) == len(want), "buckets do not sum"
    print(f"requested {len(want)}: added {len(rows)}, "
          f"no CDR RMSD {len(no_rmsd)} {no_rmsd}, no PTI-PAE {len(no_pae)} {no_pae}")
    if not rows:
        sys.exit("nothing to add")

    new = pd.concat([merged, pd.DataFrame(rows)], ignore_index=True)

    # Compare values exactly, but not dtypes: appending rows that leave tv/js/n empty
    # promotes the integer 'n' column to float, which is a storage change and not a
    # moved result.
    after = new[new["mhc_class"] == "I"].reset_index(drop=True)
    try:
        pd.testing.assert_frame_equal(before.reset_index(drop=True), after,
                                      check_dtype=False, check_exact=True)
    except AssertionError as e:
        sys.exit(f"class I rows changed — refusing to write\n{e}")
    print(f"class I unchanged ({len(after)} rows); wrote {out} with {len(new)} total")
    new.to_csv(out, index=False)

    for r in rows:
        print(f"  {r['pdb']}  cdr_rmsd={r['cdr_rmsd']:.3f}  "
              f"pti_pae={r['pti_pae']:.3f}  (31-pae={31 - r['pti_pae']:.2f})")


if __name__ == "__main__":
    main()
