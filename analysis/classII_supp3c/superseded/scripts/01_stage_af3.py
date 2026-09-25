#!/usr/bin/env python3
"""01_stage_af3.py — turn John's class II TSV into AF3 jobs + a pipeline-schema manifest.

Two outputs, because the TSV serves two masters:

  af3_inputs/<pdb>.json    one AlphaFold 3 job per triad
  data/manifest.csv        the same rows translated into the column names the rest of this
                           repo uses (tcr3d_triad_sequences.csv), with RCSB release dates
                           merged in from 00_fetch_metadata.py

**Run configuration is fixed by the manuscript, not chosen here.** Methods, "AlphaFold 3
configuration": *"AlphaFold 3 was run with MSAs and templates turned on. No custom templates or
alignments were used. All triads were run with 5 seeds (1, 2, 3, 4, 5) and 1 diffusion sample...
The default number of recycles was used (10)."* So each JSON carries `modelSeeds` 1-5 and
deliberately OMITS `unpairedMsa`/`pairedMsa`/`templates` — their absence is what makes AF3 run its
own data pipeline and build real MSAs and template hits. The submit script must therefore NOT pass
`--norun_data_pipeline`; the JSON and the flag are a matched pair.

**Chain convention: A = MHC alpha, B = MHC beta, C = peptide, D = TCR alpha, E = TCR beta.**
C/D/E match the class I convention used everywhere else in this project (peptide, TCRa, TCRb) so
downstream chain mapping stays uniform; only A/B differ, because class II has two MHC chains where
class I has heavy chain + beta-2 microglobulin.

**A caution for whoever does the RMSD step.** Three of these deposited structures (8vd0, 8vq8,
9aud) tether the peptide to the MHC beta chain, so the crystal has FOUR polymer entities where the
prediction will have five, and several entries have two copies in the asymmetric unit. Chain
correspondence there cannot be positional -- see data/pdb_metadata.csv (`polymer_entity_count`)
and the note in docs/PROVENANCE.md.

    python analysis/classII_supp3c/scripts/01_stage_af3.py
    python analysis/classII_supp3c/scripts/01_stage_af3.py --include 27eb   # if John confirms

No dependencies beyond the standard library.
"""
import argparse, csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BASE = os.path.join(ROOT, "analysis", "classII_supp3c")
TSVS = [os.path.join(ROOT, "data", "hla2_9_triads.tsv"),
        os.path.join(ROOT, "data", "hla2_derived_extra.tsv")]
META = os.path.join(BASE, "data", "pdb_metadata.csv")

SEEDS = [1, 2, 3, 4, 5]              # manuscript Methods
AA = set("ACDEFGHIKLMNPQRSTVWY")
# (chain id, TSV column) in the canonical order
CHAINS = [("A", "mhc_1_seq"), ("B", "mhc_2_seq"), ("C", "peptide"),
          ("D", "tcr_1_seq"), ("E", "tcr_2_seq")]


def check(row):
    """Reject anything AF3 would silently mangle, loudly and before the job is written."""
    problems = []
    for cid, col in CHAINS:
        s = (row.get(col) or "").strip().upper()
        if not s:
            problems.append(f"chain {cid}: {col} is empty")
            continue
        bad = sorted(set(s) - AA)
        if bad:
            problems.append(f"chain {cid} ({col}): non-standard residues {bad} — AF3 needs a CCD "
                            f"code for these, not a bare letter")
    pep = (row.get("peptide") or "").strip().upper()
    for cid, col in (("A", "mhc_1_seq"), ("B", "mhc_2_seq")):
        if pep and pep in (row.get(col) or "").upper():
            problems.append(f"peptide also occurs inside {col} — tethered construct not split; "
                            f"submitting both would fold the peptide twice")
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", nargs="+", default=TSVS,
                    help="input tables; John's original and any derived rows")
    ap.add_argument("--out", default=os.path.join(BASE, "af3_inputs"))
    ap.add_argument("--manifest", default=os.path.join(BASE, "data", "manifest.csv"))
    a = ap.parse_args()

    rows, seen = [], set()
    for t in a.tsv:
        if not os.path.exists(t):
            print(f"note: {t} absent — skipping"); continue
        for r in csv.DictReader(open(t), delimiter="\t"):
            k = r["job_name"].strip().lower()
            if k in seen:
                print(f"note: {k} already staged from an earlier table — skipping dup")
                continue
            seen.add(k); rows.append(r)
    meta = {}
    if os.path.exists(META):
        meta = {m["pdb_id"]: m for m in csv.DictReader(open(META))}
    else:
        print(f"NOTE: {META} absent — run 00_fetch_metadata.py first for release dates")

    os.makedirs(a.out, exist_ok=True)
    os.makedirs(os.path.dirname(a.manifest), exist_ok=True)

    fatal, written, man = [], [], []
    for r in rows:
        pdb = r["job_name"].strip().lower()
        probs = check(r)
        if probs:
            fatal += [f"{pdb}: {p}" for p in probs]
            continue
        job = {"name": pdb, "modelSeeds": SEEDS,
               "sequences": [{"protein": {"id": cid,
                                          "sequence": r[col].strip().upper()}}
                             for cid, col in CHAINS],
               "dialect": "alphafold3", "version": 1}
        path = os.path.join(a.out, f"{pdb}.json")
        with open(path, "w") as f:
            json.dump(job, f, indent=2)
        written.append(pdb)

        m = meta.get(pdb, {})
        man.append(dict(
            pdb_id=pdb, mhc_class=r.get("mhc_class", "II"),
            release_date=m.get("release_date", ""),
            pre_af3_cutoff=m.get("pre_af3_cutoff", ""),
            resolution=m.get("resolution", ""),
            polymer_entity_count=m.get("polymer_entity_count", ""),
            peptide_seq=r["peptide"].strip().upper(),
            mhc_a_seq=r["mhc_1_seq"].strip().upper(),
            mhc_b_seq=r["mhc_2_seq"].strip().upper(),
            tcr_a_seq=r["tcr_1_seq"].strip().upper(),
            tcr_b_seq=r["tcr_2_seq"].strip().upper(),
            mhc_allele=f"{r.get('mhc_1_name','')}/{r.get('mhc_2_name','')}",
            mhc_species=r.get("mhc_1_species", ""),
            source=r.get("source", ""), cognate=r.get("cognate", ""),
        ))

    if fatal:
        print("REFUSING to stage — fix these first:")
        for p in fatal:
            print("   " + p)
        sys.exit(1)

    with open(a.manifest, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(man[0]))
        w.writeheader(); w.writerows(man)

    # count before concluding: a staged job that never became a file is the classic silent hole
    n_json = len([x for x in os.listdir(a.out) if x.endswith(".json")])
    print(f"staged {len(written)} job(s) -> {a.out}  ({n_json} json on disk)")
    print(f"manifest -> {a.manifest}")
    if n_json != len(rows):
        sys.exit(f"FATAL: {len(rows)} input rows but {n_json} JSON files written")
    missing_dates = [m["pdb_id"] for m in man if not m["release_date"]]
    if missing_dates:
        print(f"WARNING: no release date for {missing_dates} — post-cutoff status unverified")
    print(f"seeds {SEEDS}; MSA/template fields omitted so AF3 runs its data pipeline "
          f"(do NOT pass --norun_data_pipeline)")


if __name__ == "__main__":
    main()
