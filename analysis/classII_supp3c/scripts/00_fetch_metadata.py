#!/usr/bin/env python3
"""00_fetch_metadata.py — verify the class II triads are post-AF3-cutoff, from RCSB.

`data/hla2_9_triads.tsv` was hand-assembled and carries no release dates, so "post-training-cutoff"
is an assertion until it is checked. That claim is the entire point of Supp Fig 3c -- the reviewer's
comment is precisely about it -- so it gets verified against RCSB rather than trusted, and the
verification is written to a file the figure can cite.

Also records resolution: two of the gliadin entries sit at 3.2 A, where CDR loop coordinates carry
real uncertainty. If one of them lands as an outlier, that is the first thing to check.

    python analysis/classII_supp3c/scripts/00_fetch_metadata.py

Writes analysis/classII_supp3c/data/pdb_metadata.csv. Read-only network access to RCSB.
"""
import csv, json, os, sys, time, urllib.request

CUTOFF = "2023-01-12"          # AF3 training cutoff, per the manuscript
# every class II input table: John's original plus any row we derived ourselves
TSVS = ["data/hla2_9_triads.tsv", "data/hla2_derived_extra.tsv"]
OUT = "analysis/classII_supp3c/data/pdb_metadata.csv"


def entry(pdb):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb.upper()}"
    with urllib.request.urlopen(url, timeout=30) as f:
        return json.load(f)


def main():
    ids = []
    for t in TSVS:
        if os.path.exists(t):
            ids += [r["job_name"] for r in csv.DictReader(open(t), delimiter="\t")]
    ids = list(dict.fromkeys(ids))          # preserve order, drop any duplicate
    rows, failed = [], []
    for i in ids:
        try:
            d = entry(i)
        except Exception as e:
            failed.append((i, str(e)))
            continue
        rel = (d.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10]
        res = (d.get("rcsb_entry_info") or {}).get("resolution_combined") or [None]
        rows.append(dict(
            pdb_id=i.lower(),
            release_date=rel,
            pre_af3_cutoff=str(bool(rel and rel < CUTOFF)),
            method=((d.get("exptl") or [{}])[0]).get("method", ""),
            resolution=(f"{res[0]:.2f}" if res[0] else ""),
            polymer_entity_count=(d.get("rcsb_entry_info") or {}).get("polymer_entity_count", ""),
            title=(d.get("struct") or {}).get("title", ""),
        ))
        time.sleep(0.15)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    # A row that is NOT post-cutoff would silently poison the figure's whole claim, so this is a
    # hard failure rather than a warning -- the same rule the MD side learned the hard way.
    bad = [r for r in rows if r["pre_af3_cutoff"] != "False"]
    print(f"{len(rows)}/{len(ids)} entries fetched -> {OUT}")
    if failed:
        print(f"FETCH FAILED for {len(failed)}: " + ", ".join(i for i, _ in failed))
    if bad:
        sys.exit("FATAL: not post-AF3-cutoff: " + ", ".join(f"{r['pdb_id']} ({r['release_date']})"
                                                            for r in bad))
    print(f"all {len(rows)} released on/after {CUTOFF} — post-cutoff claim verified")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
