#!/usr/bin/env python
"""
00_refresh_stcrdab.py
---------------------
The bulk db_summary.dat download URL from STCRDab is no longer available
(returns 404). However, per-PDB summary files ARE still available at:
    https://opig.stats.ox.ac.uk/webapps/stcrdab-stcrpred/summary/{pdbid}

Each per-PDB summary uses the same TSV schema as the bulk file. This script
fetches summaries for every PDB ID in the TCR3d CSVs and concatenates them
into a fresh db_summary.dat. The resulting file is a drop-in replacement
for the one in the repo, but covers the current STCRDab content (~634
PDBs as of mid-2026, up from the ~496 in the repo's snapshot).

Usage:
    python scripts/00_refresh_stcrdab.py \\
        --tcr3d_classI  data/tcr3d_classI_complexes.csv \\
        --tcr3d_classII data/tcr3d_classII_complexes.csv \\
        --out           data/db_summary_fresh.dat

Per-PDB summaries are cached in data/stcrdab_cache/ so reruns are cheap.
Failed fetches (404 = STCRDab really doesn't have it) are logged and
excluded from the output.
"""

import argparse
import sys
import time
from pathlib import Path

import polars as pl
import requests


STCRDAB_URL = "https://opig.stats.ox.ac.uk/webapps/stcrdab-stcrpred/summary/{pdb}"


def fetch_one(pdb: str, cache_dir: Path, retries: int = 2, sleep_s: float = 0.2):
    """Fetch one per-PDB summary, returning the raw TSV text (with header)
    or None on 404."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cf = cache_dir / f"{pdb.lower()}.tsv"
    if cf.exists():
        return cf.read_text()

    url = STCRDAB_URL.format(pdb=pdb.lower())
    last_err = None
    for attempt in range(retries + 1):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            cf.write_text(r.text)
            time.sleep(sleep_s)   # politeness
            return r.text
        except requests.RequestException as e:
            last_err = e
            time.sleep(1.0)
    print(f"  {pdb}: failed after {retries+1} attempts: {last_err}", file=sys.stderr)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tcr3d_classI",  type=Path, required=True)
    ap.add_argument("--tcr3d_classII", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data" / "db_summary_fresh.dat",
    )
    ap.add_argument(
        "--cache_dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data" / "stcrdab_cache",
    )
    ap.add_argument(
        "--also_keep_old",
        type=Path,
        default=None,
        help="Path to existing db_summary.dat. If provided, its entries are "
             "merged with the fresh fetches (new takes precedence on PDB ID conflict).",
    )
    args = ap.parse_args()

    # Load TCR3d
    ci = pl.read_csv(args.tcr3d_classI,  null_values=["NA","","N/A"])
    cii = pl.read_csv(args.tcr3d_classII, null_values=["NA","","N/A"])
    pdbs = sorted(
        set(ci["PDB ID"].str.to_lowercase().to_list())
        | set(cii["PDB ID"].str.to_lowercase().to_list())
    )
    print(f"TCR3d has {len(pdbs)} unique PDB IDs")

    # Fetch per-PDB summaries
    header = None
    body_lines = []
    fetched, missing = 0, 0
    for i, pdb in enumerate(pdbs, 1):
        if i % 25 == 0 or i == len(pdbs):
            print(f"  [{i}/{len(pdbs)}] fetched={fetched} missing={missing}")
        tsv = fetch_one(pdb, args.cache_dir)
        if tsv is None:
            missing += 1
            continue
        fetched += 1
        lines = tsv.strip().split("\n")
        if not lines:
            continue
        if header is None:
            header = lines[0]
        # Append data rows (skip the per-file header)
        body_lines.extend(lines[1:])

    print()
    print(f"  Fetched {fetched} PDBs, {missing} not in STCRDab")
    print(f"  Total rows in fresh fetch: {len(body_lines)}")

    # Optionally merge with old db_summary.dat (old entries that aren't in
    # the fresh fetch are kept; conflicting PDBs use the fresh version)
    if args.also_keep_old:
        old_lines = args.also_keep_old.read_text().strip().split("\n")
        old_header = old_lines[0]
        if header and old_header != header:
            print("  WARNING: old db_summary header differs from fresh — "
                  "schema may have drifted!")
            print(f"    old:   {old_header}")
            print(f"    fresh: {header}")
            print("    Going ahead anyway, using FRESH header.")
        if header is None:
            header = old_header

        # PDB IDs already in fresh body (first column of each line)
        fresh_pdbs = {line.split("\t")[0].lower() for line in body_lines}
        old_to_keep = [
            line for line in old_lines[1:]
            if line.split("\t")[0].lower() not in fresh_pdbs
        ]
        print(f"  Merging with old db_summary.dat:"
              f" {len(old_lines)-1} old rows,"
              f" {len(old_to_keep)} kept (not superseded by fresh)")
        body_lines = body_lines + old_to_keep

    if header is None:
        print("ERROR: no fresh data fetched and no --also_keep_old provided.",
              file=sys.stderr)
        sys.exit(1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(header + "\n" + "\n".join(body_lines) + "\n")
    unique_pdbs = len({line.split("\t")[0].lower() for line in body_lines})
    print()
    print(f"Wrote {args.out}: {len(body_lines)} rows across {unique_pdbs} unique PDBs")


if __name__ == "__main__":
    main()
