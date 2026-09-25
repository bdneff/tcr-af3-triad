#!/usr/bin/env python
"""
02_apply_to_repo.py
-------------------
Patch the tcrtrifold-experiments repo with the new entries identified in
step 01. Three changes are applied (each guarded against double-application):

  (A) Append new pre-cutoff rows to
        data/pdb/raw/table_S1_structure_benchmark_complexes.csv

  (B) Inject new PDB IDs into clean_pdb.py's `post_af3_cutoff_valid_pdbs` list

  (C) For STCRDab-missing entries, append new rows to clean_pdb.py's
        `outlier_pdb` polars DataFrame

Before any change, the original file is backed up to:
    <repo>/<file>.tcr3d_backup.<timestamp>

The script is idempotent: if a PDB ID already appears in the destination, it
is skipped (with a log message). Re-running is safe.

A --dry_run flag shows what would change without writing.
"""

import argparse
import datetime as dt
import re
import shutil
from pathlib import Path

import polars as pl


# Markers inserted into clean_pdb.py so the script can find its own injection
# sites on subsequent runs (idempotency).
MARK_PDB_LIST_BEGIN = "    # === TCR3D_EXPANSION_PDB_LIST_BEGIN ==="
MARK_PDB_LIST_END   = "    # === TCR3D_EXPANSION_PDB_LIST_END ==="
MARK_OUTLIER_BEGIN  = "    # === TCR3D_EXPANSION_OUTLIER_BEGIN ==="
MARK_OUTLIER_END    = "    # === TCR3D_EXPANSION_OUTLIER_END ==="


def backup(path: Path):
    """Create a backup with a timestamp."""
    ts = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    bk = path.with_suffix(path.suffix + f".tcr3d_backup.{ts}")
    shutil.copy2(path, bk)
    return bk


def _read_text(p: Path) -> str:
    return p.read_text()


def _write_text(p: Path, text: str, dry_run: bool):
    if dry_run:
        print(f"    [DRY] would write to {p}")
    else:
        p.write_text(text)


# ----------------------------------------------------------------------------
# (A) Append rows to table_S1_structure_benchmark_complexes.csv
# ----------------------------------------------------------------------------

def patch_table_S1(repo_root: Path, pre_csv: Path, dry_run: bool):
    table_path = repo_root / "data/pdb/raw/table_S1_structure_benchmark_complexes.csv"

    existing = pl.read_csv(table_path).with_columns(pl.col("pdbid").str.to_lowercase())
    existing_ids = set(existing["pdbid"].to_list())

    new = pl.read_csv(pre_csv).with_columns(pl.col("pdbid").str.to_lowercase())

    # filter out anything already there
    to_add = new.filter(~pl.col("pdbid").is_in(list(existing_ids)))
    n_skipped = new.height - to_add.height

    print(f"[A] table_S1: {to_add.height} new rows ({n_skipped} already present, skipped)")
    if to_add.height == 0:
        return 0

    if not dry_run:
        backup(table_path)
    combined = pl.concat([existing, to_add], how="diagonal_relaxed")
    if dry_run:
        print(f"    [DRY] would write {combined.height} total rows to {table_path}")
    else:
        combined.write_csv(table_path)
        print(f"    wrote {combined.height} total rows to {table_path}")
    return to_add.height


# ----------------------------------------------------------------------------
# (B) Inject post-cutoff PDB IDs into clean_pdb.py's list
# ----------------------------------------------------------------------------

def patch_post_af3_list(repo_root: Path, post_txt: Path, dry_run: bool):
    clean_py = repo_root / "workflows/subworkflows/local/cleaning/resources/usr/bin/clean_pdb.py"
    text = _read_text(clean_py)

    new_ids = [line.strip() for line in post_txt.read_text().splitlines() if line.strip()]

    # Find the existing list
    m = re.search(
        r"(post_af3_cutoff_valid_pdbs\s*=\s*\[)([^\]]*?)(\])",
        text,
        flags=re.S,
    )
    if not m:
        raise RuntimeError(
            "Couldn't locate post_af3_cutoff_valid_pdbs in clean_pdb.py — "
            "the upstream code may have been refactored. Aborting."
        )

    body = m.group(2)
    existing_ids = set(re.findall(r"\"([^\"]+)\"", body))
    to_add = [pid for pid in new_ids if pid not in existing_ids]

    print(f"[B] post_af3_cutoff_valid_pdbs: {len(to_add)} new IDs"
          f" ({len(new_ids)-len(to_add)} already present, skipped)")
    if not to_add:
        return 0

    if not dry_run:
        backup(clean_py)

    # Build replacement: keep existing body, add new IDs with a comment marker
    new_body = body.rstrip().rstrip(",").rstrip()
    new_body += ",\n"
    new_body += f"    # --- TCR3d expansion (added {dt.date.today().isoformat()}) ---\n"
    for pid in to_add:
        new_body += f'    "{pid}",\n'

    new_text = text[:m.start(2)] + new_body + text[m.end(2):]
    _write_text(clean_py, new_text, dry_run)
    return len(to_add)


# ----------------------------------------------------------------------------
# (C) Inject new outlier rows into clean_pdb.py's outlier_pdb DataFrame
# ----------------------------------------------------------------------------

OUTLIER_TEMPLATE = '''        {{
            "pdb": "{pdb}",
            "mhc_1_segid": {mhc_1_segid!r},
            "mhc_2_segid": {mhc_2_segid!r},
            "peptide_segid": {peptide_segid!r},
            "tcr_1_segid": {tcr_1_segid!r},
            "tcr_2_segid": {tcr_2_segid!r},
            "mhc_1_species": {mhc_1_species!r},
            "mhc_2_species": {mhc_2_species!r},
            "tcr_1_species": {tcr_1_species!r},
            "tcr_2_species": {tcr_2_species!r},
            "mhc_class": "{mhc_class}",
            "mhc_1_chain": "{mhc_1_chain}",
            "mhc_2_chain": {mhc_2_chain},
            "cognate": True,
            "tcr_1_chain": "alpha",
            "tcr_2_chain": "beta",
            "pdb_date": extract_pdb_date({{"pdb": "{pdb}"}})["pdb_date"],
        }},'''


def patch_outlier_pdb(repo_root: Path, outlier_csv: Path, dry_run: bool,
                      min_confidence: str = "low"):
    """
    Append rows from outlier_csv into clean_pdb.py's `outlier_pdb` DataFrame
    construction.

    min_confidence: "high"|"medium"|"low" — rows below this are skipped (logged).
    """
    clean_py = repo_root / "workflows/subworkflows/local/cleaning/resources/usr/bin/clean_pdb.py"
    text = _read_text(clean_py)

    if not outlier_csv.exists():
        print("[C] No outlier_pdb rows to add (STCRDab covers everything new).")
        return 0

    rows = pl.read_csv(outlier_csv).to_dicts()
    if not rows:
        print("[C] outlier_csv exists but is empty.")
        return 0

    # Skip error rows
    rows = [r for r in rows if not r.get("_error")]

    # Confidence filter
    rank = {"high": 3, "medium": 2, "low": 1}
    min_rank = rank[min_confidence]
    kept, skipped = [], []
    for r in rows:
        c = r.get("_classifier_confidence", "low")
        if rank.get(c, 0) >= min_rank:
            kept.append(r)
        else:
            skipped.append(r)

    if skipped:
        print(f"[C] Skipping {len(skipped)} entries below confidence threshold"
              f" ({min_confidence})")
        for r in skipped[:5]:
            print(f"     - {r['pdb']}: conf={r.get('_classifier_confidence')}, "
                  f"notes={r.get('_classifier_notes','')[:80]}")
        if len(skipped) > 5:
            print(f"     ... ({len(skipped)-5} more — see {outlier_csv})")

    # Idempotency: locate the outlier_pdb DataFrame literal
    m = re.search(
        r"(outlier_pdb\s*=\s*pl\.DataFrame\(\s*\[)(.*?)(\s*\]\s*\)\.with_columns)",
        text,
        flags=re.S,
    )
    if not m:
        raise RuntimeError(
            "Couldn't locate outlier_pdb in clean_pdb.py. Aborting."
        )

    body = m.group(2)
    existing_pdbs = set(re.findall(r'"pdb"\s*:\s*"([^"]+)"', body))
    to_add = [r for r in kept if r["pdb"] not in existing_pdbs]

    print(f"[C] outlier_pdb: {len(to_add)} new rows"
          f" ({len(kept)-len(to_add)} already present, skipped)")
    if not to_add:
        return 0

    if not dry_run:
        backup(clean_py)

    # Build rendered rows
    rendered = []
    rendered.append(
        f"\n        # --- TCR3d expansion (added {dt.date.today().isoformat()}) ---"
    )
    for r in to_add:
        mhc_class = r.get("mhc_class") or "I"
        mhc_1_chain = "heavy" if mhc_class == "I" else "alpha"
        if mhc_class == "I":
            mhc_2_chain_literal = "None"
        else:
            mhc_2_chain_literal = '"beta"'
        rendered.append(
            OUTLIER_TEMPLATE.format(
                pdb=r["pdb"],
                mhc_1_segid=r.get("mhc_1_segid"),
                mhc_2_segid=r.get("mhc_2_segid"),
                peptide_segid=r.get("peptide_segid"),
                tcr_1_segid=r.get("tcr_1_segid"),
                tcr_2_segid=r.get("tcr_2_segid"),
                mhc_1_species=r.get("mhc_1_species"),
                mhc_2_species=r.get("mhc_2_species"),
                tcr_1_species=r.get("tcr_1_species"),
                tcr_2_species=r.get("tcr_2_species"),
                mhc_class=mhc_class,
                mhc_1_chain=mhc_1_chain,
                mhc_2_chain=mhc_2_chain_literal,
            )
        )

    new_body = body.rstrip() + "\n" + "\n".join(rendered) + "\n    "
    new_text = text[:m.start(2)] + new_body + text[m.end(2):]
    _write_text(clean_py, new_text, dry_run)
    return len(to_add)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo_root", type=Path, required=True)
    ap.add_argument(
        "--data_dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data",
        help="Where step 01 wrote its outputs.",
    )
    ap.add_argument(
        "--min_confidence",
        choices=("high", "medium", "low"),
        default="medium",
        help="Minimum classifier confidence to accept for outlier_pdb entries.",
    )
    ap.add_argument("--dry_run", action="store_true")
    args = ap.parse_args()

    pre_csv = args.data_dir / "tcr3d_new_pre_cutoff.csv"
    post_txt = args.data_dir / "tcr3d_new_post_cutoff_pdbs.txt"
    outlier_csv = args.data_dir / "tcr3d_new_outlier_rows.csv"

    print(f"Patching repo at {args.repo_root}")
    print(f"  (dry_run={args.dry_run})")
    print()
    n_a = patch_table_S1(args.repo_root, pre_csv, args.dry_run)
    n_b = patch_post_af3_list(args.repo_root, post_txt, args.dry_run)
    n_c = patch_outlier_pdb(args.repo_root, outlier_csv, args.dry_run, args.min_confidence)

    print()
    print("=" * 60)
    print(f"Summary: added {n_a} pre-cutoff rows, {n_b} post-cutoff IDs,"
          f" {n_c} outlier rows")
    print("Backups written next to each modified file (look for .tcr3d_backup.*).")


if __name__ == "__main__":
    main()
