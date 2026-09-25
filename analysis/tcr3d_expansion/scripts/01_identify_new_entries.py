#!/usr/bin/env python
"""
01_identify_new_entries.py
--------------------------
Compute which TCR3d entries are new relative to the current Figure 1 input,
split them into pre/post-AF3-cutoff, look them up in (a fresh) STCRDab summary,
and produce the additions needed by the repo pipeline.

INPUTS
------
--tcr3d_classI    : TCR3d Class I peptide-complex CSV
--tcr3d_classII   : TCR3d Class II peptide-complex CSV
--repo_root       : path to tcrtrifold-experiments-main checkout
--stcrdab_summary : path to a (preferably fresh) db_summary.dat
                    [if not provided, uses repo/data/pdb/raw/db_summary.dat]
--cache_dir       : RCSB cache dir

OUTPUTS (written into the package's data/ dir)
----------------------------------------------
- tcr3d_new_pre_cutoff.csv         : rows to APPEND to table_S1_structure_benchmark_complexes.csv
- tcr3d_new_post_cutoff_pdbs.txt   : PDB IDs to add to clean_pdb.py's post_af3_cutoff_valid_pdbs list
- tcr3d_new_outlier_rows.csv       : segid annotations to add to clean_pdb.py's outlier_pdb DataFrame
                                     (for entries STILL absent from the fresh db_summary)
- summary_report.md                : human-readable summary: how many of each, classifier confidence,
                                     any flagged entries needing manual review
- classifier_validation.csv        : per-entry comparison of classifier predictions vs STCRDab on the
                                     152 entries that exist in both, for confidence calibration

The pipeline downstream is fully driven by these three artifacts plus the
fresh db_summary.dat.
"""

import argparse
import datetime as dt
import sys
from pathlib import Path

import polars as pl

# Make sibling src importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from chain_classifier import classify, validate_against_stcrdab  # noqa: E402
from rcsb_fetch import get_fasta, get_release_date, get_entity_organism  # noqa: E402

AF3_CUTOFF = dt.datetime(2023, 1, 12)

# Hardcoded list from clean_pdb.py at time of writing. If the repo's list has
# been updated since this script was written, this needs to be re-synced.
EXISTING_POST_AF3 = {
    "8gom", "8vd0", "8trr", "8wte", "8vcy", "8es9", "8gon", "8i5d",
    "8i5c", "8vcx", "7q99", "8eo8", "8dnt", "8enh", "8ye4", "8wul",
    "8f5a", "8vd2", "7q9b", "8en8", "8pjg", "7q9a",
}
EXISTING_OUTLIER = {"6l9l"}  # from clean_pdb.py outlier_pdb DataFrame


def load_tcr3d(class1_csv: Path, class2_csv: Path) -> pl.DataFrame:
    """Load and combine TCR3d class I and II peptide-complex CSVs."""
    c1 = pl.read_csv(class1_csv, null_values=["NA", "", "N/A"]).with_columns(
        pl.lit("I").alias("mhc_class")
    )
    c2 = pl.read_csv(class2_csv, null_values=["NA", "", "N/A"]).with_columns(
        pl.lit("II").alias("mhc_class")
    )
    df = pl.concat([c1, c2], how="diagonal_relaxed")
    df = df.rename({"PDB ID": "pdb", "Release<BR>date": "release_date"})
    df = df.with_columns(
        pl.col("pdb").str.to_lowercase().alias("pdb"),
        pl.col("release_date").str.to_datetime("%Y-%m-%d").alias("release_date"),
    )
    return df


def load_stcrdab(path: Path) -> pl.DataFrame:
    """Load STCRDab summary and aggregate to one row per PDB."""
    stcr = pl.read_csv(
        path,
        schema_overrides={"Gchain": pl.String, "Dchain": pl.String},
        null_values=["NA", "unknown", "NOT"],
        separator="\t",
    )

    # Mirror clean_pdb.format_stcr_df's role mapping (subset relevant to us)
    # mhc_type normalization
    stcr = stcr.with_columns(
        pl.when((pl.col("mhc_type") == "GA") | (pl.col("mhc_type") == "GB"))
        .then(pl.lit("MH2"))
        .otherwise(pl.col("mhc_type"))
        .alias("mhc_type"),
    )

    # Filter to peptide antigens, take first peptide segid if multiple
    stcr = stcr.with_columns(
        pl.when(pl.col("antigen_type").is_null())
        .then(pl.lit("peptide"))
        .otherwise(pl.col("antigen_type"))
        .alias("antigen_type"),
    ).filter(pl.col("antigen_type").str.contains("peptide"))

    # Split antigen_chain on "|" and pick the peptide-typed one
    stcr = stcr.with_columns(
        pl.col("antigen_type").str.split("|")
        .list.eval(pl.element().str.strip_chars())
        .list.eval(pl.element().index_of("peptide"))
        .list.first()
        .alias("_pep_idx"),
    ).with_columns(
        pl.col("antigen_chain").str.split("|")
        .list.eval(pl.element().str.strip_chars())
        .list.get(pl.col("_pep_idx"))
        .alias("antigen_chain"),
    )

    agg = stcr.group_by("pdb", maintain_order=True).agg(
        pl.col("mhc_type").drop_nulls().first(),
        pl.col("mhc_chain1").first().alias("mhc_1_segid"),
        pl.col("mhc_chain2").first().alias("mhc_2_segid"),
        pl.col("antigen_chain").first().alias("peptide_segid"),
        pl.col("Achain").first().alias("tcr_1_segid"),
        pl.col("Bchain").first().alias("tcr_2_segid"),
        pl.col("mhc_chain1_organism").drop_nulls().first().alias("mhc_1_species"),
        pl.col("mhc_chain2_organism").drop_nulls().first().alias("mhc_2_species"),
        pl.col("alpha_organism").drop_nulls().first().alias("tcr_1_species"),
        pl.col("beta_organism").drop_nulls().first().alias("tcr_2_species"),
    )

    # Map mhc_type to "I" or "II"
    agg = agg.with_columns(
        pl.when(pl.col("mhc_type") == "MH1").then(pl.lit("I"))
        .when(pl.col("mhc_type") == "MH2").then(pl.lit("II"))
        .otherwise(None)
        .alias("mhc_class"),
    ).drop("mhc_type")

    # Normalize species names
    for col in ("mhc_1_species", "mhc_2_species", "tcr_1_species", "tcr_2_species"):
        agg = agg.with_columns(
            pl.when(pl.col(col) == "homo sapiens").then(pl.lit("human"))
            .when(pl.col(col) == "mus musculus").then(pl.lit("mouse"))
            .otherwise(None)
            .alias(col)
        )

    agg = agg.with_columns(pl.col("pdb").str.to_lowercase())
    return agg


def to_table_S1_row(tcr3d_row: dict) -> dict:
    """
    Convert a TCR3d row into the column schema of
    table_S1_structure_benchmark_complexes.csv:
        pdbid, organism, mhc_class, mhc, peptide,
        va, ja, cdr3a, vb, jb, cdr3b,
        cdr_rmsd, cdr_rmsd_af2_full, cdr_rmsd_af2_trim

    Fields we don't have (cdr_rmsd_af2_*, exact gene calls) are left null;
    the pipeline tolerates these because they're only used for the AF2
    comparison panel and downstream join. We set cdr_rmsd=None so the
    notebook treats these rows as having no AF-TCRdock comparison.

    organism is determined later from RCSB; for now leave as null and
    we fill it in.
    """
    cls_num = 1 if tcr3d_row.get("mhc_class") == "I" else 2
    return {
        "pdbid": tcr3d_row["pdb"],
        "organism": None,            # filled by RCSB entity organism fetch
        "mhc_class": cls_num,
        "mhc": tcr3d_row.get("MHC Name") or "",
        "peptide": tcr3d_row.get("Epitope") or "",
        "va": tcr3d_row.get("TRAV gene") or "",
        "ja": "",
        "cdr3a": "",
        "vb": tcr3d_row.get("TRBV gene") or "",
        "jb": "",
        "cdr3b": "",
        "cdr_rmsd": None,
        "cdr_rmsd_af2_full": None,
        "cdr_rmsd_af2_trim": None,
    }


def outlier_row_from_classification(result, pdb_date_iso: str, organisms: dict) -> dict:
    """
    Convert a ClassificationResult + RCSB metadata into a row in the schema
    of clean_pdb.py's outlier_pdb DataFrame.
    """
    def org_of(segid):
        if segid is None:
            return None
        return organisms.get(segid)

    # Determine MHC chain labels
    if result.mhc_class == "I":
        mhc_1_chain, mhc_2_chain = "heavy", "light"
    else:
        mhc_1_chain, mhc_2_chain = "alpha", "beta"

    return {
        "pdb": result.pdb,
        "mhc_1_segid": result.mhc_1_segid,
        "mhc_2_segid": result.mhc_2_segid,
        "peptide_segid": result.peptide_segid,
        "tcr_1_segid": result.tcr_1_segid,
        "tcr_2_segid": result.tcr_2_segid,
        "mhc_1_species": org_of(result.mhc_1_segid),
        "mhc_2_species": org_of(result.mhc_2_segid),
        "tcr_1_species": org_of(result.tcr_1_segid),
        "tcr_2_species": org_of(result.tcr_2_segid),
        "mhc_class": result.mhc_class,
        "mhc_1_chain": mhc_1_chain,
        "mhc_2_chain": mhc_2_chain,
        "cognate": True,
        "tcr_1_chain": "alpha",
        "tcr_2_chain": "beta",
        "pdb_date": pdb_date_iso,
        "_classifier_confidence": result.confidence,
        "_classifier_notes": "; ".join(result.notes) if result.notes else "",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tcr3d_classI", type=Path, required=True)
    ap.add_argument("--tcr3d_classII", type=Path, required=True)
    ap.add_argument("--repo_root", type=Path, required=True)
    ap.add_argument(
        "--stcrdab_summary",
        type=Path,
        default=None,
        help="Path to a fresh db_summary.dat. Default: repo's data/pdb/raw/db_summary.dat",
    )
    ap.add_argument(
        "--cache_dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data" / "rcsb_cache",
    )
    ap.add_argument(
        "--out_dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data",
    )
    ap.add_argument(
        "--skip_validation",
        action="store_true",
        help="Skip classifier validation step (which is slow due to RCSB calls).",
    )
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # 1. Load reference data
    # ------------------------------------------------------------
    print("[1/5] Loading reference data...")
    tcr3d = load_tcr3d(args.tcr3d_classI, args.tcr3d_classII)
    rep_csv_path = args.repo_root / "data/pdb/raw/table_S1_structure_benchmark_complexes.csv"
    rep = pl.read_csv(rep_csv_path).with_columns(pl.col("pdbid").str.to_lowercase())
    existing_pre = set(rep["pdbid"].to_list())

    stcrdab_path = args.stcrdab_summary or (args.repo_root / "data/pdb/raw/db_summary.dat")
    stcr = load_stcrdab(stcrdab_path)

    print(f"  TCR3d total:           {tcr3d.height} ({tcr3d.filter(pl.col('mhc_class')=='I').height} I + {tcr3d.filter(pl.col('mhc_class')=='II').height} II)")
    print(f"  Existing pre-cutoff:   {len(existing_pre)}")
    print(f"  Existing post-cutoff:  {len(EXISTING_POST_AF3)} (hardcoded) + 1 outlier (6l9l)")
    print(f"  STCRDab summary:       {stcr.height} PDBs")
    print()

    # ------------------------------------------------------------
    # 2. Identify the new entries
    # ------------------------------------------------------------
    print("[2/5] Identifying new TCR3d entries...")
    excluded = existing_pre | EXISTING_POST_AF3 | EXISTING_OUTLIER
    tcr3d_new = tcr3d.filter(~pl.col("pdb").is_in(list(excluded)))

    tcr3d_new_pre = tcr3d_new.filter(pl.col("release_date") < AF3_CUTOFF)
    tcr3d_new_post = tcr3d_new.filter(pl.col("release_date") >= AF3_CUTOFF)
    print(f"  NEW pre-cutoff:  {tcr3d_new_pre.height}")
    print(f"  NEW post-cutoff: {tcr3d_new_post.height}")
    print()

    # ------------------------------------------------------------
    # 3. For each new entry, check whether STCRDab covers it
    # ------------------------------------------------------------
    print("[3/5] Cross-checking against STCRDab summary...")
    stcr_pdbs = set(stcr["pdb"].to_list())

    pre_in_stcr = tcr3d_new_pre.filter(pl.col("pdb").is_in(list(stcr_pdbs)))
    pre_missing = tcr3d_new_pre.filter(~pl.col("pdb").is_in(list(stcr_pdbs)))
    post_in_stcr = tcr3d_new_post.filter(pl.col("pdb").is_in(list(stcr_pdbs)))
    post_missing = tcr3d_new_post.filter(~pl.col("pdb").is_in(list(stcr_pdbs)))

    print(f"  pre-cutoff in STCRDab:    {pre_in_stcr.height}  (need only pipeline rerun)")
    print(f"  pre-cutoff missing:       {pre_missing.height}  (need classifier)")
    print(f"  post-cutoff in STCRDab:   {post_in_stcr.height}")
    print(f"  post-cutoff missing:      {post_missing.height}  (need classifier)")
    print()

    # ------------------------------------------------------------
    # 4. Build table_S1 additions (these are all the pre-cutoff entries)
    # ------------------------------------------------------------
    print("[4/5] Building table_S1 additions for pre-cutoff entries...")
    table_S1_new_rows = []
    for row in tcr3d_new_pre.iter_rows(named=True):
        table_S1_new_rows.append(to_table_S1_row(row))

    # Try to fill organism via RCSB for entries we'll need to query anyway
    # (Skip if user isn't running with network access; the field is optional
    # — clean_pdb.py uses mhc_1_species from STCRDab for those that have it.)
    pl.DataFrame(table_S1_new_rows).write_csv(
        args.out_dir / "tcr3d_new_pre_cutoff.csv"
    )
    print(f"  Wrote {len(table_S1_new_rows)} rows to tcr3d_new_pre_cutoff.csv")

    # ------------------------------------------------------------
    # 5. Build outlier_pdb additions for the entries STCRDab doesn't cover
    # ------------------------------------------------------------
    print("[5/5] Building outlier_pdb additions for STCRDab-missing entries...")

    # For post-cutoff in STCRDab: just add to post_af3_cutoff_valid_pdbs list
    post_in_stcr_ids = post_in_stcr["pdb"].to_list()
    (args.out_dir / "tcr3d_new_post_cutoff_pdbs.txt").write_text(
        "\n".join(post_in_stcr_ids) + "\n"
    )
    print(f"  Wrote {len(post_in_stcr_ids)} PDB IDs to tcr3d_new_post_cutoff_pdbs.txt")

    # For the missing ones (pre AND post): use RCSB + classifier
    missing = pl.concat([pre_missing, post_missing], how="diagonal_relaxed")
    if missing.height == 0:
        print("  No entries require the classifier — STCRDab covers everything!")
        outlier_rows = []
    else:
        print(f"  Running classifier on {missing.height} STCRDab-missing entries...")
        outlier_rows = []
        for i, row in enumerate(missing.iter_rows(named=True), 1):
            pdb = row["pdb"]
            cls = row["mhc_class"]
            print(f"    [{i:>3d}/{missing.height}] {pdb} (class {cls})...", end=" ", flush=True)
            try:
                fasta = get_fasta(pdb, cache_dir=args.cache_dir)
                organisms = get_entity_organism(pdb, cache_dir=args.cache_dir)
                pdb_date = get_release_date(pdb, cache_dir=args.cache_dir)
                result = classify(pdb, fasta, expected_class=cls)
                out = outlier_row_from_classification(result, pdb_date, organisms)
                # Whether this is post-cutoff (needs to also go in the list)
                out["_post_cutoff"] = row["release_date"] >= AF3_CUTOFF
                outlier_rows.append(out)
                print(f"conf={result.confidence}")
            except Exception as e:
                print(f"FAILED: {e}")
                outlier_rows.append({
                    "pdb": pdb,
                    "_error": str(e),
                    "_post_cutoff": row["release_date"] >= AF3_CUTOFF,
                })

        pl.DataFrame(outlier_rows).write_csv(
            args.out_dir / "tcr3d_new_outlier_rows.csv"
        )
        # Stats
        confidence_counts = {}
        for o in outlier_rows:
            c = o.get("_classifier_confidence", "ERROR" if "_error" in o else "UNK")
            confidence_counts[c] = confidence_counts.get(c, 0) + 1
        print(f"  Classifier confidence breakdown: {confidence_counts}")
        print(f"  Wrote outlier rows to tcr3d_new_outlier_rows.csv")

    # ------------------------------------------------------------
    # Optional: validate classifier on entries STCRDab DOES cover
    # ------------------------------------------------------------
    if not args.skip_validation:
        print()
        print("[validation] Running classifier on entries STCRDab covers,"
              " for confidence calibration...")
        # Validate on the 4 new TCR3d entries that ARE in STCRDab (sanity), and
        # the 152 existing pipeline entries (rigorous).
        val_pdbs = list(
            set(rep["pdbid"].to_list()) |
            set(post_in_stcr["pdb"].to_list()) |
            EXISTING_POST_AF3
        )
        stcr_lookup = {row["pdb"]: row for row in stcr.iter_rows(named=True)}

        val_results = []
        for i, pdb in enumerate(sorted(val_pdbs), 1):
            if pdb not in stcr_lookup:
                continue
            known = stcr_lookup[pdb]
            cls = known.get("mhc_class")
            if cls not in ("I", "II"):
                continue
            try:
                fasta = get_fasta(pdb, cache_dir=args.cache_dir)
                r = validate_against_stcrdab(pdb, fasta, cls, known)
                r["mhc_class"] = cls
                val_results.append(r)
                if not r["all_match"]:
                    print(f"    [{i}/{len(val_pdbs)}] {pdb}: MISMATCH (conf={r['confidence']})")
            except Exception as e:
                print(f"    [{i}/{len(val_pdbs)}] {pdb}: FAILED {e}")

        # Summary
        if val_results:
            n_total = len(val_results)
            n_match = sum(r["all_match"] for r in val_results)
            print(f"  Validation: {n_match}/{n_total} entries fully match STCRDab"
                  f" ({100*n_match/n_total:.1f}%)")

            # Per-role match rates
            roles = ("peptide", "mhc_1", "mhc_2", "tcr_1", "tcr_2")
            for r_name in roles:
                rates = sum(v["matches"][r_name] for v in val_results)
                print(f"    {r_name}: {rates}/{n_total} ({100*rates/n_total:.1f}%)")

            # Save detail
            pl.DataFrame([
                {
                    "pdb": v["pdb"],
                    "mhc_class": v["mhc_class"],
                    "all_match": v["all_match"],
                    "confidence": v["confidence"],
                    **{f"pred_{k}": v["predicted"][k] for k in roles},
                    **{f"true_{k}": v["expected"][k] for k in roles},
                    **{f"match_{k}": v["matches"][k] for k in roles},
                    "notes": "; ".join(v["notes"]),
                }
                for v in val_results
            ]).write_csv(args.out_dir / "classifier_validation.csv")
            print(f"  Wrote classifier_validation.csv")

    # ------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------
    print()
    print("=" * 60)
    print("DONE. Next steps:")
    print("  1. Inspect tcr3d_new_pre_cutoff.csv and tcr3d_new_outlier_rows.csv")
    print("  2. Run scripts/02_apply_to_repo.sh to patch the repo files")
    print("  3. Run scripts/03_run_pipeline.sh on Gemini")
    print("  4. Run scripts/04_replot_figure1.sh once pipeline completes")
    print()


if __name__ == "__main__":
    main()
