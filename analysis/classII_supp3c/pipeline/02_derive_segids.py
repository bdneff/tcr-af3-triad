#!/usr/bin/env python3
"""02_derive_segids.py — attach deposited chain IDs to the staged triad parquet.

`format_true_pdbs.py` cuts the crystal apart with `select_atoms("segid <X> and name CA")`
and renames the result to canonical A/B/C/D/E, so it needs the **deposited chain ID** for
each of the five roles. John's TSV stores role *labels* ("alpha", "beta"), not chain IDs,
so they have to be derived from RCSB by matching sequences to polymer entities.

    python 02_derive_segids.py --in  hla2_triad.staged.parquet \
                               --out hla2_triad.segid.parquet

Conventions taken from the published table rather than invented — `8vd0`, a tethered
class II deposit that IS in the paper's set, reads:

    pdb   peptide_segid  mhc_1_segid  mhc_2_segid  tcr_1_segid  tcr_2_segid
    8vd0        C             A            C            D            E      <- tethered
    8pjg        C             A            B            D            E
    8trl        C             A            B            I            J      <- not D/E!

Two things follow. **A tethered peptide shares its chain with MHC beta** (`peptide_segid ==
mhc_2_segid`); `format_true_pdbs.py` splits them by residue range, which is what its
"multiple-chains-under-one-segid" branch is for. And **chain letters vary per deposit**
(`8trl` uses I/J for the TCR), so every one must be derived, never assumed.

Matching reuses the classification `02_xcheck_sequences.py` already validated on this set
(42/45 exact): a staged sequence is `exact` if it equals a deposited entity, or
`substring` if contained in one — the expected relation for the MHC beta chain of a
tethered construct, where the deposit is peptide+linker+beta.

**Multi-copy deposits are settled by geometry, never by index order.** When an entity
appears as several chains the crystal holds more than one complex (`9yaf` reads
`A,F  B,G  C,H  D,I  E,J`). Taking the first of each looks right, but `8trl` above pairs
MHC `A`/`B` with TCR `I`/`J`, so index order is not reliably copy order — and an MHC
paired with a TCR from the other copy is a structure that never existed, whose RMSD would
look perfectly reasonable. Each candidate index set is therefore checked against the
deposited coordinates: every chain must sit within 20 Å of the MHC α chain. The first set
that forms a real complex wins; if none does, the triad is refused rather than guessed.
"""
import argparse
import json
import sys
import time
import urllib.request

import polars as pl

ROLES = [("peptide_segid", "peptide"), ("mhc_1_segid", "mhc_1_seq"),
         ("mhc_2_segid", "mhc_2_seq"), ("tcr_1_segid", "tcr_1_seq"),
         ("tcr_2_segid", "tcr_2_seq")]


def entities(pdb):
    """[(entity_id, description, sequence, [auth chain ids])] for a PDB entry."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb.upper()}"
    with urllib.request.urlopen(url, timeout=30) as f:
        entry = json.load(f)
    out = []
    for eid in entry["rcsb_entry_container_identifiers"]["polymer_entity_ids"]:
        u2 = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb.upper()}/{eid}"
        with urllib.request.urlopen(u2, timeout=30) as f:
            pe = json.load(f)
        seq = pe["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "").upper()
        desc = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", "")
        ids = pe["rcsb_polymer_entity_container_identifiers"].get("auth_asym_ids", [])
        out.append((eid, desc, seq, sorted(ids)))
        time.sleep(0.05)
    return out


def ca_coords(pdb, cache={}):
    """{auth_chain_id: [(x, y, z), ...]} of CA atoms, parsed from the deposited mmCIF.

    Hand-parsed rather than via MDAnalysis so this has no structural dependency: we need
    five columns of one loop, and a missing optional dependency must not turn into a
    silently skipped check.
    """
    if pdb in cache:
        return cache[pdb]
    url = f"https://files.rcsb.org/download/{pdb.upper()}.cif"
    with urllib.request.urlopen(url, timeout=60) as f:
        text = f.read().decode("utf-8", "replace")

    out, cols, in_loop = {}, [], False
    for line in text.splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".", 1)[1])
            in_loop = True
            continue
        if in_loop and (line.startswith("ATOM") or line.startswith("HETATM")):
            f_ = line.split()
            if len(f_) < len(cols):
                continue
            r = dict(zip(cols, f_))
            if r.get("label_atom_id") != "CA":
                continue
            ch = r.get("auth_asym_id") or r.get("label_asym_id")
            try:
                out.setdefault(ch, []).append(
                    (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])))
            except (KeyError, ValueError):
                continue
        elif in_loop and line.startswith("#"):
            if out:
                break
    cache[pdb] = out
    return out


def min_dist(a, b):
    """Minimum CA-CA distance between two chains, in angstroms."""
    if not a or not b:
        return float("inf")
    best = float("inf")
    for x1, y1, z1 in a:
        for x2, y2, z2 in b:
            d = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
            if d < best:
                best = d
    return best ** 0.5


def contacts_ok(pdb, assign, cutoff=20.0):
    """Do the chosen chains actually form one complex?

    Index order across entities is NOT reliably copy order -- the published table has 8trl
    with MHC A/B and TCR I/J -- so a candidate assignment is checked geometrically. An
    MHC paired with a TCR from the other copy in the asymmetric unit yields a structure
    that never existed and an RMSD that looks entirely reasonable.
    """
    ca = ca_coords(pdb)
    m1 = ca.get(assign["mhc_1_segid"], [])
    pairs = [("tcr_1_segid", "TCR-alpha"), ("tcr_2_segid", "TCR-beta"),
             ("peptide_segid", "peptide"), ("mhc_2_segid", "MHC-beta")]
    worst = []
    for col, label in pairs:
        ch = assign[col]
        if ch == assign["mhc_1_segid"]:
            continue
        d = min_dist(m1, ca.get(ch, []))
        worst.append((label, ch, d))
        if d > cutoff:
            return False, worst
    return True, worst


def match(staged, ents):
    """(entity, how) for the deposited entity this staged sequence belongs to."""
    if not staged:
        return None, "empty"
    s = staged.upper()
    for e in ents:
        if s == e[2]:
            return e, "exact"
    for e in ents:
        if s in e[2]:
            return e, "substring"
    return None, "MISMATCH"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--accept-first", action="store_true",
                    help="on a multi-copy deposit, take the first chain of each entity "
                         "(ONLY after checking the copies form one complex)")
    a = ap.parse_args()

    df = pl.read_parquet(a.inp)
    rows, failures = [], []

    for row in df.iter_rows(named=True):
        pdb = row["pdb"]
        try:
            ents = entities(pdb)
        except Exception as e:                      # noqa: BLE001
            failures.append(f"{pdb}: RCSB fetch failed ({e})")
            continue
        print(f"\n{pdb}: {len(ents)} polymer entities")

        picked, ncopies, bad = {}, [], False
        for col, seqcol in ROLES:
            ent, how = match(row.get(seqcol), ents)
            if ent is None:
                failures.append(f"{pdb}/{col}: {how}")
                print(f"  {col:14s} {how}")
                bad = True
                continue
            eid, desc, seq, ids = ent
            picked[col] = ids
            ncopies.append(len(ids))
            print(f"  {col:14s} -> entity {eid} ({how}, copies={','.join(ids)})  {desc[:42]}")
        if bad:
            continue

        # Single copy throughout: nothing to choose.
        if max(ncopies) == 1:
            assigned = {c: v[0] for c, v in picked.items()}
        else:
            k = min(ncopies)
            assigned, chosen = None, None
            for i in range(k):
                cand = {c: (v[i] if len(v) > i else v[0]) for c, v in picked.items()}
                ok, dists = contacts_ok(pdb, cand)
                summary = "  ".join(f"{lab}({ch}) {d:.1f}A" for lab, ch, d in dists)
                print(f"    copy {i}: {' '.join(sorted(set(cand.values())))}  "
                      f"{'CONTACT' if ok else 'apart '}  {summary}")
                if ok and assigned is None:
                    assigned, chosen = cand, i
            if assigned is None:
                failures.append(
                    f"{pdb}: {max(ncopies)} copies in the deposit and no index set forms a "
                    "complex — chains would come from different copies. Inspect manually.")
                continue
            print(f"    -> copy {chosen} selected on geometry, not index order")

        if assigned.get("peptide_segid") == assigned.get("mhc_2_segid"):
            print("  (tethered: peptide shares MHC-beta's chain, as 8vd0 does)")

        rows.append({**row, **assigned})

    print(f"\n{len(rows)} of {df.height} triads assigned")
    if failures:
        print("\nunresolved:")
        for f in failures:
            print("  " + f)
        sys.exit("refusing to write a partial table — every triad lands in exactly one "
                 "bucket, and these are in neither")

    out = pl.DataFrame(rows)
    missing = [c for c, _ in ROLES if out[c].null_count()]
    if missing:
        sys.exit(f"null segids remain in {missing}")
    out.write_parquet(a.out)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
