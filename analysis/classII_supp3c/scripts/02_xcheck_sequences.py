#!/usr/bin/env python3
"""02_xcheck_sequences.py — verify the staged sequences against the deposited structures.

`data/hla2_9_triads.tsv` is hand-assembled. The last hand-assembled sequence table in this project
carried two defects that folded the wrong antigen and were only caught by cross-checking against
the crystals (`docs/PROVENANCE.md`: 3D3V/3D39 lost the P5 anchor phenylalanine; 7BYD was a 4-mer
padded to an 8-mer). This runs the same check before any GPU time is spent.

For each triad, every chain is compared against the deposited polymer entities from RCSB:

  exact          the staged sequence is a deposited entity, verbatim
  substring      the staged sequence is contained in a deposited entity — EXPECTED for the MHC beta
                 chain of a tethered construct, where the deposit is peptide+linker+beta and the
                 staged sequence is the beta portion alone; suspicious anywhere else
  prefix/suffix  staged is a truncation of a deposited entity — flags dropped terminal residues
  MISMATCH       nothing in the deposit contains it. Investigate before folding.

A mismatch is not automatically fatal -- constructs legitimately differ from deposits (expression
tags, engineered disulfides, truncated stalks) -- but every one must be understood rather than
assumed, which is why they are printed rather than silently tolerated.

    python analysis/classII_supp3c/scripts/02_xcheck_sequences.py

Read-only network access to RCSB.
"""
import csv, json, os, sys, time, urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
TSVS = [os.path.join(ROOT, "data", "hla2_9_triads.tsv"),
        os.path.join(ROOT, "data", "hla2_derived_extra.tsv")]
COLS = [("MHC-a", "mhc_1_seq"), ("MHC-b", "mhc_2_seq"), ("peptide", "peptide"),
        ("TCR-a", "tcr_1_seq"), ("TCR-b", "tcr_2_seq")]


def entities(pdb):
    """[(entity_id, description, canonical one-letter sequence)] for a PDB entry."""
    u = f"https://data.rcsb.org/rest/v1/core/entry/{pdb.upper()}"
    with urllib.request.urlopen(u, timeout=30) as f:
        e = json.load(f)
    out = []
    for eid in e["rcsb_entry_container_identifiers"]["polymer_entity_ids"]:
        u2 = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb.upper()}/{eid}"
        with urllib.request.urlopen(u2, timeout=30) as f:
            pe = json.load(f)
        seq = pe["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "").upper()
        desc = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", "")
        out.append((eid, desc, seq))
        time.sleep(0.05)
    return out


def classify(staged, ents):
    """How does the staged sequence relate to the best-matching deposited entity?"""
    for eid, desc, seq in ents:
        if staged == seq:
            return "exact", eid, desc, ""
    for eid, desc, seq in ents:
        if staged and staged in seq:
            off = seq.index(staged)
            extra = len(seq) - len(staged)
            return "substring", eid, desc, f"offset {off}, deposit has {extra} more residue(s)"
    for eid, desc, seq in ents:
        if staged and (seq.startswith(staged[:30]) or seq.endswith(staged[-30:])):
            return "partial", eid, desc, f"staged {len(staged)} vs deposit {len(seq)}"
    return "MISMATCH", "", "", ""


def main():
    rows = []
    for t in TSVS:
        if os.path.exists(t):
            rows += list(csv.DictReader(open(t), delimiter="\t"))
    # NOTE: 27eb's sequences were themselves derived from RCSB (03_derive_27eb.py), so
    # its row here is a closed loop and will always read "exact". It is included for
    # uniformity, not as independent evidence; the real guard on that row is the
    # length/description check inside 03.
    issues = 0
    for r in rows:
        pdb = r["job_name"].lower()
        try:
            ents = entities(pdb)
        except Exception as e:
            print(f"{pdb}: RCSB fetch failed ({e})"); issues += 1; continue
        print(f"\n=== {pdb} ({len(ents)} deposited entities) ===")
        for label, col in COLS:
            staged = (r[col] or "").strip().upper()
            verdict, eid, desc, note = classify(staged, ents)
            flag = "" if verdict == "exact" else "   <<<"
            if verdict in ("MISMATCH", "partial"):
                issues += 1
            elif verdict == "substring" and label != "MHC-b":
                issues += 1          # expected only for a tethered beta chain
                flag = "   <<< unexpected for this chain"
            print(f"  {label:8} len={len(staged):4}  {verdict:9} "
                  f"{('ent ' + eid) if eid else '':7} {desc[:40]:40} {note}{flag}")
    print(f"\n{'-'*70}\n{issues} item(s) needing a look. 'substring' on MHC-b is expected for the "
          f"tethered constructs (8vd0, 8vq8, 9aud).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
