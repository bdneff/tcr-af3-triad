#!/usr/bin/env python3
"""03_derive_27eb.py — build the missing 27EB row from RCSB, in John's TSV schema.

John's write-up describes three HLA-DQ2.5 gliadin structures ("released 2026-08-19 and
2026-08-26") but `data/hla2_9_triads.tsv` contains only the two from 08-19. An RCSB search on his
stated third date returns exactly one candidate -- 27EB, "TCR-HLA-DQ2.5-glia-w2" -- so the row was
dropped from his export rather than never selected.

Rather than wait on a round trip, the sequences are taken straight from the deposited polymer
entities. That is defensible because it is demonstrably the same rule John used: 02_xcheck_sequences
found 42 of his 45 chains to be EXACT matches to deposited entities, the three exceptions being the
tethered constructs where he correctly split peptide from beta chain. 27EB has five separate
entities, so no splitting is needed.

Written to its OWN file, not appended to John's: `data/hla2_9_triads.tsv` is his artifact and
should keep saying what he sent. This row is ours, derived, and still pending his confirmation.

    python analysis/classII_supp3c/scripts/03_derive_27eb.py

Read-only network access to RCSB.
"""
import csv, json, os, sys, time, urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
OUT = os.path.join(ROOT, "data", "hla2_derived_extra.tsv")
PDB = "27EB"

# entity description -> role. Matched on substrings because depositors name chains inconsistently
# (27EB's beta chain is filed as "cDNA FLJ59004, highly similar to HLA class II histocompatibility
# antigen"), so an exact-name lookup would fail on exactly the entry we are trying to recover.
# Order matters, and so does specificity. A first attempt put the TCR roles first with a bare
# "beta chain" hint; 27EB's MHC beta is deposited as "cDNA FLJ59004, highly similar to HLA class II
# histocompatibility antigen, DQ beta chain", so "beta chain" matched the MHC and the real TCR beta
# was left over -- a silent swap that the length check would NOT have caught, since both chains are
# ~200 aa. MHC roles are therefore claimed first, and every TCR hint now demands an explicit
# receptor token.
ROLE_HINTS = [
    ("peptide", ("peptide",)),
    ("mhc_a",   ("dq alpha", "dr alpha", "dq-alpha", "dr-alpha", "class ii histocompatibility antigen, dq alpha")),
    ("mhc_b",   ("dq beta", "dr beta", "dq-beta", "dr-beta", "flj59004", "class ii histocompatibility antigen, dq beta")),
    ("tcr_a",   ("tcr alpha", "tcr-alpha", "t cell receptor alpha", "t-cell receptor alpha", "t-cell-receptor, tcr alpha")),
    ("tcr_b",   ("tcr beta", "tcr-beta", "t cell receptor beta", "t-cell receptor beta", "t-cell-receptor, tcr beta")),
]

FIELDS = ["job_name", "source", "cognate", "peptide", "mhc_class",
          "mhc_1_chain", "mhc_1_species", "mhc_1_name", "mhc_1_seq",
          "mhc_2_chain", "mhc_2_species", "mhc_2_name", "mhc_2_seq",
          "tcr_1_chain", "tcr_1_species", "tcr_1_seq",
          "tcr_2_chain", "tcr_2_species", "tcr_2_seq",
          "mean_p_tcr_interface_pae"]


def fetch(pdb):
    with urllib.request.urlopen(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}", timeout=30) as f:
        e = json.load(f)
    ents = []
    for eid in e["rcsb_entry_container_identifiers"]["polymer_entity_ids"]:
        u = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{eid}"
        with urllib.request.urlopen(u, timeout=30) as f:
            pe = json.load(f)
        ents.append(dict(
            eid=eid,
            seq=pe["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "").upper(),
            desc=((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", "")),
            species=((pe.get("rcsb_entity_source_organism") or [{}])[0]
                     .get("ncbi_scientific_name", "")),
        ))
        time.sleep(0.05)
    return e, ents


def assign(ents):
    """Map entities to roles, then sanity-check by length -- names alone are not trustworthy."""
    roles = {}
    for role, hints in ROLE_HINTS:
        for en in ents:
            if en["eid"] in [v["eid"] for v in roles.values()]:
                continue
            d = en["desc"].lower()
            if any(h in d for h in hints):
                roles[role] = en
                break
    # length is the independent check: a class II peptide is short, the chains are not
    expect = {"peptide": (8, 25), "mhc_a": (150, 220), "mhc_b": (150, 230),
              "tcr_a": (180, 240), "tcr_b": (200, 270)}
    problems = [f"{r}: not assigned" for r in expect if r not in roles]
    for r, en in roles.items():
        lo, hi = expect[r]
        if not (lo <= len(en["seq"]) <= hi):
            problems.append(f"{r}: length {len(en['seq'])} outside expected {lo}-{hi} "
                            f"(entity {en['eid']}, '{en['desc'][:40]}')")
    return roles, problems


def main():
    entry, ents = fetch(PDB)
    roles, problems = assign(ents)
    if problems:
        print(f"REFUSING to write — role assignment for {PDB} is not unambiguous:")
        for p in problems:
            print("   " + p)
        print("\nentities seen:")
        for en in ents:
            print(f"   ent {en['eid']} len={len(en['seq']):4} {en['desc'][:60]}")
        sys.exit(1)

    sp = lambda e: {"Homo sapiens": "human", "Mus musculus": "mouse"}.get(e["species"],
                                                                          e["species"])
    row = {
        "job_name": PDB.lower(), "source": "pdb", "cognate": "TRUE",
        "peptide": roles["peptide"]["seq"], "mhc_class": "II",
        "mhc_1_chain": "alpha", "mhc_1_species": sp(roles["mhc_a"]),
        "mhc_1_name": "DQA1*05:01", "mhc_1_seq": roles["mhc_a"]["seq"],
        "mhc_2_chain": "beta", "mhc_2_species": sp(roles["mhc_b"]),
        "mhc_2_name": "DQB1*02:01", "mhc_2_seq": roles["mhc_b"]["seq"],
        "tcr_1_chain": "alpha", "tcr_1_species": sp(roles["tcr_a"]),
        "tcr_1_seq": roles["tcr_a"]["seq"],
        "tcr_2_chain": "beta", "tcr_2_species": sp(roles["tcr_b"]),
        "tcr_2_seq": roles["tcr_b"]["seq"],
        "mean_p_tcr_interface_pae": "",
    }
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader(); w.writerow(row)

    rel = entry["rcsb_accession_info"]["initial_release_date"][:10]
    print(f"{PDB}  {entry['struct']['title']}")
    print(f"  released {rel} — post-AF3-cutoff: {rel >= '2023-01-12'}")
    for r in ("mhc_a", "mhc_b", "peptide", "tcr_a", "tcr_b"):
        print(f"  {r:8} ent {roles[r]['eid']}  len {len(roles[r]['seq']):4}  "
              f"{roles[r]['desc'][:48]}")
    print(f"\n-> {OUT}")
    print("NOTE: alleles are labelled DQA1*05:01 / DQB1*02:01 by series (DQ2.5, matching 9yaf/9yag);"
          "\n      confirm with John along with whether 27EB belongs at all.")


if __name__ == "__main__":
    main()
