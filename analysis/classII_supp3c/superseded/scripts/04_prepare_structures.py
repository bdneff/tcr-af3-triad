#!/usr/bin/env python3
"""04_prepare_structures.py — fetch each deposited structure and write it in canonical 5-chain form.

CDR RMSD compares predicted TCR loops against the deposited ones after superposing on the MHC, so
it needs the two structures to agree about which chain is which. For class I that is free: every
entry is five separate chains. For class II it is not, and this is the one genuinely new piece of
work in this analysis.

Three things have to be fixed per entry:

**Tethered peptides.** `8vd0`, `8vq8` and `9aud` fuse the peptide to the MHC beta chain through a
short linker, so the deposit has FOUR polymer entities where the prediction has five. Left alone,
any chain-by-position pairing silently mismatches -- and would still return a plausible number,
which is the worst kind of failure. The fused entity is cut at the offsets measured against the
staged sequences (`02_xcheck_sequences.py`): peptide occupies residues [0:len(peptide)], the linker
runs to the offset where the beta chain begins, and the linker is DISCARDED because the prediction
has no linker to compare it with.

**Multiple copies.** Several entries carry two complexes in the asymmetric unit (chains A,F / C,H
/ ...). One is chosen -- the copy whose chains are all present -- since the prediction models one.

**Chain naming.** Output is relabelled A = MHC alpha, B = MHC beta, C = peptide, D = TCR alpha,
E = TCR beta, matching the AF3 JSONs staged by `01_stage_af3.py` and the `segid D`/`segid E`
selections in `analysis/triad_analysis/rmsd.py`.

**Chains are assigned from RCSB's entity records, not by matching coordinates.** A first attempt
matched the staged sequence against the sequence observed in the ATOM records and failed on every
long chain: crystals omit disordered residues, so the observed sequence has gaps the entity
sequence does not. Peptides matched (they are fully ordered) and nothing else did. The deposit's
own entity -> auth_asym_id mapping is authoritative and carries no such problem. Coverage is then
reported as an identity fraction rather than demanded exact -- gaps are expected, and
`rmsd.py` aligns sequences itself (`seq_align_residue_groups`) before comparing.

    python analysis/classII_supp3c/scripts/04_prepare_structures.py

Writes analysis/classII_supp3c/structures/<pdb>.pdb (+ a chain_map.csv recording every decision).
Needs biopython. Read-only network access to RCSB.
"""
import csv, gzip, io, json, os, sys, time, urllib.request, warnings

from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

warnings.filterwarnings("ignore")

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
BASE = os.path.join(ROOT, "analysis", "classII_supp3c")
TSVS = [os.path.join(ROOT, "data", "hla2_9_triads.tsv"),
        os.path.join(ROOT, "data", "hla2_derived_extra.tsv")]
OUTD = os.path.join(BASE, "structures")
CIFD = os.path.join(BASE, "structures", "_cif")
MAP = os.path.join(BASE, "data", "chain_map.csv")

ROLE_ORDER = [("A", "mhc_1_seq"), ("B", "mhc_2_seq"), ("C", "peptide"),
              ("D", "tcr_1_seq"), ("E", "tcr_2_seq")]
THREE2ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}


def seq_of(chain):
    """one-letter sequence and the residue objects behind it, standard amino acids only."""
    res, s = [], []
    for r in chain:
        if not is_aa(r, standard=True):
            continue
        one = THREE2ONE.get(r.get_resname().upper())
        if one:
            res.append(r); s.append(one)
    return "".join(s), res


def fetch_cif(pdb):
    os.makedirs(CIFD, exist_ok=True)
    path = os.path.join(CIFD, f"{pdb}.cif")
    if os.path.exists(path) and os.path.getsize(path):
        return path
    url = f"https://files.rcsb.org/download/{pdb.upper()}.cif.gz"
    with urllib.request.urlopen(url, timeout=60) as r:
        data = gzip.decompress(r.read())
    with open(path, "wb") as f:
        f.write(data)
    time.sleep(0.2)
    return path


class Keep(Select):
    """write only the chosen residues, and only their standard-amino-acid atoms."""
    def __init__(self, keep):
        self.keep = {id(r) for r in keep}
    def accept_residue(self, residue):
        return id(residue) in self.keep
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() in (" ", "A")


def entity_map(pdb):
    """{entity_id: (description, canonical sequence, [auth chain ids])} from RCSB."""
    u = f"https://data.rcsb.org/rest/v1/core/entry/{pdb.upper()}"
    with urllib.request.urlopen(u, timeout=30) as f:
        e = json.load(f)
    out = {}
    for eid in e["rcsb_entry_container_identifiers"]["polymer_entity_ids"]:
        u2 = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb.upper()}/{eid}"
        with urllib.request.urlopen(u2, timeout=30) as f:
            pe = json.load(f)
        out[eid] = (
            (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", ""),
            pe["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "").upper(),
            pe["rcsb_polymer_entity_container_identifiers"].get("auth_asym_ids", []),
        )
        time.sleep(0.05)
    return out


def identity(a, b):
    """fraction of the shorter sequence recovered, gap-tolerant (LCS-free approximation)."""
    if not a or not b:
        return 0.0
    import difflib
    m = difflib.SequenceMatcher(None, a, b).find_longest_match(0, len(a), 0, len(b))
    return m.size / min(len(a), len(b))


def main():
    rows = []
    for t in TSVS:
        if os.path.exists(t):
            rows += list(csv.DictReader(open(t), delimiter="\t"))

    os.makedirs(OUTD, exist_ok=True)
    os.makedirs(os.path.dirname(MAP), exist_ok=True)
    parser, io_ = MMCIFParser(QUIET=True), PDBIO()
    records, failures = [], []

    for r in rows:
        pdb = r["job_name"].strip().lower()
        try:
            cif = fetch_cif(pdb)
            model = next(iter(parser.get_structure(pdb, cif)))
            ents = entity_map(pdb)
        except Exception as e:
            failures.append(f"{pdb}: could not read structure/entities ({e})"); continue

        # role -> entity, by matching the staged sequence against ENTITY sequences (the check
        # 02_xcheck already validated: 42/45 exact, the rest substrings of tethered entities)
        role_ent, notes = {}, []
        for out_id, col in ROLE_ORDER:
            staged = (r[col] or "").strip().upper()
            hit = None
            for eid, (desc, eseq, chs) in ents.items():
                if staged == eseq:
                    hit = (eid, 0, "exact"); break
            if hit is None:
                for eid, (desc, eseq, chs) in ents.items():
                    if staged and staged in eseq:
                        hit = (eid, eseq.index(staged), "in-entity"); break
            if hit is None:
                notes.append(f"{out_id}: NO ENTITY"); continue
            role_ent[out_id] = hit
            notes.append(f"{out_id}<-ent{hit[0]}{'*' if hit[2] != 'exact' else ''}")

        missing = [o for o, _ in ROLE_ORDER if o not in role_ent]
        if missing:
            failures.append(f"{pdb}: no entity for {missing} — " + " ".join(notes)); continue

        # pick ONE biological copy: the auth chain set that covers every role
        by_chain = {ch.id: seq_of(ch) for ch in model}
        copies = []
        for eid, (_, _, chs) in ents.items():
            copies.append(chs)
        n_copies = min(len(c) for c in copies) if copies else 1
        chosen, ok = {}, True
        for ci in range(n_copies):
            trial = {}
            for out_id, (eid, off, mode) in role_ent.items():
                chs = ents[eid][2]
                if ci >= len(chs) or chs[ci] not in by_chain:
                    trial = None; break
                trial[out_id] = chs[ci]
            if trial:
                chosen = trial; break
        if not chosen:
            failures.append(f"{pdb}: no complete copy in the asymmetric unit — " + " ".join(notes))
            continue

        # cut residues per role. A tethered entity backs BOTH the peptide and the MHC beta, so the
        # two are sliced out of one chain by locating each staged sequence in the OBSERVED residues.
        picked, cov = {}, {}
        for out_id, col in ROLE_ORDER:
            staged = (r[col] or "").strip().upper()
            cid = chosen[out_id]
            oseq, ores = by_chain[cid]
            eid, off, mode = role_ent[out_id]
            if mode == "exact" and len([o for o in chosen.values() if o == cid]) == 1:
                res = ores                      # whole chain is this role
            else:
                # tethered: find where this role's sequence sits inside the observed chain
                import difflib
                m = difflib.SequenceMatcher(None, staged, oseq).find_longest_match(
                    0, len(staged), 0, len(oseq))
                if m.size < 8:
                    failures.append(f"{pdb}: cannot locate {out_id} inside chain {cid}")
                    res = None
                else:
                    start = m.b - m.a
                    res = ores[max(0, start):max(0, start) + len(staged)]
            if not res:
                picked = {}; break
            picked[out_id] = res
            got = "".join(THREE2ONE.get(x.get_resname().upper(), "X") for x in res)
            cov[out_id] = identity(staged, got)

        if not picked or len(picked) != 5:
            continue

        weak = [f"{k} {v:.2f}" for k, v in cov.items() if v < 0.90]
        if weak:
            failures.append(f"{pdb}: low sequence coverage — " + ", ".join(weak)); continue

        for out_id, _ in ROLE_ORDER:
            for res in picked[out_id]:
                res.parent.id = f"_{out_id}"
        for ch in list(model):
            if ch.id.startswith("_"):
                ch.id = ch.id[1:]

        keep = [x for out_id, _ in ROLE_ORDER for x in picked[out_id]]
        io_.set_structure(model)
        io_.save(os.path.join(OUTD, f"{pdb}.pdb"), select=Keep(keep))
        records.append(dict(pdb_id=pdb,
                            copy=" ".join(f"{k}={v}" for k, v in chosen.items()),
                            tethered=str(any("*" in n for n in notes)),
                            coverage=" ".join(f"{k}:{v:.2f}" for k, v in cov.items()),
                            entities=" ".join(notes)))
        print(f"{pdb}: " + " ".join(notes) + "  | copy " +
              ",".join(chosen[o] for o, _ in ROLE_ORDER) +
              "  | cov " + " ".join(f"{cov[o]:.2f}" for o, _ in ROLE_ORDER))

    if records:
        with open(MAP, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(records[0]))
            w.writeheader(); w.writerows(records)
    print(f"\nwrote {len(records)}/{len(rows)} structures -> {OUTD}")
    print(f"chain map -> {MAP}")
    if failures:
        print("\nFAILED:")
        for x in failures:
            print("   " + x)
        sys.exit(1)


if __name__ == "__main__":
    main()
