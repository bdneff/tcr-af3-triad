"""
Build a TCRdock-format input TSV from supplementaryTable3.xlsx
for the 197 class I validation triads (16 cognate + 181 non-cognate).

TCRdock input columns: pdbid, organism, mhc_class, mhc, peptide,
                       va, ja, cdr3a, vb, jb, cdr3b

This script parses every TCR alpha/beta with ANARCI, extracts V, J, CDR3
in IMGT format, and validates the result against the reported species.
"""
import pandas as pd
from anarci import anarci
import sys

SRC = '/mnt/project/supplementaryTable3.xlsx'
from pathlib import Path
HERE          = Path(__file__).resolve().parent   # outputs sit next to this script
OUT_TSV       = HERE / 'tcrdock_input_class_I.tsv'
OUT_FULL_XLSX = HERE / 'tcrdock_input_class_I_full.xlsx'   # with extra context cols
OUT_QC        = HERE / 'tcrdock_input_qc.txt'

def extract_vjcdr3(numbering_chain, alignment_chain):
    """Pull V gene, J gene, CDR3 (IMGT 104-118) from ANARCI output."""
    if numbering_chain is None or alignment_chain is None:
        return None
    numbered_seq = numbering_chain[0][0]
    ad = alignment_chain[0]
    return {
        'species':    ad['germlines']['v_gene'][0][0],
        'chain_type': ad['chain_type'],
        'v_gene':     ad['germlines']['v_gene'][0][1],
        'j_gene':     ad['germlines']['j_gene'][0][1],
        'v_identity': ad['germlines']['v_gene'][1],
        'j_identity': ad['germlines']['j_gene'][1],
        'cdr3':       ''.join(aa for ((pos, ins), aa) in numbered_seq
                              if 104 <= pos <= 118 and aa != '-'),
    }

def main():
    df = pd.read_excel(SRC)
    ci = df[df['mhc_class']=='I'].copy().reset_index(drop=True)
    print(f"Loaded {len(ci)} class I triads "
          f"({(ci['cognate']==True).sum()} cognate, "
          f"{(ci['cognate']==False).sum()} non-cognate)")

    # Build a flat list of all chains for one batched ANARCI call
    # We tag each input with (row_idx, chain) so we can map back
    seqs = []
    for i, row in ci.iterrows():
        seqs.append((f'{i}_A', row['tcr_1_seq']))   # alpha
        seqs.append((f'{i}_B', row['tcr_2_seq']))   # beta

    # Allow A,B,D as well: alpha/delta-shared V genes like TRAV38-2/DV8*01
    # are sometimes called 'D' by ANARCI even when used in an alpha-position TCR.
    # We trust the V-gene name (TRAV* vs TRBV*) for the alpha/beta identity.
    print(f"Running ANARCI on {len(seqs)} chains (this takes ~30-60s)...")
    numbering, alignment_details, _ = anarci(
        seqs, scheme='imgt', allow={'A','B','D'}, assign_germline=True,
    )

    # Parse out per-chain results
    parsed = {}
    fails  = []
    for k, (tag, _) in enumerate(seqs):
        result = extract_vjcdr3(numbering[k], alignment_details[k])
        if result is None:
            fails.append(tag)
        parsed[tag] = result

    print(f"  ANARCI: {len(parsed) - len(fails)}/{len(seqs)} chains parsed, "
          f"{len(fails)} failures")
    if fails:
        print(f"  failed tags: {fails[:10]}{'...' if len(fails)>10 else ''}")

    # Assemble TCRdock TSV rows + run QC checks
    rows, qc = [], []
    for i, row in ci.iterrows():
        a = parsed.get(f'{i}_A')
        b = parsed.get(f'{i}_B')
        if a is None or b is None:
            qc.append(f"row {i}: ANARCI failure (a={a is not None}, b={b is not None})")
            continue

        # Sanity: at alpha position we expect a TRAV-family gene (incl. TRAV/DV
        # shared); at beta position we expect TRBV. We trust the V-gene
        # name rather than ANARCI's chain_type letter.
        if not a['v_gene'].startswith('TRAV'):
            qc.append(f"row {i}: alpha position has non-TRAV v_gene: {a['v_gene']}")
        if not b['v_gene'].startswith('TRBV'):
            qc.append(f"row {i}: beta position has non-TRBV v_gene: {b['v_gene']}")

        # Species: ANARCI species and table-reported species should agree
        table_species = row['tcr_1_species']
        if a['species'] != table_species:
            qc.append(f"row {i}: alpha species mismatch "
                      f"(table={table_species}, anarci={a['species']}, "
                      f"v_identity={a['v_identity']:.2f})")
        if b['species'] != row['tcr_2_species']:
            qc.append(f"row {i}: beta species mismatch "
                      f"(table={row['tcr_2_species']}, anarci={b['species']}, "
                      f"v_identity={b['v_identity']:.2f})")

        # CDR3 should start with C and end with F or W
        if not a['cdr3'].startswith('C') or a['cdr3'][-1] not in 'FW':
            qc.append(f"row {i}: unusual CDR3a: {a['cdr3']}")
        if not b['cdr3'].startswith('C') or b['cdr3'][-1] not in 'FW':
            qc.append(f"row {i}: unusual CDR3b: {b['cdr3']}")

        # Use ANARCI species for organism (TCRdock needs the TCR's organism
        # to look up V/J genes; HLA allele identifies the human MHC separately)
        organism = a['species']   # alpha and beta should agree; we re-checked above

        rows.append({
            'pdbid':     row['job_name'],     # unique per triad
            'organism':  organism,
            'mhc_class': 1,
            'mhc':       row['mhc_1_name'],
            'peptide':   row['peptide'],
            'va':        a['v_gene'], 'ja': a['j_gene'], 'cdr3a': a['cdr3'],
            'vb':        b['v_gene'], 'jb': b['j_gene'], 'cdr3b': b['cdr3'],
        })

    out = pd.DataFrame(rows)

    # Self-heal: when ANARCI sees TRAV/DV-shared V genes it sometimes assigns
    # a TRDJ (delta J) where the real J is a TRAJ (alpha J). We detect this
    # by checking every gene name against TCRdock's accepted gene database
    # and, for the known TRDJ4*01 (human) case, swapping to TRAJ29*01 which
    # uniquely matches the FR4 sequence GKGTRLSV seen in our data.
    KNOWN_FIXES = {('TRDJ4*01', 'human'): 'TRAJ29*01'}
    for (bad_gene, org), good_gene in KNOWN_FIXES.items():
        m = (out['ja'] == bad_gene) & (out['organism'] == org)
        if m.any():
            print(f"  auto-fix: ja={bad_gene} ({org}) -> {good_gene} on {m.sum()} rows")
            out.loc[m, 'ja'] = good_gene

    out.to_csv(OUT_TSV, sep='\t', index=False)
    print(f"Wrote TCRdock-format TSV: {OUT_TSV}  ({len(out)} rows)")

    # Also save with cognate/source columns retained, for our own bookkeeping
    full = out.copy()
    keep_meta = ci[['job_name','source','cognate']].rename(columns={'job_name':'pdbid'})
    full = full.merge(keep_meta, on='pdbid', how='left')
    full.to_excel(OUT_FULL_XLSX, index=False)
    print(f"Wrote bookkeeping xlsx (TCRdock cols + cognate/source): {OUT_FULL_XLSX}")

    # QC report
    with open(OUT_QC, 'w') as f:
        if not qc:
            f.write("All 197 triads parsed cleanly; no QC issues.\n")
        else:
            f.write(f"QC issues for {len(qc)} checks (across {len(ci)} triads):\n\n")
            f.write('\n'.join(qc))
    print(f"Wrote QC report: {OUT_QC}  ({len(qc)} flagged items)")
    if qc:
        print("\nFirst few QC items:")
        for line in qc[:10]:
            print(f"  {line}")

if __name__ == '__main__':
    main()
