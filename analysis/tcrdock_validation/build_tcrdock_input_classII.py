"""
Build a TCRdock-format input TSV from supplementaryTable3_updated.xlsx for the
1,433 class II validation triads (205 cognate + 1228 non-cognate from CRESTA).

Differences from the class I builder:
  - mhc column is "{alpha},{beta}" (e.g. "DRA*01:01,DRB1*15:03")
  - DRA1*01:02 -> DRA*01:01 normalization (TCRdock uses canonical IMGT, and DRA
    is functionally monomorphic so the substitution is safe)
  - mhc_class field is 2 (vs 1)

Run on Gemini:
  conda activate /scratch/bneff/tcrdock_validation/env
  pip install anarci          # if not already installed
  python build_tcrdock_input_classII.py
"""
import pandas as pd
from anarci import anarci

SRC          = 'supplementaryTable3_updated.xlsx'
OUT_TSV      = 'tcrdock_input_class_II.tsv'
OUT_FULL_XLSX= 'tcrdock_input_class_II_full.xlsx'
OUT_QC       = 'tcrdock_input_class_II_qc.txt'

# DRA chain naming: TCRdock uses standard IMGT 'DRA*01:01'; the CRESTA table
# has the non-standard 'DRA1*01:02'. DRA is near-monomorphic, so normalize.
DRA_MAP = {'DRA1*01:02': 'DRA*01:01',
           'DRA1*01:01': 'DRA*01:01'}

def normalize_mhc_alpha(name):
    return DRA_MAP.get(name, name)

def extract_vjcdr3(numbering_chain, alignment_chain):
    if numbering_chain is None or alignment_chain is None:
        return None
    numbered_seq = numbering_chain[0][0]
    ad = alignment_chain[0]
    return {
        'species':    ad['germlines']['v_gene'][0][0],
        'chain_type': ad['chain_type'],
        'v_gene':     ad['germlines']['v_gene'][0][1],
        'j_gene':     ad['germlines']['j_gene'][0][1],
        'cdr3':       ''.join(aa for ((pos, ins), aa) in numbered_seq
                              if 104 <= pos <= 118 and aa != '-'),
    }

def main():
    df = pd.read_excel(SRC)
    cii = df[df['mhc_class']=='II'].copy().reset_index(drop=True)
    print(f"Loaded {len(cii)} class II triads "
          f"({(cii['cognate']==True).sum()} cog, {(cii['cognate']==False).sum()} non-cog)")

    # All chains for one batched ANARCI call
    seqs = []
    for i, row in cii.iterrows():
        seqs.append((f'{i}_A', row['tcr_1_seq']))
        seqs.append((f'{i}_B', row['tcr_2_seq']))

    print(f"Running ANARCI on {len(seqs)} chains (~2-4 min for class II)...")
    numbering, alignment_details, _ = anarci(
        seqs, scheme='imgt', allow={'A','B','D'}, assign_germline=True,
    )

    parsed, fails = {}, []
    for k, (tag, _) in enumerate(seqs):
        result = extract_vjcdr3(numbering[k], alignment_details[k])
        if result is None:
            fails.append(tag)
        parsed[tag] = result
    print(f"  ANARCI: {len(parsed) - len(fails)}/{len(seqs)} chains parsed, "
          f"{len(fails)} failures")

    rows, qc = [], []
    for i, row in cii.iterrows():
        a = parsed.get(f'{i}_A')
        b = parsed.get(f'{i}_B')
        if a is None or b is None:
            qc.append(f"row {i}: ANARCI failure (a={a is not None}, b={b is not None})")
            continue
        if not a['v_gene'].startswith('TRAV'):
            qc.append(f"row {i}: alpha position has non-TRAV v_gene: {a['v_gene']}")
        if not b['v_gene'].startswith('TRBV'):
            qc.append(f"row {i}: beta position has non-TRBV v_gene: {b['v_gene']}")
        if not a['cdr3'].startswith('C') or a['cdr3'][-1] not in 'FW':
            qc.append(f"row {i}: unusual CDR3a: {a['cdr3']}")
        if not b['cdr3'].startswith('C') or b['cdr3'][-1] not in 'FW':
            qc.append(f"row {i}: unusual CDR3b: {b['cdr3']}")

        organism = a['species']
        # Class II MHC: '{alpha},{beta}' with DRA normalized
        mhc = f"{normalize_mhc_alpha(row['mhc_1_name'])},{row['mhc_2_name']}"

        rows.append({
            'pdbid':     row['job_name'],
            'organism':  organism,
            'mhc_class': 2,
            'mhc':       mhc,
            'peptide':   row['peptide'],
            'va':        a['v_gene'], 'ja': a['j_gene'], 'cdr3a': a['cdr3'],
            'vb':        b['v_gene'], 'jb': b['j_gene'], 'cdr3b': b['cdr3'],
        })

    out = pd.DataFrame(rows)

    # Same auto-fix as class I (TRDJ4*01 -> TRAJ29*01 for alpha/delta-shared V genes)
    KNOWN_FIXES = {('TRDJ4*01', 'human'): 'TRAJ29*01'}
    for (bad, org), good in KNOWN_FIXES.items():
        m = (out['ja'] == bad) & (out['organism'] == org)
        if m.any():
            print(f"  auto-fix: ja={bad} ({org}) -> {good} on {m.sum()} rows")
            out.loc[m, 'ja'] = good

    out.to_csv(OUT_TSV, sep='\t', index=False)
    print(f"Wrote: {OUT_TSV}  ({len(out)} rows)")

    # Bookkeeping xlsx with cognate/source columns retained
    meta = cii[['job_name','source','cognate']].rename(columns={'job_name':'pdbid'})
    full = out.merge(meta, on='pdbid', how='left')
    full.to_excel(OUT_FULL_XLSX, index=False)
    print(f"Wrote: {OUT_FULL_XLSX}")

    with open(OUT_QC, 'w') as f:
        if not qc:
            f.write("All class II triads parsed cleanly; no QC issues.\n")
        else:
            f.write(f"QC issues ({len(qc)}):\n\n" + '\n'.join(qc))
    print(f"Wrote: {OUT_QC} ({len(qc)} flagged items)")

    # Quick gene-validation against TCRdock's accepted gene db (if present)
    import os
    gene_db_path = 'TCRdock/tcrdock/tcrdist/db/combo_xcr_w_mouse_and_human.tsv'
    if os.path.exists(gene_db_path):
        gene_db = pd.read_csv(gene_db_path, sep='\t')
        ok = set(zip(gene_db['id'], gene_db['organism']))
        issues = [(c, r[c], r['organism']) for c in ['va','ja','vb','jb']
                  for _, r in out.iterrows() if (r[c], r['organism']) not in ok]
        print(f"\nGene validation: {len(issues)} of {4*len(out)} (gene, organism) pairs not in TCRdock db")
        if issues:
            unique = set(issues)
            for col, g, o in sorted(unique):
                n = sum(1 for x in issues if x == (col, g, o))
                print(f"  {col}={g} ({o}): {n} rows")

    # Per-antigen sanity
    out_ = out.assign(pmhc=out['mhc']+':'+out['peptide']).merge(meta, on='pdbid')
    print(f"\nAntigen breakdown ({out_['pmhc'].nunique()} unique antigens):")
    for pmhc, g in out_.groupby('pmhc'):
        print(f"  {pmhc}: {(g['cognate']==True).sum()} cog, {(g['cognate']==False).sum()} non-cog")

if __name__ == '__main__':
    main()
