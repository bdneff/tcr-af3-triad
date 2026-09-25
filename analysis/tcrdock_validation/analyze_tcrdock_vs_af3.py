"""
Post-process TCRdock predictions to per-antigen AUCs for direct comparison
with AF3 PTI-PAE in Figure 4a.

Inputs:
  --tcrdock_tsv  TSV output from add_pmhc_tcr_pae_to_tsvfile.py
                 (must contain at minimum: pdbid, pmhc_tcr_pae)
  --meta_xlsx    bookkeeping xlsx written by build_tcrdock_input.py
                 (pdbid, mhc, peptide, cognate, ...) - default uses
                 tcrdock_input_class_I_full.xlsx in cwd

Outputs:
  --out_per_antigen   CSV of per-antigen AUC for TCRdock and AF3
  --out_long          CSV in long form for plotting (one row per triad with
                      both scores, cognate label, antigen)
  Console: median, range, paired comparison

Scoring convention
------------------
  TCRdock pmhc_tcr_pae: lower = better cognate, so the discriminative
                        direction is -pmhc_tcr_pae (or equivalently we
                        pass roc_auc_score with the negative score).
  AF3 PTI-PAE: in the manuscript this is 31 - mean_p_tcr_interface_pae
               (higher = better). The supplementary table provides
               mean_p_tcr_interface_pae raw (lower = better), so we
               flip its sign here too. Both methods then use the same
               "higher score = more likely cognate" convention.
"""
import argparse
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from pathlib import Path

DEFAULT_META = str(Path(__file__).resolve().parent / 'tcrdock_input_class_I_full.xlsx')
DEFAULT_VAL  = '/mnt/project/supplementaryTable3.xlsx'

def per_antigen_aucs(df, score_col, score_higher_is_cognate):
    """Return dict {pmhc: AUC}. df must have 'pmhc' and 'cognate' columns.
    Each antigen needs >=1 cognate and >=1 non-cognate to score."""
    out = {}
    for pmhc, g in df.groupby('pmhc'):
        if g['cognate'].nunique() < 2:
            continue
        y_true = g['cognate'].astype(int).values     # 1 = cognate, 0 = not
        s = g[score_col].astype(float).values
        if not score_higher_is_cognate:
            s = -s
        out[pmhc] = roc_auc_score(y_true, s)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tcrdock_tsv', required=True,
                    help='Output of add_pmhc_tcr_pae_to_tsvfile.py')
    ap.add_argument('--meta_xlsx', default=DEFAULT_META,
                    help='bookkeeping file with pdbid -> mhc/peptide/cognate')
    ap.add_argument('--af3_val_xlsx', default=DEFAULT_VAL,
                    help='supplementaryTable3.xlsx (for AF3 PAE values)')
    ap.add_argument('--out_per_antigen', default='per_antigen_auc_TCRdock_vs_AF3.csv')
    ap.add_argument('--out_long', default='per_triad_scores_long.csv')
    args = ap.parse_args()

    # Load TCRdock results (one row per triad, with pmhc_tcr_pae column)
    td = pd.read_csv(args.tcrdock_tsv, sep='\t')
    if 'pmhc_tcr_pae' not in td.columns:
        raise SystemExit(f'pmhc_tcr_pae column not found in {args.tcrdock_tsv}. '
                         f'Did you run add_pmhc_tcr_pae_to_tsvfile.py?')
    td = td[['pdbid','pmhc_tcr_pae']].copy()

    # Merge with our bookkeeping (cognate label, antigen)
    meta = pd.read_excel(args.meta_xlsx)[['pdbid','mhc','peptide','cognate']]
    merged = meta.merge(td, on='pdbid', how='left')
    n_missing = merged['pmhc_tcr_pae'].isna().sum()
    if n_missing:
        print(f"WARNING: {n_missing}/{len(merged)} triads have no TCRdock score "
              f"(skipped during scoring)")

    # Merge with AF3 PAE from the manuscript's supplementary table
    af3 = pd.read_excel(args.af3_val_xlsx)
    af3 = af3[af3['mhc_class']=='I'][['job_name','mean_p_tcr_interface_pae']].rename(
        columns={'job_name':'pdbid'})
    merged = merged.merge(af3, on='pdbid', how='left')

    merged['pmhc'] = merged['mhc'] + ':' + merged['peptide']
    merged.dropna(subset=['pmhc_tcr_pae','mean_p_tcr_interface_pae'], inplace=True)

    # Per-antigen AUCs (both metrics: lower raw value = better cognate, so flip)
    td_aucs  = per_antigen_aucs(merged, 'pmhc_tcr_pae',           score_higher_is_cognate=False)
    af3_aucs = per_antigen_aucs(merged, 'mean_p_tcr_interface_pae', score_higher_is_cognate=False)

    pa = pd.DataFrame({
        'pmhc':         list(td_aucs.keys()),
        'AUC_TCRdock':  [td_aucs[k]  for k in td_aucs],
        'AUC_AF3_PTI_PAE':[af3_aucs.get(k, np.nan) for k in td_aucs],
    }).sort_values('AUC_AF3_PTI_PAE', ascending=False)
    pa.to_csv(args.out_per_antigen, index=False)

    print(f"\n=== Per-antigen AUC (n={len(pa)} antigens) ===")
    print(pa.to_string(index=False))
    print(f"\nMedian AUC TCRdock pmhc_tcr_pae : {pa['AUC_TCRdock'].median():.3f} "
          f"(range {pa['AUC_TCRdock'].min():.2f}-{pa['AUC_TCRdock'].max():.2f})")
    print(f"Median AUC AF3 PTI-PAE          : {pa['AUC_AF3_PTI_PAE'].median():.3f} "
          f"(range {pa['AUC_AF3_PTI_PAE'].min():.2f}-{pa['AUC_AF3_PTI_PAE'].max():.2f})")

    # Wilcoxon signed-rank paired test (same antigens, both methods)
    from scipy.stats import wilcoxon
    valid = pa.dropna(subset=['AUC_TCRdock','AUC_AF3_PTI_PAE'])
    if len(valid) >= 5:
        w = wilcoxon(valid['AUC_AF3_PTI_PAE'], valid['AUC_TCRdock'])
        print(f"\nPaired Wilcoxon (AF3 vs TCRdock, per-antigen AUCs): "
              f"statistic={w.statistic:.1f}, p={w.pvalue:.3g}")

    merged.to_csv(args.out_long, index=False)
    print(f"\nWrote: {args.out_per_antigen}, {args.out_long}")

if __name__ == '__main__':
    main()
