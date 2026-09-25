"""
TCRdock vs AF3 PTI-PAE comparison, with BOTH orientations:
  - antigen-centric: per-antigen AUC, AF3 vs TCRdock
  - TCR-centric:     per-cognate-TCR AUC (rank of cognate score among
                     non-cognates sharing the same antigen)

Works for class I and class II. Pass --mhc_class to control which subset of
supplementaryTable3 to pull AF3 PAE values from.

Usage:
  python analyze_tcrdock_vs_af3_v2.py \
      --tcrdock_tsv   user_output_w_pae.tsv \
      --meta_xlsx     tcrdock_input_class_II_full.xlsx \
      --af3_val_xlsx  supplementaryTable3_updated.xlsx \
      --mhc_class     II \
      --prefix        classII
"""
import argparse
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon

def per_antigen_aucs(df, score_col, score_higher_is_cognate):
    out = {}
    for pmhc, g in df.groupby('pmhc'):
        if g['cognate'].nunique() < 2:
            continue
        y = g['cognate'].astype(int).values
        s = g[score_col].astype(float).values
        if not score_higher_is_cognate: s = -s
        out[pmhc] = roc_auc_score(y, s)
    return out

def per_tcr_aucs(df, score_col, score_higher_is_cognate):
    """For each cognate triad, AUC = rank of its score vs the non-cognates
    sharing the same antigen. Returns dict {pdbid_of_cognate: AUC}."""
    out = {}
    for pmhc, g in df.groupby('pmhc'):
        cogs = g[g['cognate'] == True]
        noncogs = g[g['cognate'] == False]
        if len(cogs) == 0 or len(noncogs) == 0:
            continue
        n_scores = noncogs[score_col].astype(float).values
        for _, c in cogs.iterrows():
            c_s = float(c[score_col])
            if score_higher_is_cognate:
                # AUC = fraction of non-cognates the cognate beats
                wins = (c_s > n_scores).sum() + 0.5 * (c_s == n_scores).sum()
            else:
                wins = (c_s < n_scores).sum() + 0.5 * (c_s == n_scores).sum()
            out[c['pdbid']] = wins / len(n_scores)
    return out

def summarize(name, vals):
    arr = np.array(list(vals.values()))
    print(f"  {name:<30s} n={len(arr):4d}  "
          f"median={np.median(arr):.3f}  "
          f"range={arr.min():.2f}-{arr.max():.2f}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tcrdock_tsv',   required=True)
    ap.add_argument('--meta_xlsx',     required=True)
    ap.add_argument('--af3_val_xlsx',  required=True)
    ap.add_argument('--mhc_class',     default='I', choices=['I','II'])
    ap.add_argument('--prefix',        default='comparison')
    args = ap.parse_args()

    # 1. TCRdock output -> pdbid + pmhc_tcr_pae
    td = pd.read_csv(args.tcrdock_tsv, sep='\t')
    if 'pmhc_tcr_pae' not in td.columns:
        raise SystemExit(f"pmhc_tcr_pae column missing in {args.tcrdock_tsv}")
    td = td[['pdbid','pmhc_tcr_pae']]

    # 2. Bookkeeping (cognate, mhc, peptide)
    meta = pd.read_excel(args.meta_xlsx)[['pdbid','mhc','peptide','cognate']]

    # 3. AF3 PAE values from sup table 3 (filtered to the right class)
    af3 = pd.read_excel(args.af3_val_xlsx)
    af3 = af3[af3['mhc_class'] == args.mhc_class][['job_name','mean_p_tcr_interface_pae']]
    af3 = af3.rename(columns={'job_name':'pdbid'})

    merged = meta.merge(td, on='pdbid', how='left') \
                 .merge(af3, on='pdbid', how='left')
    n_missing_td  = merged['pmhc_tcr_pae'].isna().sum()
    n_missing_af3 = merged['mean_p_tcr_interface_pae'].isna().sum()
    if n_missing_td:  print(f"WARNING: {n_missing_td} triads missing TCRdock score (excluded)")
    if n_missing_af3: print(f"WARNING: {n_missing_af3} triads missing AF3 score   (excluded)")
    merged['pmhc'] = merged['mhc'] + ':' + merged['peptide']
    merged = merged.dropna(subset=['pmhc_tcr_pae','mean_p_tcr_interface_pae'])
    print(f"Comparing on {len(merged)} triads across {merged['pmhc'].nunique()} antigens "
          f"({(merged['cognate']==True).sum()} cog, {(merged['cognate']==False).sum()} non-cog)")

    # --- Antigen-centric AUCs
    td_ag  = per_antigen_aucs(merged, 'pmhc_tcr_pae',           False)
    af3_ag = per_antigen_aucs(merged, 'mean_p_tcr_interface_pae', False)
    # --- TCR-centric AUCs
    td_tcr = per_tcr_aucs(merged, 'pmhc_tcr_pae',           False)
    af3_tcr= per_tcr_aucs(merged, 'mean_p_tcr_interface_pae', False)

    print("\n=== Summary ===")
    print("antigen-centric:")
    summarize("AF3 PTI-PAE",      af3_ag)
    summarize("TCRdock pmhc_pae", td_ag)
    print("TCR-centric:")
    summarize("AF3 PTI-PAE",      af3_tcr)
    summarize("TCRdock pmhc_pae", td_tcr)

    # Paired Wilcoxon for both views
    def paired(d_af3, d_td, label):
        ks = sorted(set(d_af3) & set(d_td))
        if len(ks) < 6:
            print(f"  paired Wilcoxon ({label}): n={len(ks)}, too few for test")
            return
        a = np.array([d_af3[k] for k in ks])
        t = np.array([d_td[k]  for k in ks])
        try:
            w = wilcoxon(a, t)
            print(f"  paired Wilcoxon ({label}, n={len(ks)}): stat={w.statistic:.1f}, p={w.pvalue:.3g}")
        except Exception as e:
            print(f"  paired Wilcoxon ({label}): {e}")
    print()
    paired(af3_ag, td_ag,  "antigen-centric")
    paired(af3_tcr, td_tcr, "TCR-centric")

    # Write CSVs
    pd.DataFrame({'pmhc': list(af3_ag.keys()),
                  'AUC_AF3_PTI_PAE': list(af3_ag.values()),
                  'AUC_TCRdock_pmhc_tcr_pae': [td_ag.get(k, np.nan) for k in af3_ag],
                 }).to_csv(f'{args.prefix}_per_antigen_AUC.csv', index=False)
    pd.DataFrame({'pdbid': list(af3_tcr.keys()),
                  'AUC_AF3_PTI_PAE': list(af3_tcr.values()),
                  'AUC_TCRdock_pmhc_tcr_pae': [td_tcr.get(k, np.nan) for k in af3_tcr],
                 }).to_csv(f'{args.prefix}_per_TCR_AUC.csv', index=False)
    merged.to_csv(f'{args.prefix}_per_triad_long.csv', index=False)
    print(f"\nWrote:\n  {args.prefix}_per_antigen_AUC.csv\n  {args.prefix}_per_TCR_AUC.csv\n  {args.prefix}_per_triad_long.csv")

if __name__ == '__main__':
    main()
