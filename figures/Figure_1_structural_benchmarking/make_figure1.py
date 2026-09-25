"""
make_figure1.py
================
Regenerate Figure 1 of Woods et al. (TCR3d-expanded version).

Inputs (expected in the same directory as this script — see data/ folder):
    pdb_triad.af3_rmsd.parquet       - AF3 + AF-TCRdock RMSDs (300 rows, 269 unique PDBs)
    pdb_triad.boltz_rmsd.parquet     - Boltz-2 RMSDs (matched to AF3)
    pdb_triad.conf_af3.parquet       - AF3 confidence metrics (152 rows, Lawson's original)
    table_S1_structure_benchmark_complexes.csv  - Lawson's input PDB list (130 entries)

Outputs (PNG + PDF, in ./figures/):
    figure1_panel_b{_lawson130,}.{pdf,png}  - 3-method CDR RMSD boxplots
    figure1_panel_c{_lawson130,}.{pdf,png}  - Per-component RMSDs for AF3
    figure1_panel_d{_lawson130,}.{pdf,png}  - AF3 ipTM vs CDR RMSD scatter

The `_lawson130` versions subset to Lawson's original 130 PDBs plus the
post-AF3-cutoff STCRDab additions, reproducing the manuscript Figure 1.
The non-suffixed versions use the full TCR3d-expanded dataset (~282 PDBs
once the clean_pdb fix is applied; currently 269 entries pending pipeline rerun).

Run:
    python3 make_figure1.py
    # → writes figures to ./figures/
"""

from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

# ----------------------------- CONFIG -----------------------------

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / 'data'
FIG_DIR  = HERE / 'figures'
FIG_DIR.mkdir(exist_ok=True)

AF3_PARQUET   = DATA_DIR / 'pdb_triad.af3_rmsd.parquet'
BOLTZ_PARQUET = DATA_DIR / 'pdb_triad.boltz_rmsd.parquet'
CONF_PARQUET  = DATA_DIR / 'pdb_triad.conf_af3.parquet'
TABLE_S1      = DATA_DIR / 'table_S1_structure_benchmark_complexes.csv'

# Per-method training cutoffs (each model trained at a different time).
# Cutoffs verified against Lawson's notebooks/fig_1/3_method_comparison_rmsd.ipynb:
#   af2_cutoff = dt.datetime(2018, 5, 1) ; boltz_cutoff = dt.datetime(2023, 6, 1)
# AF3 cutoff is from the AF3 paper (Abramson et al. 2024).
METHOD_CUTOFFS = {
    'AF-TCRdock': pd.Timestamp('2018-05-01', tz='UTC'),
    'AF3':        pd.Timestamp('2023-01-12', tz='UTC'),
    'Boltz-2':    pd.Timestamp('2023-06-01', tz='UTC'),
}

# ----------------------------- STYLE -----------------------------

mpl.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':         10,
    'axes.titlesize':    11,
    'axes.labelsize':    10,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.linewidth':    0.8,
    'xtick.labelsize':   9,
    'ytick.labelsize':   9,
    'xtick.major.size':  3,
    'ytick.major.size':  3,
    'legend.fontsize':   9,
    'legend.frameon':    False,
    'pdf.fonttype':      42,    # editable text in PDF
    'ps.fonttype':       42,
})

# Class color: blue for class I, red for class II (used for panel C/D)
COLOR_CI       = '#185FA5'
COLOR_CI_LIGHT = '#85B7EB'
COLOR_CII      = '#A32D2D'
COLOR_CII_LIGHT = '#F09595'

# Method colors for panel B
METHOD_COLORS = {
    'AF-TCRdock': '#888780',
    'AF3':        '#534AB7',
    'Boltz-2':    '#1D9E75',
}
METHOD_FILLS = {
    'AF-TCRdock': '#D3D1C7',
    'AF3':        '#CECBF6',
    'Boltz-2':    '#9FE1CB',
}


# ----------------------------- LOAD -----------------------------

def load_data():
    """Load and lightly clean the three parquets."""
    af3   = pd.read_parquet(AF3_PARQUET)
    boltz = pd.read_parquet(BOLTZ_PARQUET)
    conf  = pd.read_parquet(CONF_PARQUET)

    # Dedupe multi-chain crystal duplicates (e.g. 6l9l × 8, others × 4).
    # All duplicate rows share identical RMSDs, so first-row keep is lossless.
    af3   = af3.drop_duplicates('pdb', keep='first').reset_index(drop=True)
    boltz = boltz.drop_duplicates('pdb', keep='first').reset_index(drop=True)

    # Coerce types
    af3['pdb_date']   = pd.to_datetime(af3['pdb_date'], utc=True)
    boltz['pdb_date'] = pd.to_datetime(boltz['pdb_date'], utc=True)
    af3['replication']   = af3['replication'].astype(bool)
    boltz['replication'] = boltz['replication'].astype(bool)
    conf['replication']  = conf['replication'].astype(bool)

    print(f"  af3:   {len(af3):>4} rows, {af3['pdb'].nunique()} unique PDBs")
    print(f"  boltz: {len(boltz):>4} rows, {boltz['pdb'].nunique()} unique PDBs")
    print(f"  conf:  {len(conf):>4} rows, {conf['pdb'].nunique()} unique PDBs (Lawson's original)")
    return af3, boltz, conf


def subset_to_lawson_original(af3, boltz, conf):
    """
    Subset to Lawson's original Figure 1 set (130 in table_S1 + post-cutoff STCRDab additions).
    This reproduces the manuscript's Figure 1 from our pipeline output.
    """
    if not TABLE_S1.exists():
        raise FileNotFoundError(f"Missing {TABLE_S1}; needed for the Lawson-130 reproducibility subset.")
    orig = pd.read_csv(TABLE_S1)
    lawson_130 = set(orig['pdbid'].str.lower())

    # Take union of:
    #   (a) the 130 PDBs in Bradley/Lawson's table_S1 (the "pre-cutoff" benchmark)
    #   (b) the post-cutoff STCRDab additions Lawson used in his manuscript Figure 1.
    #
    # The post-cutoff set is hardcoded here based on Lawson's original
    # pdb_triad.conf_af3.parquet (which had replication=False for these specific
    # PDBs). Our current conf parquet has 3 extras (8es8, 8shi for class I; 8trl
    # for class II) that we exclude from the reproducibility subset so it matches
    # the manuscript's Figure 1b numbers exactly.
    LAWSON_POST_CI = {
        '7q99','7q9a','7q9b','8dnt','8en8','8enh','8eo8','8es9','8f5a','8gom',
        '8gon','8i5c','8i5d','8wte','8wul','8ye4',
    }
    LAWSON_POST_CII = {
        '8pjg','8trr','8vcx','8vcy','8vd0','8vd2',
    }
    lawson_post = LAWSON_POST_CI | LAWSON_POST_CII
    keep_pdbs = lawson_130 | lawson_post

    af3_s   = af3  [af3  ['pdb'].str.lower().isin(keep_pdbs)].reset_index(drop=True)
    boltz_s = boltz[boltz['pdb'].str.lower().isin(keep_pdbs)].reset_index(drop=True)

    print(f"  Lawson-130 subset: af3 {len(af3_s)}, boltz {len(boltz_s)}, "
          f"target was {len(lawson_130) + len(lawson_post)}")
    return af3_s, boltz_s


# ----------------------------- PANEL B -----------------------------

def panel_b(af3, boltz, suffix=''):
    """3-method comparison: AF-TCRdock vs AF3 vs Boltz-2, CDR RMSD, by class and pre/post.
    Now draws Mann-Whitney pre-vs-post p-value brackets at the top of each subpanel."""
    from scipy.stats import mannwhitneyu

    records = []
    for _, row in af3.iterrows():
        if pd.notna(row['cdr_rmsd']):
            records.append({
                'pdb': row['pdb'], 'mhc_class': row['mhc_class'], 'method': 'AF-TCRdock',
                'pre_training': row['pdb_date'] < METHOD_CUTOFFS['AF-TCRdock'],
                'cdr_rmsd': row['cdr_rmsd'],
            })
        if pd.notna(row['cdr_rmsd_af3_4']):
            records.append({
                'pdb': row['pdb'], 'mhc_class': row['mhc_class'], 'method': 'AF3',
                'pre_training': row['pdb_date'] < METHOD_CUTOFFS['AF3'],
                'cdr_rmsd': row['cdr_rmsd_af3_4'],
            })
    for _, row in boltz.iterrows():
        if pd.notna(row['cdr_rmsd_boltz_4']):
            records.append({
                'pdb': row['pdb'], 'mhc_class': row['mhc_class'], 'method': 'Boltz-2',
                'pre_training': row['pdb_date'] < METHOD_CUTOFFS['Boltz-2'],
                'cdr_rmsd': row['cdr_rmsd_boltz_4'],
            })
    tidy = pd.DataFrame(records)

    print('  per-method, per-class sample sizes (pre / post):')
    for m in ['AF-TCRdock', 'AF3', 'Boltz-2']:
        for c in ['I', 'II']:
            sub = tidy[(tidy['method']==m) & (tidy['mhc_class']==c)]
            pre  = (sub['pre_training'] == True).sum()
            post = (sub['pre_training'] == False).sum()
            print(f'    {m:<11} class {c}:  {pre:>4} pre / {post:>4} post')

    def _fmt_p(p):
        if pd.isna(p): return ''
        if p < 1e-3:  return f"p={p:.0e}"
        if p < 0.01:  return f"p={p:.3f}"
        return f"p={p:.2f}"

    def _draw_bracket(ax, x1, x2, y, label, h=0.55, lw=1.0, color='0.25', fontsize=10):
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], color=color, lw=lw, clip_on=False)
        ax.text((x1+x2)/2, y+h+0.25, label, ha='center', va='bottom',
                fontsize=fontsize, color=color, clip_on=False)

    # narrower figure: compress horizontally, slightly taller for bigger fonts
    fig, axes = plt.subplots(2, 1, figsize=(5.8, 7.4), sharex=False)
    methods = ['AF-TCRdock', 'AF3', 'Boltz-2']
    GROUP_DX = 2.2          # compressed from 3.0
    BOX_W = 0.85            # wider from 0.7
    OFFSET = 0.50           # pre/post centers at -OFFSET / +OFFSET (gap stays small)

    for ax_i, (cls, cls_title, cls_color) in enumerate([
        ('I',  'MHC class I',  COLOR_CI),
        ('II', 'MHC class II', COLOR_CII),
    ]):
        ax = axes[ax_i]
        ticks_x, ticks_labels = [], []
        method_positions = {}
        for m_i, method in enumerate(methods):
            group_x = m_i * GROUP_DX
            pos_pre  = group_x - OFFSET
            pos_post = group_x + OFFSET
            method_positions[method] = (pos_pre, pos_post)
            for pos, pre_flag in [(pos_pre, True), (pos_post, False)]:
                sub = tidy[
                    (tidy['method']==method) &
                    (tidy['mhc_class']==cls) &
                    (tidy['pre_training']==pre_flag)
                ]['cdr_rmsd'].dropna().values
                if len(sub) == 0:
                    continue
                bp = ax.boxplot([sub], positions=[pos], widths=BOX_W,
                                patch_artist=True, showfliers=False,
                                medianprops=dict(color='black', linewidth=1.6),
                                boxprops=dict(linewidth=1.0),
                                whiskerprops=dict(linewidth=1.0),
                                capprops=dict(linewidth=1.0),
                                whis=(5, 95))
                fill = METHOD_FILLS[method] if pre_flag else 'white'
                edge = METHOD_COLORS[method]
                bp['boxes'][0].set(facecolor=fill, edgecolor=edge)
                for w in bp['whiskers'] + bp['caps']:
                    w.set(color=edge)
                rng = np.random.default_rng(seed=int(m_i*100 + (0 if pre_flag else 10) + (1 if cls=='I' else 2)))
                jitter = rng.uniform(-0.22, 0.22, size=len(sub))
                ax.scatter(np.full(len(sub), pos) + jitter, sub,
                           s=7, alpha=0.50, color=edge, edgecolor='none', zorder=3)
                ax.text(pos, -2.9, f'n={len(sub)}', ha='center', va='top',
                        fontsize=10, color=edge, fontweight='bold')

            ticks_x.append(group_x)
            cutoff_str = METHOD_CUTOFFS[method].strftime('%Y-%m')
            ticks_labels.append(f"{method}\n(cutoff {cutoff_str})")

        # ---- Mann-Whitney cross-method brackets ----
        # 4 comparisons per subpanel:
        #   AF-TCRdock-pre  vs AF3-pre        ;  AF-TCRdock-post vs AF3-post
        #   AF3-pre         vs Boltz-2-pre    ;  AF3-post        vs Boltz-2-post
        def _vals(method, pre_flag):
            return tidy[(tidy['method']==method) & (tidy['mhc_class']==cls) &
                        (tidy['pre_training']==pre_flag)]['cdr_rmsd'].dropna().values

        comparisons = [
            ('AF-TCRdock', 'AF3',     True,  25.0, 'pre'),
            ('AF-TCRdock', 'AF3',     False, 28.0, 'post'),
            ('AF3',        'Boltz-2', True,  31.0, 'pre'),
            ('AF3',        'Boltz-2', False, 34.0, 'post'),
        ]
        for m1, m2, pre_flag, y_bracket, tag in comparisons:
            a = _vals(m1, pre_flag); b = _vals(m2, pre_flag)
            if len(a) < 2 or len(b) < 2:
                continue
            # one-tailed: H1 is that the older model (a) has higher RMSD than the newer (b)
            _, p = mannwhitneyu(a, b, alternative='greater')
            x1 = method_positions[m1][0 if pre_flag else 1]
            x2 = method_positions[m2][0 if pre_flag else 1]
            _draw_bracket(ax, x1, x2, y_bracket, f"{tag}: {_fmt_p(p)}", fontsize=10)

        ax.set_xticks(ticks_x)
        ax.set_xticklabels(ticks_labels, fontsize=11)
        ax.set_ylabel('CDR RMSD (Å)', fontsize=13)
        ax.set_title(cls_title, loc='left', color=cls_color, fontweight='bold', fontsize=13)
        ax.set_ylim(-4.8, 36)
        ax.set_yticks([0, 5, 10, 15, 20])
        ax.tick_params(axis='y', labelsize=11)
        ax.axhline(0, color='#888888', linewidth=0.7, zorder=1)
        ax.yaxis.grid(True, linewidth=0.3, alpha=0.4)
        ax.set_axisbelow(True)
        ax.set_xlim(-1.3, methods.__len__() * GROUP_DX - 1.3 + 0.7)

    handles = [
        plt.Rectangle((0,0),1,1, facecolor='#bbbbbb', edgecolor='#444444',
                      label='Pre-training (deposited before model cutoff)'),
        plt.Rectangle((0,0),1,1, facecolor='white',   edgecolor='#444444',
                      label='Post-training (deposited after model cutoff)'),
    ]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.00),
               ncol=2, frameon=False, fontsize=10.5)
    fig.suptitle('AF-TCRdock vs AlphaFold 3 vs Boltz-2: CDR RMSD on TCR3d structures',
                 y=1.03, fontsize=12.5, fontweight='bold')
    fig.tight_layout()
    save(fig, f'figure1_panel_b{suffix}')

    # also print the p-values (sanity check)
    print('\n  Mann-Whitney U-test (one-tailed, H1: a > b), pre vs post:')
    for m in ['AF-TCRdock', 'AF3', 'Boltz-2']:
        for c in ['I', 'II']:
            pre  = tidy[(tidy['method']==m) & (tidy['mhc_class']==c) &
                        (tidy['pre_training']==True)]['cdr_rmsd'].dropna().values
            post = tidy[(tidy['method']==m) & (tidy['mhc_class']==c) &
                        (tidy['pre_training']==False)]['cdr_rmsd'].dropna().values
            if len(pre) > 0 and len(post) > 0:
                u, p = mannwhitneyu(post, pre, alternative='greater')
                med_pre, med_post = np.median(pre), np.median(post)
                print(f'    {m:<11} class {c}: n_pre={len(pre):>3} med={med_pre:5.2f} | '
                      f'n_post={len(post):>3} med={med_post:5.2f} | p={p:.2e}')



# ----------------------------- PANEL C -----------------------------

def panel_c(af3, suffix=''):
    """Per-component AF3 RMSDs: peptide, MHC, TCR. Class I left, class II right.
    Now with wider boxes, Mann-Whitney pre-vs-post brackets, bigger fonts."""
    from scipy.stats import mannwhitneyu

    components = [
        ('peptide_rmsd_af3_4', 'Peptide'),
        ('mhc_rmsd_af3_4',     'MHC'),
        ('tcr_rmsd_af3_4',     'TCR'),
    ]

    def _fmt_p(p):
        if pd.isna(p): return ''
        if p < 1e-3:  return f"p={p:.0e}"
        if p < 0.01:  return f"p={p:.3f}"
        return f"p={p:.2f}"

    def _draw_bracket(ax, x1, x2, y, label, h=0.18, lw=1.0, color='0.25', fontsize=10):
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], color=color, lw=lw, clip_on=False)
        ax.text((x1+x2)/2, y+h+0.1, label, ha='center', va='bottom',
                fontsize=fontsize, color=color, clip_on=False)

    GROUP_DX = 2.2
    BOX_W    = 0.85
    OFFSET   = 0.50

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 5.0), sharey=True)
    for ax_i, (cls, cls_title, cls_color, fill_color) in enumerate([
        ('I',  'MHC class I',  COLOR_CI,  COLOR_CI_LIGHT),
        ('II', 'MHC class II', COLOR_CII, COLOR_CII_LIGHT),
    ]):
        ax = axes[ax_i]
        ticks_x, ticks_labels = [], []
        comp_positions = {}
        for ci, (col, label) in enumerate(components):
            group_x = ci * GROUP_DX
            pos_pre  = group_x - OFFSET
            pos_post = group_x + OFFSET
            comp_positions[col] = (pos_pre, pos_post)
            for pos, pre_flag in [(pos_pre, True), (pos_post, False)]:
                sub = af3[
                    (af3['mhc_class']==cls) &
                    ((af3['pdb_date'] < METHOD_CUTOFFS['AF3']) == pre_flag)
                ][col].dropna().values
                if len(sub) == 0:
                    continue
                fill = fill_color if pre_flag else 'white'
                edge = cls_color
                bp = ax.boxplot([sub], positions=[pos], widths=BOX_W,
                                patch_artist=True, showfliers=False,
                                medianprops=dict(color='black', linewidth=1.6),
                                boxprops=dict(linewidth=1.0),
                                whiskerprops=dict(linewidth=1.0),
                                capprops=dict(linewidth=1.0),
                                whis=(5, 95))
                bp['boxes'][0].set(facecolor=fill, edgecolor=edge)
                for w in bp['whiskers'] + bp['caps']:
                    w.set(color=edge)
                rng = np.random.default_rng(seed=int(ci*100 + (0 if pre_flag else 10) + (1 if cls=='I' else 2)))
                jitter = rng.uniform(-0.22, 0.22, size=len(sub))
                ax.scatter(np.full(len(sub), pos) + jitter, sub,
                           s=6, alpha=0.45, color=edge, edgecolor='none', zorder=3)

            ticks_x.append(group_x)
            ticks_labels.append(label)

        # ---- pre-vs-post Mann-Whitney brackets per component (one-tailed: pre < post) ----
        bracket_y = 7.5
        for col, label in components:
            pre  = af3[(af3['mhc_class']==cls) & ((af3['pdb_date'] < METHOD_CUTOFFS['AF3']) == True)][col].dropna().values
            post = af3[(af3['mhc_class']==cls) & ((af3['pdb_date'] < METHOD_CUTOFFS['AF3']) == False)][col].dropna().values
            if len(pre) > 1 and len(post) > 1:
                # H1: pre RMSD < post RMSD  (model expected to be worse on unseen data)
                _, p = mannwhitneyu(pre, post, alternative='less')
                xp, xq = comp_positions[col]
                _draw_bracket(ax, xp, xq, bracket_y, _fmt_p(p), fontsize=10)

        ax.set_xticks(ticks_x)
        ax.set_xticklabels(ticks_labels, fontsize=12)
        ax.set_title(cls_title, loc='left', color=cls_color, fontweight='bold', fontsize=13)
        if ax_i == 0:
            ax.set_ylabel('RMSD (Å, AF3 best of 5)', fontsize=13)
        ax.set_ylim(-0.6, 9.0)
        ax.set_yticks([0, 2, 4, 6, 8])
        ax.tick_params(axis='y', labelsize=11)
        ax.axhline(0, color='#888888', linewidth=0.7, zorder=1)
        ax.yaxis.grid(True, linewidth=0.3, alpha=0.4)
        ax.set_axisbelow(True)
        ax.set_xlim(-1.3, len(components) * GROUP_DX - 1.3 + 0.7)

    handles = [
        plt.Rectangle((0,0),1,1, facecolor=COLOR_CI_LIGHT, edgecolor=COLOR_CI, label='Pre-training'),
        plt.Rectangle((0,0),1,1, facecolor='white',         edgecolor=COLOR_CI, label='Post-training'),
    ]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.02),
               ncol=2, frameon=False, fontsize=10.5)
    fig.suptitle('AlphaFold 3: per-component RMSDs',
                 y=1.07, fontsize=12.5, fontweight='bold')
    fig.tight_layout()
    save(fig, f'figure1_panel_c{suffix}')



# ----------------------------- PANEL D -----------------------------

def _dgeom_to_ndarr(dgeom_series):
    """Convert a series of dgeom dicts into a (N, 7) ndarray (cos-sin embedding on torsion).
    Matches the embedding used by Lawson's tcrdock_utils.cossin_embed.
    """
    arr = np.array([
        [d['d'], d['torsion'], d['tcr_unit_y'], d['tcr_unit_z'], d['mhc_unit_y'], d['mhc_unit_z']]
        for d in dgeom_series
    ])
    # cos-sin embed on the torsion (col index 1)
    out = np.empty((arr.shape[0], 7))
    out[:, 0] = arr[:, 0]
    out[:, 1] = np.sin(arr[:, 1])
    out[:, 2] = np.cos(arr[:, 1])
    out[:, 3:] = arr[:, 2:]
    return out


def _mn_distance(dgeom_ndarr, mu, invcov):
    """Mahalanobis distance from the consensus (mu, invcov)."""
    diff = dgeom_ndarr - mu
    dm2 = np.sum((diff @ invcov) * diff, axis=1)
    return np.sqrt(dm2)


def panel_d(af3, conf, suffix=''):
    """4-subpanel: AF3 confidence vs CDR RMSD; docking Z-score; peptide length; species.
    Mirrors Figure 1d of the manuscript. Bigger text and wider species boxes.
    """
    conf_dedup = conf.drop_duplicates('job_name', keep='first')
    merged = af3[['pdb','mhc_class','mhc_1_species','peptide','job_name','cdr_rmsd_af3_4',
                  'true_dgeom','pred_dgeom_4']].merge(
        conf_dedup[['job_name', 'iptm', 'ptm', 'ranking_score']],
        on='job_name', how='inner',
    )
    print(f"  merged for panel d: {len(merged)} rows ({merged['pdb'].nunique()} unique PDBs)")

    merged = merged.copy()
    merged['dgeom_zscore'] = np.nan
    valid_dgeom = merged['true_dgeom'].notna() & merged['pred_dgeom_4'].notna()
    for cls in ['I', 'II']:
        cls_mask = (merged['mhc_class'] == cls) & valid_dgeom
        if cls_mask.sum() < 7:
            continue
        true_arr = _dgeom_to_ndarr(merged.loc[cls_mask, 'true_dgeom'])
        pred_arr = _dgeom_to_ndarr(merged.loc[cls_mask, 'pred_dgeom_4'])
        mu = true_arr.mean(axis=0)
        cov = np.cov(true_arr, rowvar=False)
        invcov = np.linalg.inv(cov)
        merged.loc[cls_mask, 'dgeom_zscore'] = _mn_distance(pred_arr, mu, invcov)

    merged['peptide_len'] = merged['peptide'].str.len()

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0))
    A_LBL=16; A_TICK=13; LEG=11; TITLE=16; STAT=12

    # (a) ipTM vs CDR RMSD
    ax = axes[0, 0]
    for cls, cls_title, color in [('I','MHC class I',COLOR_CI),('II','MHC class II',COLOR_CII)]:
        sub = merged[merged['mhc_class']==cls].dropna(subset=['iptm','cdr_rmsd_af3_4'])
        ax.scatter(sub['iptm'], sub['cdr_rmsd_af3_4'], s=22, alpha=0.55,
                   color=color, edgecolor='none', label=f"{cls_title} (n={len(sub)})")
    valid = merged.dropna(subset=['iptm','cdr_rmsd_af3_4'])
    if len(valid) > 5:
        r, p = stats.spearmanr(valid['iptm'], valid['cdr_rmsd_af3_4'])
        ax.text(0.03, 0.97, f"Spearman r = {r:+.2f}\np = {p:.1e}",
                transform=ax.transAxes, va='top', ha='left', fontsize=STAT, color='#444441',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', linewidth=0.5))
    ax.set_xlabel('ipTM (AF3 confidence)', fontsize=A_LBL)
    ax.set_ylabel('CDR RMSD (\u00c5)', fontsize=A_LBL)
    ax.tick_params(labelsize=A_TICK)
    ax.set_xlim(0.4, 1.0); ax.set_ylim(-0.5, 25)
    ax.yaxis.grid(True, linewidth=0.3, alpha=0.4); ax.set_axisbelow(True)
    ax.legend(loc='upper right', frameon=False, fontsize=LEG)
    ax.set_title('AF3 confidence vs accuracy', loc='left', fontsize=TITLE, fontweight='bold')

    # (b) docking geometry Z-score vs CDR RMSD
    ax = axes[0, 1]
    for cls, cls_title, color in [('I','MHC class I',COLOR_CI),('II','MHC class II',COLOR_CII)]:
        sub = merged[merged['mhc_class']==cls].dropna(subset=['dgeom_zscore','cdr_rmsd_af3_4'])
        ax.scatter(sub['dgeom_zscore'], sub['cdr_rmsd_af3_4'], s=22, alpha=0.55,
                   color=color, edgecolor='none', label=f"{cls_title} (n={len(sub)})")
    valid = merged.dropna(subset=['dgeom_zscore','cdr_rmsd_af3_4'])
    if len(valid) > 5:
        r, p = stats.spearmanr(valid['dgeom_zscore'], valid['cdr_rmsd_af3_4'])
        ax.text(0.97, 0.97, f"Spearman r = {r:+.2f}\np = {p:.1e}",
                transform=ax.transAxes, va='top', ha='right', fontsize=STAT, color='#444441',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', linewidth=0.5))
    ax.set_xlabel('Z-score of predicted docking geometry', fontsize=A_LBL)
    ax.set_ylabel('CDR RMSD (\u00c5)', fontsize=A_LBL)
    ax.tick_params(labelsize=A_TICK)
    ax.set_ylim(-0.5, 25)
    ax.yaxis.grid(True, linewidth=0.3, alpha=0.4); ax.set_axisbelow(True)
    ax.legend(loc='upper left', frameon=False, fontsize=LEG)
    ax.set_title('Docking geometry vs accuracy', loc='left', fontsize=TITLE, fontweight='bold')

    # (c) Peptide length vs CDR RMSD
    ax = axes[1, 0]
    for cls, cls_title, color in [('I','MHC class I',COLOR_CI),('II','MHC class II',COLOR_CII)]:
        sub = merged[merged['mhc_class']==cls].dropna(subset=['peptide_len','cdr_rmsd_af3_4'])
        rng = np.random.default_rng(42 if cls=='I' else 1337)
        jit = rng.uniform(-0.2, 0.2, size=len(sub))
        ax.scatter(sub['peptide_len']+jit, sub['cdr_rmsd_af3_4'], s=22, alpha=0.55,
                   color=color, edgecolor='none', label=f"{cls_title} (n={len(sub)})")
    valid = merged.dropna(subset=['peptide_len','cdr_rmsd_af3_4'])
    if len(valid) > 5:
        r, p = stats.spearmanr(valid['peptide_len'], valid['cdr_rmsd_af3_4'])
        ax.text(0.03, 0.97, f"Spearman r = {r:+.2f}\np = {p:.2f}",
                transform=ax.transAxes, va='top', ha='left', fontsize=STAT, color='#444441',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', linewidth=0.5))
    ax.set_xlabel('Peptide length (residues)', fontsize=A_LBL)
    ax.set_ylabel('CDR RMSD (\u00c5)', fontsize=A_LBL)
    ax.tick_params(labelsize=A_TICK)
    ax.set_ylim(-0.5, 25)
    ax.yaxis.grid(True, linewidth=0.3, alpha=0.4); ax.set_axisbelow(True)
    ax.legend(loc='upper right', frameon=False, fontsize=LEG)
    ax.set_title('Peptide length vs accuracy', loc='left', fontsize=TITLE, fontweight='bold')

    # (d) Species boxplot -- WIDER boxes (0.55 -> 0.85)
    ax = axes[1, 1]
    species_order = ['human', 'mouse']
    species_colors = {'human': '#5680A8', 'mouse': '#A86E56'}
    box_positions, box_labels = [], []
    for i, sp in enumerate(species_order):
        sub = merged[merged['mhc_1_species']==sp].dropna(subset=['cdr_rmsd_af3_4'])
        if len(sub) == 0:
            continue
        bp = ax.boxplot([sub['cdr_rmsd_af3_4'].values],
                        positions=[i], widths=0.85,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color='black', linewidth=1.6),
                        boxprops=dict(linewidth=1.0, facecolor=species_colors[sp], alpha=0.5,
                                      edgecolor=species_colors[sp]),
                        whiskerprops=dict(linewidth=1.0, color=species_colors[sp]),
                        capprops=dict(linewidth=1.0, color=species_colors[sp]),
                        whis=(5, 95))
        rng = np.random.default_rng(7 if sp=='human' else 8)
        jit = rng.uniform(-0.25, 0.25, size=len(sub))
        ax.scatter(np.full(len(sub), i)+jit, sub['cdr_rmsd_af3_4'].values,
                   s=14, alpha=0.5, color=species_colors[sp], edgecolor='none', zorder=3)
        box_positions.append(i)
        box_labels.append(f"{sp}\n(n={len(sub)})")
    h = merged[merged['mhc_1_species']=='human']['cdr_rmsd_af3_4'].dropna().values
    m = merged[merged['mhc_1_species']=='mouse']['cdr_rmsd_af3_4'].dropna().values
    if len(h) > 5 and len(m) > 5:
        u, p = stats.mannwhitneyu(h, m, alternative='two-sided')
        ax.text(0.97, 0.97, f"Mann-Whitney U\np = {p:.2f}",
                transform=ax.transAxes, va='top', ha='right', fontsize=STAT, color='#444441',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', linewidth=0.5))
    ax.set_xticks(box_positions)
    ax.set_xticklabels(box_labels, fontsize=14)
    ax.set_ylabel('CDR RMSD (\u00c5)', fontsize=A_LBL)
    ax.tick_params(axis='y', labelsize=A_TICK)
    ax.set_ylim(-0.5, 25); ax.set_xlim(-0.7, 1.7)
    ax.yaxis.grid(True, linewidth=0.3, alpha=0.4); ax.set_axisbelow(True)
    ax.set_title('Species vs accuracy', loc='left', fontsize=TITLE, fontweight='bold')

    fig.suptitle('AF3 prediction accuracy: confidence, docking geometry, peptide length, species',
                 y=0.99, fontsize=15, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, f'figure1_panel_d{suffix}')



# ----------------------------- HELPERS -----------------------------

def save(fig, basename):
    for ext in ['png']:
        out = FIG_DIR / f'{basename}.{ext}'
        fig.savefig(out, dpi=300, bbox_inches='tight')
        print(f"  → {out}")
    plt.close(fig)


# ----------------------------- MAIN -----------------------------

if __name__ == '__main__':
    print('=' * 60)
    print('Loading data...')
    af3, boltz, conf = load_data()

    print('\n--- FULL TCR3d-EXPANDED DATASET ---')
    print('Building panel B (full)...')
    panel_b(af3, boltz, suffix='')
    print('\nBuilding panel C (full)...')
    panel_c(af3, suffix='')
    print('\nBuilding panel D (full)...')
    panel_d(af3, conf, suffix='')

    if TABLE_S1.exists():
        print('\n--- LAWSON-ORIGINAL SUBSET (reproducibility check) ---')
        print('Subsetting to Lawson\'s 130 + post-cutoff STCRDab additions...')
        af3_lawson, boltz_lawson = subset_to_lawson_original(af3, boltz, conf)
        print('\nBuilding panel B (Lawson subset)...')
        panel_b(af3_lawson, boltz_lawson, suffix='_lawson130')
        print('\nBuilding panel C (Lawson subset)...')
        panel_c(af3_lawson, suffix='_lawson130')
        print('\nBuilding panel D (Lawson subset)...')
        panel_d(af3_lawson, conf, suffix='_lawson130')
    else:
        print(f'\n(No {TABLE_S1.name} found — skipping Lawson-original subset)')

    print('\n' + '=' * 60)
    print(f'Done. Figures in: {FIG_DIR}')
