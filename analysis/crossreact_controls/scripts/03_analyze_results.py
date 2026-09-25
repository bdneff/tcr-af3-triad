"""
analyze_results.py
------------------

Post-inference analysis script. Run this AFTER the tcrtrifold-experiments
pipeline has produced the conf_af3 parquet (or CSV) for the cross-reactivity
controls.

Inputs
------

- A results file with one row per triad and (at minimum) these columns:
      name (or job_name)               -- MD5-based unique ID from script 02
      mean_p_tcr_interface_pae         -- raw AF3 output

  Accepts either .parquet (preferred; matches pipeline output) or .csv.

- The pipeline raw CSV from script 02 (joined on `name`) which carries:
      class, bl62_distance, hamming_distance, mhc_class, peptide

Outputs
-------

- `data/cross_reactivity_controls_with_PTIPAE.csv`
      One row per triad with AF3 PTI-PAE, class label, BLOSUM62 distance,
      MHC class. Ready for plotting.

- `data/cross_reactivity_controls_plot.png`
      Two-panel scatter (class I left, class II right) of AF3-PAE-PTI
      vs. BLOSUM62 distance, with the three classes color-coded:
          green   : original
          blue    : cross-reactive
          pink    : randomized
      The dashed green line shows the PTI-PAE of the original peptide
      (the "AF3 thinks this peptide binds" threshold). This is the plot
      that goes into Figure X / Supplementary Figure X of the revision.

Notes on the PTI-PAE definition
-------------------------------

Per the manuscript's main classifier:
      PTI-PAE = 31 - mean_p_tcr_interface_pae

(31 is the maximum value of the AF3 PAE bin; this inversion flips the sign
so higher PTI-PAE = stronger predicted binding, matching cognate status.)

This script computes PTI-PAE on the fly from `mean_p_tcr_interface_pae`,
so the input file just needs the raw AF3 value.

Usage
-----

    python scripts/03_analyze_results.py                              \\
        --results data/cross_reactivity_controls.conf_af3.parquet     \\
        --raw     data/cross_reactivity_controls_pipeline.csv         \\
        --out_csv data/cross_reactivity_controls_with_PTIPAE.csv      \\
        --out_png data/cross_reactivity_controls_plot.png
"""

import argparse
import csv
import sys
from pathlib import Path


def load_results(path: Path) -> list[dict]:
    """Load the AF3 feature-extraction output (.parquet or .csv)."""
    if path.suffix == '.parquet':
        try:
            import polars as pl
        except ImportError:
            print('ERROR: polars not installed; cannot read .parquet. '
                  'Run: pip install polars  (or pass a .csv instead)', file=sys.stderr)
            sys.exit(1)
        df = pl.read_parquet(path)
        # Drop nested columns that don't round-trip to dicts cleanly
        scalar_cols = [c for c, dt in zip(df.columns, df.dtypes)
                       if not str(dt).startswith(('List', 'Struct', 'Array'))]
        return df.select(scalar_cols).to_dicts()
    elif path.suffix == '.csv':
        with open(path) as f:
            return list(csv.DictReader(f))
    else:
        raise ValueError(f'unsupported results extension: {path.suffix}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', type=Path, required=True,
                    help='AF3 feature extraction output (.parquet or .csv)')
    ap.add_argument('--raw', type=Path,
                    default=Path('data/cross_reactivity_controls_pipeline.csv'),
                    help='Pipeline raw CSV from script 02 (joined on `name`)')
    ap.add_argument('--out_csv', type=Path,
                    default=Path('data/cross_reactivity_controls_with_PTIPAE.csv'))
    ap.add_argument('--out_png', type=Path,
                    default=Path('data/cross_reactivity_controls_plot.png'))
    args = ap.parse_args()

    # --- Load and join ---
    results = load_results(args.results)
    name_col = 'name' if 'name' in results[0] else 'job_name'
    pae_col  = 'mean_p_tcr_interface_pae'
    if pae_col not in results[0]:
        raise KeyError(f'expected column {pae_col!r} in results; got {list(results[0])}')

    with open(args.raw) as f:
        raw = {r['name']: r for r in csv.DictReader(f)}

    joined = []
    missing = 0
    for r in results:
        meta = raw.get(r[name_col])
        if meta is None:
            missing += 1
            continue
        try:
            pae = float(r[pae_col])
        except (TypeError, ValueError):
            continue
        pti_pae = 31.0 - pae
        joined.append({
            'name':            r[name_col],
            'mhc_class':       meta['mhc_class'],
            'peptide':         meta['peptide'],
            'class':           meta['class'],
            'bl62_distance':   int(meta['bl62_distance']),
            'hamming_distance': meta['hamming_distance'],
            'mean_p_tcr_interface_pae': pae,
            'PTI_PAE':         pti_pae,
        })

    if missing:
        print(f'WARN: {missing} result rows had no metadata match (skipped)', file=sys.stderr)

    # --- Write joined CSV ---
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(joined[0]))
        w.writeheader()
        w.writerows(joined)
    print(f'Wrote {len(joined)} rows to {args.out_csv}', file=sys.stderr)

    # --- Plot ---
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not installed; skipping plot. Install with: pip install matplotlib',
              file=sys.stderr)
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    COLORS = {
        'original':       '#1b9e77',  # green
        'cross-reactive': '#3498db',  # blue
        'randomized':     '#e377c2',  # pink
    }
    ORDER = ['randomized', 'cross-reactive', 'original']  # draw originals on top

    for ax, klass, title in [
        (axes[0], 'I',  'Class I (B*57:03, TCR fixed)\nKAFSPEVIPMF'),
        (axes[1], 'II', 'Class II (DRB3*03:01, TCR fixed)\nPPQIAANRSQLISLV'),
    ]:
        sub = [j for j in joined if j['mhc_class'] == klass]
        for label in ORDER:
            pts = [j for j in sub if j['class'] == label]
            if not pts:
                continue
            xs = [p['bl62_distance'] for p in pts]
            ys = [p['PTI_PAE']       for p in pts]
            sz = 80 if label == 'original' else (40 if label == 'cross-reactive' else 18)
            alpha = 1.0 if label == 'original' else (0.85 if label == 'cross-reactive' else 0.55)
            ax.scatter(xs, ys, s=sz, c=COLORS[label], alpha=alpha,
                       edgecolors='k' if label == 'original' else 'none',
                       linewidths=1, label=label)

        # Reference line: PTI-PAE of the original peptide for this class
        orig_pts = [j for j in sub if j['class'] == 'original']
        if orig_pts:
            ax.axhline(orig_pts[0]['PTI_PAE'], color=COLORS['original'],
                       linestyle='--', linewidth=1, alpha=0.7)

        ax.set_xlabel('BLOSUM62 distance from original peptide')
        ax.set_title(title)
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel('AF3 PTI-PAE  (= 31 − mean_p_tcr_interface_pae)')
    axes[1].legend(loc='lower left', framealpha=0.9, title='peptide source')

    plt.tight_layout()
    plt.savefig(args.out_png, dpi=150, bbox_inches='tight')
    print(f'Wrote plot to {args.out_png}', file=sys.stderr)


if __name__ == '__main__':
    main()
