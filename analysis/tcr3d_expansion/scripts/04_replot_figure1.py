#!/usr/bin/env python
"""
04_replot_figure1.py
--------------------
Regenerate Figure 1 panels b, c, and d using the expanded RMSD parquets from
the pipeline run, with the TCR3d additions merged in automatically.

Panels regenerated:
  1b — CDR RMSD comparison: AF-TCRdock vs AF3 vs Boltz-2, stratified by
       pre/post-training-cutoff, separately for class I (top) and class II
       (bottom).
  1c — Per-chain RMSD breakdown for AF3: peptide / MHC / TCR / CDR, stratified
       by class.
  1d — AF3 confidence (ipTM) vs RMSD correlation; docking-geometry Z-score
       vs RMSD.

Color theme:
  Class I  = blue   palette  (#1f4e79, #4f81bd, #95b3d7, with darker = post-cutoff)
  Class II = red    palette  (#9c0d28, #c0504d, #d99694)
  AF-TCRdock / AF3 / Boltz-2 distinguished within a class by saturation +
  hatching (light=pre-training, solid darker=post-training).

Reads the actual pipeline output parquets (preserved from the upstream code).

INPUTS (paths default to repo-relative; override on command line)
-----------------------------------------------------------------
--af3_rmsd        : pdb_triad.af3_rmsd.parquet
--boltz_rmsd      : pdb_triad.boltz_rmsd.parquet
--out_dir         : where to save figures

OUTPUTS
-------
fig1b_cdr_rmsd_compare.png/pdf
fig1c_rmsd_by_chain.png/pdf
fig1d_confidence_correlations.png/pdf
"""

import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import polars as pl
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu, spearmanr

# -------------------------------------------------------------------------
# Color palette
# -------------------------------------------------------------------------
CLASS_COLORS = {
    "I": {
        "AF2_pre":   "#a0c8e0",
        "AF2_post":  "#1f4e79",
        "AF3_pre":   "#7faed4",
        "AF3_post":  "#0f4a78",
        "Boltz_pre": "#5896c4",
        "Boltz_post":"#003366",
    },
    "II": {
        "AF2_pre":   "#e8a39e",
        "AF2_post":  "#9c0d28",
        "AF3_pre":   "#e07e6c",
        "AF3_post":  "#7a0a1f",
        "Boltz_pre": "#c95041",
        "Boltz_post":"#5a0817",
    },
}
CLASS_TITLE_COLOR = {"I": "#1f4e79", "II": "#9c0d28"}

# Cutoff dates per model
AF2_CUTOFF   = dt.datetime(2018, 5, 1)
AF3_CUTOFF   = dt.datetime(2023, 1, 12)
BOLTZ_CUTOFF = dt.datetime(2023, 6, 1)


def _setup_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 8,
    })


def _annot_p(ax, x1, x2, y, p, text=None):
    """Draw a p-value bracket between two x positions."""
    if text is None:
        if p < 1e-3:
            text = f"p={p:.1e}"
        else:
            text = f"p={p:.3f}"
        if p >= 0.05:
            text = "n.s."
    ax.plot([x1, x1, x2, x2], [y, y * 1.02, y * 1.02, y], lw=0.7, color="black")
    ax.text((x1 + x2) / 2, y * 1.035, text, ha="center", va="bottom",
            fontsize=7.5, color="black")


# -------------------------------------------------------------------------
# Panel 1b — CDR RMSD comparison across methods
# -------------------------------------------------------------------------

def plot_fig1b(af3_df: pl.DataFrame, boltz_df: pl.DataFrame, out_dir: Path):
    """Two-row figure: class I top, class II bottom. Each row has 3 method
    groups (AF-TCRdock, AF3, Boltz-2), each with pre/post-cutoff boxes."""
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    for ax, mhc_class in zip(axes, ("I", "II")):
        af3_c = af3_df.filter(pl.col("mhc_class") == mhc_class)
        boltz_c = boltz_df.filter(pl.col("mhc_class") == mhc_class)

        # Define groups: (method, cutoff, color_key, values)
        groups = []
        # AF-TCRdock (cdr_rmsd column, cutoff at AF2_CUTOFF)
        af2_pre = af3_c.filter(pl.col("af2_pre_cutoff")).select("cdr_rmsd").drop_nulls().to_series().to_numpy()
        af2_post = af3_c.filter(~pl.col("af2_pre_cutoff")).select("cdr_rmsd").drop_nulls().to_series().to_numpy()
        af3_pre = af3_c.filter(pl.col("replication")).select("cdr_rmsd_af3_4").drop_nulls().to_series().to_numpy()
        af3_post = af3_c.filter(~pl.col("replication")).select("cdr_rmsd_af3_4").drop_nulls().to_series().to_numpy()
        boltz_pre = boltz_c.filter(pl.col("boltz_pre_cutoff")).select("cdr_rmsd_boltz_4").drop_nulls().to_series().to_numpy()
        boltz_post = boltz_c.filter(~pl.col("boltz_pre_cutoff")).select("cdr_rmsd_boltz_4").drop_nulls().to_series().to_numpy()

        positions = [0.65, 1.0, 1.95, 2.30, 3.25, 3.60]
        data = [af2_pre, af2_post, af3_pre, af3_post, boltz_pre, boltz_post]
        colors = [
            CLASS_COLORS[mhc_class]["AF2_pre"],
            CLASS_COLORS[mhc_class]["AF2_post"],
            CLASS_COLORS[mhc_class]["AF3_pre"],
            CLASS_COLORS[mhc_class]["AF3_post"],
            CLASS_COLORS[mhc_class]["Boltz_pre"],
            CLASS_COLORS[mhc_class]["Boltz_post"],
        ]
        labels = [
            f"AF-TCRdock pre (n={len(af2_pre)})",
            f"AF-TCRdock post (n={len(af2_post)})",
            f"AF3 pre (n={len(af3_pre)})",
            f"AF3 post (n={len(af3_post)})",
            f"Boltz-2 pre (n={len(boltz_pre)})",
            f"Boltz-2 post (n={len(boltz_post)})",
        ]

        for pos, vals, col in zip(positions, data, colors):
            if len(vals) == 0:
                continue
            bp = ax.boxplot(
                vals, positions=[pos], widths=0.30, patch_artist=True,
                boxprops=dict(facecolor=col, edgecolor="black", linewidth=0.7),
                whiskerprops=dict(color="black", linewidth=0.7),
                capprops=dict(color="black", linewidth=0.7),
                medianprops=dict(color="white", linewidth=1.5),
                flierprops=dict(markeredgecolor=col, markersize=3, marker="o", alpha=0.5),
            )
            # Overlay individual points (jittered)
            jitter = np.random.RandomState(42).uniform(-0.07, 0.07, size=len(vals))
            ax.scatter(np.full_like(vals, pos) + jitter, vals,
                       s=6, color=col, edgecolor="black", linewidth=0.3, alpha=0.6, zorder=3)

        # Stats: pre vs post for each method
        try:
            _, p_af2  = mannwhitneyu(af2_pre, af2_post, alternative="less") if len(af2_pre)>0 and len(af2_post)>0 else (None, None)
            _, p_af3  = mannwhitneyu(af3_pre, af3_post, alternative="less") if len(af3_pre)>0 and len(af3_post)>0 else (None, None)
            _, p_bolt = mannwhitneyu(boltz_pre, boltz_post, alternative="less") if len(boltz_pre)>0 and len(boltz_post)>0 else (None, None)

            # Pick a y for annotation just above all boxes (ignore outliers)
            ymax = max([np.percentile(v, 90) for v in data if len(v) > 0]) * 1.25
            if p_af2 is not None: _annot_p(ax, positions[0], positions[1], ymax, p_af2)
            if p_af3 is not None: _annot_p(ax, positions[2], positions[3], ymax, p_af3)
            if p_bolt is not None: _annot_p(ax, positions[4], positions[5], ymax, p_bolt)
            ax.set_ylim(0, ymax * 1.25)
        except Exception as e:
            print(f"  (skipping stat annotations for class {mhc_class}: {e})")

        ax.set_xticks([0.825, 2.125, 3.425])
        ax.set_xticklabels(["AF-TCRdock", "AlphaFold 3", "Boltz-2"])
        ax.set_ylabel("CDR RMSD (Å)")
        ax.set_title(f"MHC class {mhc_class}",
                     color=CLASS_TITLE_COLOR[mhc_class], fontweight="bold", loc="left")

        # Compact legend (sample counts)
        handles = [Patch(facecolor=c, edgecolor="black", label=l)
                   for c, l in zip(colors, labels) if len(data[colors.index(c)]) > 0]
        ax.legend(handles=handles, loc="upper right", frameon=False, ncol=3, fontsize=7.5)

    fig.suptitle("CDR RMSD: AlphaFold 3 vs prior structural-prediction methods",
                 fontsize=13, y=0.99)
    fig.tight_layout()

    for ext in ("png", "pdf"):
        out = out_dir / f"fig1b_cdr_rmsd_compare.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"  wrote {out}")
    plt.close(fig)


# -------------------------------------------------------------------------
# Panel 1c — Per-chain RMSD breakdown (AF3)
# -------------------------------------------------------------------------

def plot_fig1c(af3_df: pl.DataFrame, out_dir: Path):
    """For each MHC class (rows), show peptide / MHC / TCR RMSDs side-by-side
    with pre/post-cutoff sub-boxes."""
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    chain_types = ["peptide", "mhc", "tcr"]
    chain_labels = ["Peptide", "MHC", "TCR"]

    for ax, mhc_class in zip(axes, ("I", "II")):
        c = af3_df.filter(pl.col("mhc_class") == mhc_class)
        col_pre  = CLASS_COLORS[mhc_class]["AF3_pre"]
        col_post = CLASS_COLORS[mhc_class]["AF3_post"]

        for i, (chain, label) in enumerate(zip(chain_types, chain_labels)):
            colname = f"{chain}_rmsd_af3_4"
            if colname not in c.columns:
                # Some panels may use a different naming convention
                continue
            pre  = c.filter(pl.col("replication")).select(colname).drop_nulls().to_series().to_numpy()
            post = c.filter(~pl.col("replication")).select(colname).drop_nulls().to_series().to_numpy()
            pos_pre, pos_post = i - 0.18, i + 0.18

            for pos, vals, col in [(pos_pre, pre, col_pre), (pos_post, post, col_post)]:
                if len(vals) == 0:
                    continue
                ax.boxplot(
                    vals, positions=[pos], widths=0.30, patch_artist=True,
                    boxprops=dict(facecolor=col, edgecolor="black", linewidth=0.7),
                    whiskerprops=dict(color="black", linewidth=0.7),
                    capprops=dict(color="black", linewidth=0.7),
                    medianprops=dict(color="white", linewidth=1.5),
                    flierprops=dict(markeredgecolor=col, markersize=3, marker="o", alpha=0.5),
                )
                jitter = np.random.RandomState(42 + i).uniform(-0.07, 0.07, size=len(vals))
                ax.scatter(np.full_like(vals, pos) + jitter, vals,
                           s=5, color=col, edgecolor="black", linewidth=0.3, alpha=0.55, zorder=3)

            # Pre/post test
            try:
                if len(pre) > 0 and len(post) > 0:
                    _, p = mannwhitneyu(pre, post, alternative="less")
                    ymax = max(np.percentile(pre, 90) if len(pre) else 0,
                               np.percentile(post, 90) if len(post) else 0) * 1.15
                    _annot_p(ax, pos_pre, pos_post, ymax, p)
            except Exception:
                pass

        ax.set_xticks(range(len(chain_types)))
        ax.set_xticklabels(chain_labels)
        ax.set_ylabel("RMSD (Å)")
        ax.set_title(f"MHC class {mhc_class}",
                     color=CLASS_TITLE_COLOR[mhc_class], fontweight="bold", loc="left")

        handles = [
            Patch(facecolor=col_pre, edgecolor="black", label="Pre-training cutoff"),
            Patch(facecolor=col_post, edgecolor="black", label="Post-training cutoff"),
        ]
        ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=8)

    fig.suptitle("AlphaFold 3 per-chain RMSD breakdown",
                 fontsize=13, y=0.99)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = out_dir / f"fig1c_rmsd_by_chain.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"  wrote {out}")
    plt.close(fig)


# -------------------------------------------------------------------------
# Panel 1d — AF3 confidence vs RMSD
# -------------------------------------------------------------------------

def plot_fig1d(af3_df: pl.DataFrame, out_dir: Path):
    """Scatter: ipTM vs CDR RMSD, points colored by class. Plus a second
    panel for docking-geometry Z-score if available."""
    has_dgeom = "pred_dgeom_zscore" in af3_df.columns
    if has_dgeom:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(6, 5))

    for mhc_class in ("I", "II"):
        c = af3_df.filter(pl.col("mhc_class") == mhc_class)
        x = c.select("iptm").to_series().to_numpy()
        y = c.select("cdr_rmsd_af3_4").to_series().to_numpy()
        mask = ~np.isnan(x) & ~np.isnan(y)
        x, y = x[mask], y[mask]
        if len(x) > 0:
            col = CLASS_COLORS[mhc_class]["AF3_post"]
            ax1.scatter(x, y, s=14, c=col, edgecolor="black", linewidth=0.3, alpha=0.65,
                        label=f"Class {mhc_class} (n={len(x)})")
            rho, p = spearmanr(x, y)
            print(f"  ipTM vs CDR RMSD, class {mhc_class}: rho={rho:.3f}, p={p:.2e}")

    ax1.set_xlabel("ipTM (AlphaFold 3 confidence)")
    ax1.set_ylabel("CDR RMSD (Å)")
    ax1.set_title("Confidence vs prediction accuracy",
                  fontweight="bold", loc="left")
    ax1.legend(frameon=False, fontsize=9)
    ax1.grid(True, alpha=0.3, linewidth=0.5)

    if has_dgeom:
        for mhc_class in ("I", "II"):
            c = af3_df.filter(pl.col("mhc_class") == mhc_class)
            x = c.select("pred_dgeom_zscore").to_series().to_numpy()
            y = c.select("cdr_rmsd_af3_4").to_series().to_numpy()
            mask = ~np.isnan(x) & ~np.isnan(y)
            x, y = x[mask], y[mask]
            if len(x) > 0:
                col = CLASS_COLORS[mhc_class]["AF3_post"]
                ax2.scatter(x, y, s=14, c=col, edgecolor="black", linewidth=0.3, alpha=0.65,
                            label=f"Class {mhc_class} (n={len(x)})")
                rho, p = spearmanr(x, y)
                print(f"  dgeom Z vs CDR RMSD, class {mhc_class}: rho={rho:.3f}, p={p:.2e}")
        ax2.set_xlabel("Docking-geometry Z-score")
        ax2.set_ylabel("CDR RMSD (Å)")
        ax2.set_title("Docking geometry vs prediction accuracy",
                      fontweight="bold", loc="left")
        ax2.legend(frameon=False, fontsize=9)
        ax2.grid(True, alpha=0.3, linewidth=0.5)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = out_dir / f"fig1d_confidence_correlations.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"  wrote {out}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--af3_rmsd",  type=Path, required=True)
    ap.add_argument("--boltz_rmsd",type=Path, required=True)
    ap.add_argument("--out_dir",   type=Path, default=Path("results/fig1_v2"))
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    _setup_style()

    # Load and annotate with per-model cutoff booleans (same as original notebook)
    af3 = pl.read_parquet(args.af3_rmsd).with_columns(
        pl.when(pl.col("pdb_date") < AF2_CUTOFF)
        .then(pl.lit(True)).otherwise(pl.lit(False))
        .alias("af2_pre_cutoff")
    ).filter(pl.col("pdb") != "8trr")
    boltz = pl.read_parquet(args.boltz_rmsd).with_columns(
        pl.when(pl.col("pdb_date") < BOLTZ_CUTOFF)
        .then(pl.lit(True)).otherwise(pl.lit(False))
        .alias("boltz_pre_cutoff")
    ).filter(pl.col("pdb") != "8trr")

    print(f"Loaded AF3 RMSD: {af3.height} entries"
          f" ({af3.filter(pl.col('mhc_class')=='I').height} I + "
          f"{af3.filter(pl.col('mhc_class')=='II').height} II)")
    print(f"Loaded Boltz RMSD: {boltz.height} entries")
    print(f"  pre/post (AF3 replication):  "
          f"{af3.filter(pl.col('replication')).height} pre, "
          f"{af3.filter(~pl.col('replication')).height} post")
    print()

    print("Drawing panel b (CDR RMSD comparison)...")
    plot_fig1b(af3, boltz, args.out_dir)
    print()
    print("Drawing panel c (per-chain RMSD breakdown)...")
    plot_fig1c(af3, args.out_dir)
    print()
    print("Drawing panel d (confidence vs RMSD)...")
    plot_fig1d(af3, args.out_dir)
    print()
    print(f"Done. See figures in {args.out_dir}/")


if __name__ == "__main__":
    main()
