#!/usr/bin/env python3
"""plot_supp3c.py — Supp Fig 3c: CDR RMSD vs AF3 PTI-PAE, class I and class II.

Extends the 16-point class I panel that `06_structural.ipynb` (cell 13) produces, to answer
**reviewer round 2, minor point 1** — "if additional post-training structures are currently
available, repeat this analysis using a larger benchmark set."

Everything about the class I panel is preserved deliberately: the same input CSV, the same
`31 - pti_pae` transform, the same fixed 0-31 x axis, the same `#3a78b8`, the same stats box in the
figure margin. A reader comparing the published panel with this one should see the original sixteen
points unmoved, with new ones added — not a redrawn figure. `06_structural.ipynb` is kept as the
record of the original version rather than edited, per this repo's convention that superseded
artifacts are archived, not overwritten.

Class II points use `#A32D2D`, the class II red already used throughout Figure 1
(`figures/Figure_1_structural_benchmarking/make_figure1.py`), so the two figures agree about what
red means.

**Three correlations are printed and the figure reports all of them**, because a single combined
number would obscure the thing a reviewer will want to check: that adding class II did not move the
original class I result. Class I alone must still read r = -0.53, rho = -0.60.

    python plot_supp3c.py                    # -> plot_rmsd_vs_pae_classI_II.png
    python plot_supp3c.py --check            # regression: class I stats must match the published

Input `merged_tv_rmsd_pae.csv` needs columns: pdb, cdr_rmsd, pti_pae, mhc_class.
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.transforms
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

CSV = "merged_tv_rmsd_pae.csv"
OUT = "plot_rmsd_vs_pae_classI_II.png"
COLOR = {"I": "#3a78b8", "II": "#A32D2D"}      # I: as published; II: Figure 1's class II red
LABEL = {"I": "MHC class I", "II": "MHC class II"}
PUBLISHED_CLASS_I = {"n": 16, "pearson": -0.53, "spearman": -0.60}


# Candidate label offsets in points, tried in order: right, left, above, below, then
# progressively further out. The first that collides with nothing wins.
LABEL_OFFSETS = [(7, 3), (-25, 3), (7, -10), (-25, -10), (0, 11), (-9, -15),
                 (14, 9), (-30, 9), (14, -15), (-30, -15), (20, 2), (-36, 2),
                 (0, 20), (-9, -24), (24, 14), (-40, 14), (24, -20), (-40, -20),
                 (0, 28), (-9, -32), (30, 22), (-46, 22),
                 # Escape tier: the far-right cluster (8enh/8eo8/8en8, all within
                 # 0.6 PTI-PAE and 0.1 A of each other) has nowhere to go sideways,
                 # but the space above it is empty. These reach it, with leaders.
                 (8, 42), (-46, 42), (8, 58), (-46, 58),
                 (8, 74), (-46, 74), (8, 92), (-46, 92)]


def declutter(ax, fig, df, fontsize=8, color="0.25"):
    """Place each PDB label where it overlaps neither a marker nor an already-placed label.

    The published panel annotated every point at a fixed (6, 3) offset, which was legible
    at 16 points and is not at 23 — the bottom-right cluster overprinted into mush and
    8gom/8gon vanished under the legend. Greedy placement with a leader line when the
    label has to sit far from its point.
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()

    # Seed the occupied regions with the markers themselves, plus the legend, so labels
    # never land on data or on the key.
    taken = []
    for x, y in zip(df["pti_pae_31"], df["cdr_rmsd"]):
        px, py = ax.transData.transform((x, y))
        taken.append(matplotlib.transforms.Bbox.from_bounds(px - 7, py - 7, 14, 14))
    leg = ax.get_legend()
    if leg is not None:
        taken.append(leg.get_window_extent(renderer=rend).expanded(1.05, 1.05))

    # A label that runs off the axes is as unreadable as one that overprints — several
    # spilled below the x axis before this was enforced.
    frame = ax.get_window_extent(renderer=rend)

    # Topmost first: the crowded high-RMSD points get first choice of position.
    rows = df.sort_values("cdr_rmsd", ascending=False)
    for _, r in rows.iterrows():
        x, y, txt = r["pti_pae_31"], r["cdr_rmsd"], r["pdb"]
        for i, (dx, dy) in enumerate(LABEL_OFFSETS):
            t = ax.annotate(txt, (x, y), xytext=(dx, dy), textcoords="offset points",
                            fontsize=fontsize, color=color, zorder=4)
            bb = t.get_window_extent(renderer=rend).expanded(1.06, 1.12)
            inside = (bb.x0 >= frame.x0 and bb.x1 <= frame.x1
                      and bb.y0 >= frame.y0 and bb.y1 <= frame.y1)
            if inside and not any(bb.overlaps(o) for o in taken):
                # Far-flung labels get a hairline leader so the association is unambiguous.
                if i >= 6:
                    ax.annotate("", (x, y), xytext=(dx, dy),
                                textcoords="offset points", zorder=1,
                                arrowprops=dict(arrowstyle="-", color="0.65", lw=0.5,
                                                shrinkA=0, shrinkB=4))
                taken.append(bb)
                break
            t.remove()
        else:
            # Nowhere clean; place at the default rather than dropping the label, and
            # say so, because a silently missing label is worse than a crowded one.
            t = ax.annotate(txt, (x, y), xytext=LABEL_OFFSETS[0],
                            textcoords="offset points", fontsize=fontsize,
                            color=color, zorder=4)
            taken.append(t.get_window_extent(renderer=rend))
            print(f"  note: no clear position for {txt}; placed at default offset")


def stats(df):
    if len(df) < 3:
        return None
    rp, pp = pearsonr(df["pti_pae_31"], df["cdr_rmsd"])
    rs, ps = spearmanr(df["pti_pae_31"], df["cdr_rmsd"])
    return dict(n=len(df), rp=rp, pp=pp, rs=rs, ps=ps)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=CSV)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--check", action="store_true",
                    help="verify the class I subset still reproduces the published statistics")
    a = ap.parse_args()

    df = pd.read_csv(a.csv)
    if "mhc_class" not in df:
        sys.exit(f"{a.csv} has no mhc_class column — add one ('I' for the original 16)")
    df = df.dropna(subset=["cdr_rmsd", "pti_pae"]).copy()
    df["pti_pae_31"] = 31 - df["pti_pae"]

    groups = {c: df[df["mhc_class"] == c] for c in ("I", "II") if (df["mhc_class"] == c).any()}
    allst = stats(df)
    per = {c: stats(g) for c, g in groups.items()}

    for c, st in per.items():
        if st:
            print(f"class {c:2}  n={st['n']:3}  Pearson r={st['rp']:+.3f} p={st['pp']:.3f}   "
                  f"Spearman rho={st['rs']:+.3f} p={st['ps']:.3f}")
    if allst and len(groups) > 1:
        print(f"combined  n={allst['n']:3}  Pearson r={allst['rp']:+.3f} p={allst['pp']:.3f}   "
              f"Spearman rho={allst['rs']:+.3f} p={allst['ps']:.3f}")

    # Regression guard. The point of adding class II is that the class I claim is UNCHANGED; if
    # these drift, something upstream altered the original rows and the figure must not be trusted.
    if a.check:
        st = per.get("I")
        if not st:
            sys.exit("no class I rows to check")
        bad = []
        if st["n"] != PUBLISHED_CLASS_I["n"]:
            bad.append(f"n {st['n']} != {PUBLISHED_CLASS_I['n']}")
        if abs(round(st["rp"], 2) - PUBLISHED_CLASS_I["pearson"]) > 0.005:
            bad.append(f"Pearson {st['rp']:+.2f} != {PUBLISHED_CLASS_I['pearson']:+.2f}")
        if abs(round(st["rs"], 2) - PUBLISHED_CLASS_I["spearman"]) > 0.005:
            bad.append(f"Spearman {st['rs']:+.2f} != {PUBLISHED_CLASS_I['spearman']:+.2f}")
        if bad:
            sys.exit("REGRESSION — class I no longer matches the published figure: " + "; ".join(bad))
        print(f"check passed: class I still n={st['n']}, r={st['rp']:+.2f}, rho={st['rs']:+.2f}")

    fig, ax = plt.subplots(figsize=(8, 6))
    for c, g in groups.items():
        ax.scatter(g["pti_pae_31"], g["cdr_rmsd"], s=110, color=COLOR[c],
                   edgecolors="black", linewidths=0.7, alpha=0.85,
                   label=f"{LABEL[c]}  (n={len(g)})", zorder=3)

    # Per-class trend lines, not one pooled fit. Pooling the two classes gives a weak,
    # non-significant slope (they sit at different PTI-PAE offsets), and a single
    # confident dashed line across everything would assert exactly the claim this
    # figure cannot support.
    for c, g in groups.items():
        if len(g) < 3:
            continue
        sl, ic = np.polyfit(g["pti_pae_31"], g["cdr_rmsd"], 1)
        xs = np.linspace(g["pti_pae_31"].min() - 0.4, g["pti_pae_31"].max() + 0.4, 50)
        ax.plot(xs, sl * xs + ic, "--", color=COLOR[c], alpha=0.55, lw=1.4, zorder=2)

    ax.set_xlabel("AF3 PTI-PAE  (higher = more confident interface)", fontsize=17)
    ax.set_ylabel("CDR RMSD (Å)  —  AF3 vs crystal", fontsize=17)
    ax.set_title("Post-AF3-cutoff cognate triads", fontsize=15, pad=8)
    ax.set_xlim(0, 31)
    # Headroom so the topmost labels (8gom/8gon, the reverse-dock failures) are not
    # clipped by the top spine.
    ax.set_ylim(0, df["cdr_rmsd"].max() * 1.16)
    ax.tick_params(labelsize=15)
    if len(groups) > 1:
        # Upper LEFT: every point sits above PTI-PAE 16, so the left third is empty and
        # the legend no longer covers the two highest-RMSD points.
        ax.legend(fontsize=12, loc="upper left", framealpha=0.95)

    declutter(ax, fig, df)

    lines = []
    for c, st in per.items():
        if st:
            lines.append(f"class {c}: n={st['n']}  r={st['rp']:+.2f} (p={st['pp']:.3f})  "
                         f"ρ={st['rs']:+.2f} (p={st['ps']:.3f})")
    if allst and len(groups) > 1:
        lines.append(f"all:     n={allst['n']}  r={allst['rp']:+.2f} (p={allst['pp']:.3f})  "
                     f"ρ={allst['rs']:+.2f} (p={allst['ps']:.3f})")
    fig.text(0.5, 0.985, "\n".join(lines), ha="center", va="top", fontsize=11,
             family="monospace",
             bbox=dict(facecolor="white", edgecolor="0.6", boxstyle="round,pad=0.45", alpha=0.95))
    ax.grid(True, alpha=0.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.subplots_adjust(top=0.86 if len(lines) < 3 else 0.80)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    print(f"\nSaved: {a.out}")


if __name__ == "__main__":
    main()
