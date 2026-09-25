"""
leakage_v6_blosum_triad.py

Single, uniform similarity metric (BLOSUM62 + length-gap distance with a
sliding-offset registration) applied to all three constituents of the
triad. Three output panels, one per component.

  Panel a:  per-antigen AUC vs peptide distance (val→nearest ref peptide)
  Panel b:  per-antigen AUC vs MHC distance     (val→nearest ref MHC chain;
                                                  class II sums α+β)
  Panel c:  per-TCR    AUC vs TCR distance      (val→nearest ref TCR;
                                                  Dash et al. weighted sum
                                                  of per-CDR distances
                                                  across α and β)

Distance kernel (used everywhere):
    blosum_sliding(a, b) = min over offsets of:
        gap_penalty × |Δlen|
      + Σ_i min(4, max(0, 4 - BLOSUM62(a_i, b_i)))
    weight is applied as Dash et al. specifies (×3 for CDR3, ×1 otherwise).

For peptide and MHC: the "sliding" picks the best register between the
shorter and longer sequence (so leading methionines, trailing tags, or
peptide-length differences don't shift everything).
For CDRs: positions are already canonicalized via ANARCI IMGT numbering,
so no sliding is needed within a CDR (it's just the per-position score
plus a length-gap term).

Class I uses the single MHC α chain. Class II uses both chains, summed.

Outputs land in notebooks/leakage_v6_blosum_triad/.
"""

import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRATCH = Path("/scratch/bneff/tcrtrifold")
REF_CSV = SCRATCH / "tcr3d_triad_sequences.csv"
SUPPT3  = SCRATCH / "supplementaryTable3_20260513.xlsx"
OUTDIR  = Path("/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/"
               "notebooks/leakage_v6_blosum_triad")
OUTDIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# BLOSUM62 table
# ---------------------------------------------------------------------------
BLOSUM62 = {}
_b62_str = """
A R N D C Q E G H I L K M F P S T W Y V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
_lines = [l.split() for l in _b62_str.strip().split("\n")]
_aa_order = _lines[0]
for i, row in enumerate(_lines[1:]):
    a1 = row[0]
    for j, a2 in enumerate(_aa_order):
        BLOSUM62[(a1, a2)] = int(row[j + 1])


# ---------------------------------------------------------------------------
# Distance kernels
# ---------------------------------------------------------------------------
def _per_position_blosum_penalty(c1, c2):
    """min(4, max(0, 4 - BLOSUM62(c1, c2))).  Identical residues → 0,
    strong mismatches → 4, intermediate → between."""
    bs = BLOSUM62.get((c1, c2), -4)
    return min(4, max(0, 4 - bs))


def blosum_gapless(a, b, gap_penalty=4):
    """Per-position BLOSUM62 sum with length-gap penalty, NO sliding.
    Used for already-aligned regions like ANARCI-numbered CDRs."""
    if not a and not b:
        return 0
    la, lb = len(a), len(b)
    L = min(la, lb)
    d = gap_penalty * abs(la - lb)
    for i in range(L):
        d += _per_position_blosum_penalty(a[i], b[i])
    return d


def blosum_sliding(a, b, gap_penalty=4):
    """Slide the shorter sequence along the longer one; at each offset
    compute the per-position BLOSUM penalty sum over the overlap region.
    The final distance is the minimum (best alignment) plus a length-gap
    penalty applied to the unaligned residues at the ends.

    Returns 0 if both inputs are empty."""
    if not isinstance(a, str) or not isinstance(b, str):
        return 0
    a = a.replace(".", "").replace("-", "")
    b = b.replace(".", "").replace("-", "")
    if not a and not b:
        return 0
    if not a or not b:
        return gap_penalty * max(len(a), len(b))

    if len(a) < len(b):
        a, b = b, a  # a is now the longer of the two
    La, Lb = len(a), len(b)
    n_offsets = La - Lb + 1

    # Length difference: residues in a that are not covered at the best offset.
    gap_term = gap_penalty * (La - Lb)

    best_inner = None
    for k in range(n_offsets):
        s = 0
        for i in range(Lb):
            s += _per_position_blosum_penalty(a[k + i], b[i])
        if best_inner is None or s < best_inner:
            best_inner = s
            if best_inner == 0:
                break  # can't do better than this
    return gap_term + best_inner


# ---------------------------------------------------------------------------
# ANARCI: extract all four CDRs from each TCR chain
# ---------------------------------------------------------------------------
# IMGT canonical CDR position ranges (inclusive integer-position numbering).
CDR_RANGES = {
    "cdr1":   (27,  38),
    "cdr2":   (56,  65),
    "cdr2_5": (81,  86),
    "cdr3":   (105, 117),
}
# Dash et al. 2017 CDR weights.
CDR_WEIGHTS = {"cdr1": 1, "cdr2": 1, "cdr2_5": 1, "cdr3": 3}


def extract_cdrs_from_numbering(numbering_entry):
    cdrs = {k: "" for k in CDR_RANGES}
    for (pos_tuple, aa) in numbering_entry:
        if aa == "-":
            continue
        pos_num = pos_tuple[0]
        for cdr, (lo, hi) in CDR_RANGES.items():
            if lo <= pos_num <= hi:
                cdrs[cdr] += aa
                break
    return cdrs


def anarci_chain(seq):
    try:
        from anarci import anarci
        numbering, alignment_details, _ = anarci(
            [("q", seq)], scheme="imgt", allow=("A", "B"),
            allowed_species=("human", "mouse"),
        )
        if not alignment_details[0]:
            return None
        det = alignment_details[0][0]
        cdrs = extract_cdrs_from_numbering(numbering[0][0][0])
        return {"chain_type": det["chain_type"],
                "species":    det.get("species"),
                **cdrs}
    except Exception:
        return None


def annotate_tcr_pair(seq_a, seq_b):
    ra = anarci_chain(seq_a)
    rb = anarci_chain(seq_b)
    if ra is None or rb is None:
        return None
    if ra["chain_type"] == "A" and rb["chain_type"] == "B":
        alpha, beta = ra, rb
    elif ra["chain_type"] == "B" and rb["chain_type"] == "A":
        alpha, beta = rb, ra
    else:
        return None
    return {"alpha":   {k: alpha[k] for k in CDR_RANGES},
            "beta":    {k: beta[k]  for k in CDR_RANGES},
            "species": alpha["species"] or beta["species"] or "human"}


# ---------------------------------------------------------------------------
# Triad distances
# ---------------------------------------------------------------------------
def peptide_distance(val_pep, ref_pep):
    """BLOSUM62 sliding distance between a validation peptide and a reference
    peptide."""
    return blosum_sliding(val_pep, ref_pep)


def mhc_distance(val_a, val_b, ref_a, ref_b, mhc_class):
    """For class I: BLOSUM62 sliding distance between the val and ref α chains.
    For class II: sum of α + β BLOSUM62 sliding distances."""
    if mhc_class == "I":
        return blosum_sliding(val_a, ref_a)
    else:
        return blosum_sliding(val_a, ref_a) + blosum_sliding(val_b, ref_b)


def tcr_paired_distance(val_alpha, ref_alpha, val_beta, ref_beta):
    """Sum over CDRs (with Dash weights) of gapless BLOSUM62 distances,
    then sum across α and β chains."""
    total = 0
    for cdr, w in CDR_WEIGHTS.items():
        total += w * blosum_gapless(val_alpha[cdr], ref_alpha[cdr])
        total += w * blosum_gapless(val_beta[cdr],  ref_beta[cdr])
    return total


def longer(a, b):
    if not isinstance(a, str): return b
    if not isinstance(b, str): return a
    return a if len(a) >= len(b) else b


# ---------------------------------------------------------------------------
# Reference loading
# ---------------------------------------------------------------------------
def load_reference():
    print("[1] Loading TCR3d reference …")
    ref = pd.read_csv(REF_CSV)
    pre = ref[(ref["extraction_status"] == "ok") &
              (ref["pre_af3_cutoff"].astype(bool))].copy().reset_index(drop=True)
    print(f"  total entries in CSV: {len(ref)}")
    print(f"  pre-cutoff + extraction ok: {len(pre)}")
    print(f"    by class: {pre['mhc_class'].value_counts().to_dict()}")

    pre["ref_pep"] = [longer(a, b) for a, b in
                      zip(pre["peptide_seq"], pre["epitope_listed"])]

    print(f"  ANARCIing {len(pre)} reference TCRs …")
    rows = []
    fail = 0
    for _, r in pre.iterrows():
        ann = annotate_tcr_pair(r["tcr_a_seq"], r["tcr_b_seq"])
        rec = r.to_dict()
        if ann is None:
            fail += 1
            rec["alpha"]   = {k: "" for k in CDR_RANGES}
            rec["beta"]    = {k: "" for k in CDR_RANGES}
            rec["species"] = None
            rec["tcr_ok"]  = False
        else:
            rec["alpha"]   = ann["alpha"]
            rec["beta"]    = ann["beta"]
            rec["species"] = ann["species"]
            rec["tcr_ok"]  = True
        rows.append(rec)
    print(f"    ANARCI: {len(rows) - fail} ok, {fail} failed")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Validation loading + AUCs
# ---------------------------------------------------------------------------
def load_validation():
    print("\n[2] Loading validation set (ST3) …")
    df = pd.read_excel(SUPPT3)
    df["pti_pae_score"] = 31 - df["mean_p_tcr_interface_pae"]
    df["cog"] = df["cognate"].astype(str).str.upper() == "TRUE"
    return df


def per_antigen_aucs(val):
    rows = []
    grp_cols = ["mhc_class", "peptide", "mhc_1_name", "mhc_2_name"]
    for keys, sub in val.groupby(grp_cols, dropna=False):
        if sub["cog"].sum() in (0, len(sub)):
            continue
        try:
            auc = roc_auc_score(sub["cog"].astype(int), sub["pti_pae_score"])
        except ValueError:
            continue
        rec = dict(zip(grp_cols, keys))
        rec["auc"]          = auc
        rec["n_cognate"]    = int(sub["cog"].sum())
        rec["n_noncognate"] = int(len(sub) - sub["cog"].sum())
        rec["peptide"]      = sub["peptide"].iloc[0]
        rec["mhc_1_seq"]    = sub["mhc_1_seq"].iloc[0]
        rec["mhc_2_seq"]    = sub["mhc_2_seq"].iloc[0]
        rows.append(rec)
    return pd.DataFrame(rows)


def per_tcr_aucs(val):
    val = val.copy()
    val["tcr_id"] = val["tcr_1_seq"] + "|" + val["tcr_2_seq"]
    rows = []
    for tcr_id, sub in val.groupby("tcr_id"):
        if sub["cog"].sum() in (0, len(sub)):
            continue
        try:
            auc = roc_auc_score(sub["cog"].astype(int), sub["pti_pae_score"])
        except ValueError:
            continue
        rows.append({
            "tcr_id":     tcr_id,
            "mhc_class":  sub["mhc_class"].iloc[0],
            "tcr_1_seq":  sub["tcr_1_seq"].iloc[0],
            "tcr_2_seq":  sub["tcr_2_seq"].iloc[0],
            "auc":        auc,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Per-component val→nearest-ref distances
# ---------------------------------------------------------------------------
def antigen_distances(ag_df, ref):
    """For each validation antigen, compute (peptide_distance, mhc_distance)
    to the nearest reference *for that component independently*."""
    print(f"\n[3] Peptide and MHC distances per antigen …")
    pep_d = np.zeros(len(ag_df))
    mhc_d = np.zeros(len(ag_df))

    ref_by_class = {
        "I":  ref[ref["mhc_class"] == "I"].to_dict("records"),
        "II": ref[ref["mhc_class"] == "II"].to_dict("records"),
    }
    print(f"  {len(ag_df)} val × ({len(ref_by_class['I'])} ref I, "
          f"{len(ref_by_class['II'])} ref II)")

    for i, val in ag_df.iterrows():
        cls = val["mhc_class"]
        refs = ref_by_class.get(cls, [])
        best_pep = np.inf
        best_mhc = np.inf
        for r in refs:
            d_pep = peptide_distance(val["peptide"], r["ref_pep"])
            if d_pep < best_pep:
                best_pep = d_pep
            d_mhc = mhc_distance(val["mhc_1_seq"], val["mhc_2_seq"],
                                 r["mhc_a_seq"],   r["mhc_b_seq"],
                                 cls)
            if d_mhc < best_mhc:
                best_mhc = d_mhc
        pep_d[i] = best_pep if best_pep != np.inf else np.nan
        mhc_d[i] = best_mhc if best_mhc != np.inf else np.nan
    return pep_d, mhc_d


def annotate_validation_tcrs(tcr_df):
    print(f"\n[4] ANARCIing {len(tcr_df)} validation TCRs …")
    rows = []
    fail = 0
    for _, r in tcr_df.iterrows():
        ann = annotate_tcr_pair(r["tcr_1_seq"], r["tcr_2_seq"])
        if ann is None:
            fail += 1
            continue
        rows.append({
            "tcr_id":    r["tcr_id"],
            "mhc_class": r["mhc_class"],
            "auc":       r["auc"],
            "alpha":     ann["alpha"],
            "beta":      ann["beta"],
            "species":   ann["species"],
        })
    print(f"  success: {len(rows)}, failed: {fail}")
    return pd.DataFrame(rows)


def tcr_distances(val_tcrs, ref):
    """For each validation TCR, return min full TCRdist to any reference TCR."""
    print(f"\n[5] TCR distances (BLOSUM62 + Dash et al. weights, all CDRs) …")
    ref_ok = ref[ref["tcr_ok"]].reset_index(drop=True)
    print(f"  usable reference TCRs: {len(ref_ok)} / {len(ref)}")
    print(f"  usable validation TCRs: {len(val_tcrs)}")

    if len(ref_ok) == 0:
        return np.full(len(val_tcrs), np.nan)

    ref_alphas = list(ref_ok["alpha"])
    ref_betas  = list(ref_ok["beta"])

    out = np.full(len(val_tcrs), np.nan)
    for i, r in val_tcrs.reset_index(drop=True).iterrows():
        va, vb = r["alpha"], r["beta"]
        if not isinstance(va, dict) or not isinstance(vb, dict):
            continue
        best = np.inf
        for ra, rb in zip(ref_alphas, ref_betas):
            d = tcr_paired_distance(va, ra, vb, rb)
            if d < best:
                best = d
        out[i] = best if best != np.inf else np.nan
    return out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_three_panel(ag_df, tcr_df, out_pdf, out_png):
    CI_BLUE = "#3a78b8"
    CII_RED = "#c83737"

    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size":   10,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.labelpad":     6,
        "pdf.fonttype":      42,
        "ps.fonttype":       42,
    })

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8),
                              gridspec_kw=dict(wspace=0.30,
                                               left=0.05, right=0.99,
                                               top=0.84, bottom=0.16))

    def _panel(ax, df, xcol, xlabel, ylabel, title, big_dot_cls="I"):
        ax.axhline(0.5, ls="--", lw=0.7, color="0.65", zorder=1)
        ax.grid(True, ls=":", lw=0.5, color="0.88", zorder=0)
        ax.set_axisbelow(True)

        # Plot class II first (background), class I on top
        for cls, color, msize, alpha, edge_w in [
            ("II", CII_RED,
             70 if len(df[df["mhc_class"]=="II"]) < 30 else 38,
             0.75 if len(df[df["mhc_class"]=="II"]) < 30 else 0.55,
             0.0),
            ("I", CI_BLUE, 90, 0.95, 0.9),
        ]:
            sub = df[df["mhc_class"] == cls].dropna(subset=[xcol, "auc"])
            ax.scatter(sub[xcol], sub["auc"],
                       s=msize, color=color,
                       edgecolor="black" if edge_w > 0 else "none",
                       linewidth=edge_w, alpha=alpha, zorder=3,
                       label=f"class {cls} (n = {len(sub)})")

        y0 = 0.04
        for cls, color in [("I", CI_BLUE), ("II", CII_RED)]:
            sub = df[df["mhc_class"]==cls].dropna(subset=[xcol, "auc"])
            if len(sub) >= 3 and sub[xcol].nunique() >= 2:
                r, p = spearmanr(sub[xcol], sub["auc"])
                ax.text(0.98, y0,
                        f"class {cls}: $r$ = {r:+.2f}, $p$ = {p:.2f}",
                        transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=9, color=color, fontweight="bold")
                y0 += 0.06

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, pad=8, fontsize=10.5)
        ax.set_ylim(-0.04, 1.08)

    _panel(axes[0], ag_df,
           xcol="pep_dist",
           xlabel=r"Peptide distance $\;d_{\rm pep}$  (BLOSUM62, sliding)",
           ylabel="Per-antigen AUC",
           title="Per-antigen AUC vs peptide distance")

    _panel(axes[1], ag_df,
           xcol="mhc_dist",
           xlabel=r"MHC distance $\;d_{\rm MHC}$  (BLOSUM62, sliding)",
           ylabel="Per-antigen AUC",
           title="Per-antigen AUC vs MHC distance")

    _panel(axes[2], tcr_df,
           xcol="tcr_dist",
           xlabel=r"TCR distance $\;d_{\rm TCR}$  (BLOSUM62, Dash weights)",
           ylabel="Per-TCR AUC",
           title="Per-TCR AUC vs TCR distance")

    # Shared legend on top
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=CI_BLUE, markeredgecolor="black", markersize=9,
               label="class I"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=CII_RED, markeredgecolor="none", markersize=9,
               alpha=0.7, label="class II"),
    ]
    fig.legend(handles=legend_handles, loc="upper center",
               bbox_to_anchor=(0.5, 1.00), ncol=2, frameon=False,
               handletextpad=0.4, columnspacing=2.0, fontsize=10)

    for ax, letter in zip(axes, ["a", "b", "c"]):
        ax.text(-0.13, 1.02, letter,
                transform=ax.transAxes, fontsize=15, fontweight="bold", va="top")

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Serialize nested CDR dicts before CSV write
# ---------------------------------------------------------------------------
def flatten_cdr_dicts(df):
    df = df.copy()
    for chain in ("alpha", "beta"):
        if chain in df.columns:
            for cdr in CDR_RANGES:
                col = f"{chain}_{cdr}"
                df[col] = df[chain].apply(
                    lambda d, cdr=cdr: d.get(cdr, "") if isinstance(d, dict) else "")
            df = df.drop(columns=[chain])
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 72)
    print("Leakage analysis v6 — uniform BLOSUM62 metric, three panels")
    print("=" * 72)

    ref = load_reference()

    val   = load_validation()
    ag_df = per_antigen_aucs(val)
    tcr_df = per_tcr_aucs(val)
    print(f"\n  per-antigen AUCs: {len(ag_df)}")
    print(f"  per-TCR AUCs:     {len(tcr_df)}")

    pep_d, mhc_d = antigen_distances(ag_df, ref)
    ag_df["pep_dist"] = pep_d
    ag_df["mhc_dist"] = mhc_d

    val_tcrs = annotate_validation_tcrs(tcr_df)
    val_tcrs["tcr_dist"] = tcr_distances(val_tcrs, ref)

    ag_out  = OUTDIR / "antigen_centric.csv"
    tcr_out = OUTDIR / "tcr_centric.csv"
    ref_out = OUTDIR / "reference_annotated.csv"

    ag_df.to_csv(ag_out, index=False)
    flatten_cdr_dicts(val_tcrs).to_csv(tcr_out, index=False)
    flatten_cdr_dicts(ref[["pdb_id","mhc_class","species","tcr_ok","alpha","beta"]]
                     ).to_csv(ref_out, index=False)

    plot_three_panel(ag_df, val_tcrs,
                     OUTDIR / "leakage_fig.pdf",
                     OUTDIR / "leakage_fig.png")

    print("\n" + "=" * 72)
    print("Summary correlations (Spearman)")
    print("=" * 72)
    print("Per-antigen AUC vs peptide distance:")
    for cls in ["I","II"]:
        sub = ag_df[ag_df["mhc_class"]==cls].dropna(subset=["pep_dist","auc"])
        if len(sub) >= 3 and sub["pep_dist"].nunique() >= 2:
            r,p = spearmanr(sub["pep_dist"], sub["auc"])
            print(f"  class {cls}: n={len(sub):>3}  r={r:+.3f}  p={p:.3f}")
        else:
            print(f"  class {cls}: n={len(sub):>3}  (low variance or low n)")
    print("Per-antigen AUC vs MHC distance:")
    for cls in ["I","II"]:
        sub = ag_df[ag_df["mhc_class"]==cls].dropna(subset=["mhc_dist","auc"])
        if len(sub) >= 3 and sub["mhc_dist"].nunique() >= 2:
            r,p = spearmanr(sub["mhc_dist"], sub["auc"])
            print(f"  class {cls}: n={len(sub):>3}  r={r:+.3f}  p={p:.3f}")
        else:
            print(f"  class {cls}: n={len(sub):>3}  (low variance or low n)")
    print("Per-TCR AUC vs TCR distance:")
    for cls in ["I","II"]:
        sub = val_tcrs[val_tcrs["mhc_class"]==cls].dropna(subset=["tcr_dist","auc"])
        if len(sub) >= 3:
            r,p = spearmanr(sub["tcr_dist"], sub["auc"])
            print(f"  class {cls}: n={len(sub):>3}  r={r:+.3f}  p={p:.3f}")
        else:
            print(f"  class {cls}: n={len(sub):>3}  (too few)")
    print()
    print("Wrote:")
    print(f"  {ag_out}")
    print(f"  {tcr_out}")
    print(f"  {ref_out}")
    print(f"  {OUTDIR/'leakage_fig.png'}")
    print(f"  {OUTDIR/'leakage_fig.pdf'}")


if __name__ == "__main__":
    main()
