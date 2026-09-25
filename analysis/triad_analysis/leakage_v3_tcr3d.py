"""
leakage_v3_tcr3d.py

Leakage analysis using TCR3d as the reference set.

Reference: 308 pre-AF3-cutoff TCR-pMHC triads (226 class I + 82 class II),
           extracted from TCR3d via tcr3d_triad_sequences.csv.

Validation: updated ST3 (PDB class I + CRESTA class II validation triads).

Metrics:
  - Peptide identity: sliding-window Hamming over the overlap region.
                      Identity = matches / overlap_length, max over all offsets.
  - MHC identity:     Same sliding-window Hamming on chain sequences.
                      Class I: α chain only. Class II: min(α, β).
  - Composite pMHC distance:
        pmhc_distance = 1 - max over refs of min(I_peptide, I_mhc)
        Reference peptide = longer of (peptide_seq, epitope_listed).
  - TCR distance:     CDR3-based TCRdist (BLOSUM62 substitution scores plus
                      length-difference gap penalty, weight 3) as in
                      Dash et al. 2017.  CDR3α + CDR3β paired distance.
                      No V-gene CDR1/2/2.5 contribution (CDR3-only mode,
                      matching the manuscript's threshold-of-120 usage).

Pipeline:
  1. Load TCR3d reference, filter to pre-cutoff + extraction_ok.
  2. ANARCI all reference TCRs to extract CDR3α, CDR3β.
  3. Load validation ST3, compute per-antigen and per-TCR AUCs.
  4. ANARCI all validation TCRs.
  5. Compute pMHC distances per validation antigen using sliding Hamming.
  6. Compute CDR3 TCRdist between every validation TCR and every reference
     TCR; take val→ref minimum.
  7. Save antigen_centric.csv, tcr_centric.csv, and a quick-look plot.
"""

import os
import sys
import warnings
from pathlib import Path
from datetime import datetime
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRATCH    = Path("/scratch/bneff/tcrtrifold")
REF_CSV    = SCRATCH / "tcr3d_triad_sequences.csv"
SUPPT3     = SCRATCH / "supplementaryTable3_20260513.xlsx"
OUTDIR     = Path("/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/"
                  "notebooks/leakage_v3_tcr3d")
OUTDIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Sliding-window Hamming (peptide + MHC)
# ---------------------------------------------------------------------------
def sliding_hamming_identity(a, b):
    if not isinstance(a, str) or not isinstance(b, str):
        return 0.0
    a = a.replace(".", "").replace("-", "")
    b = b.replace(".", "").replace("-", "")
    if not a or not b:
        return 0.0
    if len(a) < len(b):
        a, b = b, a
    Lb = len(b)
    best = 0.0
    for offset in range(len(a) - Lb + 1):
        matches = sum(a[offset + i] == b[i] for i in range(Lb))
        score = matches / Lb
        if score > best:
            best = score
            if best == 1.0:
                break
    return best


def longer(a, b):
    if not isinstance(a, str):
        return b
    if not isinstance(b, str):
        return a
    return a if len(a) >= len(b) else b


# ---------------------------------------------------------------------------
# ANARCI helpers — extract CDR3 + V/J gene calls for tcrdist3 input
# ---------------------------------------------------------------------------
def anarci_chain(seq):
    """Run ANARCI on one TCR chain. Returns dict (chain_type, cdr3,
    v_gene, j_gene, species) or None if numbering failed."""
    try:
        from anarci import anarci
        numbering, alignment_details, _ = anarci(
            [("q", seq)], scheme="imgt", allow=("A", "B"),
            allowed_species=("human", "mouse"),
        )
        if not alignment_details[0]:
            return None
        det = alignment_details[0][0]
        cdr3 = "".join(
            aa for (pos, aa) in numbering[0][0][0]
            if 105 <= pos[0] <= 117 and aa != "-"
        )
        germ = det.get("germlines", {}) or {}
        def _top(gname):
            x = germ.get(gname)
            if isinstance(x, list) and len(x) > 0 and isinstance(x[0], (list, tuple)):
                return x[0][1] if len(x[0]) >= 2 else None
            return None
        return {
            "chain_type": det["chain_type"],
            "cdr3":       cdr3,
            "v_gene":     _top("v_gene"),
            "j_gene":     _top("j_gene"),
            "species":    det.get("species"),
        }
    except Exception:
        return None


def annotate_tcr_pair(seq_a, seq_b):
    """ANARCI both chains, canonicalize α/β order. Returns dict or None."""
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
    return {
        "cdr3_a_aa": alpha["cdr3"], "v_a_gene": alpha["v_gene"], "j_a_gene": alpha["j_gene"],
        "cdr3_b_aa": beta["cdr3"],  "v_b_gene": beta["v_gene"],  "j_b_gene": beta["j_gene"],
        "species":   alpha["species"] or beta["species"] or "human",
    }


# ---------------------------------------------------------------------------
# Step 1 — load reference, ANARCI its TCRs
# ---------------------------------------------------------------------------
def load_reference():
    print("[1] Loading TCR3d reference …")
    ref = pd.read_csv(REF_CSV)
    pre = ref[(ref["extraction_status"] == "ok") &
              (ref["pre_af3_cutoff"].astype(bool))].copy()
    print(f"  total entries in CSV: {len(ref)}")
    print(f"  pre-cutoff + extraction ok: {len(pre)}")
    print(f"    by class: {pre['mhc_class'].value_counts().to_dict()}")

    pre["ref_pep"] = [longer(a, b) for a, b in
                      zip(pre["peptide_seq"], pre["epitope_listed"])]

    print(f"  ANARCIing {len(pre)} reference TCRs …")
    rows = []
    failed = 0
    for _, r in pre.iterrows():
        ann = annotate_tcr_pair(r["tcr_a_seq"], r["tcr_b_seq"])
        if ann is None:
            failed += 1
            ann = {"cdr3_a_aa": "", "v_a_gene": None, "j_a_gene": None,
                   "cdr3_b_aa": "", "v_b_gene": None, "j_b_gene": None,
                   "species": None}
        rec = r.to_dict()
        rec.update(ann)
        rows.append(rec)
    pre2 = pd.DataFrame(rows)
    print(f"    ANARCI: {len(pre2) - failed} ok, {failed} failed")
    return pre2


# ---------------------------------------------------------------------------
# Step 2 — validation, AUCs
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
# Step 3 — pMHC distance
# ---------------------------------------------------------------------------
def pmhc_distance_for_validation(ag_df, ref):
    out_dist = np.zeros(len(ag_df))
    ref_by_class = {
        "I":  ref[ref["mhc_class"] == "I"].to_dict("records"),
        "II": ref[ref["mhc_class"] == "II"].to_dict("records"),
    }
    print(f"\n[3] pMHC distances:  {len(ag_df)} val × "
          f"({len(ref_by_class['I'])} ref I, {len(ref_by_class['II'])} ref II)")

    for i, val in ag_df.iterrows():
        val_class = val["mhc_class"]
        refs = ref_by_class.get(val_class, [])
        val_pep  = val["peptide"]
        va_mhc_1 = val["mhc_1_seq"]
        va_mhc_2 = val["mhc_2_seq"]

        best_comp = 0.0
        for r in refs:
            pi = sliding_hamming_identity(val_pep, r["ref_pep"])
            if pi == 0.0 or pi <= best_comp:
                continue
            if val_class == "I":
                mi = sliding_hamming_identity(va_mhc_1, r["mhc_a_seq"])
            else:
                ai = sliding_hamming_identity(va_mhc_1, r["mhc_a_seq"])
                bi = sliding_hamming_identity(va_mhc_2, r["mhc_b_seq"])
                mi = min(ai, bi)
            comp = min(pi, mi)
            if comp > best_comp:
                best_comp = comp
                if best_comp >= 0.9999:
                    break
        out_dist[i] = 1.0 - best_comp
    return out_dist


# ---------------------------------------------------------------------------
# Step 4 — ANARCI validation TCRs
# ---------------------------------------------------------------------------
def annotate_validation_tcrs(tcr_df):
    print(f"\n[4] ANARCIing {len(tcr_df)} validation TCRs …")
    rows = []
    failed = 0
    for _, r in tcr_df.iterrows():
        ann = annotate_tcr_pair(r["tcr_1_seq"], r["tcr_2_seq"])
        if ann is None:
            failed += 1
            continue
        rows.append({
            "tcr_id":    r["tcr_id"],
            "mhc_class": r["mhc_class"],
            "auc":       r["auc"],
            **ann,
        })
    print(f"  success: {len(rows)},  failed: {failed}")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 5 — paired-chain CDR3 TCRdist (custom reimplementation matching
# tcrdist3's CDR3 component: BLOSUM62 substitution scoring + length-gap
# penalty, weight=3.  No V-gene CDR1/2/2.5 contribution.)
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


def cdr3_dist(a, b, gap_penalty=4, weight=3):
    """CDR3 distance, Dash et al. 2017 style:
        weight × ( gap_penalty × |Δlen| + Σ_i min(4, 4 - BLOSUM62(a[i], b[i])) )
    """
    la, lb = len(a), len(b)
    L = min(la, lb)
    d = gap_penalty * abs(la - lb)
    for i in range(L):
        bs = BLOSUM62.get((a[i], b[i]), -4)
        d += min(4, max(0, 4 - bs))
    return weight * d


def tcr_distances_via_tcrdist3(val_tcrs, ref):
    """For each validation TCR, return min(CDR3α dist + CDR3β dist) to any
    reference TCR. Uses the BLOSUM62-based scoring above (same as tcrdist3's
    CDR3 component, same as the manuscript's threshold-of-120 metric)."""
    print(f"\n[5] CDR3 TCRdist (custom, BLOSUM62 + length penalty, weight 3) …")
    ref_ok = ref[(ref["cdr3_a_aa"].fillna("") != "") &
                 (ref["cdr3_b_aa"].fillna("") != "")].copy().reset_index(drop=True)
    print(f"  usable reference TCRs:  {len(ref_ok)} / {len(ref)}")
    print(f"  usable validation TCRs: {len(val_tcrs)}")

    if len(ref_ok) == 0:
        return np.full(len(val_tcrs), np.nan)

    ref_a = ref_ok["cdr3_a_aa"].tolist()
    ref_b = ref_ok["cdr3_b_aa"].tolist()

    out = np.full(len(val_tcrs), np.nan)
    for i, r in val_tcrs.reset_index(drop=True).iterrows():
        va, vb = r["cdr3_a_aa"], r["cdr3_b_aa"]
        if not isinstance(va, str) or not isinstance(vb, str) or not va or not vb:
            continue
        best = np.inf
        for ra, rb in zip(ref_a, ref_b):
            d = cdr3_dist(va, ra) + cdr3_dist(vb, rb)
            if d < best:
                best = d
        out[i] = best if best != np.inf else np.nan
    return out


# ---------------------------------------------------------------------------
# Step 6 — plot
# ---------------------------------------------------------------------------
def plot_panels(ag_df, tcr_df, out_pdf, out_png):
    CLASS_COLOR = {"I": "#3a78b8", "II": "#c83737"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

    ax = axes[0]
    for cls in ["I","II"]:
        sub = ag_df[ag_df["mhc_class"]==cls].dropna(subset=["pmhc_distance","auc"])
        ax.scatter(sub["pmhc_distance"], sub["auc"],
                   c=CLASS_COLOR[cls], edgecolor="white", linewidth=0.6,
                   s=70, alpha=0.85, label=f"class {cls} (n={len(sub)})")
        if len(sub) >= 3 and sub["pmhc_distance"].nunique() >= 2:
            r,p = spearmanr(sub["pmhc_distance"], sub["auc"])
            ax.text(0.97, 0.05 if cls=="II" else 0.11,
                    f"class {cls}: r={r:+.2f}, p={p:.2f}",
                    transform=ax.transAxes, ha="right",
                    color=CLASS_COLOR[cls], fontsize=10, fontweight="bold")
    ax.axhline(0.5, ls="--", color="0.6", lw=0.7)
    ax.set_xlabel("pMHC distance (sliding Hamming composite)")
    ax.set_ylabel("Per-antigen AUC")
    ax.set_title("Antigen-centric (TCR3d reference)")
    ax.set_ylim(-0.05, 1.10)
    ax.legend(loc="lower left", fontsize=10, frameon=False)

    ax = axes[1]
    for cls in ["I","II"]:
        sub = tcr_df[tcr_df["mhc_class"]==cls].dropna(subset=["tcr_dist","auc"])
        ax.scatter(sub["tcr_dist"], sub["auc"],
                   c=CLASS_COLOR[cls],
                   edgecolor="white" if cls=="II" else "black",
                   linewidth=0.6, s=40 if cls=="II" else 70,
                   alpha=0.55 if cls=="II" else 0.85,
                   label=f"class {cls} (n={len(sub)})")
        if len(sub) >= 3:
            r,p = spearmanr(sub["tcr_dist"], sub["auc"])
            ax.text(0.97, 0.05 if cls=="II" else 0.11,
                    f"class {cls}: r={r:+.2f}, p={p:.2f}",
                    transform=ax.transAxes, ha="right",
                    color=CLASS_COLOR[cls], fontsize=10, fontweight="bold")
    ax.axhline(0.5, ls="--", color="0.6", lw=0.7)
    ax.set_xlabel("Paired-chain TCRdist (tcrdist3)")
    ax.set_ylabel("Per-TCR AUC")
    ax.set_title("TCR-centric (TCR3d reference)")
    ax.set_ylim(-0.05, 1.10)
    ax.legend(loc="lower left", fontsize=10, frameon=False)

    fig.suptitle("Leakage v3 — TCR3d reference, sliding Hamming pMHC + tcrdist3 TCR",
                 y=1.02, fontsize=12)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 72)
    print("Leakage analysis v3 — TCR3d reference, tcrdist3 TCR distance")
    print("=" * 72)

    ref = load_reference()

    val   = load_validation()
    ag_df = per_antigen_aucs(val)
    tcr_df = per_tcr_aucs(val)
    print(f"\n  per-antigen AUCs: {len(ag_df)}")
    print(f"  per-TCR AUCs:     {len(tcr_df)}")

    ag_df["pmhc_distance"] = pmhc_distance_for_validation(ag_df, ref)

    val_tcrs = annotate_validation_tcrs(tcr_df)
    val_tcrs["tcr_dist"] = tcr_distances_via_tcrdist3(val_tcrs, ref)

    ag_out  = OUTDIR / "antigen_centric.csv"
    tcr_out = OUTDIR / "tcr_centric.csv"
    ag_df.to_csv(ag_out, index=False)
    val_tcrs.to_csv(tcr_out, index=False)

    plot_panels(ag_df, val_tcrs,
                OUTDIR / "leakage_fig.pdf",
                OUTDIR / "leakage_fig.png")

    print("\n" + "=" * 72)
    print("Summary correlations (Spearman)")
    print("=" * 72)
    print("Antigen-centric:")
    for cls in ["I","II"]:
        sub = ag_df[ag_df["mhc_class"]==cls].dropna(subset=["pmhc_distance","auc"])
        if len(sub) >= 3 and sub["pmhc_distance"].nunique() >= 2:
            r,p = spearmanr(sub["pmhc_distance"], sub["auc"])
            print(f"  class {cls}: n={len(sub):>3}  r={r:+.3f}  p={p:.3f}")
        else:
            print(f"  class {cls}: n={len(sub):>3}  "
                  f"(x range {sub['pmhc_distance'].min():.3f}–"
                  f"{sub['pmhc_distance'].max():.3f})")
    print("TCR-centric:")
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
    print(f"  {OUTDIR/'leakage_fig.png'}")
    print(f"  {OUTDIR/'leakage_fig.pdf'}")


if __name__ == "__main__":
    main()
