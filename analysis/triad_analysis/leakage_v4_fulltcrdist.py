"""
leakage_v4_tcr3d_fulltcrdist.py

Same as leakage_v3 except for the TCR side:

  v3:  CDR3-only TCRdist (BLOSUM62 + length penalty, weight 3 on CDR3,
       reimplemented locally — no V-gene info needed).

  v4:  Full TCRdist via the published tcrdist3 library, which uses
       CDR1+CDR2+CDR2.5+CDR3 weighted as in Dash et al. 2017.
       Requires V-gene names for every TCR, obtained via ANARCI's CLI
       germline-assignment step (a BLAST against the IMGT germline
       database).

The pMHC side is unchanged from v3.

Pipeline:
  1. Load TCR3d reference (308 pre-cutoff triads). V-gene names already
     present in the CSV (TRAV gene / TRBV gene from TCR3d's class tables).
  2. ANARCI all reference TCRs to extract CDR3α and CDR3β. Use the V-gene
     names from the TCR3d CSV directly (no germline lookup needed for ref).
  3. Load validation set (ST3), compute AUCs.
  4. ANARCI all validation TCRs to extract CDR3α and CDR3β, AND run
     ANARCI's CLI germline assignment to get V-gene names.
  5. Compute pMHC distances (same as v3).
  6. Build a tcrdist3 TCRrep with all reference + validation TCRs, compute
     paired-chain distances, read off val→ref minimums.
  7. Save outputs to notebooks/leakage_v4_fulltcrdist/

Notes:
  - ANARCI CLI is called as a subprocess. The germline-assignment step
    adds a BLAST per chain — for ~500 chains this is a few minutes.
  - tcrdist3 expects organism-specific germline databases. We run TCRrep
    once per organism (human and mouse separately) and merge results.
  - V-gene format must match tcrdist3's expectations. ANARCI returns
    names like "TRBV6-5*01"; tcrdist3 accepts these directly with its
    default `alphabeta_gammadelta_db.tsv` lookup. If a particular allele
    isn't in tcrdist3's DB we try the *01 fallback (e.g. "TRBV6-5*01").
"""

import os
import re
import sys
import subprocess
import tempfile
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
SCRATCH    = Path("/scratch/bneff/tcrtrifold")
REF_CSV    = SCRATCH / "tcr3d_triad_sequences.csv"
SUPPT3     = SCRATCH / "supplementaryTable3_20260513.xlsx"
OUTDIR     = Path("/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/"
                  "notebooks/leakage_v4_fulltcrdist")
OUTDIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Sliding-window Hamming (peptide + MHC) — unchanged from v3
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
# ANARCI CLI: extract CDR3 + V-gene from a chain
# ---------------------------------------------------------------------------
def anarci_cli_batch(seqs, scheme="imgt"):
    """Run ANARCI's CLI in batch mode with --assign_germline. Returns one dict
    per input sequence with keys: chain_type, cdr3, v_gene, j_gene, species.
    Failed entries return all-None.

    We use the CLI rather than the Python API because only the CLI exposes
    the germline assignment step (a BLAST against the IMGT germline database).
    """
    if not seqs:
        return []

    with tempfile.TemporaryDirectory() as tmpdir:
        fa_path = os.path.join(tmpdir, "in.fa")
        out_pref = os.path.join(tmpdir, "out")
        with open(fa_path, "w") as f:
            for i, s in enumerate(seqs):
                if not isinstance(s, str) or not s:
                    s = "X"  # placeholder for empty; will fail downstream
                f.write(f">q{i}\n{s}\n")

        # ANARCI CLI: -i in.fa -o out --scheme imgt --restrict TCR
        #             --assign_germline --csv
        # The --csv flag writes per-chain CSV files (one per chain type).
        cmd = [
            "ANARCI",
            "-i", fa_path,
            "-o", out_pref,
            "--scheme", scheme,
            "--restrict", "tr",
            "--assign_germline",
            "--csv",
        ]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if res.returncode != 0:
                print(f"ANARCI CLI returned {res.returncode}")
                print(f"  stderr (last 400 chars): {res.stderr[-400:]}")
        except Exception as e:
            print(f"ANARCI CLI subprocess error: {e}")
            return [{"chain_type": None, "cdr3": "", "v_gene": None,
                     "j_gene": None, "species": None}] * len(seqs)

        # Read output CSVs — ANARCI writes one per (chain_type, species) combo
        # named like out_<species>_<chain>.csv (e.g. out_A.csv or out_human_A.csv
        # depending on version). We glob the directory and merge.
        results = {i: None for i in range(len(seqs))}
        import glob
        for csv_path in glob.glob(os.path.join(tmpdir, "out*.csv")):
            df = pd.read_csv(csv_path)
            if "Id" not in df.columns or "chain_type" not in df.columns:
                continue
            for _, row in df.iterrows():
                # Id is the FASTA header (q0, q1, ...). Extract index.
                m = re.match(r"q(\d+)", str(row["Id"]))
                if not m:
                    continue
                idx = int(m.group(1))
                # CDR3: ANARCI CSV stores numbered residues as columns.
                # The CDR3 columns are at IMGT positions 105–117.
                cdr3 = ""
                for col in df.columns:
                    cm = re.match(r"^(\d+)([A-Z]?)$", col)
                    if not cm:
                        continue
                    pos = int(cm.group(1))
                    if 105 <= pos <= 117:
                        val = str(row[col]) if pd.notna(row[col]) else "-"
                        if val and val != "-":
                            cdr3 += val
                results[idx] = {
                    "chain_type": row.get("chain_type"),
                    "cdr3":       cdr3,
                    "v_gene":     row.get("v_gene"),
                    "j_gene":     row.get("j_gene"),
                    "species":    row.get("species"),
                }

        out = []
        for i in range(len(seqs)):
            if results[i] is None:
                out.append({"chain_type": None, "cdr3": "", "v_gene": None,
                            "j_gene": None, "species": None})
            else:
                out.append(results[i])
        return out


def annotate_tcr_pair_batch(seq_a_list, seq_b_list):
    """Run ANARCI on lists of α and β candidate chains. Return per-pair record
    with canonical (α, β) ordering. Mismatched chain types → None entry."""
    print(f"  ANARCI on {len(seq_a_list)} 'chain a' seqs …")
    a_results = anarci_cli_batch(seq_a_list)
    print(f"  ANARCI on {len(seq_b_list)} 'chain b' seqs …")
    b_results = anarci_cli_batch(seq_b_list)

    pairs = []
    fail_count = 0
    for ra, rb in zip(a_results, b_results):
        if ra["chain_type"] is None or rb["chain_type"] is None:
            pairs.append(None)
            fail_count += 1
            continue
        if ra["chain_type"] == "A" and rb["chain_type"] == "B":
            alpha, beta = ra, rb
        elif ra["chain_type"] == "B" and rb["chain_type"] == "A":
            alpha, beta = rb, ra
        else:
            pairs.append(None)
            fail_count += 1
            continue
        pairs.append({
            "cdr3_a_aa": alpha["cdr3"], "v_a_gene": alpha["v_gene"], "j_a_gene": alpha["j_gene"],
            "cdr3_b_aa": beta["cdr3"],  "v_b_gene": beta["v_gene"],  "j_b_gene": beta["j_gene"],
            "species":   alpha["species"] or beta["species"] or "human",
        })
    print(f"  {len(pairs) - fail_count} ok, {fail_count} failed")
    return pairs


# ---------------------------------------------------------------------------
# Step 1 — load reference, annotate
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

    print(f"  ANARCI + germline assignment on {len(pre)} reference TCR pairs …")
    pairs = annotate_tcr_pair_batch(pre["tcr_a_seq"].tolist(),
                                    pre["tcr_b_seq"].tolist())

    rows = []
    for i, p in enumerate(pairs):
        r = pre.iloc[i].to_dict()
        if p is None:
            r.update({"cdr3_a_aa":"", "v_a_gene":None, "j_a_gene":None,
                      "cdr3_b_aa":"", "v_b_gene":None, "j_b_gene":None,
                      "species":None, "tcr_ok":False})
        else:
            r.update(p)
            r["tcr_ok"] = True
        rows.append(r)
    return pd.DataFrame(rows)


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
# Step 3 — pMHC distance (unchanged from v3)
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
# Step 4 — annotate validation TCRs (CDR3 + V-gene)
# ---------------------------------------------------------------------------
def annotate_validation_tcrs(tcr_df):
    print(f"\n[4] ANARCI + germline assignment on {len(tcr_df)} validation TCR pairs …")
    pairs = annotate_tcr_pair_batch(tcr_df["tcr_1_seq"].tolist(),
                                    tcr_df["tcr_2_seq"].tolist())
    rows = []
    for i, p in enumerate(pairs):
        if p is None:
            continue
        rows.append({
            "tcr_id":    tcr_df["tcr_id"].iloc[i],
            "mhc_class": tcr_df["mhc_class"].iloc[i],
            "auc":       tcr_df["auc"].iloc[i],
            **p,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 5 — full TCRdist via tcrdist3
# ---------------------------------------------------------------------------
def normalize_v_gene(v):
    """Coerce a V-gene name into tcrdist3's expected format. tcrdist3 needs
    allele-resolved names like 'TRBV6-5*01'. ANARCI returns the gene plus
    most-likely allele; if just the gene comes back, append *01."""
    if not isinstance(v, str) or not v:
        return None
    v = v.strip()
    if "*" not in v:
        v = v + "*01"
    return v


def known_to_tcrdist3():
    """Load tcrdist3's germline DB once and return the set of recognized
    gene names (e.g. 'TRAV7-1*01'). We use this to drop any reference or
    validation TCR with V/J calls outside tcrdist3's lookup, rather than
    silently substituting alleles."""
    import os, tcrdist
    db_path = os.path.join(os.path.dirname(tcrdist.__file__),
                           "db", "alphabeta_gammadelta_db.tsv")
    g = pd.read_csv(db_path, sep="\t")
    return set(g["id"])


def filter_to_known_genes(df, known, label=""):
    """Drop rows whose V or J gene names aren't in tcrdist3's DB."""
    if not len(df):
        return df
    cols = ["v_a_gene","j_a_gene","v_b_gene","j_b_gene"]
    mask = pd.Series(True, index=df.index)
    for c in cols:
        mask &= df[c].isin(known)
    n_drop = (~mask).sum()
    if n_drop:
        bad = df.loc[~mask, cols].apply(
            lambda r: [v for v in r.values if v not in known and pd.notna(v)],
            axis=1)
        examples = {}
        for lst in bad:
            for x in lst:
                examples[x] = examples.get(x, 0) + 1
        print(f"  {label} dropping {n_drop} rows with V/J not in tcrdist3 db:")
        for x, n in sorted(examples.items(), key=lambda kv: -kv[1]):
            print(f"    {x!r}: {n} row(s)")
    return df[mask].copy()


def tcr_distances_full(val_tcrs, ref):
    """Compute full paired-chain TCRdist (all CDRs, Dash et al. weights) via
    tcrdist3. For each validation TCR, return the min distance to any
    reference TCR (within same organism)."""
    print(f"\n[5] Full TCRdist via tcrdist3 …")
    try:
        from tcrdist.repertoire import TCRrep
    except ImportError as e:
        print(f"  ERROR: could not import tcrdist3: {e}")
        print("  Install with:  pip install tcrdist3")
        return np.full(len(val_tcrs), np.nan)

    # Normalize V-gene names for both sides
    val = val_tcrs.copy()
    val["v_a_gene"] = val["v_a_gene"].map(normalize_v_gene)
    val["v_b_gene"] = val["v_b_gene"].map(normalize_v_gene)

    ref_ok = ref[ref["tcr_ok"]].copy()
    ref_ok["v_a_gene"] = ref_ok["v_a_gene"].map(normalize_v_gene)
    ref_ok["v_b_gene"] = ref_ok["v_b_gene"].map(normalize_v_gene)
    ref_ok = ref_ok[(ref_ok["cdr3_a_aa"] != "") & (ref_ok["cdr3_b_aa"] != "") &
                    ref_ok["v_a_gene"].notna() & ref_ok["v_b_gene"].notna() &
                    ref_ok["j_a_gene"].notna() & ref_ok["j_b_gene"].notna()]

    val = val[(val["cdr3_a_aa"] != "") & (val["cdr3_b_aa"] != "") &
              val["v_a_gene"].notna() & val["v_b_gene"].notna() &
              val["j_a_gene"].notna() & val["j_b_gene"].notna()]

    known = known_to_tcrdist3()
    ref_ok = filter_to_known_genes(ref_ok, known, label="reference:")
    val    = filter_to_known_genes(val,    known, label="validation:")

    print(f"  usable reference TCRs:  {len(ref_ok)} / {len(ref)}")
    print(f"  usable validation TCRs: {len(val)} / {len(val_tcrs)}")
    print(f"  ref species mix:  {ref_ok['species'].value_counts().to_dict()}")
    print(f"  val species mix:  {val['species'].value_counts().to_dict()}")

    # Map results back to the original val_tcrs row order
    out = np.full(len(val_tcrs), np.nan)
    val_idx_lookup = {tid: i for i, tid in enumerate(val_tcrs["tcr_id"])}

    for species in sorted(set(val["species"].dropna())):
        v_sub = val[val["species"] == species].reset_index(drop=True)
        r_sub = ref_ok[ref_ok["species"] == species].reset_index(drop=True)
        if len(r_sub) == 0:
            print(f"  WARNING: no {species} reference TCRs — skipping {len(v_sub)} {species} val TCRs")
            continue

        cols = ["cdr3_a_aa","v_a_gene","j_a_gene","cdr3_b_aa","v_b_gene","j_b_gene"]
        cell_df = pd.concat([v_sub[cols], r_sub[cols]], ignore_index=True)
        cell_df["count"]    = 1
        cell_df["clone_id"] = np.arange(len(cell_df))

        try:
            tr = TCRrep(cell_df=cell_df,
                        organism=species,
                        chains=["alpha","beta"],
                        db_file="alphabeta_gammadelta_db.tsv",
                        compute_distances=True,
                        deduplicate=False)
            paired = tr.pw_alpha + tr.pw_beta
        except Exception as e:
            print(f"  ERROR running TCRrep on {species}: {type(e).__name__}: {e}")
            continue

        n_val = len(v_sub)
        val_to_ref = paired[:n_val, n_val:]
        min_dists  = val_to_ref.min(axis=1)
        for local_i, tid in enumerate(v_sub["tcr_id"].values):
            orig_i = val_idx_lookup.get(tid)
            if orig_i is not None:
                out[orig_i] = float(min_dists[local_i])
        print(f"  {species}: computed {val_to_ref.shape[0]} × {val_to_ref.shape[1]} block")

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
    ax.set_xlabel("Full TCRdist (all CDRs, tcrdist3)")
    ax.set_ylabel("Per-TCR AUC")
    ax.set_title("TCR-centric (TCR3d reference)")
    ax.set_ylim(-0.05, 1.10)
    ax.legend(loc="lower left", fontsize=10, frameon=False)

    fig.suptitle("Leakage v4 — TCR3d reference, sliding Hamming pMHC + full TCRdist",
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
    print("Leakage analysis v4 — TCR3d reference, full TCRdist via tcrdist3")
    print("=" * 72)

    ref = load_reference()

    val   = load_validation()
    ag_df = per_antigen_aucs(val)
    tcr_df = per_tcr_aucs(val)
    print(f"\n  per-antigen AUCs: {len(ag_df)}")
    print(f"  per-TCR AUCs:     {len(tcr_df)}")

    ag_df["pmhc_distance"] = pmhc_distance_for_validation(ag_df, ref)

    val_tcrs = annotate_validation_tcrs(tcr_df)
    val_tcrs["tcr_dist"] = tcr_distances_full(val_tcrs, ref)

    ag_out  = OUTDIR / "antigen_centric.csv"
    tcr_out = OUTDIR / "tcr_centric.csv"
    ref_out = OUTDIR / "reference_annotated.csv"
    ag_df.to_csv(ag_out, index=False)
    val_tcrs.to_csv(tcr_out, index=False)
    ref[['pdb_id','mhc_class','tcr_ok','cdr3_a_aa','v_a_gene','j_a_gene','cdr3_b_aa','v_b_gene','j_b_gene','species']].to_csv(ref_out, index=False)

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
            print(f"  class {cls}: n={len(sub):>3}  (low variance or low n)")
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
