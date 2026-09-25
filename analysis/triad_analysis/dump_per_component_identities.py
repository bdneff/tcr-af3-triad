import sys
sys.path.insert(0, "/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/notebooks/leakage_v2")
import importlib, leakage_v2_cluster_unstrat
importlib.reload(leakage_v2_cluster_unstrat)
from leakage_v2_cluster_unstrat import (
    load_reference, load_validation, per_antigen_aucs,
    _split_class_ii_ref, _blast_identity,
)
import pandas as pd

ref = load_reference()
val = load_validation()
ag_df = per_antigen_aucs(val)
seq_lookup = (val.drop_duplicates(subset=["mhc_class","mhc_1_name","mhc_2_name","peptide"])
                 [["mhc_class","mhc_1_name","mhc_2_name","peptide","mhc_1_seq","mhc_2_seq"]])
ag_df = ag_df.merge(seq_lookup, on=["mhc_class","mhc_1_name","mhc_2_name","peptide"], how="left")

def pep_id_no_threshold(a, b):
    """Position-by-position identity over the longer peptide. No threshold."""
    if not isinstance(a, str) or not isinstance(b, str) or not a or not b:
        return 0.0
    L = max(len(a), len(b))
    return sum(x == y for x, y in zip(a, b)) / L if L else 0.0

rows = []
for _, val_row in ag_df.iterrows():
    vp = val_row["peptide"]
    val_class = val_row["mhc_class"]
    va = val_row["mhc_1_seq"]
    vb = val_row["mhc_2_seq"]

    best_pep_at_best_composite = 0.0
    best_mhc_at_best_composite = 0.0
    best_composite = 0.0

    for _, r in ref.iterrows():
        rp = r.get("pep_seq")
        if not isinstance(rp, str):
            continue
        rclass_int = r["mhc_class"]
        if val_class == "I" and rclass_int != 1:
            continue
        if val_class == "II" and rclass_int != 2:
            continue

        pep_i = pep_id_no_threshold(vp, rp)
        if pep_i == 0:
            continue

        if val_class == "I":
            mhc_i = _blast_identity(va, r["mhc_alignseq"])
        else:
            a, b = _split_class_ii_ref(r["mhc_alignseq"])
            ai = _blast_identity(va, a)
            bi = _blast_identity(vb, b)
            mhc_i = min(ai, bi)

        composite = min(pep_i, mhc_i)
        if composite > best_composite:
            best_composite = composite
            best_pep_at_best_composite = pep_i
            best_mhc_at_best_composite = mhc_i

    rows.append({
        "mhc_class": val_row["mhc_class"],
        "mhc_1_name": val_row["mhc_1_name"],
        "mhc_2_name": val_row["mhc_2_name"],
        "peptide": val_row["peptide"],
        "auc": val_row["auc"],
        "pep_identity": best_pep_at_best_composite,
        "mhc_identity": best_mhc_at_best_composite,
        "composite_identity": best_composite,
    })

out = pd.DataFrame(rows)
out.to_csv("/scratch/bneff/tcrtrifold/per_component_identities.csv", index=False)
print(out.to_string())
print(f"\nSaved to /scratch/bneff/tcrtrifold/per_component_identities.csv")
