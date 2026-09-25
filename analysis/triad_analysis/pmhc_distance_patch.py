"""
pmhc_distance_patch.py

Drop-in replacement for the pmhc_distance function in
leakage_v2_cluster_unstrat.py.

What changed:
  - MHC similarity is now BLAST percent identity on the chain sequences
    (matching the manuscript's exclusion-filter methodology), instead of
    a binary 0/1 on the field-2 allele name.

  - For class II, we BLAST both α and β chains separately and take the
    min of the two identities — mirroring the exclusion filter's
    "both chains ≥ 95%" intersection rule.

  - The composite distance keeps its intersection logic:
        pmhc_distance = 1 − max over refs of min(pep_identity, mhc_identity)

Schema assumptions (verified on Gemini 2026-05-13):

  Validation (ST3, class I):
    mhc_1_seq  → full HLA-A/B/C α-chain, ~275 res, leading "M" signal
    mhc_2_seq  → β2-microglobulin (irrelevant for similarity)

  Validation (ST3, class II):
    mhc_1_name / mhc_1_seq  → α chain (DRA1, DQA1), ~190 res
    mhc_2_name / mhc_2_seq  → β chain (DRB1, DRB5, DQB1), ~198 res

  Reference (ternary_templates_v2.tsv):
    mhc_class==1  → mhc_alignseq is α1+α2 only (~175 res, no signal, no α3)
    mhc_class==2  → mhc_alignseq is "α_seq/β_seq"
                    (single string, exactly one '/', halves ~75 and ~85 res)
                    mhc_allele is "α_allele,β_allele"

Local-mode BLAST is required for class I because ST3 chains run longer
than the reference (α3 domain is present in ST3, absent in reference).
For class II, sequence lengths are comparable on each chain.

USAGE (in leakage_v2_cluster_unstrat.py):

    # Replace the entire existing pmhc_distance function with the one below.
    # Everything else stays the same.

    def pmhc_distance(val_df, ref_df):
        ...
"""

import numpy as np
import Bio.Align


# Local-mode aligner with BLAST scoring. Local mode is important: ST3's
# class I chains have an α3 domain the reference doesn't, so a global
# alignment would force-align the unmatched tail and depress identity.
_aligner = Bio.Align.PairwiseAligner(scoring="blastp")
_aligner.mode = "local"


def _strip_gaps(s):
    """Remove '.' (MSA gaps from the reference) and '-' before scoring."""
    if not isinstance(s, str):
        return None
    return s.replace(".", "").replace("-", "")


def _blast_identity(a, b):
    """BLAST-style percent identity in [0, 1] under local alignment.

    Returns 0.0 if either input is missing/empty or alignment fails.
    Identity is matches / alignment_length over the aligned region only
    — i.e., the same definition the exclusion filter uses.
    """
    a = _strip_gaps(a)
    b = _strip_gaps(b)
    if not a or not b:
        return 0.0
    try:
        aln = sorted(_aligner.align(a, b))[0]
    except Exception:
        return 0.0
    a1, a2 = aln
    alnlen = sum(x != "-" for x in a1)
    if alnlen == 0:
        return 0.0
    return sum(x == y for x, y in zip(a1, a2)) / alnlen


def _split_class_ii_ref(ref_alignseq):
    """Class II reference alignseq is 'α_seq/β_seq'. Return (α, β) or (None, None)."""
    if not isinstance(ref_alignseq, str) or "/" not in ref_alignseq:
        return None, None
    parts = ref_alignseq.split("/")
    if len(parts) != 2:
        return None, None
    return parts[0], parts[1]


def pmhc_distance(val_df, ref_df):
    """
    Composite pMHC distance: 1 − max over ref pMHCs of min(pep_identity, mhc_identity).

    Both pep_identity and mhc_identity are continuous BLAST percent identities,
    consistent with the manuscript's exclusion filter.

    For class II, mhc_identity = min(α_identity, β_identity) on the two MHC
    chains separately, mirroring the exclusion filter's per-chain ≥95% rule.

    Inputs:
      val_df:  validation antigens — needs peptide, mhc_class,
               mhc_1_seq, mhc_2_seq
      ref_df:  reference pMHCs — needs pep_seq, mhc_class, mhc_alignseq
    """
    # Pre-process the reference once.  For class II we pre-split α/β here so
    # the inner loop doesn't repeat string work on every (val, ref) pair.
    ref_records = []
    for _, r in ref_df.dropna(subset=["pep_seq"]).iterrows():
        rec = {
            "pep_seq":  r.get("pep_seq"),
            "mhc_class": str(r.get("mhc_class", "")).strip(),
            "mhc_alignseq": r.get("mhc_alignseq"),
        }
        if rec["mhc_class"] in ("2", "II", "2.0"):
            a, b = _split_class_ii_ref(rec["mhc_alignseq"])
            rec["mhc_a"] = a
            rec["mhc_b"] = b
        else:
            rec["mhc_a"] = rec["mhc_alignseq"]
            rec["mhc_b"] = None
        ref_records.append(rec)

    out = []
    for _, val in val_df.iterrows():
        val_pep   = val.get("peptide")
        val_class = str(val.get("mhc_class", "")).strip()
        val_mhc_a = val.get("mhc_1_seq")  # class I α / class II α
        val_mhc_b = val.get("mhc_2_seq")  # class I β2m (ignore) / class II β

        best_composite = 0.0
        for r in ref_records:
            # Peptide BLAST identity first — cheapest filter.
            pep_i = _blast_identity(val_pep, r["pep_seq"])
            if pep_i == 0.0:
                continue
            # Short-circuit if pep_i can't possibly beat current best
            # (mhc_i is bounded by 1.0, so composite ≤ pep_i).
            if pep_i <= best_composite:
                continue

            # MHC identity. Cross-class never counts.
            if val_class in ("I", "1", "1.0"):
                if r["mhc_class"] not in ("1", "I", "1.0"):
                    continue
                mhc_i = _blast_identity(val_mhc_a, r["mhc_a"])
            elif val_class in ("II", "2", "2.0"):
                if r["mhc_class"] not in ("2", "II", "2.0"):
                    continue
                a_i = _blast_identity(val_mhc_a, r["mhc_a"])
                b_i = _blast_identity(val_mhc_b, r["mhc_b"])
                mhc_i = min(a_i, b_i)
            else:
                continue

            composite = min(pep_i, mhc_i)
            if composite > best_composite:
                best_composite = composite
                if best_composite >= 0.9999:
                    break
        out.append(1.0 - best_composite)
    return np.array(out)
