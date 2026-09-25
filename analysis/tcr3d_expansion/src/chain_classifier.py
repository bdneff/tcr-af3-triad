"""
chain_classifier.py
-------------------
Identify peptide / MHC / TCR chains in a TCR-pMHC complex from sequences alone.

Designed to complement STCRDab's db_summary.dat: when a PDB entry is too new
for the STCRDab release we have on disk (a real problem with ~45 of the 49
post-AF3-cutoff TCR3d additions), we fall back to sequence-based classification
to produce the same fields the rest of the pipeline expects.

Outputs match the clean_pdb.py `outlier_pdb` row schema exactly:
    pdb, mhc_1_segid, mhc_2_segid, peptide_segid, tcr_1_segid, tcr_2_segid,
    mhc_1_species, mhc_2_species, tcr_1_species, tcr_2_species,
    mhc_class, mhc_1_chain, mhc_2_chain, cognate, tcr_1_chain, tcr_2_chain,
    pdb_date

Validation strategy: this module should reproduce the STCRDab annotations on
the 152 entries that exist in both TCR3d AND db_summary.dat. Use
`validate_against_stcrdab()` to spot-check before trusting any new annotations.
"""

from dataclasses import dataclass
from typing import Optional

# ---------------------------------------------------------------------------
# Motifs used for sequence-based chain identification
# ---------------------------------------------------------------------------

# β2-microglobulin: human "IQRTPKIQVYSRHPAENGKSNFLNCYV..." very characteristic.
# Mouse β2m also starts with IQK/IQR + similar. Match a conservative core.
B2M_MOTIFS = ("IQRTPKIQVY", "IQKTPQIQVY", "IQKTPQIQTY")
B2M_LENGTH_RANGE = (95, 105)

# TCR α constant region C-terminal motif. Human: "...LTVLDMR" / "...LTVLEDLR"
# Mouse: "...VTVLDMR" / similar. Conservative: look for last ~80 residues that
# contain the "constant region" signature "WSNKSDF" (TRAC) or "WTNKSDF" (mouse).
# Easier+more robust: TCR α constant contains "WSNKSDFACAN" or close. TCR β
# constant contains "FFPDHVELSWWVNGKEVHSGV" (human TRBC1/2) or
# "VFPDHVELSWWVNGKEVHSGV" (mouse). The 'WW' double-tryptophan is unique to
# TCR β constants among the components we expect here.
TCR_A_CONST_MOTIFS = (
    "WSNKSDFACAN",      # human TRAC core
    "WANRSNFSCDA",      # mouse TRAC core
    "WSNKSDFAC",        # human truncation
    "WAKKKDFACAN",      # rare variant; permissive fallback
)
TCR_B_CONST_MOTIFS = (
    "FFPDHVELSWWVNGK",  # human TRBC1/TRBC2 core
    "VFPDHVELSWWVNGK",  # mouse TRBC core
    "DHVELSWWVNGK",     # truncation
    "WWVNG",            # last-resort: double-W signature
)
TCR_LENGTH_RANGE = (180, 320)   # V+J+C region typical length

# Class I MHC heavy chain N-terminal motif. Human HLA-A/B/C: "GSHSMRY..."
# Mouse H2: "GPHSMRY..." / "GSHSLRY...". Length ~270-360 (mature, no signal pep).
MHC_I_HEAVY_MOTIFS = ("GSHSMRY", "GPHSMRY", "GSHSLRY", "GSHSLKY", "GSHSIKY", "GPHFLRY")
MHC_I_HEAVY_LENGTH_RANGE = (260, 380)

# Class II MHC chains. α and β are harder to distinguish by motif alone;
# both are ~180-220 aa, both start near "IKADHV" / "RDSPED" / similar.
# Reasonable heuristic: HLA-DR α starts with "IKEEHV" or "IKDHVHEFCV";
# HLA-DR β starts with "GDTRPRFL" / "RFLEYS" etc. We'll use a permissive
# motif library and fall back on docking-angle/length tie-breaking if needed.
# Importantly, the TCR3d row TELLS US the mhc_class, so we don't need to
# discriminate class I MHC from class II MHC by motif — only need to call
# α vs β within class II.
MHC_II_ALPHA_MOTIFS = (
    "IKEEHV",           # HLA-DRA, HLA-DPA, HLA-DQA region
    "EDIVADHV",
    "EEHVIIQA",
    "IKADHV",
    "EDFVYQF",          # HLA-DPA
    "EDFVFQF",
)
MHC_II_BETA_MOTIFS = (
    "GDTRPRFL",         # HLA-DRB common
    "RDSPEDFV",
    "RDSPEDFVY",
    "DTRPRFLE",
    "WNSQKDLL",         # HLA-DQB
    "RDSPEDFV",
    "DTRPRFL",
    "TERVRLVT",         # HLA-DPB
)
MHC_II_LENGTH_RANGE = (160, 230)

PEPTIDE_LENGTH_RANGE = (5, 30)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _contains_any(seq: str, motifs) -> bool:
    return any(m in seq for m in motifs)


@dataclass
class ChainScore:
    chain_id: str
    sequence: str
    length: int
    is_peptide_candidate: bool
    is_b2m_candidate: bool
    is_tcr_a_candidate: bool
    is_tcr_b_candidate: bool
    is_mhc_i_heavy_candidate: bool
    is_mhc_ii_alpha_candidate: bool
    is_mhc_ii_beta_candidate: bool
    # tie-breaking: number of distinct motifs hit (specificity)
    score_per_role: dict


def _score_chain(chain_id: str, seq: str) -> ChainScore:
    L = len(seq)
    seq_u = seq.upper()
    return ChainScore(
        chain_id=chain_id,
        sequence=seq,
        length=L,
        is_peptide_candidate=PEPTIDE_LENGTH_RANGE[0] <= L <= PEPTIDE_LENGTH_RANGE[1],
        is_b2m_candidate=(
            B2M_LENGTH_RANGE[0] <= L <= B2M_LENGTH_RANGE[1]
            and _contains_any(seq_u, B2M_MOTIFS)
        ),
        is_tcr_a_candidate=(
            TCR_LENGTH_RANGE[0] <= L <= TCR_LENGTH_RANGE[1]
            and _contains_any(seq_u, TCR_A_CONST_MOTIFS)
        ),
        is_tcr_b_candidate=(
            TCR_LENGTH_RANGE[0] <= L <= TCR_LENGTH_RANGE[1]
            and _contains_any(seq_u, TCR_B_CONST_MOTIFS)
        ),
        is_mhc_i_heavy_candidate=(
            MHC_I_HEAVY_LENGTH_RANGE[0] <= L <= MHC_I_HEAVY_LENGTH_RANGE[1]
            and _contains_any(seq_u, MHC_I_HEAVY_MOTIFS)
        ),
        is_mhc_ii_alpha_candidate=(
            MHC_II_LENGTH_RANGE[0] <= L <= MHC_II_LENGTH_RANGE[1]
            and _contains_any(seq_u, MHC_II_ALPHA_MOTIFS)
        ),
        is_mhc_ii_beta_candidate=(
            MHC_II_LENGTH_RANGE[0] <= L <= MHC_II_LENGTH_RANGE[1]
            and _contains_any(seq_u, MHC_II_BETA_MOTIFS)
        ),
        score_per_role={},
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

@dataclass
class ClassificationResult:
    pdb: str
    peptide_segid: Optional[str]
    mhc_1_segid: Optional[str]
    mhc_2_segid: Optional[str]
    tcr_1_segid: Optional[str]
    tcr_2_segid: Optional[str]
    mhc_class: str                # "I" or "II"
    confidence: str               # "high" | "medium" | "low"
    notes: list                   # human-readable reasoning / flags


def classify(pdb: str, chain_seqs: dict, expected_class: str) -> ClassificationResult:
    """
    Classify chains in a PDB into peptide / MHC_1 / MHC_2 / TCR_α / TCR_β.

    Parameters
    ----------
    pdb : str
        PDB ID (lowercase preferred; only used for diagnostics).
    chain_seqs : dict
        Mapping {chain_id: sequence} as returned by an RCSB FASTA parse.
        For PDBs with multiple identical biological assemblies, deduplicate
        BEFORE calling classify (e.g., take just the first copy of each
        unique sequence).
    expected_class : str
        "I" or "II". From TCR3d. Used to constrain the role set.

    Returns
    -------
    ClassificationResult with segid assignments and a confidence label.

    Confidence
    ----------
    high   — all expected roles assigned by single non-conflicting candidate
    medium — assignments made but with fallbacks (e.g., MHC II α/β resolved
             by length only, no motif hit) or after disambiguating two
             candidates by length
    low    — at least one role had ambiguous or no candidate; manual review
             needed
    """
    notes = []
    scored = {c: _score_chain(c, s) for c, s in chain_seqs.items()}

    def pick(candidate_attr, exclude=set()):
        """Find a single chain matching candidate_attr, excluding already-used."""
        hits = [
            sc for c, sc in scored.items()
            if c not in exclude and getattr(sc, candidate_attr)
        ]
        return hits

    assigned = {}
    confidence = "high"

    # ------------------------------------------------------------------
    # 1. β2m FIRST (only for class I) — it's the most specific motif
    # ------------------------------------------------------------------
    if expected_class == "I":
        b2m_hits = pick("is_b2m_candidate")
        if len(b2m_hits) == 1:
            assigned["mhc_2"] = b2m_hits[0].chain_id
        elif len(b2m_hits) > 1:
            # Some biological assemblies have multiple β2m copies — take the first
            assigned["mhc_2"] = b2m_hits[0].chain_id
            notes.append(f"multiple β2m candidates: {[h.chain_id for h in b2m_hits]}, took first")
        else:
            notes.append("class I but no β2m candidate found")
            confidence = "low"
            assigned["mhc_2"] = None

    # ------------------------------------------------------------------
    # 2. TCR β — double-W motif is highly specific
    # ------------------------------------------------------------------
    used = set(v for v in assigned.values() if v)
    tcr_b_hits = pick("is_tcr_b_candidate", exclude=used)
    if len(tcr_b_hits) == 1:
        assigned["tcr_2"] = tcr_b_hits[0].chain_id
    elif len(tcr_b_hits) > 1:
        assigned["tcr_2"] = tcr_b_hits[0].chain_id
        notes.append(f"multiple TCR β candidates: {[h.chain_id for h in tcr_b_hits]}, took first")
        confidence = "medium" if confidence == "high" else confidence
    else:
        notes.append("no TCR β candidate matched constant-region motif")
        confidence = "low"
        assigned["tcr_2"] = None

    # ------------------------------------------------------------------
    # 3. TCR α
    # ------------------------------------------------------------------
    used = set(v for v in assigned.values() if v)
    tcr_a_hits = pick("is_tcr_a_candidate", exclude=used)
    if len(tcr_a_hits) == 1:
        assigned["tcr_1"] = tcr_a_hits[0].chain_id
    elif len(tcr_a_hits) > 1:
        assigned["tcr_1"] = tcr_a_hits[0].chain_id
        notes.append(f"multiple TCR α candidates: {[h.chain_id for h in tcr_a_hits]}, took first")
        confidence = "medium" if confidence == "high" else confidence
    else:
        notes.append("no TCR α candidate matched constant-region motif")
        confidence = "low"
        assigned["tcr_1"] = None

    # ------------------------------------------------------------------
    # 4. MHC: class I heavy OR class II α + β
    # ------------------------------------------------------------------
    used = set(v for v in assigned.values() if v)

    if expected_class == "I":
        mhc_hits = pick("is_mhc_i_heavy_candidate", exclude=used)
        if len(mhc_hits) == 1:
            assigned["mhc_1"] = mhc_hits[0].chain_id
        elif len(mhc_hits) > 1:
            # Take longest (most complete) heavy chain
            assigned["mhc_1"] = max(mhc_hits, key=lambda h: h.length).chain_id
            notes.append(f"multiple MHC I heavy candidates, took longest: {assigned['mhc_1']}")
        else:
            # Fallback: longest remaining chain in the heavy-chain length range
            fallback = [
                sc for c, sc in scored.items()
                if c not in used and MHC_I_HEAVY_LENGTH_RANGE[0] <= sc.length <= MHC_I_HEAVY_LENGTH_RANGE[1]
            ]
            if fallback:
                assigned["mhc_1"] = max(fallback, key=lambda h: h.length).chain_id
                notes.append(f"MHC I heavy fallback by length only: {assigned['mhc_1']}")
                confidence = "medium" if confidence == "high" else confidence
            else:
                notes.append("no MHC I heavy candidate")
                confidence = "low"
                assigned["mhc_1"] = None

    else:  # class II
        a_hits = pick("is_mhc_ii_alpha_candidate", exclude=used)
        b_hits = pick("is_mhc_ii_beta_candidate", exclude=used)

        # Resolve α
        if len(a_hits) == 1:
            assigned["mhc_1"] = a_hits[0].chain_id
        elif len(a_hits) > 1:
            assigned["mhc_1"] = a_hits[0].chain_id
            notes.append(f"multiple MHC II α candidates, took first: {assigned['mhc_1']}")
            confidence = "medium" if confidence == "high" else confidence
        else:
            notes.append("no MHC II α matched by motif; fallback below")

        used = set(v for v in assigned.values() if v)

        # Resolve β
        if len(b_hits) == 1 and b_hits[0].chain_id not in used:
            assigned["mhc_2"] = b_hits[0].chain_id
        elif len([h for h in b_hits if h.chain_id not in used]) > 0:
            cands = [h for h in b_hits if h.chain_id not in used]
            assigned["mhc_2"] = cands[0].chain_id
            notes.append(f"multiple MHC II β candidates, took first: {assigned['mhc_2']}")
            confidence = "medium" if confidence == "high" else confidence
        else:
            notes.append("no MHC II β matched by motif; fallback below")

        # Fallback: if α or β not resolved, pick remaining chains in length range
        used = set(v for v in assigned.values() if v)
        unresolved = [r for r in ("mhc_1", "mhc_2") if r not in assigned or assigned[r] is None]
        if unresolved:
            remaining = [
                sc for c, sc in scored.items()
                if c not in used
                and MHC_II_LENGTH_RANGE[0] <= sc.length <= MHC_II_LENGTH_RANGE[1]
            ]
            for role in unresolved:
                if remaining:
                    pick_chain = remaining.pop(0)
                    assigned[role] = pick_chain.chain_id
                    notes.append(f"{role} resolved by length-only fallback: {pick_chain.chain_id}")
                    confidence = "low"
                else:
                    assigned[role] = None
                    confidence = "low"
                    notes.append(f"{role}: no candidate available")

    # ------------------------------------------------------------------
    # 5. Peptide — shortest unassigned chain in peptide length range
    # ------------------------------------------------------------------
    used = set(v for v in assigned.values() if v)
    pep_hits = [
        sc for c, sc in scored.items()
        if c not in used and sc.is_peptide_candidate
    ]
    if pep_hits:
        # shortest
        assigned["peptide"] = min(pep_hits, key=lambda h: h.length).chain_id
    else:
        notes.append("no peptide candidate in length range; check chain composition")
        confidence = "low"
        assigned["peptide"] = None

    return ClassificationResult(
        pdb=pdb.lower(),
        peptide_segid=assigned.get("peptide"),
        mhc_1_segid=assigned.get("mhc_1"),
        mhc_2_segid=assigned.get("mhc_2"),
        tcr_1_segid=assigned.get("tcr_1"),
        tcr_2_segid=assigned.get("tcr_2"),
        mhc_class=expected_class,
        confidence=confidence,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Validation: run classifier against known db_summary annotations
# ---------------------------------------------------------------------------

def validate_against_stcrdab(
    pdb: str,
    chain_seqs: dict,
    expected_class: str,
    known: dict,
) -> dict:
    """
    Run classifier on a PDB and compare to known STCRDab annotation.

    known : dict with keys peptide_segid, mhc_1_segid, mhc_2_segid,
            tcr_1_segid, tcr_2_segid (any value may be None for class I MHC β,
            i.e. β2m is mhc_2_segid for class I in STCRDab)

    Returns dict with predicted, expected, match (bool per role), and notes.
    """
    pred = classify(pdb, chain_seqs, expected_class)
    roles = ("peptide", "mhc_1", "mhc_2", "tcr_1", "tcr_2")
    pred_dict = {
        "peptide": pred.peptide_segid,
        "mhc_1": pred.mhc_1_segid,
        "mhc_2": pred.mhc_2_segid,
        "tcr_1": pred.tcr_1_segid,
        "tcr_2": pred.tcr_2_segid,
    }
    matches = {}
    for r in roles:
        seg_key = f"{r}_segid"
        e = known.get(seg_key)
        p = pred_dict.get(r)
        matches[r] = (e == p) or (e is None and p is None)
    return {
        "pdb": pdb.lower(),
        "predicted": pred_dict,
        "expected": {r: known.get(f"{r}_segid") for r in roles},
        "matches": matches,
        "all_match": all(matches.values()),
        "confidence": pred.confidence,
        "notes": pred.notes,
    }
