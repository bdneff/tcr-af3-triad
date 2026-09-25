#!/usr/bin/env python3
"""Extract per-residue CDR phi/psi for ALL class I cognate validation triads."""
import sys, math
from pathlib import Path
import numpy as np
import pandas as pd

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.vectors import calc_dihedral
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner

try:
    from anarci import run_anarci
except ImportError:
    sys.exit("ERROR: ANARCI not installed.")
try:
    from Bio.PDB.Polypeptide import three_to_one as _t2o
    def three_to_one(n):
        try: return _t2o(n)
        except KeyError: return "X"
except ImportError:
    from Bio.PDB.Polypeptide import protein_letters_3to1
    def three_to_one(n): return protein_letters_3to1.get(n.upper(), "X")

BASE = Path("/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/data/pdb/triad")
META = BASE / "pdb_validation_triad.conf_af3.parquet"
PRED_DIR = BASE / "predictions"
CRYSTAL_DIR = BASE / "cleaned_pdb"
OUT_CSV = Path("./ramachandran_all_classI.csv")
MODEL_IDX = 0
CDR_RANGES = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr2_5": (81, 86), "cdr3": (104, 118)}


def chain_seq_residues(chain):
    seq, rs = [], []
    for r in chain:
        if is_aa(r, standard=True):
            seq.append(three_to_one(r.get_resname())); rs.append(r)
    return "".join(seq), rs

def compute_phi_psi(residues):
    n = len(residues)
    phi = np.full(n, np.nan); psi = np.full(n, np.nan)
    for i, r in enumerate(residues):
        try:
            if i > 0:
                phi[i] = math.degrees(calc_dihedral(
                    residues[i-1]["C"].get_vector(), r["N"].get_vector(),
                    r["CA"].get_vector(), r["C"].get_vector()))
        except Exception: pass
        try:
            if i < n-1:
                psi[i] = math.degrees(calc_dihedral(
                    r["N"].get_vector(), r["CA"].get_vector(),
                    r["C"].get_vector(), residues[i+1]["N"].get_vector()))
        except Exception: pass
    return phi, psi

def imgt_number(seq):
    for kwargs in [{}, {"allowed_species": ["mouse"]}, {"allowed_species": ["rat"]}]:
        try:
            res = run_anarci([("q", seq)], scheme="imgt", **kwargs)
            numbering = res[1][0]
            if numbering:
                domain, start, _ = numbering[0]
                nums = [None] * len(seq)
                si = start
                for (nt, aa) in domain:
                    if aa == "-": continue
                    if si < len(nums): nums[si] = nt[0]
                    si += 1
                if any(n is not None for n in nums): return nums
        except Exception: continue
    return [None] * len(seq)

def cdr_indices(nums):
    out = {k: [] for k in CDR_RANGES}
    for i, n in enumerate(nums):
        if n is None: continue
        for cdr, (lo, hi) in CDR_RANGES.items():
            if lo <= n <= hi: out[cdr].append(i)
    return out

def alignment_identity(s1, s2, aligner):
    try:
        aln = aligner.align(s1, s2)[0]
        sb, rb = aln.aligned
        m = 0
        for (a0, a1), (b0, b1) in zip(sb, rb):
            for k in range(a1 - a0):
                if s1[a0+k] == s2[b0+k]: m += 1
        return m / min(len(s1), len(s2))
    except Exception:
        return 0.0

def align_to_ref(s, ref, aligner):
    aln = aligner.align(s, ref)[0]
    sb, rb = aln.aligned
    m = {}
    for (s0, s1), (r0, r1) in zip(sb, rb):
        for k in range(s1 - s0): m[s0+k] = r0+k
    return m

def assign_tcr_chains(structure, t1, t2, aligner):
    cs = []
    for model in structure:
        for chain in model:
            seq, _ = chain_seq_residues(chain)
            if len(seq) >= 30: cs.append((chain, seq))
        break
    if len(cs) < 2: return None, None
    n = len(cs)
    s1 = np.array([alignment_identity(seq, t1, aligner) for _, seq in cs])
    s2 = np.array([alignment_identity(seq, t2, aligner) for _, seq in cs])
    best = (None, None, -1.0)
    for i in range(n):
        for j in range(n):
            if i == j: continue
            v = s1[i] + s2[j]
            if v > best[2]: best = (i, j, v)
    i, j, _ = best
    return cs[i][0], cs[j][0]


def main():
    meta = pd.read_parquet(META)
    cog = meta[meta["cognate"] & (meta["mhc_class"] == "I")].copy()
    print(f"[info] {len(cog)} class I cognate triads")

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2; aligner.mismatch_score = -1
    aligner.open_gap_score = -10; aligner.extend_gap_score = -1

    cif_p = MMCIFParser(QUIET=True)
    pdb_p = PDBParser(QUIET=True)
    all_rows = []

    for _, r in cog.iterrows():
        pdb_code = r["pdb"]; job = r["job_name"]
        pae = r["mean_p_tcr_interface_pae"]
        pred_f = PRED_DIR / job / f"{job}_model_{MODEL_IDX}.cif"
        cryst_f = CRYSTAL_DIR / f"{pdb_code}.pdb"
        if not pred_f.exists() or not cryst_f.exists():
            print(f"[skip] {pdb_code}: missing file"); continue
        try:
            pred = cif_p.get_structure(job, str(pred_f))
            cryst = pdb_p.get_structure(pdb_code, str(cryst_f))
        except Exception as e:
            print(f"[skip] {pdb_code}: {e}"); continue

        t1_seq, t2_seq = r["tcr_1_seq"], r["tcr_2_seq"]
        pred_t1, pred_t2 = assign_tcr_chains(pred, t1_seq, t2_seq, aligner)
        cryst_t1, cryst_t2 = assign_tcr_chains(cryst, t1_seq, t2_seq, aligner)

        n_res = 0
        for tcr_name, ref_seq, p_ch, c_ch in [("tcr_1", t1_seq, pred_t1, cryst_t1),
                                                ("tcr_2", t2_seq, pred_t2, cryst_t2)]:
            if p_ch is None or c_ch is None: continue
            nums = imgt_number(ref_seq)
            cdrs = cdr_indices(nums)
            p_seq, p_res = chain_seq_residues(p_ch)
            c_seq, c_res = chain_seq_residues(c_ch)
            p_phi, p_psi = compute_phi_psi(p_res)
            c_phi, c_psi = compute_phi_psi(c_res)
            p2r = align_to_ref(p_seq, ref_seq, aligner)
            c2r = align_to_ref(c_seq, ref_seq, aligner)
            r2p = {v:k for k,v in p2r.items()}; r2c = {v:k for k,v in c2r.items()}

            for cdr_name, ref_idxs in cdrs.items():
                for ri in ref_idxs:
                    if ri not in r2p or ri not in r2c: continue
                    pi, ci = r2p[ri], r2c[ri]
                    all_rows.append({
                        "pdb": pdb_code,
                        "peptide": r["peptide"],
                        "pti_pae": pae,
                        "tcr": tcr_name,
                        "cdr": cdr_name,
                        "imgt_num": nums[ri],
                        "aa": ref_seq[ri],
                        "phi_pred": p_phi[pi] if not np.isnan(p_phi[pi]) else None,
                        "psi_pred": p_psi[pi] if not np.isnan(p_psi[pi]) else None,
                        "phi_crystal": c_phi[ci] if not np.isnan(c_phi[ci]) else None,
                        "psi_crystal": c_psi[ci] if not np.isnan(c_psi[ci]) else None,
                    })
                    n_res += 1
        print(f"[done] {pdb_code} (PTI-PAE={pae:.2f}): {n_res} residues")

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\n[info] wrote {len(df)} residues across {df['pdb'].nunique()} triads -> {OUT_CSV}")


if __name__ == "__main__":
    main()
