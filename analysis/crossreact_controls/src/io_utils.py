"""
io_utils.py
-----------

Tiny helpers for reading the per-class cross-reactivity input CSVs (I.csv,
II.csv) and assembling the receptor/MHC context that all generated triads
share.

The input CSVs use mixed line endings ('\\r\\r\\n' in some exports), which
this module normalizes silently.
"""

import csv
from pathlib import Path
from typing import Dict, List


def read_class_csv(path: str | Path) -> List[Dict[str, str]]:
    """Read I.csv or II.csv, tolerating '\\r\\r\\n' line endings.

    Returns a list of dicts (one per row) with keys including:
        peptide, mhc_1_seq, mhc_2_seq, tcr_1_seq, tcr_2_seq,
        mean_p_tcr_interface_pae, hamming_distance, class
    """
    with open(path) as f:
        text = (
            f.read()
             .replace('\r\r\n', '\n')
             .replace('\r\n', '\n')
             .replace('\r', '\n')
        )
    return list(csv.DictReader(text.strip().split('\n')))


def extract_fixed_context(rows: List[Dict[str, str]]) -> Dict[str, str]:
    """Pull the (fixed across all rows) TCR + MHC sequences for a class.

    Each class file contains many peptides paired with one fixed TCR + MHC.
    This function asserts that assumption and returns the shared sequences.
    """
    fields = ['mhc_1_seq', 'mhc_2_seq', 'tcr_1_seq', 'tcr_2_seq']
    fixed = {}
    for key in fields:
        uniq = set(r[key] for r in rows)
        if len(uniq) != 1:
            raise ValueError(
                f"expected one unique value for {key} across rows, got {len(uniq)}"
            )
        fixed[key] = uniq.pop()
    return fixed


def get_original_peptide(rows: List[Dict[str, str]]) -> str:
    """Return the single 'original' (cognate, distance=0) peptide for the class."""
    originals = [r for r in rows if r.get('class') == 'original']
    if len(originals) != 1:
        raise ValueError(f"expected exactly 1 original peptide, found {len(originals)}")
    return originals[0]['peptide']
