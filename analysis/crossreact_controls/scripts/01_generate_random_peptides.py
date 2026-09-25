"""
generate_random_peptides.py
---------------------------

Generates the "pink" randomized peptide controls for the cross-reactivity
analysis, paired with each class's fixed TCR + MHC.

Strategy
--------

For each class (I, II), we:

  1. Read I.csv / II.csv (which already contain the original cognate peptide
     plus the published cross-reactive peptides).
  2. Hold the TCR + MHC fixed (these are constant across all rows).
  3. Generate N_PER_CLASS random peptides with substitutions designed to span
     the BLOSUM62 distance x-axis approximately evenly, from d ~= 5 (single
     conservative substitution regime) up to d_max ~= 1.5 * expected-random
     distance for that peptide (covers the fully-randomized regime).

The random sampler is accept-reject:

  - At each draw, pick `n_subs` uniformly in [1, peptide_length].
  - Substitute that many random positions with random amino acids (different
     from the original at that position).
  - Compute BLOSUM62 distance of the resulting peptide.
  - Assign it to one of N_BINS equal-width bins covering [d_min, d_max].
  - Keep the peptide only if its bin still has room (PER_BIN slots).
  - Stop when all bins are filled, or after MAX_TRIES draws.

This gives a roughly uniform distribution along the x-axis without specifying
a substitution-count -> distance mapping by hand.

Output
------

A "raw" CSV at `data/cross_reactivity_controls_raw.csv` in the same column
format as `cross_reactivity_raw.csv` in the tcrtrifold-experiments repo,
containing:

  - 1 original triad per class                 (class = "original")
  - 21-23 cross-reactive triads per class      (class = "cross-reactive")
  - N_PER_CLASS randomized triads per class    (class = "randomized")

The pipeline's downstream cleaning script will assign job names (MD5 hashes)
and add CDR/numbering columns. This script does not duplicate that work.

Usage
-----

    python scripts/01_generate_random_peptides.py \\
        --i_csv  data/I.csv                       \\
        --ii_csv data/II.csv                      \\
        --out    data/cross_reactivity_controls_raw.csv \\
        --n_per_class 400                         \\
        --seed 42

Reproducibility
---------------

The random seed (default 42) controls peptide generation. Same seed + same
input CSVs => bit-identical output CSV.
"""

import argparse
import csv
import random
import sys
from pathlib import Path

# Make src/ importable when run as `python scripts/01_...`
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent / 'src'))

from blosum62 import (
    BLOSUM62, AAS,
    bl62_dist, self_score,
    expected_uniform_random_distance, max_possible_distance,
)
from io_utils import read_class_csv, extract_fixed_context, get_original_peptide


# ---------------------------------------------------------------------------
# Per-class HLA / chain metadata
# ---------------------------------------------------------------------------
# These reproduce the constant non-sequence columns we'd need in the raw CSV.
# They match the source studies for each class:
#   class I:  HLA-B*57:03 restricted, 11mer epitope KAFSPEVIPMF
#             (Wooldridge et al. / Birnbaum-style cross-reactivity)
#   class II: HLA-DRB3*03:01 restricted, 15mer epitope PPQIAANRSQLISLV
#             (Mason 1998-style cross-reactivity, 5OB tcr)
CLASS_METADATA = {
    'I': {
        'mhc_class':     'I',
        'mhc_1_chain':   'alpha',
        'mhc_1_species': 'human',
        'mhc_2_chain':   'beta',
        'mhc_2_species': 'human',
        'tcr_1_chain':   'alpha',
        'tcr_1_species': 'human',
        'tcr_2_chain':   'beta',
        'tcr_2_species': 'human',
        'source':        'cross_reactivity_controls_class_I',
    },
    'II': {
        'mhc_class':     'II',
        'mhc_1_chain':   'alpha',
        'mhc_1_species': 'human',
        'mhc_2_chain':   'beta',
        'mhc_2_species': 'human',
        'tcr_1_chain':   'alpha',
        'tcr_1_species': 'human',
        'tcr_2_chain':   'beta',
        'tcr_2_species': 'human',
        'source':        'cross_reactivity_controls_class_II',
    },
}


# ---------------------------------------------------------------------------
# Random peptide generator
# ---------------------------------------------------------------------------
def random_substituted_peptide(p_ref: str, n_subs: int, rng: random.Random) -> str:
    """Return a peptide with `n_subs` random positions substituted (different AA each)."""
    L = len(p_ref)
    n_subs = min(n_subs, L)
    positions = rng.sample(range(L), n_subs)
    p = list(p_ref)
    for pos in positions:
        orig = p[pos]
        p[pos] = rng.choice([a for a in AAS if a != orig])
    return ''.join(p)


def sample_randomized_peptides(
    p_ref: str,
    n_target: int,
    seed: int,
    *,
    d_min: int = 5,
    d_max: int | None = None,
    n_bins: int = 20,
    max_tries: int = 1_000_000,
) -> list[tuple[str, int]]:
    """Accept-reject sample peptides to roughly span [d_min, d_max] uniformly.

    Returns a list of (peptide, blosum62_distance) tuples. Length is <= n_target;
    will be exactly n_target if max_tries is not exhausted before all bins fill.

    Defaults:
        d_max = ceil(1.5 * expected-random distance), capped at theoretical max.
        n_bins chosen so per_bin = ceil(n_target / n_bins) is a small int.
    """
    rng = random.Random(seed)
    L = len(p_ref)

    if d_max is None:
        exp_rand = expected_uniform_random_distance(p_ref)
        d_max = min(int(round(1.5 * exp_rand)), max_possible_distance(p_ref))

    per_bin = max(1, (n_target + n_bins - 1) // n_bins)  # ceil(n_target / n_bins)
    bins: dict[int, list[tuple[str, int]]] = {i: [] for i in range(n_bins)}
    seen: set[str] = {p_ref}

    bin_width = (d_max - d_min) / n_bins

    tries = 0
    while tries < max_tries:
        tries += 1
        n_subs = rng.randint(1, L)
        p = random_substituted_peptide(p_ref, n_subs, rng)
        if p in seen:
            continue
        d = bl62_dist(p, p_ref)
        if not (d_min <= d < d_max):
            continue
        bin_idx = min(int((d - d_min) / bin_width), n_bins - 1)
        if len(bins[bin_idx]) >= per_bin:
            continue
        bins[bin_idx].append((p, d))
        seen.add(p)
        if all(len(v) >= per_bin for v in bins.values()):
            break

    out = [item for bin_list in bins.values() for item in bin_list]
    out.sort(key=lambda x: x[1])  # sort by distance for readability of output CSV
    return out


# ---------------------------------------------------------------------------
# Main: assemble the raw CSV
# ---------------------------------------------------------------------------
OUTPUT_FIELDS = [
    'peptide',
    'mhc_class',
    'mhc_1_chain', 'mhc_1_species', 'mhc_1_seq',
    'mhc_2_chain', 'mhc_2_species', 'mhc_2_seq',
    'tcr_1_chain', 'tcr_1_species', 'tcr_1_seq',
    'tcr_2_chain', 'tcr_2_species', 'tcr_2_seq',
    'source',
    'class',                # original | cross-reactive | randomized
    'bl62_distance',        # 0 for original; published value for cross-reactive; sampled for randomized
    'hamming_distance',     # carried through from I/II.csv for original + cross-reactive
]


def build_rows_for_class(
    class_label: str,                # 'I' or 'II'
    class_csv_path: Path,
    n_random: int,
    seed: int,
) -> list[dict]:
    """Build all triad rows for one HLA class: original + cross-reactive + randomized."""
    src_rows = read_class_csv(class_csv_path)
    fixed = extract_fixed_context(src_rows)
    p_ref = get_original_peptide(src_rows)
    meta = CLASS_METADATA[class_label]

    def make_row(peptide: str, klass: str, hamming: str | int, bl62: int) -> dict:
        return {
            'peptide':          peptide,
            'mhc_class':        meta['mhc_class'],
            'mhc_1_chain':      meta['mhc_1_chain'],
            'mhc_1_species':    meta['mhc_1_species'],
            'mhc_1_seq':        fixed['mhc_1_seq'],
            'mhc_2_chain':      meta['mhc_2_chain'],
            'mhc_2_species':    meta['mhc_2_species'],
            'mhc_2_seq':        fixed['mhc_2_seq'],
            'tcr_1_chain':      meta['tcr_1_chain'],
            'tcr_1_species':    meta['tcr_1_species'],
            'tcr_1_seq':        fixed['tcr_1_seq'],
            'tcr_2_chain':      meta['tcr_2_chain'],
            'tcr_2_species':    meta['tcr_2_species'],
            'tcr_2_seq':        fixed['tcr_2_seq'],
            'source':           meta['source'],
            'class':            klass,
            'bl62_distance':    bl62,
            'hamming_distance': hamming,
        }

    rows_out: list[dict] = []

    # 1) original
    rows_out.append(make_row(p_ref, 'original', 0, 0))

    # 2) cross-reactive (carry through from input CSV, compute bl62 distance fresh)
    for r in src_rows:
        if r['class'] != 'cross-reactive':
            continue
        d = bl62_dist(r['peptide'], p_ref)
        rows_out.append(make_row(r['peptide'], 'cross-reactive', r['hamming_distance'], d))

    # 3) randomized
    samples = sample_randomized_peptides(p_ref, n_random, seed=seed)
    for pep, d in samples:
        # hamming for reference, not used in plot
        hd = sum(1 for a, b in zip(pep, p_ref) if a != b)
        rows_out.append(make_row(pep, 'randomized', hd, d))

    return rows_out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--i_csv',  type=Path, default=Path('data/I.csv'))
    ap.add_argument('--ii_csv', type=Path, default=Path('data/II.csv'))
    ap.add_argument('--out',    type=Path, default=Path('data/cross_reactivity_controls_raw.csv'))
    ap.add_argument('--n_per_class', type=int, default=400,
                    help='Number of randomized peptides per class (default 400; 800 total across I+II)')
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    for label, path in [('I', args.i_csv), ('II', args.ii_csv)]:
        # Use a class-specific derived seed so changing n_per_class doesn't reshuffle both classes' randoms identically.
        class_seed = args.seed + (0 if label == 'I' else 1)
        rows = build_rows_for_class(label, path, args.n_per_class, class_seed)
        n_orig    = sum(1 for r in rows if r['class'] == 'original')
        n_cr      = sum(1 for r in rows if r['class'] == 'cross-reactive')
        n_random  = sum(1 for r in rows if r['class'] == 'randomized')
        print(f'  class {label}: {n_orig} original + {n_cr} cross-reactive + {n_random} randomized '
              f'= {len(rows)} triads', file=sys.stderr)
        all_rows.extend(rows)

    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS)
        w.writeheader()
        w.writerows(all_rows)

    print(f'Wrote {len(all_rows)} triads to {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
