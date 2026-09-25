"""
blosum62.py
-----------

BLOSUM62 substitution matrix (standard NCBI values) and peptide distance
function used as the x-axis for the cross-reactivity analysis.

Distance definition
-------------------

For two equal-length peptides p and p_ref:

    d(p, p_ref) = sum_i BLOSUM62(p_ref_i, p_ref_i) - sum_i BLOSUM62(p_i, p_ref_i)

so that d = 0 when p == p_ref, and d grows as substitutions accumulate
(weighted by their BLOSUM62 dissimilarity). This is the "shift to identity = 0"
convention specified by KL.

Notes
-----

- Only the 20 standard amino acids are supported. Non-standard residues
  (B, Z, X, U, O, *) raise KeyError if they appear in input peptides.
- All peptides compared with bl62_dist must be the same length as p_ref.
  Cross-reactivity peptides in this study are equal-length to the originals
  (class I: 11mer; class II: 15mer) by design.
"""

# Standard BLOSUM62 matrix
BLOSUM62 = {
    'A': {'A': 4, 'R':-1, 'N':-2, 'D':-2, 'C': 0, 'Q':-1, 'E':-1, 'G': 0,
          'H':-2, 'I':-1, 'L':-1, 'K':-1, 'M':-1, 'F':-2, 'P':-1, 'S': 1,
          'T': 0, 'W':-3, 'Y':-2, 'V': 0},
    'R': {'A':-1, 'R': 5, 'N': 0, 'D':-2, 'C':-3, 'Q': 1, 'E': 0, 'G':-2,
          'H': 0, 'I':-3, 'L':-2, 'K': 2, 'M':-1, 'F':-3, 'P':-2, 'S':-1,
          'T':-1, 'W':-3, 'Y':-2, 'V':-3},
    'N': {'A':-2, 'R': 0, 'N': 6, 'D': 1, 'C':-3, 'Q': 0, 'E': 0, 'G': 0,
          'H': 1, 'I':-3, 'L':-3, 'K': 0, 'M':-2, 'F':-3, 'P':-2, 'S': 1,
          'T': 0, 'W':-4, 'Y':-2, 'V':-3},
    'D': {'A':-2, 'R':-2, 'N': 1, 'D': 6, 'C':-3, 'Q': 0, 'E': 2, 'G':-1,
          'H':-1, 'I':-3, 'L':-4, 'K':-1, 'M':-3, 'F':-3, 'P':-1, 'S': 0,
          'T':-1, 'W':-4, 'Y':-3, 'V':-3},
    'C': {'A': 0, 'R':-3, 'N':-3, 'D':-3, 'C': 9, 'Q':-3, 'E':-4, 'G':-3,
          'H':-3, 'I':-1, 'L':-1, 'K':-3, 'M':-1, 'F':-2, 'P':-3, 'S':-1,
          'T':-1, 'W':-2, 'Y':-2, 'V':-1},
    'Q': {'A':-1, 'R': 1, 'N': 0, 'D': 0, 'C':-3, 'Q': 5, 'E': 2, 'G':-2,
          'H': 0, 'I':-3, 'L':-2, 'K': 1, 'M': 0, 'F':-3, 'P':-1, 'S': 0,
          'T':-1, 'W':-2, 'Y':-1, 'V':-2},
    'E': {'A':-1, 'R': 0, 'N': 0, 'D': 2, 'C':-4, 'Q': 2, 'E': 5, 'G':-2,
          'H': 0, 'I':-3, 'L':-3, 'K': 1, 'M':-2, 'F':-3, 'P':-1, 'S': 0,
          'T':-1, 'W':-3, 'Y':-2, 'V':-2},
    'G': {'A': 0, 'R':-2, 'N': 0, 'D':-1, 'C':-3, 'Q':-2, 'E':-2, 'G': 6,
          'H':-2, 'I':-4, 'L':-4, 'K':-2, 'M':-3, 'F':-3, 'P':-2, 'S': 0,
          'T':-2, 'W':-2, 'Y':-3, 'V':-3},
    'H': {'A':-2, 'R': 0, 'N': 1, 'D':-1, 'C':-3, 'Q': 0, 'E': 0, 'G':-2,
          'H': 8, 'I':-3, 'L':-3, 'K':-1, 'M':-2, 'F':-1, 'P':-2, 'S':-1,
          'T':-2, 'W':-2, 'Y': 2, 'V':-3},
    'I': {'A':-1, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-3, 'E':-3, 'G':-4,
          'H':-3, 'I': 4, 'L': 2, 'K':-3, 'M': 1, 'F': 0, 'P':-3, 'S':-2,
          'T':-1, 'W':-3, 'Y':-1, 'V': 3},
    'L': {'A':-1, 'R':-2, 'N':-3, 'D':-4, 'C':-1, 'Q':-2, 'E':-3, 'G':-4,
          'H':-3, 'I': 2, 'L': 4, 'K':-2, 'M': 2, 'F': 0, 'P':-3, 'S':-2,
          'T':-1, 'W':-2, 'Y':-1, 'V': 1},
    'K': {'A':-1, 'R': 2, 'N': 0, 'D':-1, 'C':-3, 'Q': 1, 'E': 1, 'G':-2,
          'H':-1, 'I':-3, 'L':-2, 'K': 5, 'M':-1, 'F':-3, 'P':-1, 'S': 0,
          'T':-1, 'W':-3, 'Y':-2, 'V':-2},
    'M': {'A':-1, 'R':-1, 'N':-2, 'D':-3, 'C':-1, 'Q': 0, 'E':-2, 'G':-3,
          'H':-2, 'I': 1, 'L': 2, 'K':-1, 'M': 5, 'F': 0, 'P':-2, 'S':-1,
          'T':-1, 'W':-1, 'Y':-1, 'V': 1},
    'F': {'A':-2, 'R':-3, 'N':-3, 'D':-3, 'C':-2, 'Q':-3, 'E':-3, 'G':-3,
          'H':-1, 'I': 0, 'L': 0, 'K':-3, 'M': 0, 'F': 6, 'P':-4, 'S':-2,
          'T':-2, 'W': 1, 'Y': 3, 'V':-1},
    'P': {'A':-1, 'R':-2, 'N':-2, 'D':-1, 'C':-3, 'Q':-1, 'E':-1, 'G':-2,
          'H':-2, 'I':-3, 'L':-3, 'K':-1, 'M':-2, 'F':-4, 'P': 7, 'S':-1,
          'T':-1, 'W':-4, 'Y':-3, 'V':-2},
    'S': {'A': 1, 'R':-1, 'N': 1, 'D': 0, 'C':-1, 'Q': 0, 'E': 0, 'G': 0,
          'H':-1, 'I':-2, 'L':-2, 'K': 0, 'M':-1, 'F':-2, 'P':-1, 'S': 4,
          'T': 1, 'W':-3, 'Y':-2, 'V':-2},
    'T': {'A': 0, 'R':-1, 'N': 0, 'D':-1, 'C':-1, 'Q':-1, 'E':-1, 'G':-2,
          'H':-2, 'I':-1, 'L':-1, 'K':-1, 'M':-1, 'F':-2, 'P':-1, 'S': 1,
          'T': 5, 'W':-2, 'Y':-2, 'V': 0},
    'W': {'A':-3, 'R':-3, 'N':-4, 'D':-4, 'C':-2, 'Q':-2, 'E':-3, 'G':-2,
          'H':-2, 'I':-3, 'L':-2, 'K':-3, 'M':-1, 'F': 1, 'P':-4, 'S':-3,
          'T':-2, 'W':11, 'Y': 2, 'V':-3},
    'Y': {'A':-2, 'R':-2, 'N':-2, 'D':-3, 'C':-2, 'Q':-1, 'E':-2, 'G':-3,
          'H': 2, 'I':-1, 'L':-1, 'K':-2, 'M':-1, 'F': 3, 'P':-3, 'S':-2,
          'T':-2, 'W': 2, 'Y': 7, 'V':-1},
    'V': {'A': 0, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-2, 'E':-2, 'G':-3,
          'H':-3, 'I': 3, 'L': 1, 'K':-2, 'M': 1, 'F':-1, 'P':-2, 'S':-2,
          'T': 0, 'W':-3, 'Y':-1, 'V': 4},
}

AAS = list(BLOSUM62.keys())  # 20 standard amino acids


def bl62_dist(p: str, p_ref: str) -> int:
    """BLOSUM62 distance: identity-shifted negative sum of BLOSUM62 scores.

    d(p, p_ref) = sum_i BLOSUM62(p_ref_i, p_ref_i) - sum_i BLOSUM62(p_i, p_ref_i)

    Identity (p == p_ref) gives d = 0. Higher = more dissimilar.
    """
    if len(p) != len(p_ref):
        raise ValueError(f"length mismatch: {len(p)} vs {len(p_ref)}")
    self_score = sum(BLOSUM62[a][a]      for a    in p_ref)
    obs_score  = sum(BLOSUM62[a][b] for a, b in zip(p, p_ref))
    return self_score - obs_score


def self_score(p_ref: str) -> int:
    """Sum of BLOSUM62 self-substitution scores for p_ref. Used as identity baseline."""
    return sum(BLOSUM62[a][a] for a in p_ref)


def max_possible_distance(p_ref: str) -> int:
    """Theoretical maximum BLOSUM62 distance for any peptide of equal length."""
    self_s = self_score(p_ref)
    min_obs = sum(min(BLOSUM62[aa][a] for aa in AAS) for a in p_ref)
    return self_s - min_obs


def expected_uniform_random_distance(p_ref: str) -> float:
    """Expected BLOSUM62 distance for a fully-random equal-length peptide (uniform AAs).

    This is the BLOSUM62 distance regime we expect a "non-cognate-by-design" peptide
    drawn uniformly from amino acid space to occupy.
    """
    self_s = self_score(p_ref)
    exp_obs = sum(sum(BLOSUM62[aa][a] for aa in AAS) / 20.0 for a in p_ref)
    return self_s - exp_obs
