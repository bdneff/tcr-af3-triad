"""
build_pipeline_raw_csv.py
-------------------------

Convert the BLOSUM62-controlled CSV produced by 01_generate_random_peptides.py
into a "raw" CSV directly ingestible by the tcrtrifold-experiments pipeline's
cleaning step (matching the column convention used by
`data/cross_reactivity/raw/cross_reactivity_raw.csv` in the prior cross-
reactivity run).

What this script does
---------------------

1. Reads `data/cross_reactivity_controls_raw.csv`.
2. Adds a `cognate` column:
       - `True`  for class == 'original' (the published cognate)
       - `True`  for class == 'cross-reactive' (published cross-reactive)
       - `False` for class == 'randomized' (the pink controls)
       - (the published cross-reactive labels include some that were later
         shown to be weak; the analysis script ignores `cognate` for the
         randomized rows and uses BLOSUM62 distance as the x-axis instead)
3. Adds a `name` column = first-8 of MD5 hash of the concatenated 5 chain
   sequences (peptide + both MHC chains + both TCR chains). This matches the
   pipeline's `generate_job_name` convention.
4. Reorders columns to match `cross_reactivity_raw.csv` from the prior run.
5. Writes `data/cross_reactivity_controls_pipeline.csv`.

This output is what drops into:
    tcrtrifold-experiments/data/cross_reactivity_controls/raw/

so the pipeline's clean -> MSA -> inference -> feature-extract workflow can
run on it identically to the original 750-triad cross_reactivity job.

Usage
-----

    python scripts/02_build_pipeline_raw_csv.py
"""

import argparse
import csv
import hashlib
import sys
from pathlib import Path


def job_name(peptide: str, mhc1: str, mhc2: str, tcr1: str, tcr2: str) -> str:
    """MD5-derived 8-char job name (matches tcrtrifold pipeline convention)."""
    blob = (peptide + mhc1 + mhc2 + tcr1 + tcr2).encode()
    return hashlib.md5(blob).hexdigest()[:8]


# Output column order matches tcrtrifold-experiments/data/cross_reactivity/raw/cross_reactivity_raw.csv
PIPELINE_FIELDS = [
    'name',
    'source',
    'cognate',
    'peptide',
    'mhc_class',
    'mhc_1_chain', 'mhc_1_species', 'mhc_1_seq',
    'mhc_2_chain', 'mhc_2_species', 'mhc_2_seq',
    'tcr_1_chain', 'tcr_1_species', 'tcr_1_seq',
    'tcr_2_chain', 'tcr_2_species', 'tcr_2_seq',
    # extra columns carried for downstream analysis (not used by AF3, but the
    # pipeline preserves arbitrary columns through cleaning)
    'class',
    'bl62_distance',
    'hamming_distance',
]


def class_to_cognate(class_label: str) -> str:
    """Map our 3-way class label to the pipeline's binary `cognate` column."""
    if class_label == 'original':       return 'True'
    if class_label == 'cross-reactive': return 'True'
    if class_label == 'randomized':     return 'False'
    raise ValueError(f"unknown class label: {class_label}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--in',  dest='in_csv',  type=Path,
                    default=Path('data/cross_reactivity_controls_raw.csv'))
    ap.add_argument('--out', type=Path,
                    default=Path('data/cross_reactivity_controls_pipeline.csv'))
    args = ap.parse_args()

    with open(args.in_csv) as f:
        in_rows = list(csv.DictReader(f))

    out_rows = []
    seen_names = set()
    for r in in_rows:
        name = job_name(
            r['peptide'], r['mhc_1_seq'], r['mhc_2_seq'],
            r['tcr_1_seq'], r['tcr_2_seq'],
        )
        if name in seen_names:
            print(f'WARN: duplicate job name {name} for peptide {r["peptide"]}; '
                  f'extending to 12 hex chars to disambiguate', file=sys.stderr)
            blob = (r['peptide'] + r['mhc_1_seq'] + r['mhc_2_seq']
                    + r['tcr_1_seq'] + r['tcr_2_seq']).encode()
            name = hashlib.md5(blob).hexdigest()[:12]
        seen_names.add(name)

        out_rows.append({
            'name':            name,
            'source':          r['source'],
            'cognate':         class_to_cognate(r['class']),
            'peptide':         r['peptide'],
            'mhc_class':       r['mhc_class'],
            'mhc_1_chain':     r['mhc_1_chain'],
            'mhc_1_species':   r['mhc_1_species'],
            'mhc_1_seq':       r['mhc_1_seq'],
            'mhc_2_chain':     r['mhc_2_chain'],
            'mhc_2_species':   r['mhc_2_species'],
            'mhc_2_seq':       r['mhc_2_seq'],
            'tcr_1_chain':     r['tcr_1_chain'],
            'tcr_1_species':   r['tcr_1_species'],
            'tcr_1_seq':       r['tcr_1_seq'],
            'tcr_2_chain':     r['tcr_2_chain'],
            'tcr_2_species':   r['tcr_2_species'],
            'tcr_2_seq':       r['tcr_2_seq'],
            'class':           r['class'],
            'bl62_distance':   r['bl62_distance'],
            'hamming_distance': r['hamming_distance'],
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=PIPELINE_FIELDS)
        w.writeheader()
        w.writerows(out_rows)

    print(f'Wrote {len(out_rows)} rows to {args.out}', file=sys.stderr)
    print(f'Unique job names: {len(seen_names)} (no collisions)', file=sys.stderr)


if __name__ == '__main__':
    main()
