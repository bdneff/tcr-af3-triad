"""
rcsb_fetch.py
-------------
Small wrappers around RCSB endpoints. Designed to run on a machine with
network access to data.rcsb.org and www.rcsb.org.

- get_fasta(pdb_id): {chain_id: sequence}
- get_release_date(pdb_id): ISO datetime string

Both functions cache to a JSON file to avoid repeated network calls. The
cache directory should be tracked alongside the package so reruns are
deterministic across machines.
"""

import json
from io import StringIO
from pathlib import Path
from typing import Dict, Optional

import requests
from Bio import SeqIO


def _parse_chain_token(chain_token: str):
    """Parse 'Chain A' or 'Chains A, B' or 'Chain B[auth K]' from a FASTA description."""
    def _strip_one(s: str) -> str:
        s = s.strip()
        if "[" in s:
            return s.split("[auth ")[1].split("]")[0].strip()
        return s.replace(" ", "")

    if chain_token.startswith("Chain "):
        return [_strip_one(chain_token.split("Chain ")[1])]
    elif chain_token.startswith("Chains "):
        return [_strip_one(c) for c in chain_token.split("Chains ")[1].split(",")]
    else:
        return [_strip_one(chain_token)]


def _cache_path(cache_dir: Optional[Path], kind: str, pdb_id: str) -> Optional[Path]:
    if cache_dir is None:
        return None
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"{pdb_id.lower()}.{kind}.json"


def get_fasta(pdb_id: str, cache_dir: Optional[Path] = None) -> Dict[str, str]:
    """
    Return a dict mapping {chain_id: sequence} for every chain in pdb_id's
    deposited FASTA. Cached on disk.
    """
    pdb_id = pdb_id.lower()
    cp = _cache_path(cache_dir, "fasta", pdb_id)
    if cp and cp.exists():
        return json.loads(cp.read_text())

    r = requests.get(f"https://www.rcsb.org/fasta/entry/{pdb_id}", timeout=30)
    r.raise_for_status()

    seq_dict: Dict[str, str] = {}
    for rec in SeqIO.parse(StringIO(r.text), "fasta"):
        # description looks like e.g. "6L9L_1|Chains A|MHC class I antigen|..."
        try:
            chain_token = rec.description.split("|")[1]
            chains = _parse_chain_token(chain_token)
        except (IndexError, ValueError):
            continue
        for c in chains:
            seq_dict[c] = str(rec.seq)

    if cp:
        cp.write_text(json.dumps(seq_dict, indent=2))
    return seq_dict


def get_release_date(pdb_id: str, cache_dir: Optional[Path] = None) -> Optional[str]:
    """Return ISO string for initial_release_date, or None on failure."""
    pdb_id = pdb_id.lower()
    cp = _cache_path(cache_dir, "release", pdb_id)
    if cp and cp.exists():
        return json.loads(cp.read_text()).get("release_date")

    r = requests.get(
        f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=30
    )
    r.raise_for_status()
    d = r.json().get("rcsb_accession_info", {}).get("initial_release_date")
    if cp:
        cp.write_text(json.dumps({"release_date": d}, indent=2))
    return d


def get_entity_organism(pdb_id: str, cache_dir: Optional[Path] = None) -> Dict[str, str]:
    """
    Return a dict mapping {chain_id: organism} for polymer entities.
    Organism is 'human', 'mouse', or None if unknown.

    Uses the polymer_entity endpoint plus the entity-to-chain mapping
    in the entry payload.
    """
    pdb_id = pdb_id.lower()
    cp = _cache_path(cache_dir, "organism", pdb_id)
    if cp and cp.exists():
        return json.loads(cp.read_text())

    # First, get the list of polymer entity IDs and their chain assignments
    r = requests.get(
        f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=30
    )
    r.raise_for_status()
    entry = r.json()
    polymer_entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "polymer_entity_ids", []
    )

    result: Dict[str, str] = {}
    for ent_id in polymer_entity_ids:
        r = requests.get(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{ent_id}",
            timeout=30,
        )
        if r.status_code != 200:
            continue
        ent = r.json()
        chains = ent.get("rcsb_polymer_entity_container_identifiers", {}).get(
            "auth_asym_ids", []
        )
        organisms = []
        for src in ent.get("rcsb_entity_source_organism", []) or []:
            sci = (src.get("scientific_name") or "").lower()
            if "homo sapiens" in sci:
                organisms.append("human")
            elif "mus musculus" in sci:
                organisms.append("mouse")
        org = organisms[0] if organisms else None
        for c in chains:
            result[c] = org

    if cp:
        cp.write_text(json.dumps(result, indent=2))
    return result
