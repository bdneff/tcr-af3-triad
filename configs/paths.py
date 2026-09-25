"""configs/paths.py — the ONLY file in this repo holding absolute cluster paths.

Code lives in git; bulk data, model weights, and run outputs live on the cluster.
Everything that is not tracked is reachable from here, so porting the repo to a new
machine means editing this file and nothing else.

    from configs.paths import TRIAD_PARQUETS, AF3_SIF
    df = pd.read_parquet(TRIAD_PARQUETS / "pdb_triad.af3_rmsd.parquet")

Anything that resolves to None is simply absent on this host; callers should say so
rather than failing with a confusing FileNotFoundError deeper in.
"""
from pathlib import Path

# --- this repo --------------------------------------------------------------
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
ANALYSIS = REPO / "analysis"
FIGURES = REPO / "figures"
STRUCTURES = REPO / "structures"
MD = REPO / "md"

# --- cluster (gemini) -------------------------------------------------------
# $WS is the durable lab workspace; /scratch is EPHEMERAL (swept ~30 days), so a
# checkout there is disposable but anything long-lived must be copied to $WS.
WS = Path("/tgen_labs/altin/alphafold3/workspace")
SCRATCH = Path("/scratch/bneff")

# The upstream AF3 + Boltz-2 Nextflow pipeline that produced the manuscript's
# predictions: github.com/AltinLab/tcrtrifold-experiments
TRIFOLD = WS / "tcrtrifold-experiments"
TRIAD_PARQUETS = TRIFOLD / "data/pdb/triad/staged"

# --- AlphaFold 3 (lab container) --------------------------------------------
AF3_SIF = Path("/tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif")
AF3_MODELS = Path("/ref_genomes/alphafold/alphafold3/models")
AF3_DBS = Path("/ref_genomes/alphafold/alphafold3/")
AF3_BINDS = "/home,/scratch,/tgen_labs,/ref_genomes"

# --- AF-TCRdock baseline ----------------------------------------------------
# Built by analysis/tcrdock_validation/setup_env.sh; see docs/TCRDOCK.md. Pinned to
# CUDA 11.8 / cuDNN 8.9 because TCRdock runs an AlphaFold-2.3.2-era JAX.
TCRDOCK_WORK = SCRATCH / "tcrdock_work"
TCRDOCK_ENV = TCRDOCK_WORK / "env"
TCRDOCK_REPO = TCRDOCK_WORK / "TCRdock"
TCRDOCK_PARAMS = TCRDOCK_WORK / "alphafold_params"
# The fine-tuned TCR:pMHC weights -- the irreplaceable file if the Dropbox links rot.
TCRDOCK_FT_WEIGHTS = TCRDOCK_PARAMS / "params/tcrpmhc_run4_af_mhc_params_891.pkl"

# --- IEDB -------------------------------------------------------------------
# 13.5 GB on the cluster; a 90 MB pruned copy is at data/iedb_public.db (gitignored).
IEDB_DB = Path("/tgen_labs/altin/alphafold3/IEDB/2025-04-15/iedb_public.db")
IEDB_LOCAL = DATA / "iedb_public.db"



def first_existing(*candidates):
    """Return the first path that exists, else None. Use for cluster-or-local data."""
    for c in candidates:
        if c is not None and Path(c).exists():
            return Path(c)
    return None


def iedb() -> Path | None:
    """Prefer the full cluster DB, fall back to the local pruned copy."""
    return first_existing(IEDB_DB, IEDB_LOCAL)
