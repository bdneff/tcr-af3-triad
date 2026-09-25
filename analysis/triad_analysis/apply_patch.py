"""Replace pmhc_distance + its helpers in leakage_v2_cluster_unstrat.py.

Backs up the original. Idempotent: bails out cleanly if already patched.
Run from anywhere — paths are absolute.
"""
import re, shutil
from datetime import datetime
from pathlib import Path

TARGET = Path("/tgen_labs/altin/alphafold3/workspace/tcrtrifold-experiments/notebooks/leakage_v2/leakage_v2_cluster_unstrat.py")
PATCH  = Path("/scratch/bneff/tcrtrifold/pmhc_distance_patch.py")
SENTINEL = "# === PATCHED: composite pMHC distance with BLAST MHC identity ==="

assert TARGET.exists(), f"missing: {TARGET}"
assert PATCH.exists(),  f"missing: {PATCH}"

src = TARGET.read_text()

if SENTINEL in src:
    print("Already patched — no changes made.")
    raise SystemExit(0)

# Find the start of the old function block. The function we're replacing
# is `def pmhc_distance(`. Everything from that line through the end of
# the function (next top-level `def ` at column 0) gets removed.
m = re.search(r"^def pmhc_distance\(", src, flags=re.MULTILINE)
assert m, "could not find `def pmhc_distance(` in target"
start = m.start()

# Find the next top-level def AFTER our function start
nxt = re.search(r"^def \w+\(", src[m.end():], flags=re.MULTILINE)
assert nxt, "could not find a following top-level def to bound the replacement"
end = m.end() + nxt.start()

# Compose the replacement: marker + patch body, sandwiched between the
# code before and after.
patch_body = PATCH.read_text()
new_src = src[:start] + SENTINEL + "\n" + patch_body + "\n\n" + src[end:]

# Backup + write
stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
backup = TARGET.with_suffix(f".py.bak.{stamp}")
shutil.copy(TARGET, backup)
TARGET.write_text(new_src)

# Validate: file must still parse as Python, must contain the sentinel,
# must still define pmhc_distance and main.
import ast
try:
    ast.parse(new_src)
except SyntaxError as e:
    shutil.copy(backup, TARGET)
    raise SystemExit(f"PATCH PRODUCED INVALID PYTHON, rolled back: {e}")

assert SENTINEL in new_src
assert "def pmhc_distance(" in new_src
assert "def main(" in new_src

print(f"OK. Backup at: {backup}")
print(f"Original size: {len(src):>7} chars")
print(f"Patched size:  {len(new_src):>7} chars")
print(f"Sentinel + pmhc_distance + main all present.")
