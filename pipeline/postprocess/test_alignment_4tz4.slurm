#!/bin/bash
#SBATCH --job-name=merge4_audit
#SBATCH -p aichem
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/merge4_audit_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/merge4_audit_%j.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
OUT="${ROOT}/analysis/alignment_test_4tz4"
PY="/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python"
mkdir -p "${OUT}"

"${PY}" "${ROOT}/src/analysis/merge_ligands_docking_from_dir.py" 4tz4 \
  --base-dir "${ROOT}" \
  --output-dir "${OUT}" \
  --results-dir "${ROOT}/results" \
  --boltz-results-dir "${ROOT}/results" \
  --dynamicbind-subdir dynamicbind_new

"${PY}" - <<'PY'
from pathlib import Path

import pandas as pd

root = Path("/mnt/tank/scratch/okonovalova/AffiTox")
merged = pd.read_csv(
    root / "analysis/alignment_test_4tz4/merged_ligands_docking_4tz4.csv"
)
source = pd.read_csv(
    root
    / "results/4tz4/docking/dynamicbind_new/"
    "4tz4_CRBN_Ki_WT_ChEMBL_127_nodubl/affinity_prediction.csv"
)
actual = float(merged.iloc[0]["dynamicbind_affinity_bestpose"])
expected = float(source.loc[source["name"] == "idx_0", "affinity"].iloc[0])
assert actual == expected, (actual, expected)
print("CORRECT_4TZ4_SOURCE", actual, expected)
PY
