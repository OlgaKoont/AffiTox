#!/bin/bash
#SBATCH --job-name=merge7_audit
#SBATCH -p aichem
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/merge7_audit_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/merge7_audit_%j.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
OUT="${ROOT}/analysis/alignment_test_7kk3"
PY="/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python"
mkdir -p "${OUT}"

"${PY}" "${ROOT}/src/analysis/merge_ligands_docking_from_dir.py" 7kk3 \
  --base-dir "${ROOT}" \
  --output-dir "${OUT}" \
  --results-dir "${ROOT}/results" \
  --boltz-results-dir "${ROOT}/results" \
  --dynamicbind-subdir dynamicbind_new

"${PY}" "${ROOT}/src/analysis/add_pvalue_column.py" \
  --input-dir "${OUT}" \
  --output-dir "${OUT}"

"${PY}" - <<'PY'
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr

path = Path(
    "/mnt/tank/scratch/okonovalova/AffiTox/analysis/alignment_test_7kk3/"
    "merged_ligands_docking_7kk3.csv"
)
data = pd.read_csv(path)
paired = data[["pValue", "dynamicbind_affinity_bestpose"]].dropna()
print(
    "CORRECTED_7KK3",
    len(paired),
    pearsonr(paired.iloc[:, 0], paired.iloc[:, 1]).statistic,
    spearmanr(paired.iloc[:, 0], paired.iloc[:, 1]).statistic,
)
PY
