#!/bin/bash
#SBATCH --job-name=affitox_final
#SBATCH -p aichem
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-final-publication_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-final-publication_%j.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
RESULTS="${ROOT}/results"
OUT="${ROOT}/analysis/excluding_2z5x_3mjg"
TABLES="${OUT}/tables"
PB="${TABLES}/posebuster"
TARGETS_SPACE="1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 5lf3 7kk3"

source "/mnt/tank/scratch/okonovalova/miniconda3/etc/profile.d/conda.sh"
set +u
conda activate docking
set -u

export TOXAFFINITY_ROOT="${ROOT}"
export RESULTS_DIR="${RESULTS}"
export ANALYSIS_ROOT="${OUT}"
export DATA_DIR="${TABLES}"
export MERGED_DATA_DIR="${TABLES}"
export POSEBUSTERS_DIR="${PB}"
export POSEBUSTERS_BOLTZ2_DIR="${PB}"
export BOLTZ_RESULTS_DIR="${RESULTS}"
export DYNAMICBIND_POSEBUSTERS_LABEL="dynamicbind_new"
export TARGETS="${TARGETS_SPACE}"
export METHODS="boltz2 dynamicbind gnina plapt qvina"
export PYTHON="$(command -v python)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export COLOR_LOW="#3558C5"
export COLOR_MID="#FFF9EE"
export COLOR_HIGH="#EF4938"
export FIGURE_DPI=300
export FONT_FAMILY="DejaVu Sans"

echo "started=$(date -Iseconds)"
echo "analysis_root=${OUT}"

echo "=== Audit all regenerated PoseBusters outputs ==="
"${PYTHON}" "${ROOT}/src/analysis/consolidate_posebusters_results.py" \
  --root "${ROOT}" \
  --analysis-root "${OUT}"

"${PYTHON}" - <<'PY'
from pathlib import Path
import pandas as pd

path = Path("/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/tables/posebusters_input_output_audit.csv")
audit = pd.read_csv(path)
bad = audit.loc[audit["status"] != "complete"]
if len(bad):
    raise SystemExit("Incomplete PoseBusters combinations:\n" + bad.to_string(index=False))
print(f"PoseBusters audit complete: {len(audit)}/{len(audit)} target-method combinations")
PY

echo "=== Rebuild merged target tables from canonical results ==="
for target in ${TARGETS}; do
  "${PYTHON}" "${ROOT}/src/analysis/merge_ligands_docking_from_dir.py" \
    "${target}" \
    --base-dir "${ROOT}" \
    --output-dir "${TABLES}" \
    --results-dir "${RESULTS}" \
    --boltz-results-dir "${RESULTS}" \
    --dynamicbind-subdir dynamicbind_new
done

"${PYTHON}" "${ROOT}/src/analysis/add_pvalue_column.py" \
  --input-dir "${TABLES}" \
  --output-dir "${TABLES}"

echo "=== Full 10k/10k statistics and PDB-labelled figures ==="
export TARGET_LABEL_MODE=pdb
export FIGURE_VARIANT=by_pdb
bash "${ROOT}/pipeline/postprocess/run_si_full.sh"

echo "=== Protein-labelled figure set from the same final statistics ==="
export TARGET_LABEL_MODE=protein
export FIGURE_VARIANT=by_protein
export RUN_CORRELATIONS=0
export RUN_ENRICHMENT=0
export RUN_POSEBUSTERS=0
export RUN_INFERENTIAL=0
export RUN_FIGURES=1
bash "${ROOT}/pipeline/postprocess/run_article_analysis.sh"

echo "finished=$(date -Iseconds)"
echo "statistics=${TABLES}"
echo "pdb_figures=${OUT}/figures/by_pdb"
echo "protein_figures=${OUT}/figures/by_protein"
