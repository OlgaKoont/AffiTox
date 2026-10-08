#!/bin/bash
#SBATCH -p aichem
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pb-%x_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pb-%x_%j.err

set -euo pipefail

: "${PB_TARGET:?PB_TARGET is required}"
: "${PB_METHODS:?PB_METHODS is required}"
: "${PB_RUN_BOLTZ2:?PB_RUN_BOLTZ2 is required}"

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
RESULTS="${ROOT}/results"
PB="${ROOT}/analysis/excluding_2z5x_3mjg/tables/posebuster"
BUST="/mnt/tank/scratch/okonovalova/miniconda3/envs/posebuster/bin/bust"

source "/mnt/tank/scratch/okonovalova/miniconda3/etc/profile.d/conda.sh"
set +u
conda activate docking
set -u

PYTHON="$(command -v python)"
METHODS_CSV="${PB_METHODS//:/,}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "target=${PB_TARGET} methods=${METHODS_CSV} boltz2=${PB_RUN_BOLTZ2}"
echo "cpus=${SLURM_CPUS_PER_TASK:-16} started=$(date -Iseconds)"

"${PYTHON}" "${ROOT}/src/analysis/prepare_and_run_posebusters.py" \
  --poses-dir "${RESULTS}" \
  --proteins-dir "${ROOT}/processed/proteins" \
  --output-dir "${PB}" \
  --methods "${METHODS_CSV}" \
  --proteins "${PB_TARGET}" \
  --dynamicbind-subdir dynamicbind_new \
  --dynamicbind-method-label dynamicbind_new \
  --bust-bin "${BUST}" \
  --max-workers "${SLURM_CPUS_PER_TASK:-16}"

if [[ "${PB_RUN_BOLTZ2}" == "1" ]]; then
  "${PYTHON}" "${ROOT}/src/analysis/prepare_and_run_posebusters_boltz2.py" \
    --boltz-results-dir "${RESULTS}" \
    --output-dir "${PB}" \
    --proteins "${PB_TARGET}" \
    --bust-bin "${BUST}" \
    --max-workers "${SLURM_CPUS_PER_TASK:-16}"
fi

if [[ -d "${PB}/results_by_method" ]]; then
  cp -f "${PB}/results_by_method/posebusters_results_${PB_TARGET}_"*.csv "${PB}/"
fi

echo "finished=$(date -Iseconds)"
