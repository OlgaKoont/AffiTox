#!/bin/bash
#SBATCH --job-name=pb_all
#SBATCH -p aichem
#SBATCH --array=0-63%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pb-all_%A_%a.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pb-all_%A_%a.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
RESULTS="${ROOT}/results"
PB="${ROOT}/analysis/excluding_2z5x_3mjg/tables/posebuster"
BUST="/mnt/tank/scratch/okonovalova/miniconda3/envs/posebuster/bin/bust"

TARGETS=(
  1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65
  4tz4 4zau 5jkv 5mo4 6gqj 6jok 7awe 7kk3
)
METHODS=(gnina qvina dynamicbind boltz2)

task_id="${SLURM_ARRAY_TASK_ID}"
target_index=$((task_id / ${#METHODS[@]}))
method_index=$((task_id % ${#METHODS[@]}))
target="${TARGETS[$target_index]}"
method="${METHODS[$method_index]}"

if [[ "${method}" == "dynamicbind" ]]; then
  result_tag="dynamicbind_new"
else
  result_tag="${method}"
fi

source "/mnt/tank/scratch/okonovalova/miniconda3/etc/profile.d/conda.sh"
set +u
conda activate docking
set -u

PYTHON="$(command -v python)"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p "${PB}/results_by_method" "${PB}/tables_by_method" "${PB}/sdf"

result_name="posebusters_results_${target}_${result_tag}.csv"
result_path="${PB}/results_by_method/${result_name}"
parts_path="${PB}/results_by_method/.posebusters_results_${target}_${result_tag}_parts"

# This is a chemistry-corrected full rebuild. Old result chunks must not be
# resumed because they may have been produced from Open Babel-inferred bonds.
rm -f "${result_path}" "${PB}/${result_name}"
rm -rf "${parts_path}"

echo "array_job=${SLURM_ARRAY_JOB_ID} task=${task_id}"
echo "target=${target} method=${method} cpus=${SLURM_CPUS_PER_TASK:-16}"
echo "started=$(date -Iseconds)"

if [[ "${method}" == "boltz2" ]]; then
  "${PYTHON}" "${ROOT}/src/analysis/prepare_and_run_posebusters_boltz2.py" \
    --boltz-results-dir "${RESULTS}" \
    --output-dir "${PB}" \
    --proteins "${target}" \
    --bust-bin "${BUST}" \
    --max-workers "${SLURM_CPUS_PER_TASK:-16}" \
    --chunk-size 100 \
    --force
else
  "${PYTHON}" "${ROOT}/src/analysis/prepare_and_run_posebusters.py" \
    --poses-dir "${RESULTS}" \
    --proteins-dir "${ROOT}/processed/proteins" \
    --output-dir "${PB}" \
    --methods "${method}" \
    --proteins "${target}" \
    --dynamicbind-subdir dynamicbind_new \
    --dynamicbind-method-label dynamicbind_new \
    --bust-bin "${BUST}" \
    --max-workers "${SLURM_CPUS_PER_TASK:-16}" \
    --chunk-size 100 \
    --force
fi

if [[ ! -s "${result_path}" ]]; then
  echo "ERROR: expected PoseBusters result was not created: ${result_path}" >&2
  exit 1
fi

cp -f "${result_path}" "${PB}/${result_name}"

row_count=$(($(wc -l < "${result_path}") - 1))
echo "result=${PB}/${result_name}"
echo "rows=${row_count}"
echo "finished=$(date -Iseconds)"
