#!/bin/bash
#SBATCH --job-name=affitox_pb_prep
#SBATCH -p aichem
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pb-prep_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pb-prep_%j.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
RESULTS="${ROOT}/results"
OUT="${ROOT}/analysis/excluding_2z5x_3mjg"
TABLES="${OUT}/tables"
PB="${TABLES}/posebuster"
CANONICAL_PB="${ROOT}/analysis/tables/posebuster"
TARGETS="1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 7awe 7kk3"

source "/mnt/tank/scratch/okonovalova/miniconda3/etc/profile.d/conda.sh"
set +u
conda activate docking
set -u
PYTHON="$(command -v python)"

mkdir -p "${TABLES}" "${PB}"

"${PYTHON}" "${ROOT}/src/analysis/audit_result_alignment.py" \
  --root "${ROOT}" \
  --results-dir "${RESULTS}" \
  --targets "${TARGETS}" \
  --output-csv "${TABLES}/result_alignment_audit.csv" \
  --output-json "${TABLES}/result_alignment_audit_details.json"

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

# Keep the already validated PoseBusters outputs for unchanged targets.
cp -f "${CANONICAL_PB}"/posebusters_results_*.csv "${PB}/"

# These outputs will be generated afresh by the four worker jobs.
for target in 2v5z 4tz4 11ue; do
  rm -f "${PB}/posebusters_results_${target}_"*.csv
  rm -f "${PB}/results_by_method/posebusters_results_${target}_"*.csv
done
rm -f "${PB}/posebusters_results_7kk3_gnina.csv"
rm -f "${PB}/results_by_method/posebusters_results_7kk3_gnina.csv"
