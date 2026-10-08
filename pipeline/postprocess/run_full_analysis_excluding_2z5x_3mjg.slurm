#!/bin/bash
#SBATCH --job-name=affitox_full
#SBATCH -p aichem
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-full-analysis_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-full-analysis_%j.err

# Full documented AffiTox analysis:
# merge raw results -> recompute PoseBusters -> 10k/10k SI statistics and figures.
# Uses the corrected MAO-B/PDGFRB structures (2v5z, 11ue), excluding 2z5x/3mjg.

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
RESULTS="${ROOT}/results"
OUT="${ROOT}/analysis/excluding_2z5x_3mjg"
TABLES="${OUT}/tables"
PB="${TABLES}/posebuster"
CONDA_SH="/mnt/tank/scratch/okonovalova/miniconda3/etc/profile.d/conda.sh"
BUST="/mnt/tank/scratch/okonovalova/miniconda3/envs/posebuster/bin/bust"

TARGETS_SPACE="1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 5lf3 7kk3"
TARGETS_CSV="${TARGETS_SPACE// /,}"
METHODS="boltz2 dynamicbind gnina plapt qvina"

mkdir -p "${OUT}/logs" "${TABLES}" "${PB}"

source "${CONDA_SH}"
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
export METHODS
export PYTHON="$(command -v python)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "================================================================"
echo "AffiTox corrected-structure full benchmark"
echo "job=${SLURM_JOB_ID:-n/a} node=${SLURM_NODELIST:-n/a}"
echo "started=$(date -Iseconds)"
echo "results=${RESULTS}"
echo "analysis=${OUT}"
echo "targets=${TARGETS}"
echo "excluded=2z5x 3mjg"
echo "python=${PYTHON}"
echo "================================================================"

command -v obabel
[[ -x "${BUST}" ]] || { echo "ERROR: PoseBusters CLI missing: ${BUST}" >&2; exit 1; }
echo "${BUST}"

echo "=== [0/4] Audit result-to-ligand identifiers ==="
"${PYTHON}" "${ROOT}/src/analysis/audit_result_alignment.py" \
  --root "${ROOT}" \
  --results-dir "${RESULTS}" \
  --targets "${TARGETS}" \
  --output-csv "${TABLES}/result_alignment_audit.csv" \
  --output-json "${TABLES}/result_alignment_audit_details.json"

echo "=== [1/4] Merge raw method results ==="
for target in ${TARGETS}; do
  echo "--- merge ${target} ---"
  "${PYTHON}" "${ROOT}/src/analysis/merge_ligands_docking_from_dir.py" \
    "${target}" \
    --base-dir "${ROOT}" \
    --output-dir "${TABLES}" \
    --results-dir "${RESULTS}" \
    --boltz-results-dir "${RESULTS}" \
    --dynamicbind-subdir dynamicbind_new
done

echo "=== [2/4] Add pKi column ==="
"${PYTHON}" "${ROOT}/src/analysis/add_pvalue_column.py" \
  --input-dir "${TABLES}" \
  --output-dir "${TABLES}"

echo "=== [3/4] Recompute PoseBusters from result poses ==="
"${PYTHON}" "${ROOT}/src/analysis/prepare_and_run_posebusters.py" \
  --poses-dir "${RESULTS}" \
  --proteins-dir "${ROOT}/processed/proteins" \
  --output-dir "${PB}" \
  --methods gnina,qvina,dynamicbind \
  --proteins "${TARGETS_CSV}" \
  --dynamicbind-subdir dynamicbind_new \
  --dynamicbind-method-label dynamicbind_new \
  --bust-bin "${BUST}" \
  --max-workers "${SLURM_CPUS_PER_TASK:-16}"

if [[ -d "${PB}/results_by_method" ]]; then
  cp -f "${PB}/results_by_method/"*.csv "${PB}/"
fi

"${PYTHON}" "${ROOT}/src/analysis/prepare_and_run_posebusters_boltz2.py" \
  --boltz-results-dir "${RESULTS}" \
  --output-dir "${PB}" \
  --proteins "${TARGETS_CSV}" \
  --bust-bin "${BUST}" \
  --max-workers "${SLURM_CPUS_PER_TASK:-16}"

if [[ -d "${PB}/results_by_method" ]]; then
  cp -f "${PB}/results_by_method/"*.csv "${PB}/"
fi

echo "=== [4/4] Full SI statistics and figures (10k/10k) ==="
bash "${ROOT}/pipeline/postprocess/run_si_full.sh"

echo "================================================================"
echo "FINISHED $(date -Iseconds)"
echo "Tables:  ${TABLES}"
echo "Figures: ${OUT}/figures"
echo "================================================================"
