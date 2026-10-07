#!/usr/bin/env bash
# Stage 3: merge ligand tables on (pdb_id, ligand_id), then add pKi / is_active.
#
# EXAMPLE_N=5  # first five ligand_id rows (ligand_0001..ligand_0005); same five
#              # for all methods. Default: all ligand_id rows for each protein.
# Cluster: pipeline/hpc/run_merge_ligand_id.slurm is an aichem example only.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/project.env.sh"

export PYTHONPATH="${TOXAFFINITY_ROOT}/src:${PYTHONPATH:-}"
export BOLTZ_RESULTS_DIR="${BOLTZ_RESULTS_DIR:-${RESULTS_DIR}}"

POST_PY="${TOXAFFINITY_ROOT}/src/analysis"
mkdir -p "${MERGED_DATA_DIR}" "${POSEBUSTERS_DIR}" "${POSEBUSTERS_BOLTZ2_DIR}"

cd "${TOXAFFINITY_ROOT}"

echo "=== Merge docking metrics -> ${MERGED_DATA_DIR} ==="
"${PYTHON}" "${POST_PY}/merge_ligands_docking_from_dir.py" \
  --base-dir "${TOXAFFINITY_ROOT}" \
  --output-dir "${MERGED_DATA_DIR}" \
  --results-dir "${RESULTS_DIR}" \
  --dynamicbind-subdir dynamicbind_new \
  ${BOLTZ_RESULTS_DIR:+--boltz-results-dir "${BOLTZ_RESULTS_DIR}"} \
  "$@"

echo "=== Add pValue column ==="
"${PYTHON}" "${POST_PY}/add_pvalue_column.py" \
  --input-dir "${MERGED_DATA_DIR}" \
  --output-dir "${MERGED_DATA_DIR}"

echo "Postprocess merge complete."
