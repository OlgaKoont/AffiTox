#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/project.env.sh"

POST_PY="${TOXAFFINITY_ROOT}/src/analysis"
mkdir -p "${POSEBUSTERS_DIR}" "${POSEBUSTERS_BOLTZ2_DIR}"

cd "${TOXAFFINITY_ROOT}"

deposited_csv=0
shopt -s nullglob
for f in "${TOXAFFINITY_ROOT}/analysis/tables/posebuster"/posebusters_results_*.csv; do
  deposited_csv=1
  break
done
shopt -u nullglob

if ! command -v bust >/dev/null 2>&1; then
  if [[ "${deposited_csv}" -eq 1 ]]; then
    echo "PoseBusters CLI 'bust' is not on PATH."
    echo "Keeping deposited CSVs in analysis/tables/posebuster/ (analysis can use them as-is)."
    echo "Install bust only if you want to recompute pass rates from poses."
    exit 0
  fi
  echo "ERROR: PoseBusters CLI 'bust' is not on PATH, and no deposited posebusters_results_*.csv files were found." >&2
  echo "Install bust, or clone AffiTox with analysis/tables/posebuster already present." >&2
  exit 1
fi

if [[ ! -d "${RESULTS_DIR}" ]]; then
  if [[ "${deposited_csv}" -eq 1 ]]; then
    echo "No ${RESULTS_DIR} tree (raw poses). Keeping deposited PoseBusters CSVs."
    echo "Unpack results_<pdb>.zip from Zenodo to recompute bust."
    exit 0
  fi
  echo "ERROR: ${RESULTS_DIR} does not exist. Unpack Zenodo results_<pdb>.zip or run dock." >&2
  exit 1
fi

echo "=== PoseBusters (gnina, qvina, dynamicbind) ==="
"${PYTHON}" "${POST_PY}/prepare_and_run_posebusters.py" \
  --poses-dir "${RESULTS_DIR}" \
  --proteins-dir "${PROCESSED_DIR}/proteins" \
  --output-dir "${ANALYSIS_ROOT}/tables/posebuster" \
  "$@"
if [[ -d "${POSEBUSTERS_DIR}/results_by_method" ]]; then
  mv "${POSEBUSTERS_DIR}/results_by_method/"*.csv "${POSEBUSTERS_DIR}/" 2>/dev/null || true
fi

if [[ -n "${BOLTZ_RESULTS_DIR}" ]]; then
  echo "=== PoseBusters (Boltz-2) ==="
  "${PYTHON}" "${POST_PY}/prepare_and_run_posebusters_boltz2.py" \
    --boltz-results-dir "${BOLTZ_RESULTS_DIR}" \
    --output-dir "${ANALYSIS_ROOT}/tables/posebuster" \
    "$@"
  if [[ -d "${POSEBUSTERS_BOLTZ2_DIR}/results_by_method" ]]; then
    mv "${POSEBUSTERS_BOLTZ2_DIR}/results_by_method/"*.csv "${POSEBUSTERS_BOLTZ2_DIR}/" 2>/dev/null || true
  fi
fi
