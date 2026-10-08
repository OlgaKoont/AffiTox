#!/usr/bin/env bash
# Package Zenodo-ready archives from AffiTox workspace.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

OUT_DIR="${OUT_DIR:-${ROOT}/zenodo_packages}"
RESULTS_MODE="${RESULTS_MODE:-per_target}"   # per_target | single | skip
ZIP_LEVEL="${ZIP_LEVEL:-0}"                  # 0 = store (fast), 1-9 = compressed

mkdir -p "${OUT_DIR}" "${OUT_DIR}/results_by_target" "${OUT_DIR}/logs"

echo "=== Zenodo packaging ==="
echo "root        : ${ROOT}"
echo "out_dir     : ${OUT_DIR}"
echo "results_mode: ${RESULTS_MODE}"
echo "zip_level   : ${ZIP_LEVEL}"
echo "started     : $(date -Iseconds)"

cd "${ROOT}"

# Always refresh core archives.
rm -f "${OUT_DIR}/toxdock-input-raw.zip" "${OUT_DIR}/toxdock-analysis-artifacts.zip"
zip -r -"${ZIP_LEVEL}" "${OUT_DIR}/toxdock-input-raw.zip" "input/bindingdb"
zip -r -"${ZIP_LEVEL}" "${OUT_DIR}/toxdock-analysis-artifacts.zip" "analysis/tables" "analysis/figures"

case "${RESULTS_MODE}" in
  per_target)
    rm -f "${OUT_DIR}/results_by_target/"*.zip
    for d in results/*; do
      [[ -d "${d}" ]] || continue
      t="$(basename "${d}")"
      zip -r -"${ZIP_LEVEL}" "${OUT_DIR}/results_by_target/results_${t}.zip" "${d}"
    done
    ;;
  single)
    rm -f "${OUT_DIR}/toxdock-results-raw.zip"
    zip -r -"${ZIP_LEVEL}" "${OUT_DIR}/toxdock-results-raw.zip" "results"
    ;;
  skip)
    echo "Skipping results archiving (RESULTS_MODE=skip)"
    ;;
  *)
    echo "Unknown RESULTS_MODE='${RESULTS_MODE}'. Use: per_target | single | skip" >&2
    exit 2
    ;;
esac

cat > "${OUT_DIR}/README_ZENODO_PACKAGES.txt" <<'EOF'
Zenodo package manifest
=======================

1) toxdock-input-raw.zip
   - contains: input/bindingdb/
   - purpose: raw upstream snapshot / provenance rebuilds

2) toxdock-analysis-artifacts.zip
   - contains: analysis/tables/ and analysis/figures/
   - purpose: manuscript/SI-ready tables and figures

3) results archives
   - RESULTS_MODE=per_target: results_by_target/results_<target>.zip
   - RESULTS_MODE=single: toxdock-results-raw.zip
   - purpose: full raw docking outputs by target/method
EOF

echo "finished    : $(date -Iseconds)"
echo "Artifacts:"
ls -lh "${OUT_DIR}"
