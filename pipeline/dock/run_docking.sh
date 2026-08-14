#!/usr/bin/env bash
# Stage 2: run docking methods (qvina, gnina, plapt, dynamicbind).
# Boltz-2 is run externally; set BOLTZ_RESULTS_DIR before merge/postprocess.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/project.env.sh"

cd "${TOXAFFINITY_ROOT}"
export PYTHONPATH="${TOXAFFINITY_ROOT}/src:${PYTHONPATH:-}"

echo "=== AffiTox: docking ==="
"${PYTHON}" -m docking_benchmark2.cli.run_benchmark \
  --config "${TOXDOCK_CONFIG}" \
  --methods-config "${METHODS_CONFIG}" \
  --stage docking "$@"
