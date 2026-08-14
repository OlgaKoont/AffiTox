#!/usr/bin/env bash
# Stage 1: prepare proteins, ligands, and docking boxes.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/project.env.sh"

cd "${TOXAFFINITY_ROOT}"
export PYTHONPATH="${TOXAFFINITY_ROOT}/src:${PYTHONPATH:-}"

echo "=== AffiTox: preparation ==="
"${PYTHON}" -m docking_benchmark2.cli.run_benchmark \
  --config "${TOXDOCK_CONFIG}" \
  --methods-config "${METHODS_CONFIG}" \
  --stage preparation "$@"
