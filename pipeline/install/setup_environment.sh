#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../../config/project.env.sh
source "${SCRIPT_DIR}/../../config/project.env.sh"

echo "AffiTox environment setup"
echo "  root   : ${TOXAFFINITY_ROOT}"
echo "  python : ${PYTHON}"

cd "${TOXAFFINITY_ROOT}"
"${PYTHON}" -m pip install -e ".[analysis]" 2>/dev/null || "${PYTHON}" -m pip install -e ".[rdkit]" 2>/dev/null || "${PYTHON}" -m pip install -e .

if command -v conda >/dev/null 2>&1 && [[ -f "${TOXAFFINITY_ROOT}/environment.yml" ]]; then
  echo "Optional: create conda env with: conda env create -f environment.yml"
fi

echo "Done. Activate your env and set BOLTZ_RESULTS_DIR if using Boltz-2."
