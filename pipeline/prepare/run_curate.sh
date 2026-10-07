#!/usr/bin/env bash
# Layer 1 → 2: ChEMBL snapshot → curated tables with ligand_id.
#
# Default --mode frozen-nodubl keeps the row order that was actually docked
# (ligand_0001 = first row of input/ligands_nodubl). Use --mode from-raw to
# apply Ki / unit / MW / SMILES filters to input/ligands/.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/project.env.sh"

cd "${TOXAFFINITY_ROOT}"
export PYTHONPATH="${TOXAFFINITY_ROOT}/src:${PYTHONPATH:-}"

# EXAMPLE_N=5  # first five ligand_id rows only (see ligand_identity.py --example-n)
# default: all ligand_id rows for each protein

echo "=== AffiTox: curate ligands (layer 1 → 2) ==="
"${PYTHON}" -m docking_benchmark2.preprocessing.curate_chembl \
  --base-dir "${TOXAFFINITY_ROOT}" \
  "$@"
