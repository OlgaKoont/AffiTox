#!/usr/bin/env bash
# Zip the six AffiTox information blocks for Zenodo.
# Record 03 exceeds the 50 GiB / record limit (~193 GiB raw), so docking zips
# are grouped into sibling 03a–03e records. Do not reuse 10.5281/zenodo.20825057–067.
set -euo pipefail
ROOT="${TOXAFFINITY_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
PYTHON="${PYTHON:-/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python}"
export TOXAFFINITY_ROOT="${ROOT}"
"${PYTHON}" "${ROOT}/pipeline/postprocess/pack_zenodo_blocks.py"
