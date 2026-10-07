#!/usr/bin/env bash
#SBATCH --job-name=affitox_zenodo_up
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/release/zenodo/payloads/logs/upload_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/release/zenodo/payloads/logs/upload_%j.err
#SBATCH --partition=aichem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=72:00:00

set -euo pipefail
ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
PYTHON="${PYTHON:-/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python}"
mkdir -p "${ROOT}/release/zenodo/payloads/logs"
if [[ -f "${HOME}/.zenodo/token" ]]; then
  export ZENODO_ACCESS_TOKEN
  ZENODO_ACCESS_TOKEN="$(tr -d '[:space:]' < "${HOME}/.zenodo/token")"
fi
"${PYTHON}" "${ROOT}/pipeline/postprocess/upload_zenodo_blocks.py" "$@"
