#!/usr/bin/env bash
#SBATCH --job-name=affitox_zenodo_pack
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/release/zenodo/payloads/logs/pack_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/release/zenodo/payloads/logs/pack_%j.err
#SBATCH --partition=aichem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00

set -euo pipefail
ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
PYTHON="${PYTHON:-/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python}"
mkdir -p "${ROOT}/release/zenodo/payloads/logs"
export TOXAFFINITY_ROOT="${ROOT}"
"${PYTHON}" "${ROOT}/pipeline/postprocess/pack_zenodo_blocks.py"
