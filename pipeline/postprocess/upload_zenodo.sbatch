#!/usr/bin/env bash
#SBATCH --job-name=zenodo-upload
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/zenodo_packages/logs/zenodo_upload_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/zenodo_packages/logs/zenodo_upload_%j.err
#SBATCH --account=ailab
#SBATCH --partition=aichem
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=72:00:00

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
PYTHON="${PYTHON:-/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python}"

mkdir -p "${ROOT}/zenodo_packages/logs"

echo "Zenodo upload started: $(date -Iseconds)"
"${PYTHON}" "${ROOT}/pipeline/postprocess/upload_zenodo.py" "$@"
echo "Zenodo upload finished: $(date -Iseconds)"
