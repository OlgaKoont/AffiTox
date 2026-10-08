#!/usr/bin/env bash
#SBATCH --job-name=toxdock-zenodo-pack
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/zenodo_packages/logs/zenodo_pack_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/zenodo_packages/logs/zenodo_pack_%j.err
#SBATCH --account=ailab
#SBATCH --partition=aichem
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
cd "${ROOT}"

# Optional overrides at submit time:
#   sbatch --export=ALL,RESULTS_MODE=single,ZIP_LEVEL=0 pipeline/postprocess/package_zenodo.sbatch
#   sbatch --export=ALL,RESULTS_MODE=per_target,ZIP_LEVEL=0 pipeline/postprocess/package_zenodo.sbatch

bash "${ROOT}/pipeline/postprocess/package_zenodo.sh"
