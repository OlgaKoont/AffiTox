#!/bin/bash
#SBATCH --job-name=pb_convert_4tz4
#SBATCH -p aichem
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pb-convert-4tz4_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pb-convert-4tz4_%j.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
PYTHON="/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python"

cd "${ROOT}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${PYTHON}" "${ROOT}/pipeline/postprocess/benchmark_pose_conversion_4tz4.py"
