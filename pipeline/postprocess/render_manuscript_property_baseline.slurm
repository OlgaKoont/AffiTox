#!/bin/bash
#SBATCH --job-name=pearson_baseline_pdf
#SBATCH -p aichem
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pearson-baseline_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-pearson-baseline_%j.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
PYTHON="/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python"

cd "${ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

"${PYTHON}" "${ROOT}/pipeline/postprocess/render_manuscript_property_baseline.py"
