#!/bin/bash
#SBATCH --job-name=affitox_bmc_pdf
#SBATCH -p aichem
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-bmc-pdf_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-bmc-pdf_%j.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
PYTHON="/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python"

"${PYTHON}" "${ROOT}/pipeline/postprocess/export_bmc_vector_pdfs.py" \
  --figures-dir "${ROOT}/analysis/figures" \
  --workers "${SLURM_CPUS_PER_TASK:-16}"
