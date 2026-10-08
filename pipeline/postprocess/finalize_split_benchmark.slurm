#!/bin/bash
#SBATCH --job-name=affitox_finalize
#SBATCH -p aichem
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-finalize_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/slurm-finalize_%j.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
OUT="${ROOT}/analysis/excluding_2z5x_3mjg"
TABLES="${OUT}/tables"
PB="${TABLES}/posebuster"

source "/mnt/tank/scratch/okonovalova/miniconda3/etc/profile.d/conda.sh"
set +u
conda activate docking
set -u

export TOXAFFINITY_ROOT="${ROOT}"
export RESULTS_DIR="${ROOT}/results"
export ANALYSIS_ROOT="${OUT}"
export DATA_DIR="${TABLES}"
export MERGED_DATA_DIR="${TABLES}"
export POSEBUSTERS_DIR="${PB}"
export POSEBUSTERS_BOLTZ2_DIR="${PB}"
export BOLTZ_RESULTS_DIR="${ROOT}/results"
export DYNAMICBIND_POSEBUSTERS_LABEL="dynamicbind_new"
export TARGETS="1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 7awe 7kk3"
export METHODS="boltz2 dynamicbind gnina plapt qvina"
export PYTHON="$(command -v python)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

for target in 2v5z 4tz4 11ue; do
  for method in gnina qvina dynamicbind_new boltz2; do
    test -s "${PB}/posebusters_results_${target}_${method}.csv"
  done
done
test -s "${PB}/posebusters_results_7kk3_gnina.csv"

bash "${ROOT}/pipeline/postprocess/run_si_full.sh"
