#!/bin/bash
#SBATCH --job-name=affitox_fig_text
#SBATCH -p aichem
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/rebuild_by_protein_vector_pdfs_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/rebuild_by_protein_vector_pdfs_%j.err

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
OUT="${ROOT}/analysis/excluding_2z5x_3mjg"
TABLES="${OUT}/tables"
PB="${TABLES}/posebuster"

source "/mnt/tank/scratch/okonovalova/miniconda3/etc/profile.d/conda.sh"
set +u
conda activate docking
set -u

export MPLBACKEND=Agg
export TOXAFFINITY_ROOT="${ROOT}"
export RESULTS_DIR="${ROOT}/results"
export ANALYSIS_ROOT="${OUT}"
export DATA_DIR="${TABLES}"
export MERGED_DATA_DIR="${TABLES}"
export POSEBUSTERS_DIR="${PB}"
export POSEBUSTERS_BOLTZ2_DIR="${PB}"
export DYNAMICBIND_POSEBUSTERS_LABEL="dynamicbind_new"
export TARGETS="1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 5lf3 7kk3"
export METHODS="boltz2 dynamicbind gnina plapt qvina"
export PYTHON="$(command -v python)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export COLOR_LOW="#3558C5"
export COLOR_MID="#FFF9EE"
export COLOR_HIGH="#EF4938"
export FIGURE_DPI=300
export FONT_FAMILY="DejaVu Sans"
export TARGET_LABEL_MODE=protein
export FIGURE_VARIANT=by_protein
export RUN_CORRELATIONS=0
export RUN_ENRICHMENT=0
export RUN_POSEBUSTERS=0
export RUN_INFERENTIAL=0
export RUN_FIGURES=1

echo "started=$(date -Iseconds)"
echo "figures=${OUT}/figures/by_protein"
bash "${ROOT}/pipeline/postprocess/run_article_analysis.sh"
echo "finished=$(date -Iseconds)"
echo "pdf_count=$(find "${OUT}/figures/by_protein" -name '*.pdf' | wc -l)"
