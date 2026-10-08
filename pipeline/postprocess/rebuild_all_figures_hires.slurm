#!/bin/bash
#SBATCH --job-name=affitox_fig_hires
#SBATCH -p aichem
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/rebuild_all_figures_hires_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/rebuild_all_figures_hires_%j.err

# Rebuild all excluding_2z5x_3mjg figures at print PNG resolution (no stats recompute).

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
export RUN_CORRELATIONS=0
export RUN_ENRICHMENT=0
export RUN_POSEBUSTERS=0
export RUN_INFERENTIAL=0
export RUN_FIGURES=1

echo "started=$(date -Iseconds)"
echo "node=${SLURM_NODELIST:-n/a} job=${SLURM_JOB_ID:-n/a}"

export TARGET_LABEL_MODE=pdb
export FIGURE_VARIANT=by_pdb
echo "=== by_pdb ==="
bash "${ROOT}/pipeline/postprocess/run_article_analysis.sh"

export TARGET_LABEL_MODE=protein
export FIGURE_VARIANT=by_protein
echo "=== by_protein ==="
bash "${ROOT}/pipeline/postprocess/run_article_analysis.sh"

echo "=== manuscript heatmaps ==="
"${PYTHON}" "${ROOT}/pipeline/postprocess/render_corrected_ultramarine_heatmap.py"
"${PYTHON}" "${ROOT}/pipeline/postprocess/render_manuscript_property_baseline.py"
cp -f "${OUT}/figures/correlations/heatmaps/"pearson_heatmap_ultramarine_ivory_vermilion_square_cells.{png,pdf,svg} \
  "${OUT}/figures/by_protein/correlations/heatmaps/"
cp -f "${OUT}/figures/correlations/heatmaps/"pearson_heatmap_ultramarine_ivory_vermilion_square_cells.{png,pdf,svg} \
  "${OUT}/figures/by_pdb/correlations/heatmaps/"
cp -f "${OUT}/figures/correlations/heatmaps/"pearson_heatmap_with_property_baselines.{png,pdf,svg} \
  "${OUT}/figures/by_protein/correlations/heatmaps/"
cp -f "${OUT}/figures/correlations/heatmaps/"pearson_heatmap_with_property_baselines.{png,pdf,svg} \
  "${OUT}/figures/by_pdb/correlations/heatmaps/"

echo "=== png sizes ==="
python - <<'PY'
from pathlib import Path
from PIL import Image
root = Path("/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/figures")
for rel in (
    "by_pdb/correlations/heatmaps/heatmap_correlations_combined.png",
    "by_pdb/correlations/heatmaps/heatmap_pearson_r.png",
    "by_protein/correlations/heatmaps/heatmap_correlations_combined.png",
    "correlations/heatmaps/pearson_heatmap_ultramarine_ivory_vermilion_square_cells.png",
):
    p = root / rel
    im = Image.open(p)
    print(f"{im.size[0]}x{im.size[1]}  {rel}")
PY

echo "finished=$(date -Iseconds)"
echo "png_count=$(find "${OUT}/figures" -name '*.png' | wc -l)"
