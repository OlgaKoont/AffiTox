#!/bin/bash
#SBATCH --job-name=affitox_5lf3
#SBATCH -p aichem
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/slurm-refresh-5lf3_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/slurm-refresh-5lf3_%j.err

# Replace 7AWE with 5LF3 (constitutive 20S PSMB5) in excluding_2z5x_3mjg:
# merge raw 5lf3 docking, PoseBusters, 10k/10k stats, by_pdb + by_protein figures.

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
RESULTS="${ROOT}/results"
OUT="${ROOT}/analysis/excluding_2z5x_3mjg"
TABLES="${OUT}/tables"
PB="${TABLES}/posebuster"
CANONICAL_PB="${ROOT}/analysis/tables/posebuster"
CONDA_SH="/mnt/tank/scratch/okonovalova/miniconda3/etc/profile.d/conda.sh"
BUST="/mnt/tank/scratch/okonovalova/miniconda3/envs/posebuster/bin/bust"
BOLTZ_SRC="/mnt/tank/scratch/okonovalova/boltz/data/results/5lf3/docking/boltz2_s80_d10_seed42_aff_s80_d1"

TARGETS_SPACE="1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 5lf3 7kk3"
METHODS="boltz2 dynamicbind gnina plapt qvina"

mkdir -p "${OUT}/logs" "${TABLES}" "${PB}" "${RESULTS}/5lf3/docking"

source "${CONDA_SH}"
set +u
conda activate docking
set -u

export TOXAFFINITY_ROOT="${ROOT}"
export RESULTS_DIR="${RESULTS}"
export ANALYSIS_ROOT="${OUT}"
export DATA_DIR="${TABLES}"
export MERGED_DATA_DIR="${TABLES}"
export POSEBUSTERS_DIR="${PB}"
export POSEBUSTERS_BOLTZ2_DIR="${PB}"
export BOLTZ_RESULTS_DIR="${RESULTS}"
export DYNAMICBIND_POSEBUSTERS_LABEL="dynamicbind_new"
export TARGETS="${TARGETS_SPACE}"
export METHODS
export PYTHON="$(command -v python)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export COLOR_LOW="#3558C5"
export COLOR_MID="#FFF9EE"
export COLOR_HIGH="#EF4938"
export FIGURE_DPI=300
export FONT_FAMILY="DejaVu Sans"
export MPLBACKEND=Agg

echo "================================================================"
echo "AffiTox refresh: 5LF3 replaces 7AWE"
echo "job=${SLURM_JOB_ID:-n/a} node=${SLURM_NODELIST:-n/a}"
echo "started=$(date -Iseconds)"
echo "================================================================"

command -v obabel
[[ -x "${BUST}" ]] || { echo "ERROR: PoseBusters CLI missing: ${BUST}" >&2; exit 1; }
[[ -d "${BOLTZ_SRC}" ]] || { echo "ERROR: missing Boltz results ${BOLTZ_SRC}" >&2; exit 1; }

cp -f "${ROOT}/processed/proteins/cleaned/5lf3_chainK_protein_only.pdb" \
  "${ROOT}/processed/proteins/5lf3_chainK.pdb"
ln -sfn "${BOLTZ_SRC}" "${RESULTS}/5lf3/docking/boltz2_s80_d10_seed42_aff_s80_d1"

echo "=== [1/5] Merge 5lf3 docking into analysis tables ==="
"${PYTHON}" "${ROOT}/src/analysis/merge_ligands_docking_from_dir.py" \
  5lf3 \
  --base-dir "${ROOT}" \
  --output-dir "${TABLES}" \
  --results-dir "${RESULTS}" \
  --boltz-results-dir "${RESULTS}" \
  --dynamicbind-subdir dynamicbind_new

"${PYTHON}" "${ROOT}/src/analysis/add_pvalue_column.py" \
  --input-dir "${TABLES}" \
  --output-dir "${TABLES}"

echo "=== [2/5] PoseBusters for 5lf3 (keep other targets) ==="
rm -f "${PB}/posebusters_results_5lf3_"*.csv
rm -f "${PB}/results_by_method/posebusters_results_5lf3_"*.csv

"${PYTHON}" "${ROOT}/src/analysis/prepare_and_run_posebusters.py" \
  --poses-dir "${RESULTS}" \
  --proteins-dir "${ROOT}/processed/proteins" \
  --output-dir "${PB}" \
  --methods gnina,qvina,dynamicbind \
  --proteins 5lf3 \
  --dynamicbind-subdir dynamicbind_new \
  --dynamicbind-method-label dynamicbind_new \
  --bust-bin "${BUST}" \
  --max-workers "${SLURM_CPUS_PER_TASK:-32}"

if [[ -d "${PB}/results_by_method" ]]; then
  cp -f "${PB}/results_by_method/"*.csv "${PB}/"
fi

"${PYTHON}" "${ROOT}/src/analysis/prepare_and_run_posebusters_boltz2.py" \
  --boltz-results-dir "${RESULTS}" \
  --output-dir "${PB}" \
  --proteins 5lf3 \
  --bust-bin "${BUST}" \
  --max-workers "${SLURM_CPUS_PER_TASK:-32}"

if [[ -d "${PB}/results_by_method" ]]; then
  cp -f "${PB}/results_by_method/"*.csv "${PB}/"
fi

echo "=== [3/5] Full SI statistics and PDB-labelled figures ==="
export TARGET_LABEL_MODE=pdb
export FIGURE_VARIANT=by_pdb
bash "${ROOT}/pipeline/postprocess/run_si_full.sh"

echo "=== [4/5] Protein-labelled figure set ==="
export TARGET_LABEL_MODE=protein
export FIGURE_VARIANT=by_protein
export RUN_CORRELATIONS=0
export RUN_ENRICHMENT=0
export RUN_POSEBUSTERS=0
export RUN_INFERENTIAL=0
export RUN_FIGURES=1
bash "${ROOT}/pipeline/postprocess/run_article_analysis.sh"

echo "=== [5/5] Drop stale 7AWE per-target figures ==="
find "${OUT}/figures" -iname '*7awe*' -delete
rm -f "${TABLES}/merged_ligands_docking_7awe.csv" || true

echo "================================================================"
echo "finished=$(date -Iseconds)"
echo "merged=${TABLES}/merged_ligands_docking_5lf3.csv"
echo "protein_figures=${OUT}/figures/by_protein"
echo "================================================================"
