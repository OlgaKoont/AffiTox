#!/bin/bash
#SBATCH --job-name=affitox_merge_id
#SBATCH -p aichem
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/merge_ligand_id_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/excluding_2z5x_3mjg/merge_ligand_id_%j.err
#
# Cluster example. Portable entrypoint is:
#   bash pipeline/postprocess/run_merge.sh
# Replace this wrapper on another site (partition, account, modules).

set -euo pipefail
ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
CONDA_SH="/mnt/tank/scratch/okonovalova/miniconda3/etc/profile.d/conda.sh"

source "${CONDA_SH}"
set +u
conda activate docking
set -u

export TOXAFFINITY_ROOT="${ROOT}"
export PYTHON="$(command -v python)"
export PYTHONPATH="${ROOT}/src:${PYTHONPATH:-}"
export ANALYSIS_ROOT="${ROOT}/analysis/excluding_2z5x_3mjg"
export MERGED_DATA_DIR="${ANALYSIS_ROOT}/tables"
export RESULTS_DIR="${ROOT}/results"
export BOLTZ_RESULTS_DIR="${ROOT}/results"
export TARGETS="1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 5lf3 7kk3"

mkdir -p "${MERGED_DATA_DIR}" "${ANALYSIS_ROOT}"

echo "=== maps + curated tables ==="
"${PYTHON}" "${ROOT}/src/analysis/ligand_identity.py" --base-dir "${ROOT}"
"${PYTHON}" -m docking_benchmark2.preprocessing.curate_chembl --base-dir "${ROOT}" --mode frozen-nodubl

echo "=== full-panel merge on ligand_id ==="
bash "${ROOT}/pipeline/postprocess/run_merge.sh"

echo "=== example 16x5 merge ==="
# EXAMPLE_N=5  # first five ligand_id rows; same five for all methods
EXAMPLE_N=5 bash "${ROOT}/pipeline/postprocess/run_merge.sh" \
  --output-dir "${ROOT}/example/04_merged"

echo "=== copy method weights into Zenodo record 02 ==="
bash "${ROOT}/release/zenodo/02_prepared/weights/pack_weights.sh" || true

echo "=== refresh Zenodo 04/06 + example slices ==="
"${PYTHON}" "${ROOT}/src/analysis/package_release.py" --base-dir "${ROOT}"

echo "finished=$(date -Iseconds) N=$(cat "${ROOT}/input/ligands_curated/N_panel.txt")"
