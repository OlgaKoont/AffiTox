#!/usr/bin/env bash
#SBATCH --job-name=toxdock-analysis
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/analysis/logs/article_analysis_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/analysis/logs/article_analysis_%j.err
#SBATCH --account=ailab
#SBATCH --partition=aichem
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
cd "${ROOT}"

export PYTHON="/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python"

# Full-data statistical run (canonical 10k bootstrap / 10k permutation).
# Override at submit time if needed, e.g.:
#   sbatch --export=ALL,RUN_FIGURES=1 pipeline/postprocess/run_article_analysis.sbatch

export N_BOOTSTRAP="${N_BOOTSTRAP:-10000}"
export N_PERMUTATION="${N_PERMUTATION:-10000}"
export RANDOM_SEED="${RANDOM_SEED:-42}"

export RUN_CORRELATIONS="${RUN_CORRELATIONS:-1}"
export RUN_ENRICHMENT="${RUN_ENRICHMENT:-1}"
export RUN_POSEBUSTERS="${RUN_POSEBUSTERS:-1}"
export RUN_INFERENTIAL="${RUN_INFERENTIAL:-1}"
export RUN_FIGURES="${RUN_FIGURES:-0}"

echo "================================================================"
echo "AffiTox article analysis (SLURM)"
echo "  started : $(date -Iseconds)"
echo "  host    : $(hostname)"
echo "  N_BOOTSTRAP=${N_BOOTSTRAP}  N_PERMUTATION=${N_PERMUTATION}"
echo "  RUN_CORRELATIONS=${RUN_CORRELATIONS} RUN_ENRICHMENT=${RUN_ENRICHMENT}"
echo "  RUN_POSEBUSTERS=${RUN_POSEBUSTERS} RUN_INFERENTIAL=${RUN_INFERENTIAL}"
echo "  RUN_FIGURES=${RUN_FIGURES}"
echo "================================================================"

bash "${ROOT}/pipeline/postprocess/run_article_analysis.sh"

echo "================================================================"
echo "Finished: $(date -Iseconds)"
echo "================================================================"
