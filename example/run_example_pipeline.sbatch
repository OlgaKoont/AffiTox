#!/usr/bin/env bash
#SBATCH --job-name=toxaff-example
#SBATCH --output=/mnt/tank/scratch/okonovalova/AffiTox/example/logs/run_example_%j.out
#SBATCH --error=/mnt/tank/scratch/okonovalova/AffiTox/example/logs/run_example_%j.err
#SBATCH --account=ailab
#SBATCH --partition=aichem
#SBATCH --qos=high_gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

ROOT="/mnt/tank/scratch/okonovalova/AffiTox"
cd "${ROOT}"
mkdir -p "${ROOT}/example/logs"

# Toolchain: per-method conda envs (see config/example_toolchain.env.sh via project.env.example.sh).
export CONDA_BASE="${CONDA_BASE:-/mnt/tank/scratch/okonovalova/miniconda3}"
export SKIP_INSTALL=1

# Fast stats for smoke test (override at submit time for heavier run).
export N_BOOTSTRAP="${N_BOOTSTRAP:-200}"
export N_PERMUTATION="${N_PERMUTATION:-200}"
export FIGURE_DPI="${FIGURE_DPI:-150}"

# Boltz: set BOLTZ_SKIP_RUN=0 to run predict on GPU inside this job.
export BOLTZ_SKIP_RUN="${BOLTZ_SKIP_RUN:-1}"
export BOLTZ_SYNC_MODE="${BOLTZ_SYNC_MODE:-symlink}"

# Optional: only subset of targets, e.g. TARGETS_ONLY="5mo4 4f65"
export TARGETS_ONLY="${TARGETS_ONLY:-}"

echo "================================================================"
echo "AffiTox example pipeline (SLURM)"
echo "  started : $(date -Iseconds)"
echo "  host    : $(hostname)"
echo "  BOLTZ_SKIP_RUN=${BOLTZ_SKIP_RUN}"
echo "  N_BOOTSTRAP=${N_BOOTSTRAP}  N_PERMUTATION=${N_PERMUTATION}"
echo "================================================================"

# Stages passed to run_example_pipeline.sh (default: full GPU+postprocess tail).
export EXAMPLE_STAGES="${EXAMPLE_STAGES:-prepare-inputs boltz-prepare prepare dock boltz-sync merge posebusters analysis}"

echo "  stages   : ${EXAMPLE_STAGES}"
echo "================================================================"

bash "${ROOT}/example/run_example_pipeline.sh" ${EXAMPLE_STAGES}

if [[ "${BOLTZ_SKIP_RUN}" == "0" ]] && [[ " ${EXAMPLE_STAGES} " != *" boltz-run "* ]]; then
  echo "=== boltz-run (extra GPU step) ==="
  if [[ -n "${TARGETS_ONLY}" ]]; then
    TARGETS_ONLY="${TARGETS_ONLY}" bash "${ROOT}/example/scripts/run_boltz2_example.sh"
  else
    bash "${ROOT}/example/scripts/run_boltz2_example.sh"
  fi
  bash "${ROOT}/example/run_example_pipeline.sh" boltz-sync merge posebusters analysis
fi

echo "================================================================"
echo "Finished: $(date -Iseconds)"
echo "Outputs : ${ROOT}/example/{results,analysis}"
echo "================================================================"
