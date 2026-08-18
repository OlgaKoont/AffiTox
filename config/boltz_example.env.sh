#!/usr/bin/env bash
# Boltz-2 path defaults for the AffiTox example pipeline.
# Source after config/project.env.example.sh

# External Boltz workspace (same layout as run_boltz2_docking_bindingdb_slurm.sh).
# Set BOLTZ_ROOT and BOLTZ_REPO to your Boltz checkout/data tree before running.
export BOLTZ_ROOT="${BOLTZ_ROOT:-}"
export BOLTZ_REPO="${BOLTZ_REPO:-}"

# Inputs (defaults: mini example ligands + precomputed assets from BOLTZ_ROOT).
export BOLTZ_PROTEIN_CIF_DIR="${BOLTZ_PROTEIN_CIF_DIR:-${BOLTZ_ROOT}/input/proteins}"
export BOLTZ_LIGAND_DIR="${BOLTZ_LIGAND_DIR:-${EXAMPLE_ROOT}/input/ligands_nodubl}"
export BOLTZ_MSA_DIR="${BOLTZ_MSA_DIR:-${BOLTZ_ROOT}/msa_cache/precomputed/a3m}"

# Where `boltz predict` writes (external tree by default).
export BOLTZ_WRITE_DIR="${BOLTZ_WRITE_DIR:-${BOLTZ_ROOT}/results}"

# Run tag subdirectory under {pdb}/docking/ (must contain "aff" for merge).
export BOLTZ_RUN_TAG="${BOLTZ_RUN_TAG:-boltz2_s80_d10_seed42_aff_s80_d1}"

# Boltz CLI / env (override on compute node).
export BOLTZ_CONDA_ENV="${BOLTZ_CONDA_ENV:-boltz-env}"
export BOLTZ_SEED="${BOLTZ_SEED:-42}"
export BOLTZ_SAMPLING_STEPS="${BOLTZ_SAMPLING_STEPS:-80}"
export BOLTZ_DIFFUSION_SAMPLES="${BOLTZ_DIFFUSION_SAMPLES:-10}"
export BOLTZ_SAMPLING_STEPS_AFFINITY="${BOLTZ_SAMPLING_STEPS_AFFINITY:-80}"
export BOLTZ_DIFFUSION_SAMPLES_AFFINITY="${BOLTZ_DIFFUSION_SAMPLES_AFFINITY:-1}"

# Sync external Boltz output into example/results/ (symlink | copy).
export BOLTZ_SYNC_MODE="${BOLTZ_SYNC_MODE:-symlink}"
export BOLTZ_SKIP_RUN="${BOLTZ_SKIP_RUN:-1}"
export BOLTZ_SKIP_SYNC="${BOLTZ_SKIP_SYNC:-0}"

# merge/posebusters read path:
#   - after sync: example/results (default once sync ran)
#   - direct external: BOLTZ_USE_EXTERNAL=1 BOLTZ_RESULTS_DIR=${BOLTZ_WRITE_DIR}
export BOLTZ_USE_EXTERNAL="${BOLTZ_USE_EXTERNAL:-0}"
if [[ "${BOLTZ_USE_EXTERNAL}" == "1" ]]; then
  export BOLTZ_RESULTS_DIR="${BOLTZ_RESULTS_DIR:-${BOLTZ_WRITE_DIR}}"
else
  export BOLTZ_RESULTS_DIR="${BOLTZ_RESULTS_DIR:-${RESULTS_DIR}}"
fi
