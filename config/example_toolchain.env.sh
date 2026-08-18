#!/usr/bin/env bash
# Conda/tool paths for the example pipeline on HPC (login + GPU nodes).
# Source from config/project.env.example.sh
# Export CONDA_BASE / GNINA_BIN / HF_HOME for your machine; no cluster paths are assumed.

export CONDA_BASE="${CONDA_BASE:-}"

export CONDA_ENV_DOCKING="${CONDA_ENV_DOCKING:-docking}"
export CONDA_ENV_QVINA="${CONDA_ENV_QVINA:-docking}"
export CONDA_ENV_GNINA="${CONDA_ENV_GNINA:-gnina}"
export CONDA_ENV_PLAPT="${CONDA_ENV_PLAPT:-plapt}"
export CONDA_ENV_DYNAMICBIND="${CONDA_ENV_DYNAMICBIND:-dynamicbind}"
export CONDA_ENV_POSEBUSTER="${CONDA_ENV_POSEBUSTER:-posebuster}"

_default_py="$(command -v python3 || echo python3)"
if [[ -n "${CONDA_BASE}" ]]; then
  export PYTHON_DOCKING="${PYTHON_DOCKING:-${CONDA_BASE}/envs/${CONDA_ENV_DOCKING}/bin/python}"
  export PYTHON_GNINA="${PYTHON_GNINA:-${CONDA_BASE}/envs/${CONDA_ENV_GNINA}/bin/python}"
  export PYTHON_PLAPT="${PYTHON_PLAPT:-${CONDA_BASE}/envs/${CONDA_ENV_PLAPT}/bin/python}"
  export PYTHON_DYNAMICBIND="${PYTHON_DYNAMICBIND:-${CONDA_BASE}/envs/${CONDA_ENV_DYNAMICBIND}/bin/python}"
  export PYTHON_POSEBUSTER="${PYTHON_POSEBUSTER:-${CONDA_BASE}/envs/${CONDA_ENV_POSEBUSTER}/bin/python}"
  export OBABEL_BIN="${OBABEL_BIN:-${CONDA_BASE}/envs/${CONDA_ENV_DOCKING}/bin/obabel}"
  export BUST_BIN="${BUST_BIN:-${CONDA_BASE}/envs/${CONDA_ENV_POSEBUSTER}/bin/bust}"
else
  export PYTHON_DOCKING="${PYTHON_DOCKING:-${_default_py}}"
  export PYTHON_GNINA="${PYTHON_GNINA:-${_default_py}}"
  export PYTHON_PLAPT="${PYTHON_PLAPT:-${_default_py}}"
  export PYTHON_DYNAMICBIND="${PYTHON_DYNAMICBIND:-${_default_py}}"
  export PYTHON_POSEBUSTER="${PYTHON_POSEBUSTER:-${_default_py}}"
  export OBABEL_BIN="${OBABEL_BIN:-obabel}"
  export BUST_BIN="${BUST_BIN:-bust}"
fi
unset _default_py
export PYTHON_QVINA="${PYTHON_QVINA:-${PYTHON_DOCKING}}"
# Standalone GNINA binary, or `gnina` on PATH.
export GNINA_BIN="${GNINA_BIN:-gnina}"

# HuggingFace cache (PLAPT needs Rostlab/prot_bert).
export HF_HOME="${HF_HOME:-${HOME}/.cache/huggingface}"
export TRANSFORMERS_CACHE="${TRANSFORMERS_CACHE:-${HF_HOME}}"
export HF_HUB_CACHE="${HF_HUB_CACHE:-${HF_HOME}/hub}"
mkdir -p "${HF_HUB_CACHE}" 2>/dev/null || true

# Default driver python for prepare/merge/analysis (docking env).
export PYTHON="${PYTHON:-${PYTHON_DOCKING}}"

init_example_conda() {
  if [[ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
  fi
}

run_in_conda_env() {
  local env_name="$1"
  shift
  init_example_conda
  if [[ -n "${env_name}" ]]; then
    conda run -n "${env_name}" --no-capture-output "$@"
  else
    "$@"
  fi
}
