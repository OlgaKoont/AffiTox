#!/usr/bin/env bash
# AffiTox — canonical path contract (JCIM reproducibility).
# Source from any pipeline script:
#   source "$(dirname "${BASH_SOURCE[0]}")/project.env.sh"

if [[ -z "${TOXAFFINITY_ROOT:-}" ]]; then
  _cfg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  TOXAFFINITY_ROOT="$(cd "${_cfg_dir}/.." && pwd)"
  unset _cfg_dir
fi
export TOXAFFINITY_ROOT

# --- layout (matches published artifact tree) ---
export INPUT_DIR="${INPUT_DIR:-${TOXAFFINITY_ROOT}/input}"
export PROCESSED_DIR="${PROCESSED_DIR:-${TOXAFFINITY_ROOT}/processed}"
export RESULTS_DIR="${RESULTS_DIR:-${TOXAFFINITY_ROOT}/results}"
export ANALYSIS_ROOT="${ANALYSIS_ROOT:-${TOXAFFINITY_ROOT}/analysis}"
export DATA_DIR="${DATA_DIR:-${ANALYSIS_ROOT}/tables}"
export MERGED_DATA_DIR="${MERGED_DATA_DIR:-${DATA_DIR}}"
export POSEBUSTERS_DIR="${POSEBUSTERS_DIR:-${DATA_DIR}/posebuster}"
export POSEBUSTERS_BOLTZ2_DIR="${POSEBUSTERS_BOLTZ2_DIR:-${DATA_DIR}/posebuster}"

# --- external tools (override on HPC / local install) ---
export PYTHON="${PYTHON:-python3}"
export BOLTZ_RESULTS_DIR="${BOLTZ_RESULTS_DIR:-}"
export DYNAMICBIND_POSEBUSTERS_LABEL="${DYNAMICBIND_POSEBUSTERS_LABEL:-dynamicbind_new}"

# --- curated 16-target BindingDB panel ---
export TARGETS="${TARGETS:-1g5m 2z5x 3eyg 3jy9 3lxk 3mjg 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 7awe 7kk3}"
export METHODS="${METHODS:-boltz2 dynamicbind gnina plapt qvina}"

# --- pipeline config ---
export TOXDOCK_CONFIG="${TOXDOCK_CONFIG:-${TOXAFFINITY_ROOT}/config/toxdock_config.yaml}"
export METHODS_CONFIG="${METHODS_CONFIG:-${TOXAFFINITY_ROOT}/config/methods_config.yaml}"
