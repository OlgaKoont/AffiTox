#!/usr/bin/env bash
# Path contract for the example/ smoke-test pipeline.
# Source from example/run_example_pipeline.sh (not from production scripts).

_cfg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOXAFFINITY_ROOT="$(cd "${_cfg_dir}/.." && pwd)"
export TOXAFFINITY_ROOT
unset _cfg_dir

EXAMPLE_ROOT="${TOXAFFINITY_ROOT}/example"
export EXAMPLE_ROOT

# --- data tree (mirrors repo layout, rooted at example/) ---
export INPUT_DIR="${EXAMPLE_ROOT}/input"
export PROCESSED_DIR="${EXAMPLE_ROOT}/processed"
export RESULTS_DIR="${EXAMPLE_ROOT}/results"
export ANALYSIS_ROOT="${EXAMPLE_ROOT}/analysis"
export DATA_DIR="${ANALYSIS_ROOT}/tables"
export MERGED_DATA_DIR="${DATA_DIR}"
export POSEBUSTERS_DIR="${DATA_DIR}/posebuster"
export POSEBUSTERS_BOLTZ2_DIR="${POSEBUSTERS_DIR}"

# Merge reads ligands from <MERGE_BASE_DIR>/input/ligands_nodubl/
export MERGE_BASE_DIR="${EXAMPLE_ROOT}"

# --- pipeline configs (code stays in repo root) ---
export TOXDOCK_CONFIG="${TOXAFFINITY_ROOT}/config/toxdock_config.example.yaml"
export METHODS_CONFIG="${METHODS_CONFIG:-${TOXAFFINITY_ROOT}/config/methods_config.yaml}"

# Example dock writes to docking/dynamicbind/ (production uses dynamicbind_new).
export DYNAMICBIND_POSEBUSTERS_LABEL="${DYNAMICBIND_POSEBUSTERS_LABEL:-dynamicbind}"

# Per-method conda envs and tool binaries (see example/run_example_pipeline.sh).
# shellcheck source=example_toolchain.env.sh
source "${TOXAFFINITY_ROOT}/config/example_toolchain.env.sh"

# Boltz-2 bridge (see config/boltz_example.env.sh + example/scripts/run_boltz2_example.sh).
# shellcheck source=boltz_example.env.sh
source "${TOXAFFINITY_ROOT}/config/boltz_example.env.sh"

# --- 16-target panel ---
export TARGETS="${TARGETS:-1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 5lf3 7kk3}"
export METHODS="${METHODS:-boltz2 dynamicbind gnina plapt qvina}"

# --- fast statistics for smoke test (override for full SI run) ---
export N_BOOTSTRAP="${N_BOOTSTRAP:-200}"
export N_PERMUTATION="${N_PERMUTATION:-200}"
export RANDOM_SEED="${RANDOM_SEED:-42}"
export FIGURE_DPI="${FIGURE_DPI:-150}"

