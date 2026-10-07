#!/usr/bin/env bash
# Smoke-test driver: full AffiTox pipeline into example/ (mini dataset).
#
# Usage:
#   bash example/run_example_pipeline.sh [stage ...]
#
# Stages (default: all):
#   prepare-inputs   2 ligands/target + protein links + dir scaffold
#   install          pip install -e . (skipped when SKIP_INSTALL=1)
#   prepare          protein/ligand/box preparation
#   dock             qvina, gnina, plapt, dynamicbind
#   boltz-prepare    stage mini ligands + CIF + MSA for Boltz-2
#   boltz-run        call external boltz predict (GPU; skipped when BOLTZ_SKIP_RUN=1)
#   boltz-sync       link/copy boltz/data/results -> example/results
#   merge            merged tables + pValue
#   posebusters      PoseBusters pass-rate tables
#   analysis         correlations, enrichment, inferential tests, figures
#   all              all of the above (except install unless SKIP_INSTALL=0)
#
# Examples:
#   bash example/run_example_pipeline.sh prepare-inputs
#   SKIP_INSTALL=1 bash example/run_example_pipeline.sh all
#   CLEAN_EXAMPLE=1 bash example/run_example_pipeline.sh prepare-inputs all
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=config/project.env.example.sh
source "${REPO_ROOT}/config/project.env.example.sh"

POST_PY="${TOXAFFINITY_ROOT}/src/analysis"
ANALYSIS_DEFAULTS="${TOXAFFINITY_ROOT}/analysis/config/defaults.sh"

usage() {
  sed -n '3,20p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
}

run_prepare_inputs() {
  local extra=()
  [[ "${CLEAN_EXAMPLE:-0}" == "1" ]] && extra+=(--clean)
  [[ "${COPY_PROTEINS:-0}" == "1" ]] && extra+=(--copy-proteins)
  "${PYTHON}" "${EXAMPLE_ROOT}/scripts/prepare_inputs.py" \
    --repo-root "${REPO_ROOT}" \
    --ligands-per-target "${EXAMPLE_LIGANDS_PER_TARGET:-2}" \
    "${extra[@]}"
}

run_install() {
  if [[ "${SKIP_INSTALL:-1}" == "1" ]]; then
    echo "=== install (skipped, SKIP_INSTALL=1) ==="
    return 0
  fi
  bash "${REPO_ROOT}/pipeline/install/setup_environment.sh"
}

run_prepare() {
  cd "${REPO_ROOT}"
  export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"
  echo "=== prepare -> ${PROCESSED_DIR} (env: ${CONDA_ENV_DOCKING}) ==="
  run_in_conda_env "${CONDA_ENV_DOCKING}" \
    env PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}" \
    "${PYTHON_DOCKING}" -m docking_benchmark2.cli.run_benchmark \
    --config "${TOXDOCK_CONFIG}" \
    --methods-config "${METHODS_CONFIG}" \
    --stage preparation \
    "$@"
}

_example_dock_methods() {
  local m
  for m in ${METHODS}; do
    [[ "${m}" == "boltz2" ]] && continue
    echo "${m}"
  done
}

_run_dock_method() {
  local method="$1"
  shift
  local env_name="${CONDA_ENV_DOCKING}"
  local py="${PYTHON_DOCKING}"
  case "${method}" in
    qvina)
      env_name="${CONDA_ENV_QVINA}"
      py="${PYTHON_QVINA}"
      ;;
    gnina)
      env_name="${CONDA_ENV_GNINA}"
      py="${PYTHON_GNINA}"
      ;;
    plapt)
      env_name="${CONDA_ENV_PLAPT}"
      py="${PYTHON_PLAPT}"
      ;;
    dynamicbind)
      env_name="${CONDA_ENV_DYNAMICBIND}"
      py="${PYTHON_DYNAMICBIND}"
      ;;
  esac

  local extra_args=()
  if [[ "${method}" == "gnina" && -n "${GNINA_BIN:-}" ]]; then
    extra_args+=(--gnina-binary "${GNINA_BIN}")
  fi

  echo "=== dock ${method} -> ${RESULTS_DIR} (env: ${env_name}) ==="
  run_in_conda_env "${env_name}" \
    env PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}" \
    "${py}" -m docking_benchmark2.cli.run_benchmark \
    --config "${TOXDOCK_CONFIG}" \
    --methods-config "${METHODS_CONFIG}" \
    --stage docking \
    --methods "${method}" \
    "${extra_args[@]}" \
    "$@"
}

run_dock() {
  cd "${REPO_ROOT}"
  export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"
  init_example_conda
  local method
  while IFS= read -r method; do
    [[ -n "${method}" ]] || continue
    _run_dock_method "${method}" "$@"
  done < <(_example_dock_methods)
}

run_boltz_prepare() {
  echo "=== boltz-prepare (staging under ${EXAMPLE_ROOT}/boltz) ==="
  "${PYTHON}" "${EXAMPLE_ROOT}/scripts/prepare_boltz_inputs.py" \
    --example-root "${EXAMPLE_ROOT}" \
    --boltz-root "${BOLTZ_ROOT}" \
    --boltz-protein-cif-dir "${BOLTZ_PROTEIN_CIF_DIR}" \
    --boltz-msa-dir "${BOLTZ_MSA_DIR}"
}

run_boltz_run() {
  if [[ "${BOLTZ_SKIP_RUN:-1}" == "1" ]]; then
    echo "=== boltz-run (skipped, BOLTZ_SKIP_RUN=1) ==="
    echo "    Run manually on GPU: bash example/scripts/run_boltz2_example.sh"
    return 0
  fi
  echo "=== boltz-run -> ${BOLTZ_WRITE_DIR} ==="
  bash "${EXAMPLE_ROOT}/scripts/run_boltz2_example.sh"
}

run_boltz_sync() {
  if [[ "${BOLTZ_SKIP_SYNC:-0}" == "1" ]]; then
    echo "=== boltz-sync (skipped, BOLTZ_SKIP_SYNC=1) ==="
    return 0
  fi
  if [[ "${BOLTZ_USE_EXTERNAL:-0}" == "1" ]]; then
    echo "=== boltz-sync skipped (BOLTZ_USE_EXTERNAL=1, merge reads ${BOLTZ_RESULTS_DIR}) ==="
    return 0
  fi
  echo "=== boltz-sync ${BOLTZ_WRITE_DIR} -> ${RESULTS_DIR} (${BOLTZ_SYNC_MODE}) ==="
  "${PYTHON}" "${EXAMPLE_ROOT}/scripts/boltz_sync_results.py" \
    --source-root "${BOLTZ_WRITE_DIR}" \
    --dest-root "${RESULTS_DIR}" \
    --run-tag "${BOLTZ_RUN_TAG}" \
    --mode "${BOLTZ_SYNC_MODE}" \
    --manifest "${EXAMPLE_ROOT}/input/manifest.json" \
    --targets ${TARGETS}
  export BOLTZ_RESULTS_DIR="${RESULTS_DIR}"
}

run_merge() {
  mkdir -p "${MERGED_DATA_DIR}" "${POSEBUSTERS_DIR}" "${POSEBUSTERS_BOLTZ2_DIR}"
  cd "${REPO_ROOT}"
  echo "=== merge -> ${MERGED_DATA_DIR} (base: ${MERGE_BASE_DIR}) ==="
  boltz_args=()
  if [[ -n "${BOLTZ_RESULTS_DIR}" ]]; then
    boltz_args+=(--boltz-results-dir "${BOLTZ_RESULTS_DIR}")
  fi
  "${PYTHON}" "${POST_PY}/merge_ligands_docking_from_dir.py" \
    --base-dir "${MERGE_BASE_DIR}" \
    --output-dir "${MERGED_DATA_DIR}" \
    --results-dir "${RESULTS_DIR}" \
    --dynamicbind-subdir "${DYNAMICBIND_POSEBUSTERS_LABEL}" \
    "${boltz_args[@]}" \
    "$@"

  echo "=== add pValue -> ${MERGED_DATA_DIR} ==="
  "${PYTHON}" "${POST_PY}/add_pvalue_column.py" \
    --input-dir "${MERGED_DATA_DIR}" \
    --output-dir "${MERGED_DATA_DIR}"
}

run_posebusters() {
  mkdir -p "${POSEBUSTERS_DIR}" "${POSEBUSTERS_BOLTZ2_DIR}"
  cd "${REPO_ROOT}"
  local pb_args=(
    --poses-dir "${RESULTS_DIR}"
    --proteins-dir "${PROCESSED_DIR}/proteins"
    --output-dir "${ANALYSIS_ROOT}/tables/posebuster"
    --obabel-bin "${OBABEL_BIN}"
    --bust-bin "${BUST_BIN}"
  )
  echo "=== PoseBusters (gnina, qvina, dynamicbind) -> ${POSEBUSTERS_DIR} ==="
  echo "    bust=${BUST_BIN}  obabel=${OBABEL_BIN}"
  run_in_conda_env "${CONDA_ENV_POSEBUSTER}" \
    env PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}" \
    "${PYTHON_POSEBUSTER}" "${POST_PY}/prepare_and_run_posebusters.py" \
    "${pb_args[@]}" \
    "$@"
  if [[ -d "${POSEBUSTERS_DIR}/results_by_method" ]]; then
    mv "${POSEBUSTERS_DIR}/results_by_method/"*.csv "${POSEBUSTERS_DIR}/" 2>/dev/null || true
  fi

  if [[ -n "${BOLTZ_RESULTS_DIR}" ]]; then
    echo "=== PoseBusters (Boltz-2) -> ${POSEBUSTERS_BOLTZ2_DIR} ==="
    run_in_conda_env "${CONDA_ENV_POSEBUSTER}" \
      env PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}" \
      "${PYTHON_POSEBUSTER}" "${POST_PY}/prepare_and_run_posebusters_boltz2.py" \
      --boltz-results-dir "${BOLTZ_RESULTS_DIR}" \
      --output-dir "${ANALYSIS_ROOT}/tables/posebuster" \
      --obabel-bin "${OBABEL_BIN}" \
      --bust-bin "${BUST_BIN}" \
      "$@"
    if [[ -d "${POSEBUSTERS_BOLTZ2_DIR}/results_by_method" ]]; then
      mv "${POSEBUSTERS_BOLTZ2_DIR}/results_by_method/"*.csv "${POSEBUSTERS_BOLTZ2_DIR}/" 2>/dev/null || true
    fi
  else
    echo "=== PoseBusters Boltz-2 skipped (BOLTZ_RESULTS_DIR empty) ==="
  fi
}

run_analysis() {
  # shellcheck source=analysis/config/defaults.sh
  source "${ANALYSIS_DEFAULTS}"

  local pipeline="${TOXAFFINITY_ROOT}/src/analysis/run_pipeline.py"
  local args=(
    --analysis-root "${ANALYSIS_ROOT}"
    --merged-dir "${MERGED_DATA_DIR}"
    --posebusters-dir "${POSEBUSTERS_DIR}"
    --posebusters-boltz2-dir "${POSEBUSTERS_BOLTZ2_DIR}"
    --dynamicbind-posebusters-label "${DYNAMICBIND_POSEBUSTERS_LABEL}"
    --targets "${TARGETS}"
    --methods "${METHODS}"
    --exp-col "${EXP_COL}"
    --color-low "${COLOR_LOW}"
    --color-mid "${COLOR_MID}"
    --color-high "${COLOR_HIGH}"
    --method-colors "${METHOD_COLORS}"
    --figure-dpi "${FIGURE_DPI}"
    --font-family "${FONT_FAMILY}"
    --n-bootstrap "${N_BOOTSTRAP}"
    --n-permutation "${N_PERMUTATION}"
    --random-seed "${RANDOM_SEED}"
  )

  [[ "${RUN_CORRELATIONS:-1}" == "0" ]] && args+=(--skip-correlations)
  [[ "${RUN_ENRICHMENT:-1}" == "0" ]] && args+=(--skip-enrichment)
  [[ "${RUN_POSEBUSTERS:-1}" == "0" ]] && args+=(--skip-posebusters)
  [[ "${RUN_INFERENTIAL:-1}" == "0" ]] && args+=(--skip-inferential)
  [[ "${RUN_FIGURES:-1}" == "0" ]] && args+=(--skip-figures)

  echo "=== analysis -> ${ANALYSIS_ROOT} ==="
  "${PYTHON}" "${pipeline}" "${args[@]}"
}

run_stage() {
  case "$1" in
    prepare-inputs) run_prepare_inputs ;;
    install)        run_install ;;
    prepare)        run_prepare ;;
    dock)           run_dock ;;
    boltz-prepare)  run_boltz_prepare ;;
    boltz-run)      run_boltz_run ;;
    boltz-sync)     run_boltz_sync ;;
    merge)          run_merge ;;
    posebusters)    run_posebusters ;;
    analysis)       run_analysis ;;
    all)
      run_prepare_inputs
      run_install
      run_prepare
      run_dock
      run_boltz_prepare
      run_boltz_run
      run_boltz_sync
      run_merge
      run_posebusters
      run_analysis
      ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown stage: $1"; usage; exit 1 ;;
  esac
}

STAGES=("$@")
if [[ ${#STAGES[@]} -eq 0 ]]; then
  STAGES=(all)
fi

echo "AffiTox example pipeline"
echo "  repo root     : ${REPO_ROOT}"
echo "  example root  : ${EXAMPLE_ROOT}"
echo "  stages        : ${STAGES[*]}"
echo "  targets       : ${TARGETS}"
echo "  methods       : ${METHODS}"
echo "  conda base    : ${CONDA_BASE}"
echo "  dock envs     : qvina=${CONDA_ENV_QVINA} gnina=${CONDA_ENV_GNINA} plapt=${CONDA_ENV_PLAPT} dynamicbind=${CONDA_ENV_DYNAMICBIND}"
echo "  posebuster    : ${BUST_BIN}"
echo "  ligands/target: ${EXAMPLE_LIGANDS_PER_TARGET:-2}"
echo "  boltz root     : ${BOLTZ_ROOT}"
echo "  boltz write    : ${BOLTZ_WRITE_DIR}"
echo "  boltz merge dir: ${BOLTZ_RESULTS_DIR}"
echo "  boltz skip run : ${BOLTZ_SKIP_RUN:-1}"
echo ""

for stage in "${STAGES[@]}"; do
  run_stage "${stage}"
done

echo ""
echo "Example pipeline finished."
