#!/usr/bin/env bash
# AffiTox article analysis — main entry point.
#
# Usage:
#   bash run_article_analysis.sh
#
# Override defaults (examples):
#   TARGETS="1g5m 3eyg 6gqj" \
#   METHODS="gnina qvina" \
#   COLOR_LOW="#1a237e" COLOR_MID="#fff9c4" COLOR_HIGH="#b71c1c" \
#   METHOD_COLORS="gnina:#2E7D32,qvina:#EF6C00" \
#   bash run_article_analysis.sh
#
# Skip stages:
#   RUN_FIGURES=0 bash run_article_analysis.sh
#
# Full SI run (10k bootstrap + 10k permutation, all stages, ~30–60 min):
#   bash run_si_full.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=analysis/config/defaults.sh
source "${ROOT}/analysis/config/defaults.sh"

PIPELINE="${TOXAFFINITY_ROOT}/src/analysis/run_pipeline.py"

ARGS=(
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
  --target-label-mode "${TARGET_LABEL_MODE}"
  --figure-variant "${FIGURE_VARIANT}"
  --n-bootstrap "${N_BOOTSTRAP}"
  --n-permutation "${N_PERMUTATION}"
  --random-seed "${RANDOM_SEED}"
)

[[ "${RUN_CORRELATIONS}" == "0" ]] && ARGS+=(--skip-correlations)
[[ "${RUN_ENRICHMENT}" == "0" ]] && ARGS+=(--skip-enrichment)
[[ "${RUN_POSEBUSTERS}" == "0" ]] && ARGS+=(--skip-posebusters)
[[ "${RUN_INFERENTIAL}" == "0" ]] && ARGS+=(--skip-inferential)
[[ "${RUN_FIGURES}" == "0" ]] && ARGS+=(--skip-figures)

echo "Starting analysis -> ${ANALYSIS_ROOT}"
"${PYTHON}" "${PIPELINE}" "${ARGS[@]}"
