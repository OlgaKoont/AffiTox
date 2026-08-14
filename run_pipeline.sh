#!/usr/bin/env bash
# AffiTox end-to-end driver (no dependency on docking-benchmark-2/).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT}/config/project.env.sh"

usage() {
  cat <<EOF
Usage: bash run_pipeline.sh [stage ...]

Stages (default: analysis):
  install      pip install -e .
  prepare      protein/ligand/box preparation
  dock         run docking methods
  merge        merge tables + pValue
  posebusters  PoseBusters pass-rate tables
  analysis     article statistics + figures (src/analysis via pipeline/postprocess)
  all          prepare -> dock -> merge -> posebusters -> analysis

Environment: see config/project.env.sh
EOF
}

STAGES=("$@")
if [[ ${#STAGES[@]} -eq 0 ]]; then
  STAGES=(analysis)
fi

run_stage() {
  case "$1" in
    install)     bash "${ROOT}/pipeline/install/setup_environment.sh" ;;
    prepare)     bash "${ROOT}/pipeline/prepare/run_prepare.sh" ;;
    dock)        bash "${ROOT}/pipeline/dock/run_docking.sh" ;;
    merge)       bash "${ROOT}/pipeline/postprocess/run_merge.sh" ;;
    posebusters) bash "${ROOT}/pipeline/postprocess/run_posebusters.sh" ;;
    analysis)    bash "${ROOT}/pipeline/postprocess/run_article_analysis.sh" ;;
    all)
      bash "${ROOT}/pipeline/prepare/run_prepare.sh"
      bash "${ROOT}/pipeline/dock/run_docking.sh"
      bash "${ROOT}/pipeline/postprocess/run_merge.sh"
      bash "${ROOT}/pipeline/postprocess/run_posebusters.sh"
      bash "${ROOT}/pipeline/postprocess/run_article_analysis.sh"
      ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown stage: $1"; usage; exit 1 ;;
  esac
}

for s in "${STAGES[@]}"; do
  run_stage "$s"
done
