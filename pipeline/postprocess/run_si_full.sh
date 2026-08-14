#!/usr/bin/env bash
# Full SI reproducibility run for AffiTox analysis pipeline.
#
# Computes all tables + figures with canonical 10k bootstrap / 10k permutation
# resampling for correlation uncertainty (correlations_with_ci_perm.csv).
# Expected runtime: ~30–60 min (dominated by bootstrap/permutation on 80 pairs).
#
# Usage:
#   bash run_si_full.sh
#   bash run_si_full.sh 2>&1 | tee logs/si_full_$(date +%Y%m%d_%H%M%S).log
#
# Quick dev run (do NOT use for SI):
#   N_BOOTSTRAP=200 N_PERMUTATION=200 bash run_article_analysis.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=analysis/config/defaults.sh
source "${ROOT}/analysis/config/defaults.sh"

LOG_DIR="${ANALYSIS_ROOT}/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/si_full_$(date +%Y%m%d_%H%M%S).log"

export N_BOOTSTRAP=10000
export N_PERMUTATION=10000
export RANDOM_SEED="${RANDOM_SEED:-42}"

export RUN_CORRELATIONS=1
export RUN_ENRICHMENT=1
export RUN_POSEBUSTERS=1
export RUN_INFERENTIAL=1
export RUN_FIGURES=1

{
  echo "================================================================"
  echo "AffiTox SI full run"
  echo "  started : $(date -Iseconds)"
  echo "  root    : ${ANALYSIS_ROOT}"
  echo "  merged  : ${MERGED_DATA_DIR}"
  echo "  N_BOOTSTRAP=${N_BOOTSTRAP}  N_PERMUTATION=${N_PERMUTATION}"
  echo "================================================================"
} | tee "${LOG_FILE}"

bash "${SCRIPT_DIR}/run_article_analysis.sh" 2>&1 | tee -a "${LOG_FILE}"

{
  echo "================================================================"
  echo "Finished: $(date -Iseconds)"
  echo "Log: ${LOG_FILE}"
  echo "Key outputs:"
  echo "  tables/correlations/correlations_with_ci_perm.csv"
  echo "  tables/correlations/summary_all_proteins.csv"
  echo "  tables/enrichment/ef_summary_all_proteins.csv"
  echo "  tables/posebusters/pass_rates_by_target_method.csv"
  echo "  tables/inferential/pairwise_method_tests.csv"
  echo "  figures/"
  echo "================================================================"
} | tee -a "${LOG_FILE}"
