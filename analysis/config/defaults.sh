#!/usr/bin/env bash
# Defaults for AffiTox article analysis pipeline.
# Override any variable before calling run_article_analysis.sh, e.g.:
#   TARGETS="1g5m 3eyg" METHODS="gnina qvina" bash run_article_analysis.sh

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../../config/project.env.sh
source "${_SCRIPT_DIR}/../../config/project.env.sh"
unset _SCRIPT_DIR

# Path migration guard: if stale env vars still point to removed data/, fall back to analysis/tables.
if [[ ! -d "${MERGED_DATA_DIR}" && -d "${ANALYSIS_ROOT}/tables" ]]; then
  MERGED_DATA_DIR="${ANALYSIS_ROOT}/tables"
fi
if [[ ! -d "${POSEBUSTERS_DIR}" && -d "${ANALYSIS_ROOT}/tables/posebuster" ]]; then
  POSEBUSTERS_DIR="${ANALYSIS_ROOT}/tables/posebuster"
fi
if [[ ! -d "${POSEBUSTERS_BOLTZ2_DIR}" && -d "${ANALYSIS_ROOT}/tables/posebuster" ]]; then
  POSEBUSTERS_BOLTZ2_DIR="${ANALYSIS_ROOT}/tables/posebuster"
fi

# --- analysis-specific overrides ---
# project.env.sh sets PYTHON=python3; analysis always needs the docking env unless overridden at submit time.
if [[ -z "${PYTHON:-}" || "${PYTHON}" == "python3" ]]; then
  PYTHON="/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python"
fi
export PYTHON

# --- experimental axis ---
EXP_COL="${EXP_COL:-pValue}"

# --- article colormap (salmon positive / blue negative; see plots/style.py) ---
COLOR_LOW="${COLOR_LOW:-#6584E1}"
COLOR_MID="${COLOR_MID:-#DEDAD7}"
COLOR_HIGH="${COLOR_HIGH:-#E5885F}"

METHOD_COLORS="${METHOD_COLORS:-boltz2:#E5885F,dynamicbind:#FDBBAA,gnina:#AFC6F2,plapt:#DEDAD7,qvina:#6584E1}"

FIGURE_DPI="${FIGURE_DPI:-700}"
FONT_FAMILY="${FONT_FAMILY:-DejaVu Sans}"

N_BOOTSTRAP="${N_BOOTSTRAP:-10000}"
N_PERMUTATION="${N_PERMUTATION:-10000}"
RANDOM_SEED="${RANDOM_SEED:-42}"

RUN_CORRELATIONS="${RUN_CORRELATIONS:-1}"
RUN_ENRICHMENT="${RUN_ENRICHMENT:-1}"
RUN_POSEBUSTERS="${RUN_POSEBUSTERS:-1}"
RUN_INFERENTIAL="${RUN_INFERENTIAL:-1}"
RUN_FIGURES="${RUN_FIGURES:-1}"
