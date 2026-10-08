#!/usr/bin/env bash
# Copy pinned checkpoints into this folder before Zenodo upload.
set -euo pipefail
DEST="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BOLTZ_DIR="${BOLTZ_WEIGHTS_DIR:-$HOME/.boltz}"
DB_DIR="${DYNAMICBIND_WEIGHTS_DIR:-/mnt/tank/scratch/okonovalova/DynamicBind/workdir/big_score_model_sanyueqi_with_time}"
PLAPT="${PLAPT_ONNX:-/mnt/tank/scratch/okonovalova/WELP-PLAPT/models/affinity_predictor.onnx}"

cp -n "${BOLTZ_DIR}/boltz2_conf.ckpt" "${DEST}/" || true
cp -n "${BOLTZ_DIR}/boltz2_aff.ckpt" "${DEST}/" || true
cp -n "${DB_DIR}/ema_inference_epoch314_model.pt" "${DEST}/" || true
cp -n "${DB_DIR}/pro_ema_inference_epoch138_model.pt" "${DEST}/" || true
cp -n "${PLAPT}" "${DEST}/" || true
( cd "${DEST}" && sha256sum *.ckpt *.pt *.onnx README.md pack_weights.sh > SHA256SUMS.txt 2>/dev/null || true )
ls -lh "${DEST}"
