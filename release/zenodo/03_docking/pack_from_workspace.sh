#!/usr/bin/env bash
# Zip full docking dumps from the workspace results/ tree.
# GPU re-runs of Boltz-2 / DynamicBind are not bit-identical; this zip is the record.
set -euo pipefail
ROOT="${TOXAFFINITY_ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}"
DEST="${1:-$ROOT/release/zenodo/03_docking/payload}"
mkdir -p "$DEST"
for pdb in 1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 5lf3 7kk3; do
  zip -r -0 "$DEST/results_${pdb}.zip" "$ROOT/results/$pdb" || true
done
