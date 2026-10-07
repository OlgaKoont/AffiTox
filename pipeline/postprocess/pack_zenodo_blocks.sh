#!/usr/bin/env bash
# Zip the six AffiTox information blocks for Zenodo.
# Record 03 exceeds the 50 GiB / record limit (~193 GiB raw), so docking zips
# are grouped into sibling 03a–03e records. Do not reuse 10.5281/zenodo.20825057–067.
set -euo pipefail

ROOT="${TOXAFFINITY_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
ZENODO="${ROOT}/release/zenodo"
OUT="${OUT_DIR:-${ZENODO}/payloads}"
ZIP_LEVEL="${ZIP_LEVEL:-0}"
mkdir -p "${OUT}/03" "${OUT}/logs"

echo "=== pack Zenodo blocks ==="
echo "root ${ROOT}"
echo "out  ${OUT}"
echo "started $(date -Iseconds)"
cd "${ROOT}"

zip_tree() {
  local name="$1"
  local src="$2"
  local dest="${OUT}/${name}.zip"
  rm -f "${dest}"
  (cd "${src}/.." && zip -r -"${ZIP_LEVEL}" "${dest}" "$(basename "${src}")")
  ls -lh "${dest}"
}

zip_tree "01_raw" "${ZENODO}/01_raw"
zip_tree "02_prepared" "${ZENODO}/02_prepared"
zip_tree "04_merged" "${ZENODO}/04_merged"
zip_tree "05_posebusters" "${ZENODO}/05_posebusters"
zip_tree "06_metrics_figures" "${ZENODO}/06_metrics_figures"

# Per-target docking dumps (store, no compression — already binary-heavy).
for pdb in 1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 5lf3 7kk3; do
  dest="${OUT}/03/results_${pdb}.zip"
  rm -f "${dest}"
  if [[ -d "${ROOT}/results/${pdb}" ]]; then
    (cd "${ROOT}" && zip -r -"${ZIP_LEVEL}" "${dest}" "results/${pdb}")
  else
    echo "WARNING: missing results/${pdb}" >&2
  fi
  ls -lh "${dest}" 2>/dev/null || true
done

python3 - "${OUT}" <<'PY'
import json, sys
from pathlib import Path
out = Path(sys.argv[1])
# Pre-measured working sizes (GiB). Re-group if a zip is larger than 45 GiB.
buckets = [
    ("03a", ["3eyg", "4f65", "5lf3", "2v5z"]),
    ("03b", ["3jy9", "3lxk"]),
    ("03c", ["5mo4", "6gqj", "4zau"]),
    ("03d", ["7kk3", "5jkv", "4tz4"]),
    ("03e", ["4ase", "6jok", "11ue", "1g5m"]),
]
manifest = {"records": []}
limit = 45 * 1024**3
for name, pdbs in buckets:
    files = []
    total = 0
    for pdb in pdbs:
        z = out / "03" / f"results_{pdb}.zip"
        size = z.stat().st_size if z.is_file() else 0
        files.append({"pdb": pdb, "zip": str(z.relative_to(out)), "bytes": size})
        total += size
    if total > limit:
        raise SystemExit(f"{name} is {total/1024**3:.1f} GiB > 45 GiB; regroup targets")
    manifest["records"].append({"id": name, "files": files, "bytes": total})
(out / "03_buckets.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps({r["id"]: round(r["bytes"] / 1024**3, 2) for r in manifest["records"]}, indent=2))
PY

echo "finished $(date -Iseconds)"
ls -lh "${OUT}"/*.zip "${OUT}/03"/*.zip | head -40
