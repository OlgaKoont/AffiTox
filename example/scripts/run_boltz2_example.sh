#!/usr/bin/env bash
# Run Boltz-2 on the example mini panel.
#
# Writes to ${BOLTZ_WRITE_DIR}/{pdb}/docking/${BOLTZ_RUN_TAG}/ (same layout as production).
# Requires GPU node + conda env ${BOLTZ_CONDA_ENV}.
#
# Typical usage (after prepare-inputs + boltz-prepare):
#   source config/project.env.example.sh
#   source config/boltz_example.env.sh
#   bash example/scripts/run_boltz2_example.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=../../config/project.env.example.sh
source "${REPO_ROOT}/config/project.env.example.sh"
# shellcheck source=../../config/boltz_example.env.sh
source "${REPO_ROOT}/config/boltz_example.env.sh"

# Prefer staged example/boltz inputs when present.
STAGING_META="${EXAMPLE_ROOT}/boltz/boltz_prepare.json"
read_boltz_meta() {
  python3 - "$1" "$2" <<'PY'
import json
import sys

print(json.load(open(sys.argv[1], encoding="utf-8"))[sys.argv[2]])
PY
}
if [[ -f "${STAGING_META}" ]]; then
  BOLTZ_LIGAND_DIR="$(read_boltz_meta "${STAGING_META}" ligand_dir)"
  BOLTZ_PROTEIN_CIF_DIR="$(read_boltz_meta "${STAGING_META}" protein_cif_dir)"
  BOLTZ_MSA_DIR="$(read_boltz_meta "${STAGING_META}" msa_dir)"
fi

PROTEINS=(1g5m 2z5x 3eyg 3jy9 3lxk 3mjg 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 7awe 7kk3)
LIGANDS=(
  BCL2_Ki_WT_ChEMBL_252_nodubl.csv
  MAO-B_Ki_WT_ChEMBL_246_nodubl.csv
  JAK1_Ki_WT_ChEMBL_2255_nodubl.csv
  JAK2_Ki_WT_ChEMBL_2027_nodubl.csv
  JAK3_Ki_WT_ChEMBL_786_nodubl.csv
  PDGFRB_Ki_WT_ChEMBL_275_nodubl.csv
  VEGFR2_Ki_WT_ChEMBL_875_nodubl.csv
  FGFR1_Ki_WT_ChEMBL_134_nodubl.csv
  CRBN_Ki_WT_ChEMBL_127_nodubl.csv
  EGFR_Ki_WT_curated_251_nodubl.csv
  CYP19A1_Aromatase_Ki_WT_ChEMBL_548_nodubl.csv
  ABL1_BCR-ABL_Ki_WT_ChEMBL_693_nodubl.csv
  KIT_Ki_WT_curated_1298_nodubl.csv
  PDGFRA_Ki_WT_curated_250_nodubl.csv
  PSMB5_Ki_WT_ChEMBL_88_nodubl.csv
  PARP1_Ki_WT_ChEMBL_1075_nodubl.csv
)
SAFE_CHAINS=(A A A A A B A A C A A A A A L C)

if [[ "${#PROTEINS[@]}" -ne "${#LIGANDS[@]}" || "${#PROTEINS[@]}" -ne "${#SAFE_CHAINS[@]}" ]]; then
  echo "ERROR: protein/ligand/chain array length mismatch" >&2
  exit 1
fi

echo "Boltz-2 example run"
echo "  ligands : ${BOLTZ_LIGAND_DIR}"
echo "  proteins: ${BOLTZ_PROTEIN_CIF_DIR}"
echo "  msa     : ${BOLTZ_MSA_DIR}"
echo "  write   : ${BOLTZ_WRITE_DIR}"
echo "  run tag : ${BOLTZ_RUN_TAG}"
echo "  targets : ${TARGETS_ONLY:-all 16}"

CONDA_SH="${CONDA_BASE:-}/etc/profile.d/conda.sh"
if [[ -f "${CONDA_SH}" ]]; then
  # shellcheck source=/dev/null
  source "${CONDA_SH}"
  conda activate "${BOLTZ_CONDA_ENV}"
fi

export PYTHONDONTWRITEBYTECODE=1
export PYTHONHASHSEED="${BOLTZ_SEED}"

find_cif() {
  local protein="$1"
  for candidate in \
    "${BOLTZ_PROTEIN_CIF_DIR}/${protein}.cif" \
    "${BOLTZ_PROTEIN_CIF_DIR}/${protein}.CIF" \
    "${BOLTZ_PROTEIN_CIF_DIR}/${protein^^}.cif"
  do
    if [[ -f "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

run_one_pair() {
  local protein="$1"
  local chain="$2"
  local ligand_csv_name="$3"

  local ligand_csv="${BOLTZ_LIGAND_DIR}/${ligand_csv_name}"
  local cif_path
  cif_path="$(find_cif "${protein}")" || {
    echo "ERROR: CIF not found for ${protein}" >&2
    return 1
  }

  local msa_path="${BOLTZ_MSA_DIR}/${protein}_chain${chain}.a3m"
  [[ -s "${msa_path}" ]] || { echo "ERROR: MSA missing: ${msa_path}" >&2; return 1; }
  [[ -s "${ligand_csv}" ]] || { echo "ERROR: ligand CSV missing: ${ligand_csv}" >&2; return 1; }

  local out_dir="${BOLTZ_WRITE_DIR}/${protein}/docking/${BOLTZ_RUN_TAG}"
  local yaml_dir="${out_dir}/yaml_inputs/${ligand_csv_name%.csv}"
  local pred_dir="${out_dir}/predictions/${ligand_csv_name%.csv}"
  mkdir -p "${yaml_dir}" "${pred_dir}"

  echo ""
  echo "------------------------------------------"
  echo "Protein: ${protein} chain ${chain}"
  echo "Ligands: ${ligand_csv_name}"
  echo "YAML   : ${yaml_dir}"
  echo "OUT    : ${pred_dir}"
  echo "------------------------------------------"

  python3 - "${protein}" "${chain}" "${cif_path}" "${msa_path}" "${ligand_csv}" "${yaml_dir}" <<'PY'
from pathlib import Path
import csv
import re
import sys

import gemmi

protein, chain_id, cif_path, msa_path, ligand_csv, yaml_dir = map(str, sys.argv[1:])
cif_path = Path(cif_path)
msa_path = Path(msa_path)
ligand_csv = Path(ligand_csv)
yaml_dir = Path(yaml_dir)
yaml_dir.mkdir(parents=True, exist_ok=True)

st = gemmi.read_structure(str(cif_path))
st.setup_entities()
model = st[0]
chains = [c.name for c in model]
if chain_id not in chains:
    raise RuntimeError(f"Chain '{chain_id}' not found in {cif_path}. Available: {chains}")
chain = model[chain_id]
resnames = [res.name for res in chain.get_polymer()]
seq = gemmi.one_letter_code(resnames).replace("?", "X")
if len(seq) < 20:
    raise RuntimeError(f"Suspicious sequence length for {protein}: {len(seq)}")

with ligand_csv.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
    reader = csv.DictReader(handle, delimiter=";")
    if reader.fieldnames is None:
        raise RuntimeError(f"No header in {ligand_csv}")
    smiles_col = next(
        (c for c in reader.fieldnames if c.lower() in {"canonical_smiles", "smiles"}),
        None,
    )
    id_col = next(
        (c for c in reader.fieldnames if c.lower() in {"molecule_chembl_id", "chembl_id"}),
        None,
    )
    if smiles_col is None:
        raise RuntimeError(f"No SMILES column in {ligand_csv}")

    n_written = 0
    for row in reader:
        smi = (row.get(smiles_col) or "").strip()
        if not smi:
            continue
        n_written += 1
        lig_id = str(n_written)
        raw_name = (row.get(id_col) or "").strip() if id_col else ""
        stem = raw_name if raw_name else f"lig{n_written}"
        stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", stem)[:80]
        yml = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {seq}
      msa: {msa_path.as_posix()}
  - ligand:
      id: "{lig_id}"
      smiles: '{smi}'
templates:
  - cif: {cif_path.as_posix()}
    chain_id: {chain_id}
properties:
  - affinity:
      binder: "{lig_id}"
"""
        (yaml_dir / f"{stem}.yaml").write_text(yml)

if n_written == 0:
    raise RuntimeError(f"No SMILES parsed from {ligand_csv}")
print(f"Generated {n_written} YAML files in {yaml_dir}")
PY

  if ! boltz predict "${yaml_dir}" \
    --model boltz2 \
    --out_dir "${pred_dir}" \
    --accelerator gpu \
    --devices 1 \
    --seed "${BOLTZ_SEED}" \
    --sampling_steps "${BOLTZ_SAMPLING_STEPS}" \
    --diffusion_samples "${BOLTZ_DIFFUSION_SAMPLES}" \
    --sampling_steps_affinity "${BOLTZ_SAMPLING_STEPS_AFFINITY}" \
    --diffusion_samples_affinity "${BOLTZ_DIFFUSION_SAMPLES_AFFINITY}" \
    --affinity_mw_correction \
    --max_parallel_samples 1 \
    --num_workers 2
  then
    echo "WARN: boltz predict failed for ${protein}/${ligand_csv_name}" >&2
    return 0
  fi

  echo "OK: ${protein}/${ligand_csv_name}"
}

TARGETS_ONLY="${TARGETS_ONLY:-}"
for i in "${!PROTEINS[@]}"; do
  protein="${PROTEINS[$i]}"
  if [[ -n "${TARGETS_ONLY}" ]] && ! echo " ${TARGETS_ONLY} " | grep -qF " ${protein} "; then
    echo "Skip (not in TARGETS_ONLY): ${protein}"
    continue
  fi
  run_one_pair "${PROTEINS[$i]}" "${SAFE_CHAINS[$i]}" "${LIGANDS[$i]}"
done

echo ""
echo "Boltz-2 example run finished."
echo "Results: ${BOLTZ_WRITE_DIR}"
echo "Next: bash example/run_example_pipeline.sh boltz-sync"
