#!/usr/bin/env python3
"""Prepare Boltz-2 inputs for the example mini panel."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

PDB_TO_CSV = {
    "1g5m": "BCL2_Ki_WT_ChEMBL_252",
    "2z5x": "MAO-B_Ki_WT_ChEMBL_246",
    "3eyg": "JAK1_Ki_WT_ChEMBL_2255",
    "3jy9": "JAK2_Ki_WT_ChEMBL_2027",
    "3lxk": "JAK3_Ki_WT_ChEMBL_786",
    "3mjg": "PDGFRB_Ki_WT_ChEMBL_275",
    "4ase": "VEGFR2_Ki_WT_ChEMBL_875",
    "4f65": "FGFR1_Ki_WT_ChEMBL_134",
    "4tz4": "CRBN_Ki_WT_ChEMBL_127",
    "4zau": "EGFR_Ki_WT_curated_251",
    "5jkv": "CYP19A1_Aromatase_Ki_WT_ChEMBL_548",
    "5mo4": "ABL1_BCR-ABL_Ki_WT_ChEMBL_693",
    "6gqj": "KIT_Ki_WT_curated_1298",
    "6jok": "PDGFRA_Ki_WT_curated_250",
    "7awe": "PSMB5_Ki_WT_ChEMBL_88",
    "7kk3": "PARP1_Ki_WT_ChEMBL_1075",
}

PDB_TO_CHAIN = {
    "1g5m": "A",
    "2z5x": "A",
    "3eyg": "A",
    "3jy9": "A",
    "3lxk": "A",
    "3mjg": "B",
    "4ase": "A",
    "4f65": "A",
    "4tz4": "C",
    "4zau": "A",
    "5jkv": "A",
    "5mo4": "A",
    "6gqj": "A",
    "6jok": "A",
    "7awe": "L",
    "7kk3": "C",
}


def _link_or_copy(src: Path, dst: Path, copy: bool) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    if copy:
        shutil.copy2(src, dst)
    else:
        os.symlink(src.resolve(), dst)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--example-root", type=Path, required=True)
    parser.add_argument("--boltz-root", type=Path, required=True)
    parser.add_argument("--boltz-protein-cif-dir", type=Path, required=True)
    parser.add_argument("--boltz-msa-dir", type=Path, required=True)
    parser.add_argument(
        "--staging-dir",
        type=Path,
        default=None,
        help="Optional example/boltz staging tree (mirrors boltz/data layout)",
    )
    parser.add_argument(
        "--copy-proteins",
        action="store_true",
        help="Copy CIF files instead of symlinking",
    )
    args = parser.parse_args()

    example_root = args.example_root.resolve()
    example_ligands = example_root / "input" / "ligands_nodubl"
    staging = (args.staging_dir or (example_root / "boltz")).resolve()
    staging_input = staging / "input"
    staging_ligands = staging_input / "ligands_nodubl"
    staging_proteins = staging_input / "proteins"
    staging_msa = staging / "msa_cache" / "precomputed" / "a3m"

    for path in (staging_ligands, staging_proteins, staging_msa):
        path.mkdir(parents=True, exist_ok=True)

    manifest_path = example_root / "input" / "manifest.json"
    targets = list(PDB_TO_CSV.keys())
    if manifest_path.is_file():
        with manifest_path.open(encoding="utf-8") as handle:
            targets = list(json.load(handle).get("targets", PDB_TO_CSV).keys())

    print(f"Staging Boltz inputs under {staging}")
    for pdb in targets:
        csv_stem = PDB_TO_CSV[pdb]
        src_csv = example_ligands / f"{csv_stem}_nodubl.csv"
        dst_csv = staging_ligands / f"{csv_stem}_nodubl.csv"
        if not src_csv.is_file():
            raise FileNotFoundError(f"Missing example ligand CSV: {src_csv}")
        _link_or_copy(src_csv, dst_csv, copy=False)
        print(f"  ligand: {dst_csv.name}")

        chain = PDB_TO_CHAIN[pdb]
        for ext in (".cif", ".CIF"):
            src_cif = args.boltz_protein_cif_dir / f"{pdb}{ext}"
            if src_cif.is_file():
                dst_cif = staging_proteins / f"{pdb}.cif"
                _link_or_copy(src_cif, dst_cif, copy=args.copy_proteins)
                print(f"  protein: {dst_cif.name}")
                break
        else:
            raise FileNotFoundError(f"No CIF for {pdb} in {args.boltz_protein_cif_dir}")

        src_msa = args.boltz_msa_dir / f"{pdb}_chain{chain}.a3m"
        dst_msa = staging_msa / f"{pdb}_chain{chain}.a3m"
        if not src_msa.is_file():
            raise FileNotFoundError(f"Missing precomputed MSA: {src_msa}")
        _link_or_copy(src_msa, dst_msa, copy=False)
        print(f"  msa: {dst_msa.name}")

    meta = {
        "staging_root": str(staging),
        "ligand_dir": str(staging_ligands),
        "protein_cif_dir": str(staging_proteins),
        "msa_dir": str(staging_msa),
        "targets": targets,
        "pdb_to_chain": {p: PDB_TO_CHAIN[p] for p in targets},
    }
    meta_path = staging / "boltz_prepare.json"
    with meta_path.open("w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(f"Wrote {meta_path}")


if __name__ == "__main__":
    main()
