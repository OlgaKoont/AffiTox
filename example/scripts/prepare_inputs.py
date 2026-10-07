#!/usr/bin/env python3
"""Build minimal example/input dataset (2 ligands per target) and scaffold dirs."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

import pandas as pd

# 16-target BindingDB panel (same mapping as merge / training_overlap_audit).
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


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _write_mini_csv(src: Path, dst: Path, n_rows: int) -> int:
    df = pd.read_csv(src, sep=";")
    mini = df.head(n_rows)
    mini.to_csv(dst, sep=";", index=False)
    return len(mini)


def _write_mini_grouped_json(src: Path, dst: Path, chembl_ids: list[str]) -> None:
    if not src.exists():
        return
    with src.open(encoding="utf-8") as handle:
        data = json.load(handle)
    kept = []
    for group in data.get("groups", []):
        ligands = [
            lig
            for lig in group.get("ligands", [])
            if lig.get("molecule_chembl_id") in chembl_ids
        ]
        if ligands:
            kept.append({**group, "ligands": ligands, "group_size": len(ligands)})
    out = {
        "total_groups": len(kept),
        "total_ligands": len(chembl_ids),
        "groups": kept,
    }
    with dst.open("w", encoding="utf-8") as handle:
        json.dump(out, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def _link_or_copy(src: Path, dst: Path, use_symlinks: bool) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    if use_symlinks:
        os.symlink(src.resolve(), dst)
    else:
        shutil.copy2(src, dst)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=_repo_root(),
        help="AffiTox repository root",
    )
    parser.add_argument(
        "--ligands-per-target",
        type=int,
        default=int(os.environ.get("EXAMPLE_LIGANDS_PER_TARGET", "2")),
        help="Number of first rows to keep from each ligand CSV",
    )
    parser.add_argument(
        "--copy-proteins",
        action="store_true",
        help="Copy protein PDBs instead of symlinking",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove example/processed, example/results, example/analysis before scaffold",
    )
    args = parser.parse_args()

    repo = args.repo_root.resolve()
    example = repo / "example"
    src_ligands = repo / "input" / "ligands_nodubl"
    src_proteins = repo / "input" / "proteins"
    dst_ligands = example / "input" / "ligands_nodubl"
    dst_proteins = example / "input" / "proteins"

    if args.clean:
        for sub in ("processed", "results", "analysis"):
            target = example / sub
            if target.exists():
                shutil.rmtree(target)

    for sub in (
        "input/ligands_nodubl",
        "input/proteins",
        "processed",
        "results",
        "analysis/tables",
        "analysis/figures",
        "logs",
    ):
        (example / sub).mkdir(parents=True, exist_ok=True)

    manifest: dict[str, object] = {
        "ligands_per_target": args.ligands_per_target,
        "targets": {},
    }

    print(f"Writing mini ligands -> {dst_ligands}")
    for pdb, csv_stem in PDB_TO_CSV.items():
        src_csv = src_ligands / f"{csv_stem}_nodubl.csv"
        dst_csv = dst_ligands / f"{csv_stem}_nodubl.csv"
        if not src_csv.exists():
            raise FileNotFoundError(f"Missing source ligand CSV: {src_csv}")

        n_written = _write_mini_csv(src_csv, dst_csv, args.ligands_per_target)
        mini_df = pd.read_csv(dst_csv, sep=";")
        chembl_ids = mini_df["molecule_chembl_id"].astype(str).tolist()

        src_json = src_ligands / f"{csv_stem}_nodubl_grouped.json"
        dst_json = dst_ligands / f"{csv_stem}_nodubl_grouped.json"
        _write_mini_grouped_json(src_json, dst_json, chembl_ids)

        manifest["targets"][pdb] = {
            "csv": dst_csv.name,
            "n_ligands": n_written,
            "chembl_ids": chembl_ids,
        }
        print(f"  {pdb}: {n_written} ligands -> {dst_csv.name}")

    print(f"Linking proteins -> {dst_proteins}")
    for pdb in PDB_TO_CSV:
        src_pdb = src_proteins / f"{pdb}.pdb"
        if not src_pdb.exists():
            raise FileNotFoundError(f"Missing protein PDB: {src_pdb}")
        dst_pdb = dst_proteins / f"{pdb}.pdb"
        _link_or_copy(src_pdb, dst_pdb, use_symlinks=not args.copy_proteins)
        print(f"  {pdb}.pdb")

    manifest_path = example / "input" / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    print(f"Done. Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
