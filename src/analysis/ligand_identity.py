#!/usr/bin/env python3
"""Sidecar maps: ligand_id ↔ engine-native identifiers."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_SRC = Path(__file__).resolve().parent.parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from analysis.panel import (  # noqa: E402
    CANONICAL_TARGETS,
    EXAMPLE_LIGAND_COUNT,
    nodubl_csv_name,
    project_root,
)
from docking_benchmark2.ligand_ids import (  # noqa: E402
    dynamicbind_legacy_id,
    format_ligand_id,
    vina_legacy_id,
)

MAP_COLUMNS = [
    "pdb_id",
    "ligand_id",
    "molecule_chembl_id",
    "canonical_smiles",
    "csv_row_1based",
    "vina_legacy_id",
    "dynamicbind_legacy_id",
    "boltz_native_id",
    "engine_native_id",
]


def map_rows_from_nodubl(pdb_id: str, nodubl_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(nodubl_csv, sep=";", low_memory=False)
    rows = []
    for i, row in df.iterrows():
        one_based = int(i) + 1
        chembl = row.get("molecule_chembl_id")
        chembl_s = (
            str(chembl).strip()
            if pd.notna(chembl) and str(chembl).strip() not in ("", "nan")
            else ""
        )
        ligand_id = format_ligand_id(one_based)
        smiles = row.get("canonical_smiles", "")
        rows.append(
            {
                "pdb_id": pdb_id.lower(),
                "ligand_id": ligand_id,
                "molecule_chembl_id": chembl_s,
                "canonical_smiles": smiles if pd.notna(smiles) else "",
                "csv_row_1based": one_based,
                "vina_legacy_id": vina_legacy_id(one_based),
                "dynamicbind_legacy_id": dynamicbind_legacy_id(one_based),
                "boltz_native_id": chembl_s,
                "engine_native_id": ligand_id,
            }
        )
    return pd.DataFrame(rows, columns=MAP_COLUMNS)


def chembl_to_ligand_ids(map_df: pd.DataFrame) -> dict[str, list[str]]:
    """One Boltz CHEMBL folder may belong to several ligand_id rows."""
    out: dict[str, list[str]] = {}
    for _, row in map_df.iterrows():
        chembl = row.get("boltz_native_id")
        ligand_id = str(row["ligand_id"])
        if pd.isna(chembl):
            continue
        token = str(chembl).strip()
        if not token:
            continue
        out.setdefault(token, [])
        out.setdefault(token.upper(), [])
        if ligand_id not in out[token]:
            out[token].append(ligand_id)
        if ligand_id not in out[token.upper()]:
            out[token.upper()].append(ligand_id)
    return out


def native_to_ligand_id(map_df: pd.DataFrame) -> dict[str, str]:
    """Lookup from any historical or canonical token to ligand_id."""
    out: dict[str, str] = {}
    for _, row in map_df.iterrows():
        ligand_id = str(row["ligand_id"])
        for col in (
            "ligand_id",
            "vina_legacy_id",
            "dynamicbind_legacy_id",
            "boltz_native_id",
            "engine_native_id",
            "molecule_chembl_id",
        ):
            value = row.get(col)
            if pd.isna(value):
                continue
            token = str(value).strip()
            if not token:
                continue
            out[token] = ligand_id
            out[token.upper()] = ligand_id
    return out


def write_maps(
    base_dir: Path,
    targets: list[str] | None = None,
    example_n: int = EXAMPLE_LIGAND_COUNT,
) -> pd.DataFrame:
    targets = [t.lower() for t in (targets or list(CANONICAL_TARGETS))]
    nodubl_dir = base_dir / "input" / "ligands_nodubl"
    map_dir = base_dir / "processed" / "id_maps"
    map_dir.mkdir(parents=True, exist_ok=True)
    parts = []
    for pdb in targets:
        src = nodubl_dir / nodubl_csv_name(pdb)
        if not src.is_file():
            raise FileNotFoundError(src)
        part = map_rows_from_nodubl(pdb, src)
        part.to_csv(map_dir / f"{pdb}_ligand_id_map.csv", index=False)
        example = part.head(example_n)
        example.to_csv(map_dir / f"{pdb}_example_ligand_id_map.csv", index=False)
        parts.append(part)
        print(f"{pdb}\t{len(part)}\texample={len(example)}")
    all_maps = pd.concat(parts, ignore_index=True)
    all_maps.to_csv(map_dir / "panel_ligand_id_map.csv", index=False)
    print(f"N_panel={len(all_maps)} -> {map_dir / 'panel_ligand_id_map.csv'}")
    return all_maps


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=None)
    parser.add_argument("--targets", nargs="*", default=None)
    parser.add_argument("--example-n", type=int, default=EXAMPLE_LIGAND_COUNT)
    args = parser.parse_args()
    write_maps(args.base_dir or project_root(), args.targets, args.example_n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
