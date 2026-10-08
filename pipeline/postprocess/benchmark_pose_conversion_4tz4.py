#!/usr/bin/env python3
"""Benchmark template-based SDF vs MOL conversion on five 4tz4 poses/method."""

from __future__ import annotations

import csv
import statistics
import subprocess
import sys
from pathlib import Path
from time import perf_counter

import pandas as pd


ROOT = Path("/mnt/tank/scratch/okonovalova/AffiTox")
sys.path.insert(0, str(ROOT / "src"))

from analysis.pose_chemistry import (
    connectivity_key,
    load_pose_molecule,
    pdb_to_template_mol,
    pdbqt_to_template_mol,
    write_pose_molecule,
)
from analysis.prepare_and_run_posebusters_boltz2 import (
    _extract_ligand_from_name,
    _split_cif_to_protein_ligand_pdb,
)


TARGET = "4tz4"
N_PER_METHOD = 5
DOCKING = ROOT / "results" / TARGET / "docking"
PROTEIN_PDB = ROOT / "processed" / "proteins" / f"{TARGET}.pdb"
LIGAND_CSV = (
    ROOT / "input" / "ligands_nodubl" / "CRBN_Ki_WT_ChEMBL_127_nodubl.csv"
)
OUT = ROOT / "analysis" / "pose_conversion_benchmark_4tz4"
BUST = Path(
    "/mnt/tank/scratch/okonovalova/miniconda3/envs/posebuster/bin/bust"
)


def _load_templates() -> tuple[dict[str, str], dict[str, str]]:
    with LIGAND_CSV.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    by_ligand = {
        f"ligand_{index}": row["canonical_smiles"]
        for index, row in enumerate(rows, start=1)
    }
    by_chembl = {
        row["molecule_chembl_id"]: row["canonical_smiles"] for row in rows
    }
    return by_ligand, by_chembl


def _coordinate_max_delta(mol_a: object, mol_b: object) -> float:
    import numpy as np

    coords_a = mol_a.GetConformer().GetPositions()
    coords_b = mol_b.GetConformer().GetPositions()
    if coords_a.shape != coords_b.shape:
        return float("inf")
    return float(np.max(np.abs(coords_a - coords_b)))


def _write_and_validate(
    *,
    method: str,
    ligand: str,
    fixed_mol: object,
    graph_seconds: float,
    mol_cond: Path,
) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    records: list[dict[str, object]] = []
    table_rows: list[dict[str, str]] = []
    written: dict[str, object] = {}
    for suffix in (".sdf", ".mol"):
        output = OUT / "poses" / method / f"{ligand}{suffix}"
        started = perf_counter()
        write_pose_molecule(fixed_mol, output)
        write_seconds = perf_counter() - started

        started = perf_counter()
        loaded = load_pose_molecule(output)
        load_seconds = perf_counter() - started
        written[suffix] = loaded
        graph_ok = connectivity_key(loaded) == connectivity_key(fixed_mol)
        records.append(
            {
                "method": method,
                "ligand": ligand,
                "format": suffix.removeprefix("."),
                "graph_seconds": graph_seconds,
                "write_seconds": write_seconds,
                "end_to_end_seconds": graph_seconds + write_seconds,
                "load_seconds": load_seconds,
                "bytes": output.stat().st_size,
                "connectivity_matches": graph_ok,
                "coordinate_max_delta_angstrom": _coordinate_max_delta(
                    fixed_mol, loaded
                ),
            }
        )
        table_rows.append(
            {
                "mol_pred": str(output),
                "mol_cond": str(mol_cond),
                "protein": TARGET.upper(),
                "method": method,
                "ligand": ligand,
                "format": suffix.removeprefix("."),
            }
        )

    sdf_vs_mol_delta = _coordinate_max_delta(written[".sdf"], written[".mol"])
    if sdf_vs_mol_delta > 1e-3:
        raise ValueError(
            f"{method}/{ligand}: SDF/MOL coordinate delta {sdf_vs_mol_delta}"
        )
    return records, table_rows


def _collect_pdbqt(
    method: str, smiles_by_ligand: dict[str, str]
) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    records: list[dict[str, object]] = []
    table_rows: list[dict[str, str]] = []
    candidates = sorted(
        (DOCKING / method).glob("ligand_*_out.pdbqt"),
        key=lambda path: int(path.stem.split("_")[1]),
    )
    for source in candidates:
        ligand = source.stem.removesuffix("_out")
        smiles = smiles_by_ligand.get(ligand)
        if not smiles:
            continue
        try:
            started = perf_counter()
            fixed = pdbqt_to_template_mol(source, smiles)
            graph_seconds = perf_counter() - started
            current_records, current_rows = _write_and_validate(
                method=method,
                ligand=ligand,
                fixed_mol=fixed,
                graph_seconds=graph_seconds,
                mol_cond=PROTEIN_PDB,
            )
        except Exception as exc:
            print(f"SKIP {method}/{ligand}: {exc}", flush=True)
            continue
        records.extend(current_records)
        table_rows.extend(current_rows)
        if len({row["ligand"] for row in table_rows}) == N_PER_METHOD:
            break
    return records, table_rows


def _collect_boltz(
    smiles_by_chembl: dict[str, str],
) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    records: list[dict[str, object]] = []
    table_rows: list[dict[str, str]] = []
    cifs = sorted(
        (DOCKING / "boltz2_s80_d10_seed42_aff_s80_d1").glob(
            "**/*_model_0.cif"
        )
    )
    for cif_path in cifs:
        ligand = _extract_ligand_from_name(cif_path)
        smiles = smiles_by_chembl.get(ligand)
        if not smiles:
            continue
        split_dir = OUT / "boltz_split" / ligand
        protein_pdb = split_dir / "protein.pdb"
        ligand_pdb = split_dir / "ligand.pdb"
        if not _split_cif_to_protein_ligand_pdb(cif_path, protein_pdb, ligand_pdb):
            print(f"SKIP boltz2/{ligand}: CIF split failed", flush=True)
            continue
        try:
            started = perf_counter()
            fixed = pdb_to_template_mol(ligand_pdb, smiles)
            graph_seconds = perf_counter() - started
            current_records, current_rows = _write_and_validate(
                method="boltz2",
                ligand=ligand,
                fixed_mol=fixed,
                graph_seconds=graph_seconds,
                mol_cond=protein_pdb,
            )
        except Exception as exc:
            print(f"SKIP boltz2/{ligand}: {exc}", flush=True)
            continue
        records.extend(current_records)
        table_rows.extend(current_rows)
        if len({row["ligand"] for row in table_rows}) == N_PER_METHOD:
            break
    return records, table_rows


def _run_posebusters(table_rows: list[dict[str, str]]) -> pd.DataFrame:
    result_frames: list[pd.DataFrame] = []
    for method in ("gnina", "qvina", "boltz2"):
        for output_format in ("sdf", "mol"):
            rows = [
                row
                for row in table_rows
                if row["method"] == method and row["format"] == output_format
            ]
            if len(rows) != N_PER_METHOD:
                raise RuntimeError(
                    f"Expected {N_PER_METHOD} rows for {method}/{output_format}, "
                    f"found {len(rows)}"
                )
            table = OUT / "tables" / f"input_{method}_{output_format}.csv"
            result = OUT / "results" / f"posebusters_{method}_{output_format}.csv"
            table.parent.mkdir(parents=True, exist_ok=True)
            result.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(rows).drop(columns=["format"]).to_csv(table, index=False)
            started = perf_counter()
            subprocess.run(
                [
                    str(BUST),
                    "-t",
                    str(table),
                    "--outfmt",
                    "csv",
                    "--output",
                    str(result),
                    "--max-workers",
                    "5",
                ],
                check=True,
            )
            wall_seconds = perf_counter() - started
            frame = pd.read_csv(result)
            frame.insert(0, "format", output_format)
            frame.insert(0, "method_group", method)
            frame["posebusters_wall_seconds"] = wall_seconds
            result_frames.append(frame)
    return pd.concat(result_frames, ignore_index=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    smiles_by_ligand, smiles_by_chembl = _load_templates()
    conversion_records: list[dict[str, object]] = []
    table_rows: list[dict[str, str]] = []

    for method in ("gnina", "qvina"):
        records, rows = _collect_pdbqt(method, smiles_by_ligand)
        conversion_records.extend(records)
        table_rows.extend(rows)
    records, rows = _collect_boltz(smiles_by_chembl)
    conversion_records.extend(records)
    table_rows.extend(rows)

    conversion = pd.DataFrame(conversion_records)
    conversion.to_csv(OUT / "conversion_records.csv", index=False)
    summary = (
        conversion.groupby(["method", "format"], as_index=False)
        .agg(
            n=("ligand", "count"),
            graph_seconds_mean=("graph_seconds", "mean"),
            write_seconds_mean=("write_seconds", "mean"),
            write_seconds_median=("write_seconds", "median"),
            end_to_end_seconds_mean=("end_to_end_seconds", "mean"),
            load_seconds_mean=("load_seconds", "mean"),
            bytes_mean=("bytes", "mean"),
            max_coordinate_delta=("coordinate_max_delta_angstrom", "max"),
            all_connectivity_match=("connectivity_matches", "all"),
        )
    )
    summary.to_csv(OUT / "conversion_summary.csv", index=False)

    posebusters = _run_posebusters(table_rows)
    posebusters.to_csv(OUT / "posebusters_combined.csv", index=False)
    pb_summary_rows: list[dict[str, object]] = []
    for (method, output_format), frame in posebusters.groupby(
        ["method_group", "format"]
    ):
        pb_summary_rows.append(
            {
                "method": method,
                "format": output_format,
                "n": len(frame),
                "wall_seconds": float(frame["posebusters_wall_seconds"].iloc[0]),
                "inchi_convertible_true": int(
                    pd.to_numeric(
                        frame.get("inchi_convertible", pd.Series(dtype=float)),
                        errors="coerce",
                    )
                    .fillna(0)
                    .gt(0)
                    .sum()
                ),
                "internal_energy_nonnull": int(
                    frame.get("internal_energy", pd.Series(dtype=float))
                    .notna()
                    .sum()
                ),
            }
        )
    pd.DataFrame(pb_summary_rows).to_csv(
        OUT / "posebusters_summary.csv", index=False
    )

    print(summary.to_string(index=False))
    print(pd.DataFrame(pb_summary_rows).to_string(index=False))


if __name__ == "__main__":
    main()
