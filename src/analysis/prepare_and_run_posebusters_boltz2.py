#!/usr/bin/env python3
"""
Prepare PoseBusters inputs from Boltz-2 CIF predictions and optionally run checks.

Pipeline:
1) Find `*_model_0.cif` files in Boltz-2 results tree.
2) Split each CIF complex into:
   - protein PDB (polymer residues only)
   - ligand SDF (non-polymer residues only)
3) Build PoseBusters tables and (optionally) run `bust`.

Output layout (protein + method = one CSV):
- <output-dir>/tables_by_method/posebusters_input_<protein>_boltz2.csv
- <output-dir>/results_by_method/posebusters_results_<protein>_boltz2.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_BOLTZ_RESULTS_DIR = Path(os.environ.get("BOLTZ_RESULTS_DIR", ""))
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "analysis" / "tables" / "posebuster"
DEFAULT_POSEBUSTERS_CONFIG = os.environ.get("POSEBUSTERS_CONFIG", "")


def _short_chain_name(idx: int) -> str:
    """Return a PDB-compatible short chain identifier."""
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    return alphabet[idx % len(alphabet)]


def _extract_protein_from_path(cif_path: Path) -> Optional[str]:
    """Extract protein/pdb id from path segment .../results/<protein>/docking/..."""
    m = re.search(r"/results/([^/]+)/docking/", str(cif_path))
    if not m:
        return None
    return m.group(1).lower()


def _extract_ligand_from_name(cif_path: Path) -> str:
    """CHEMBL535_model_0.cif -> CHEMBL535"""
    return cif_path.stem.replace("_model_0", "")


def _is_valid_file(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def _split_cif_to_protein_ligand_pdb(
    cif_path: Path, protein_pdb_path: Path, ligand_pdb_path: Path
) -> bool:
    """Split CIF into polymer-only protein PDB and non-polymer ligand PDB via gemmi."""
    try:
        import gemmi  # type: ignore[import]
    except Exception:
        return False

    try:
        st = gemmi.read_structure(str(cif_path))
        model = st[0]
    except Exception:
        return False

    protein = gemmi.Structure()
    protein.name = st.name
    protein.cell = st.cell
    protein.spacegroup_hm = st.spacegroup_hm
    protein_model = gemmi.Model("1")

    ligand = gemmi.Structure()
    ligand.name = st.name
    ligand.cell = st.cell
    ligand.spacegroup_hm = st.spacegroup_hm
    ligand_model = gemmi.Model("1")

    for idx, chain in enumerate(model):
        chain_name = _short_chain_name(idx)
        p_chain = gemmi.Chain(chain_name)
        l_chain = gemmi.Chain(chain_name)
        for res in chain:
            if res.entity_type == gemmi.EntityType.Polymer:
                p_chain.add_residue(res.clone())
            elif res.entity_type == gemmi.EntityType.NonPolymer:
                l_chain.add_residue(res.clone())
        if len(p_chain) > 0:
            protein_model.add_chain(p_chain)
        if len(l_chain) > 0:
            ligand_model.add_chain(l_chain)

    if len(protein_model) == 0 or len(ligand_model) == 0:
        return False

    protein.add_model(protein_model)
    ligand.add_model(ligand_model)
    protein_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    ligand_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    protein.write_pdb(str(protein_pdb_path))
    ligand.write_pdb(str(ligand_pdb_path))
    return _is_valid_file(protein_pdb_path) and _is_valid_file(ligand_pdb_path)


def _convert_pdb_to_sdf(pdb_path: Path, sdf_path: Path, obabel_bin: str) -> bool:
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [obabel_bin, "-ipdb", str(pdb_path), "-osdf", "-O", str(sdf_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        return False
    return _is_valid_file(sdf_path)


def _run_posebusters(
    table_csv: Path,
    output_csv: Path,
    bust_bin: str,
    max_workers: int,
    posebusters_config: Optional[str],
) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        bust_bin,
        "-t",
        str(table_csv),
        "--outfmt",
        "csv",
        "--output",
        str(output_csv),
        "--max-workers",
        str(max_workers),
    ]
    if posebusters_config:
        cmd.extend(["--config", posebusters_config])
    subprocess.run(cmd, check=True)


def _write_table(rows: List[Dict[str, str]], table_path: Path) -> None:
    table_path.parent.mkdir(parents=True, exist_ok=True)
    with table_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["mol_pred", "mol_cond", "protein", "method", "ligand", "source_cif"]
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare/run PoseBusters for Boltz-2 CIF predictions."
    )
    parser.add_argument(
        "--boltz-results-dir",
        type=str,
        default=str(DEFAULT_BOLTZ_RESULTS_DIR),
        help="Root folder with Boltz-2 results (default: /mnt/.../boltz/data/results).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(DEFAULT_OUTPUT_DIR),
        help="Output directory for split structures, tables and PoseBusters results.",
    )
    parser.add_argument(
        "--proteins",
        type=str,
        default="",
        help="Optional comma-separated proteins (e.g. 1g5m,3mjg).",
    )
    parser.add_argument(
        "--obabel-bin",
        type=str,
        default=os.environ.get("OBABEL_BIN", "obabel"),
        help="Path/name of Open Babel executable.",
    )
    parser.add_argument(
        "--bust-bin",
        type=str,
        default=os.environ.get("BUST_BIN", "bust"),
        help="Path/name of PoseBusters executable.",
    )
    parser.add_argument(
        "--posebusters-config",
        type=str,
        default=DEFAULT_POSEBUSTERS_CONFIG,
        help="PoseBusters config file passed to --config (dock_fast recommended).",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=16,
        help="PoseBusters max workers.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Only prepare split files/tables; skip PoseBusters run.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate split/converted files even if they already exist.",
    )
    args = parser.parse_args()

    boltz_root = Path(args.boltz_results_dir)
    output_dir = Path(args.output_dir)
    split_dir = output_dir / "split_structures"
    tables_dir = output_dir / "tables_by_method"
    results_dir = output_dir / "results_by_method"
    method = "boltz2"

    proteins_filter = set()
    if args.proteins.strip():
        proteins_filter = {p.strip().lower() for p in args.proteins.split(",") if p.strip()}

    cif_paths = list(boltz_root.rglob("*_model_0.cif"))
    rows_by_protein: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    failures = 0

    print(f"Scanning CIF files under: {boltz_root}")
    print(f"Found model_0 CIF files: {len(cif_paths)}")

    for cif_path in cif_paths:
        protein = _extract_protein_from_path(cif_path)
        if protein is None:
            failures += 1
            continue
        if proteins_filter and protein not in proteins_filter:
            continue
        ligand = _extract_ligand_from_name(cif_path)

        protein_pdb = split_dir / protein / method / "protein_pdb" / f"{ligand}.pdb"
        ligand_pdb = split_dir / protein / method / "ligand_pdb" / f"{ligand}.pdb"
        ligand_sdf = split_dir / protein / method / "ligand_sdf" / f"{ligand}.sdf"

        if args.force or not (_is_valid_file(protein_pdb) and _is_valid_file(ligand_pdb)):
            ok_split = _split_cif_to_protein_ligand_pdb(cif_path, protein_pdb, ligand_pdb)
            if not ok_split:
                failures += 1
                continue

        if args.force or not _is_valid_file(ligand_sdf):
            ok_sdf = _convert_pdb_to_sdf(ligand_pdb, ligand_sdf, args.obabel_bin)
            if not ok_sdf:
                failures += 1
                continue

        rows_by_protein[protein].append(
            {
                "mol_pred": str(ligand_sdf),
                "mol_cond": str(protein_pdb),
                "protein": protein.upper(),
                "method": method,
                "ligand": ligand,
                "source_cif": str(cif_path),
            }
        )

    proteins = sorted(rows_by_protein.keys())
    total_rows = sum(len(v) for v in rows_by_protein.values())
    print(f"Prepared rows: {total_rows}")
    print(f"Proteins with rows: {len(proteins)}")
    print(f"Failed entries: {failures}")

    created = 0
    for protein in proteins:
        rows = rows_by_protein[protein]
        if not rows:
            continue
        created += 1
        table_path = tables_dir / f"posebusters_input_{protein}_{method}.csv"
        result_path = results_dir / f"posebusters_results_{protein}_{method}.csv"
        _write_table(rows, table_path)
        print(f"Table [{protein.upper()} {method.upper()}]: {table_path} (rows={len(rows)})")
        if args.prepare_only:
            continue
        _run_posebusters(
            table_csv=table_path,
            output_csv=result_path,
            bust_bin=args.bust_bin,
            max_workers=args.max_workers,
            posebusters_config=(args.posebusters_config.strip() or None),
        )
        print(f"PoseBusters [{protein.upper()} {method.upper()}]: {result_path}")

    if args.prepare_only:
        print(f"Prepare-only complete. Protein+method tables created: {created}")
    else:
        print(f"PoseBusters complete. Protein+method result CSVs created: {created}")


if __name__ == "__main__":
    main()

