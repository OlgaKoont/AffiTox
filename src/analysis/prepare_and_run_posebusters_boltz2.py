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
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from .posebusters_runner import run_posebusters_chunked
    from .prepare_and_run_posebusters import PDB_TO_LIGAND_CSV
    from .pose_chemistry import pdb_to_template_mol, write_pose_molecule
except ImportError:
    from posebusters_runner import run_posebusters_chunked
    from prepare_and_run_posebusters import PDB_TO_LIGAND_CSV
    from pose_chemistry import pdb_to_template_mol, write_pose_molecule


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_BOLTZ_RESULTS_DIR = Path(os.environ.get("BOLTZ_RESULTS_DIR", ""))
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "analysis" / "tables" / "posebuster"
DEFAULT_POSEBUSTERS_CONFIG = os.environ.get("POSEBUSTERS_CONFIG", "")


def _short_chain_name(idx: int) -> str:
    """Return a PDB-compatible short chain identifier."""
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    return alphabet[idx % len(alphabet)]


def _extract_protein_from_path(
    cif_path: Path, boltz_root: Optional[Path] = None
) -> Optional[str]:
    """Extract protein/pdb id from path segment .../results/<protein>/docking/..."""
    m = re.search(r"/results/([^/]+)/docking/", str(cif_path))
    if m:
        return m.group(1).lower()
    if boltz_root is not None:
        try:
            relative = cif_path.relative_to(boltz_root)
        except ValueError:
            return None
        parts = relative.parts
        if len(parts) >= 2 and parts[1] == "docking":
            return parts[0].lower()
    return None


def _extract_ligand_from_name(cif_path: Path) -> str:
    """CHEMBL535_model_0.cif -> CHEMBL535"""
    return cif_path.stem.replace("_model_0", "")


def _is_valid_file(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def _discover_model_cifs(boltz_root: Path) -> List[Path]:
    """Discover CIFs, including method directories linked into results/."""
    candidates = set(boltz_root.rglob("*_model_0.cif"))
    for docking_dir in boltz_root.glob("*/docking"):
        for method_dir in docking_dir.glob("boltz2*"):
            if method_dir.is_dir():
                candidates.update(method_dir.rglob("*_model_0.cif"))
    return sorted(candidates)


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


def _load_smiles_by_molecule_id(protein: str) -> Dict[str, str]:
    csv_name = PDB_TO_LIGAND_CSV.get(protein.lower())
    if not csv_name:
        return {}
    csv_path = PROJECT_ROOT / "input" / "ligands_nodubl" / csv_name
    if not csv_path.is_file():
        return {}
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle, delimiter=";")
        return {
            row["molecule_chembl_id"]: row["canonical_smiles"]
            for row in rows
            if row.get("molecule_chembl_id") and row.get("canonical_smiles")
        }


def _convert_pdb_to_sdf(
    pdb_path: Path, sdf_path: Path, canonical_smiles: str
) -> bool:
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    sdf_path.unlink(missing_ok=True)
    try:
        fixed = pdb_to_template_mol(pdb_path, canonical_smiles)
        write_pose_molecule(fixed, sdf_path)
    except Exception as exc:
        print(f"Template chemistry conversion failed [{pdb_path}]: {exc}")
        sdf_path.unlink(missing_ok=True)
        return False
    return _is_valid_file(sdf_path)


def _run_posebusters(
    table_csv: Path,
    output_csv: Path,
    bust_bin: str,
    max_workers: int,
    posebusters_config: Optional[str],
    chunk_size: int,
) -> None:
    run_posebusters_chunked(
        table_csv=table_csv,
        output_csv=output_csv,
        bust_bin=bust_bin,
        max_workers=max_workers,
        posebusters_config=posebusters_config,
        chunk_size=chunk_size,
    )


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
        "--chunk-size",
        type=int,
        default=100,
        help="Rows per resumable PoseBusters chunk; use 0 to disable chunking.",
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

    cif_paths = _discover_model_cifs(boltz_root)
    rows_by_protein: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    failures = 0

    print(f"Scanning CIF files under: {boltz_root}")
    print(f"Found model_0 CIF files: {len(cif_paths)}")

    # Keep one prediction per protein+ligand if copied/nested result trees overlap.
    cif_candidates: Dict[Tuple[str, str], Path] = {}
    for cif_path in cif_paths:
        protein = _extract_protein_from_path(cif_path, boltz_root)
        if protein is None:
            failures += 1
            continue
        if proteins_filter and protein not in proteins_filter:
            continue
        ligand = _extract_ligand_from_name(cif_path)
        key = (protein, ligand)
        previous = cif_candidates.get(key)
        if previous is None or len(cif_path.relative_to(boltz_root).parts) < len(
            previous.relative_to(boltz_root).parts
        ):
            cif_candidates[key] = cif_path

    smiles_cache: Dict[str, Dict[str, str]] = {}
    for (protein, ligand), cif_path in sorted(cif_candidates.items()):
        protein_pdb = split_dir / protein / method / "protein_pdb" / f"{ligand}.pdb"
        ligand_pdb = split_dir / protein / method / "ligand_pdb" / f"{ligand}.pdb"
        ligand_sdf = split_dir / protein / method / "ligand_sdf" / f"{ligand}.sdf"

        if args.force or not (_is_valid_file(protein_pdb) and _is_valid_file(ligand_pdb)):
            ok_split = _split_cif_to_protein_ligand_pdb(cif_path, protein_pdb, ligand_pdb)
            if not ok_split:
                failures += 1
                continue

        if protein not in smiles_cache:
            smiles_cache[protein] = _load_smiles_by_molecule_id(protein)
        canonical_smiles = smiles_cache[protein].get(ligand)
        if not canonical_smiles:
            failures += 1
            continue
        if args.force or not _is_valid_file(ligand_sdf):
            ok_sdf = _convert_pdb_to_sdf(
                ligand_pdb, ligand_sdf, canonical_smiles
            )
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
            chunk_size=args.chunk_size,
        )
        print(f"PoseBusters [{protein.upper()} {method.upper()}]: {result_path}")

    if args.prepare_only:
        print(f"Prepare-only complete. Protein+method tables created: {created}")
    else:
        print(f"PoseBusters complete. Protein+method result CSVs created: {created}")


if __name__ == "__main__":
    main()

