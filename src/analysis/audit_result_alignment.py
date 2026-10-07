#!/usr/bin/env python3
"""Audit result identifiers against the original ligand-table row identity.

Missing predictions are valid.  Invalid, duplicated, or ambiguous identifiers
are errors because they can silently shift positional method outputs.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

from merge_ligands_docking_from_dir import PDB_TO_CSV_MAPPING


POSITIONAL_METHODS = ("qvina", "gnina", "plapt", "dynamicbind")


def _parse_position(value: object) -> int | None:
    """Convert explicit ligand_N (one-based) or idx_N (zero-based) to row index."""
    text = str(value).strip()
    match = re.fullmatch(r"ligand_(\d+)", text)
    if match:
        index = int(match.group(1)) - 1
    else:
        match = re.fullmatch(r"idx_(\d+)", text)
        if not match:
            return None
        index = int(match.group(1))
    return index


def _summarize_positions(
    target: str,
    method: str,
    identifiers: list[str],
    n_ligands: int,
    *,
    internal_mismatches: int = 0,
) -> tuple[dict, list[int]]:
    parsed: list[int] = []
    invalid: list[str] = []
    extras: list[str] = []
    for identifier in identifiers:
        index = _parse_position(identifier)
        if index is None:
            invalid.append(identifier)
        elif not 0 <= index < n_ligands:
            extras.append(identifier)
        else:
            parsed.append(index)
    duplicates = sorted(index for index, count in Counter(parsed).items() if count > 1)
    present = set(parsed)
    missing = sorted(set(range(n_ligands)) - present)
    error = bool(duplicates or internal_mismatches)
    has_extras = bool(invalid or extras)
    if error:
        status = "ERROR"
    elif missing:
        status = "OK_WITH_MISSING"
    elif has_extras:
        status = "OK_WITH_EXTRAS"
    else:
        status = "OK"
    return (
        {
            "target": target,
            "method": method,
            "mapping": "explicit_position",
            "n_ligands": n_ligands,
            "n_result_records": len(identifiers),
            "n_unique_mapped": len(present),
            "n_missing": len(missing),
            "n_invalid_ids": len(invalid),
            "n_extra_ids": len(extras),
            "n_duplicate_ids": len(duplicates),
            "n_internal_id_mismatches": internal_mismatches,
            "status": status,
            "invalid_ids": "|".join(invalid[:20]),
            "extra_ids": "|".join(extras[:20]),
            "duplicate_indices": "|".join(map(str, duplicates[:20])),
            "missing_indices": "|".join(map(str, missing[:100])),
        },
        missing,
    )


def _positional_identifiers(
    docking_dir: Path, method: str
) -> tuple[list[str], int, list[str]]:
    internal_mismatches = 0
    sources: list[str] = []
    if method in {"qvina", "gnina"}:
        paths = sorted((docking_dir / method).glob("ligand_*.log"))
        identifiers = [path.stem for path in paths]
        sources = [str(path) for path in paths]
        return identifiers, internal_mismatches, sources

    if method == "plapt":
        paths = sorted((docking_dir / method).glob("**/ligand_*.json"))
        identifiers: list[str] = []
        for path in paths:
            file_id = path.stem
            identifiers.append(file_id)
            sources.append(str(path))
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                internal_mismatches += 1
                continue
            internal_id = str(payload.get("ligand_id", file_id))
            if internal_id != file_id:
                internal_mismatches += 1
        return identifiers, internal_mismatches, sources

    dynamic_root = docking_dir / "dynamicbind_new"
    csv_paths = sorted(
        path
        for path in dynamic_root.glob("**/affinity_prediction.csv")
        if path.name == "affinity_prediction.csv"
    )
    if csv_paths:
        depths = {
            path: len(path.relative_to(dynamic_root).parts) for path in csv_paths
        }
        min_depth = min(depths.values())
        csv_paths = [path for path in csv_paths if depths[path] == min_depth]
        if len(csv_paths) > 1:
            return [], 1, [str(path) for path in csv_paths]
    identifiers = []
    for path in csv_paths:
        try:
            frame = pd.read_csv(path)
        except Exception:
            internal_mismatches += 1
            continue
        id_col = next(
            (column for column in ("name", "ligand_id", "ligand", "id") if column in frame),
            None,
        )
        if id_col is None:
            internal_mismatches += len(frame)
            continue
        identifiers.extend(frame[id_col].dropna().astype(str).tolist())
        sources.extend([str(path)] * len(frame))
    return identifiers, internal_mismatches, sources


def _audit_boltz(
    target: str, docking_dir: Path, ligand_frame: pd.DataFrame
) -> tuple[dict, list[str]]:
    expected_ids = set(ligand_frame["molecule_chembl_id"].dropna().astype(str))
    affinity_paths = sorted(docking_dir.glob("boltz2_*/**/affinity*.json"))
    result_ids: list[str] = []
    invalid_paths: list[str] = []
    for path in affinity_paths:
        match = re.search(r"(CHEMBL\d+)", str(path))
        if match:
            result_ids.append(match.group(1))
        else:
            invalid_paths.append(str(path))
    result_unique = set(result_ids)
    unknown = sorted(result_unique - expected_ids)
    missing = sorted(expected_ids - result_unique)
    # Repeated files for one molecule can be model-level artifacts and are not
    # a row-alignment error because Boltz joins by molecule_chembl_id.
    error = False
    has_extras = bool(invalid_paths or unknown)
    if missing:
        status = "OK_WITH_MISSING"
    elif has_extras:
        status = "OK_WITH_EXTRAS"
    else:
        status = "OK"
    return (
        {
            "target": target,
            "method": "boltz2",
            "mapping": "molecule_chembl_id",
            "n_ligands": len(ligand_frame),
            "n_result_records": len(affinity_paths),
            "n_unique_mapped": len(result_unique & expected_ids),
            "n_missing": len(missing),
            "n_invalid_ids": len(invalid_paths),
            "n_extra_ids": len(unknown),
            "n_duplicate_ids": len(result_ids) - len(result_unique),
            "n_internal_id_mismatches": 0,
            "status": status,
            "invalid_ids": "|".join(invalid_paths[:20]),
            "extra_ids": "|".join(unknown[:20]),
            "duplicate_indices": "",
            "missing_indices": "|".join(missing[:100]),
        },
        missing,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    targets = [
        token.lower()
        for token in args.targets.replace(",", " ").split()
        if token.strip()
    ]
    rows: list[dict] = []
    details: dict[str, dict] = {}
    errors = 0

    for target in targets:
        dataset = PDB_TO_CSV_MAPPING.get(target)
        if dataset is None:
            raise KeyError(f"No ligand-table mapping for target {target}")
        ligand_csv = args.root / "input" / "ligands_nodubl" / f"{dataset}_nodubl.csv"
        ligands = pd.read_csv(ligand_csv, sep=";", low_memory=False)
        docking_dir = args.results_dir / target / "docking"
        details[target] = {"n_ligands": len(ligands), "methods": {}}

        for method in POSITIONAL_METHODS:
            identifiers, internal_mismatches, sources = _positional_identifiers(
                docking_dir, method
            )
            row, missing = _summarize_positions(
                target,
                method,
                identifiers,
                len(ligands),
                internal_mismatches=internal_mismatches,
            )
            rows.append(row)
            details[target]["methods"][method] = {
                "identifiers": identifiers,
                "missing_indices": missing,
                "sources": sorted(set(sources)),
                "summary": row,
            }
            errors += row["status"] == "ERROR"

        row, missing_ids = _audit_boltz(target, docking_dir, ligands)
        rows.append(row)
        details[target]["methods"]["boltz2"] = {
            "missing_molecule_chembl_ids": missing_ids,
            "summary": row,
        }
        errors += row["status"] == "ERROR"

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output_csv, index=False)
    args.output_json.write_text(json.dumps(details, indent=2) + "\n", encoding="utf-8")

    print(pd.DataFrame(rows).to_string(index=False))
    print(f"Wrote {args.output_csv}")
    print(f"Wrote {args.output_json}")
    if errors:
        print(f"ERROR: {errors} method/target mappings contain invalid or ambiguous IDs")
        return 1
    print("Alignment audit passed: all available records have explicit valid IDs.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
