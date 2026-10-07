#!/usr/bin/env python3
"""Sync Boltz-2 outputs from external boltz/data into example/results/."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

# 16-target panel: pdb -> ligand csv stem (without _nodubl.csv)
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


def _detect_run_tag(docking_dir: Path) -> str | None:
    if not docking_dir.is_dir():
        return None
    candidates = sorted(
        p.name
        for p in docking_dir.iterdir()
        if p.is_dir() and p.name.startswith("boltz2_") and "aff" in p.name
    )
    return candidates[-1] if candidates else None


def _link_or_copy(src: Path, dst: Path, mode: str) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if dst.is_dir() and not dst.is_symlink():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    if mode == "copy":
        if src.is_dir():
            shutil.copytree(src, dst, symlinks=True)
        else:
            shutil.copy2(src, dst)
    else:
        os.symlink(src.resolve(), dst)


def _chembl_prediction_dirs(src_run: Path, ligand_dataset: str, chembl_id: str) -> list[Path]:
    """Return existing source directories for one CHEMBL id (old + new layouts)."""
    found: list[Path] = []
    old = (
        src_run
        / "predictions"
        / ligand_dataset
        / f"boltz_results_{ligand_dataset}"
        / "predictions"
        / chembl_id
    )
    if old.is_dir():
        found.append(old)
    new = src_run / f"boltz_results_{chembl_id}" / "predictions" / chembl_id
    if new.is_dir():
        found.append(new)
    return found


def sync_target(
    protein: str,
    source_root: Path,
    dest_root: Path,
    run_tag: str | None,
    mode: str,
    ligand_dataset: str | None,
    chembl_ids: list[str] | None,
) -> dict[str, str | bool | int]:
    src_docking = source_root / protein / "docking"
    tag = run_tag or _detect_run_tag(src_docking)
    if not tag:
        return {"protein": protein, "synced": False, "reason": "no boltz2_*_aff_* run found"}

    src_run = src_docking / tag
    if not src_run.is_dir():
        return {"protein": protein, "synced": False, "reason": f"missing {src_run}"}

    dst_run = dest_root / protein / "docking" / tag
    if dst_run.exists() or dst_run.is_symlink():
        if dst_run.is_dir() and not dst_run.is_symlink():
            shutil.rmtree(dst_run)
        else:
            dst_run.unlink()

    info: dict[str, str | bool | int] = {
        "protein": protein,
        "run_tag": tag,
        "source": str(src_run),
        "dest": str(dst_run),
    }

    if not chembl_ids:
        _link_or_copy(src_run, dst_run, mode)
        info.update({"synced": True, "mode": "full"})
        if ligand_dataset:
            src_pred = src_run / "predictions" / ligand_dataset
            if src_pred.is_dir():
                info["predictions"] = str(src_pred)
            else:
                info["predictions_missing"] = str(src_pred)
        return info

    linked = 0
    missing: list[str] = []
    for chembl_id in chembl_ids:
        src_dirs = _chembl_prediction_dirs(src_run, ligand_dataset or "", chembl_id)
        if not src_dirs:
            missing.append(chembl_id)
            continue
        for src_chembl in src_dirs:
            rel = src_chembl.relative_to(src_run)
            dst_chembl = dst_run / rel
            _link_or_copy(src_chembl, dst_chembl, mode)
            linked += 1

    if linked == 0:
        return {
            "protein": protein,
            "synced": False,
            "reason": f"no manifest ligands found under {src_run}",
            "missing": ",".join(missing),
        }

    info.update(
        {
            "synced": True,
            "mode": "manifest",
            "linked_paths": linked,
            "n_chembl": len(chembl_ids),
            "missing_chembl": ",".join(missing),
        }
    )
    if ligand_dataset:
        info["ligand_dataset"] = ligand_dataset
    return info


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        type=Path,
        default=Path(os.environ.get("BOLTZ_WRITE_DIR", "")),
        help="External Boltz results root (e.g. /path/boltz/data/results)",
    )
    parser.add_argument(
        "--dest-root",
        type=Path,
        default=Path(os.environ.get("RESULTS_DIR", "")),
        help="Unified AffiTox results root (e.g. example/results)",
    )
    parser.add_argument(
        "--run-tag",
        type=str,
        default=os.environ.get("BOLTZ_RUN_TAG", ""),
        help="Boltz run subdirectory (default: auto-detect latest *aff*)",
    )
    parser.add_argument(
        "--mode",
        choices=("symlink", "copy"),
        default=os.environ.get("BOLTZ_SYNC_MODE", "symlink"),
    )
    parser.add_argument(
        "--targets",
        nargs="*",
        default=os.environ.get("TARGETS", "").split(),
        help="PDB ids to sync (default: 16-target panel)",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="example/input/manifest.json — sync only chembl_ids listed per target",
    )
    parser.add_argument(
        "--no-filter",
        action="store_true",
        help="Symlink/copy full boltz run dir (ignore manifest chembl_ids)",
    )
    args = parser.parse_args()

    if not args.source_root:
        raise SystemExit("--source-root or BOLTZ_WRITE_DIR is required")
    if not args.dest_root:
        raise SystemExit("--dest-root or RESULTS_DIR is required")

    manifest_targets: dict = {}
    if args.manifest and args.manifest.is_file():
        with args.manifest.open(encoding="utf-8") as handle:
            manifest_targets = json.load(handle).get("targets", {})

    targets = args.targets or list(PDB_TO_CSV.keys())
    run_tag = args.run_tag or None

    print(f"Sync Boltz-2: {args.source_root} -> {args.dest_root} ({args.mode})")
    if args.manifest and args.manifest.is_file() and not args.no_filter:
        print(f"  filter: manifest chembl_ids from {args.manifest}")
    elif args.no_filter:
        print("  filter: disabled (--no-filter)")
    if run_tag:
        print(f"  run tag: {run_tag}")
    else:
        print("  run tag: auto-detect per target")

    synced = 0
    for protein in targets:
        pdb = protein.lower()
        csv_stem = PDB_TO_CSV.get(pdb)
        ligand_dataset = f"{csv_stem}_nodubl" if csv_stem else None
        chembl_ids: list[str] | None = None
        if pdb in manifest_targets:
            entry = manifest_targets[pdb]
            ligand_dataset = entry.get("csv", "").replace(".csv", "") or ligand_dataset
            if not args.no_filter:
                chembl_ids = entry.get("chembl_ids") or None

        result = sync_target(
            pdb,
            args.source_root,
            args.dest_root,
            run_tag,
            args.mode,
            ligand_dataset,
            chembl_ids,
        )
        if result.get("synced"):
            synced += 1
            extra = ""
            if result.get("mode") == "manifest":
                extra = f" ({result.get('linked_paths', 0)} paths, missing: {result.get('missing_chembl') or 'none'})"
            print(f"  OK {pdb}: {result['run_tag']} -> {result['dest']}{extra}")
        else:
            print(f"  SKIP {pdb}: {result.get('reason', 'unknown')}")

    print(f"Synced {synced}/{len(targets)} targets.")


if __name__ == "__main__":
    main()
