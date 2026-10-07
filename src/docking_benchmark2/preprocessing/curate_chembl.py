#!/usr/bin/env python3
"""ChEMBL Ki snapshot → curated AffiTox ligand table (layer 1 → 2).

Filter contract (public recipe):
  1. Keep experimentally measured Ki only (standard_type matching Ki).
  2. Convert M / mM / µM / nM / pM → nM; drop other units.
  3. Deduplicate by molecule_chembl_id (else SMILES); median Ki if replicates.
  4. Drop RDKit MW > 500 Da.
  5. Require parseable SMILES.
  6. Assign ligand_id in remaining row order (ligand_0001, …).
  7. pKi = 9 - log10(Ki [nM]).
  8. is_active if Ki < 1000 nM; activity_class high/medium/low as in add_pvalue.

The already-docked panel used input/ligands_nodubl row order. To keep those
poses attached to the same compounds, default ``--from-frozen-nodubl`` copies
that order and only adds ligand_id / pKi / flags. Use ``--from-raw`` to apply
the filters above to the ChEMBL API snapshot.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import pandas as pd

# Allow running as a file or as a module.
_SRC = Path(__file__).resolve().parents[2]
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from docking_benchmark2.ligand_ids import format_ligand_id  # noqa: E402

try:
    from analysis.panel import (  # noqa: E402
        CANONICAL_TARGETS,
        PDB_TO_NODUBL_STEM,
        PDB_TO_RAW_CSV,
        nodubl_csv_name,
        project_root,
    )
except ImportError:
    from analysis.panel import (  # type: ignore  # noqa: E402
        CANONICAL_TARGETS,
        PDB_TO_NODUBL_STEM,
        PDB_TO_RAW_CSV,
        nodubl_csv_name,
        project_root,
    )

KI_TYPE_RE = re.compile(r"^\s*ki(\b|[\s_\-].*)$", re.IGNORECASE)

# Multiply source value to obtain nM.
UNIT_TO_NM = {
    "nm": 1.0,
    "nanomolar": 1.0,
    "um": 1.0e3,
    "µm": 1.0e3,
    "μm": 1.0e3,
    "micromolar": 1.0e3,
    "mm": 1.0e6,
    "millimolar": 1.0e6,
    "m": 1.0e9,
    "molar": 1.0e9,
    "pm": 1.0e-3,
    "picomolar": 1.0e-3,
}


def _norm_unit(text: object) -> str:
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return ""
    raw = str(text).strip().lower().replace("μ", "µ")
    raw = raw.replace("uo_0000064", "m").replace("uo_0000063", "mm")
    raw = raw.replace("uo_0000062", "um").replace("uo_0000065", "nm")
    return raw


def ki_to_nm(value: object, units: object, standard_units: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number) or number <= 0:
        return None
    for candidate in (standard_units, units):
        key = _norm_unit(candidate)
        if key in UNIT_TO_NM:
            return number * UNIT_TO_NM[key]
    return None


def is_ki_type(text: object) -> bool:
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return False
    return KI_TYPE_RE.match(str(text).strip()) is not None


def rdkit_mw(smiles: str) -> float | None:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError as exc:
        raise RuntimeError("RDKit is required for MW filtering") from exc
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return float(Descriptors.MolWt(mol))


def activity_class(ki_nm: float) -> str:
    if ki_nm <= 100.0:
        return "high"
    if ki_nm < 1000.0:
        return "medium"
    return "low"


def _read_semicolon_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=";", low_memory=False)


def curate_from_raw(raw_df: pd.DataFrame) -> pd.DataFrame:
    df = raw_df.copy()
    type_col = "standard_type" if "standard_type" in df.columns else "type"
    df = df[df[type_col].map(is_ki_type)].copy()

    value_col = "standard_value" if "standard_value" in df.columns else "value"
    units_col = "units" if "units" in df.columns else None
    std_units_col = "standard_units" if "standard_units" in df.columns else None
    ki_nm = [
        ki_to_nm(
            row.get(value_col),
            row.get(units_col) if units_col else None,
            row.get(std_units_col) if std_units_col else None,
        )
        for _, row in df.iterrows()
    ]
    df["ki_nm"] = ki_nm
    df = df[df["ki_nm"].notna()].copy()

    smiles_col = "canonical_smiles" if "canonical_smiles" in df.columns else "smiles"
    df = df[df[smiles_col].notna() & (df[smiles_col].astype(str).str.strip() != "")].copy()
    mw = df[smiles_col].map(rdkit_mw)
    df["mol_wt"] = mw
    df = df[df["mol_wt"].notna() & (df["mol_wt"] <= 500.0)].copy()

    id_col = "molecule_chembl_id" if "molecule_chembl_id" in df.columns else None
    group_key = id_col if id_col else smiles_col
    agg = {col: "first" for col in df.columns if col not in (group_key, "ki_nm")}
    agg["ki_nm"] = "median"
    df = df.groupby(group_key, as_index=False, dropna=False).agg(agg)

    df = annotate_curated(df, ki_col="ki_nm", smiles_col=smiles_col)
    return df


def annotate_curated(
    df: pd.DataFrame,
    *,
    ki_col: str = "standard_value",
    smiles_col: str = "canonical_smiles",
) -> pd.DataFrame:
    out = df.copy().reset_index(drop=True)
    out.insert(0, "ligand_id", [format_ligand_id(i) for i in range(1, len(out) + 1)])
    ki = pd.to_numeric(out[ki_col], errors="coerce")
    out["standard_value"] = ki
    out["pKi"] = 9.0 - ki.map(lambda x: math.log10(x) if pd.notna(x) and x > 0 else float("nan"))
    out["pValue"] = out["pKi"]
    out["is_active"] = ki < 1000.0
    out["activity_class"] = ki.map(lambda x: activity_class(float(x)) if pd.notna(x) else None)
    if smiles_col != "canonical_smiles" and smiles_col in out.columns:
        out["canonical_smiles"] = out[smiles_col]
    return out


def curate_from_frozen_nodubl(path: Path) -> pd.DataFrame:
    df = _read_semicolon_csv(path)
    return annotate_curated(df, ki_col="standard_value")


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("frozen-nodubl", "from-raw"),
        default="frozen-nodubl",
        help="frozen-nodubl: keep docked table order. from-raw: apply Ki/MW filters.",
    )
    parser.add_argument("--base-dir", type=Path, default=None)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=None,
        help="Raw ChEMBL CSVs (default: <base>/input/ligands)",
    )
    parser.add_argument(
        "--nodubl-dir",
        type=Path,
        default=None,
        help="Frozen nodubl CSVs (default: <base>/input/ligands_nodubl)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Curated tables with ligand_id (default: <base>/input/ligands_curated)",
    )
    parser.add_argument(
        "--targets",
        nargs="*",
        default=list(CANONICAL_TARGETS),
    )
    args = parser.parse_args()
    base = args.base_dir or project_root()
    raw_dir = args.raw_dir or (base / "input" / "ligands")
    nodubl_dir = args.nodubl_dir or (base / "input" / "ligands_nodubl")
    out_dir = args.output_dir or (base / "input" / "ligands_curated")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for pdb in args.targets:
        pdb = pdb.lower()
        if args.mode == "frozen-nodubl":
            src = nodubl_dir / nodubl_csv_name(pdb)
            if not src.is_file():
                print(f"MISSING nodubl {pdb}: {src}", file=sys.stderr)
                continue
            curated = curate_from_frozen_nodubl(src)
            source = str(src)
        else:
            src = raw_dir / PDB_TO_RAW_CSV[pdb]
            if not src.is_file():
                print(f"MISSING raw {pdb}: {src}", file=sys.stderr)
                continue
            curated = curate_from_raw(_read_semicolon_csv(src))
            source = str(src)
        dest = out_dir / f"{pdb}_ligands.csv"
        write_csv(curated, dest)
        print(f"{pdb}\t{len(curated)}\t{dest}")
        rows.append(
            {
                "pdb_id": pdb,
                "n_ligands": len(curated),
                "source": source,
                "output": str(dest),
                "nodubl_stem": PDB_TO_NODUBL_STEM[pdb],
            }
        )
    summary = pd.DataFrame(rows)
    summary_path = out_dir / "curation_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"N_panel={int(summary['n_ligands'].sum())} -> {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
