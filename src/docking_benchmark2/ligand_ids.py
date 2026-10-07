"""Stable ligand identifiers assigned before docking.

Canonical form: ligand_0001, ligand_0002, ... (1-based, zero-padded to 4).
Historical engine names from the already-run panel are aliases, not keys:

- QVina / GNINA / PLAPT files: ligand_1 == ligand_0001
- DynamicBind folders: idx_0 == ligand_0001
- Boltz-2 folders: CHEMBL* mapped via molecule_chembl_id in the sidecar
"""

from __future__ import annotations

import re
from typing import Optional

LIGAND_ID_RE = re.compile(r"^ligand_(\d+)$")
IDX_RE = re.compile(r"^idx_(\d+)$")
CHEMBL_RE = re.compile(r"^CHEMBL\d+$", re.IGNORECASE)


def format_ligand_id(one_based: int) -> str:
    if int(one_based) < 1:
        raise ValueError(f"ligand_id is 1-based, got {one_based}")
    return f"ligand_{int(one_based):04d}"


def parse_ligand_number(token: str) -> Optional[int]:
    """Return the 1-based ligand number if token is ligand_N (padded or not)."""
    if token is None:
        return None
    match = LIGAND_ID_RE.fullmatch(str(token).strip())
    if not match:
        return None
    number = int(match.group(1))
    return number if number >= 1 else None


def normalize_ligand_id(token: str) -> Optional[str]:
    number = parse_ligand_number(token)
    return format_ligand_id(number) if number is not None else None


def vina_legacy_id(one_based: int) -> str:
    """Filename stem used by the existing QVina/GNINA/PLAPT panel run."""
    return f"ligand_{int(one_based)}"


def dynamicbind_legacy_id(one_based: int) -> str:
    """0-based idx_N used by the existing DynamicBind panel run."""
    return f"idx_{int(one_based) - 1}"


def parse_dynamicbind_idx(token: str) -> Optional[int]:
    """Return 1-based ligand number from idx_N (0-based)."""
    if token is None:
        return None
    match = IDX_RE.fullmatch(str(token).strip())
    if not match:
        return None
    return int(match.group(1)) + 1


def is_chembl_id(token: str) -> bool:
    return bool(token) and CHEMBL_RE.fullmatch(str(token).strip()) is not None
