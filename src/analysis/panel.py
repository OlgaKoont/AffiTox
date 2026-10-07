"""Canonical AffiTox 16-target panel (public release).

Do not add 7AWE, 2Z5X, or 3MJG. Canonical analysis root is
analysis/excluding_2z5x_3mjg/.
"""

from __future__ import annotations

from pathlib import Path

CANONICAL_TARGETS: tuple[str, ...] = (
    "1g5m",
    "2v5z",
    "3eyg",
    "3jy9",
    "3lxk",
    "11ue",
    "4ase",
    "4f65",
    "4tz4",
    "4zau",
    "5jkv",
    "5mo4",
    "6gqj",
    "6jok",
    "5lf3",
    "7kk3",
)

CANONICAL_METHODS: tuple[str, ...] = (
    "boltz2",
    "dynamicbind",
    "gnina",
    "plapt",
    "qvina",
)

EXAMPLE_LIGAND_COUNT = 5
RANDOM_SEED = 42

# Frozen curated tables (row order = ligand_0001...). Shared by 7awe/5lf3 etc.
PDB_TO_NODUBL_STEM: dict[str, str] = {
    "1g5m": "BCL2_Ki_WT_ChEMBL_252",
    "2v5z": "MAO-B_Ki_WT_ChEMBL_246",
    "3eyg": "JAK1_Ki_WT_ChEMBL_2255",
    "3jy9": "JAK2_Ki_WT_ChEMBL_2027",
    "3lxk": "JAK3_Ki_WT_ChEMBL_786",
    "11ue": "PDGFRB_Ki_WT_ChEMBL_275",
    "4ase": "VEGFR2_Ki_WT_ChEMBL_875",
    "4f65": "FGFR1_Ki_WT_ChEMBL_134",
    "4tz4": "CRBN_Ki_WT_ChEMBL_127",
    "4zau": "EGFR_Ki_WT_curated_251",
    "5jkv": "CYP19A1_Aromatase_Ki_WT_ChEMBL_548",
    "5mo4": "ABL1_BCR-ABL_Ki_WT_ChEMBL_693",
    "6gqj": "KIT_Ki_WT_curated_1298",
    "6jok": "PDGFRA_Ki_WT_curated_250",
    "5lf3": "PSMB5_Ki_WT_ChEMBL_88",
    "7kk3": "PARP1_Ki_WT_ChEMBL_1075",
}

# Raw ChEMBL API snapshots (semicolon CSV).
PDB_TO_RAW_CSV: dict[str, str] = {
    "1g5m": "BCL2_Ki_WT_ChEMBL.csv",
    "2v5z": "MAO-B_Ki_WT_ChEMBL.csv",
    "3eyg": "JAK1_Ki_WT_ChEMBL.csv",
    "3jy9": "JAK2_Ki_WT_ChEMBL.csv",
    "3lxk": "JAK3_Ki_WT_ChEMBL.csv",
    "11ue": "PDGFRB_Ki_WT_ChEMBL.csv",
    "4ase": "VEGFR2_Ki_WT_ChEMBL.csv",
    "4f65": "FGFR1_Ki_WT_ChEMBL.csv",
    "4tz4": "CRBN_Ki_WT_ChEMBL.csv",
    "4zau": "EGFR_Ki_WT_curated.csv",
    "5jkv": "CYP19A1_Aromatase_Ki_WT_ChEMBL.csv",
    "5mo4": "ABL1_BCR-ABL_Ki_WT_ChEMBL.csv",
    "6gqj": "KIT_Ki_WT_curated.csv",
    "6jok": "PDGFRA_Ki_WT_curated.csv",
    "5lf3": "PSMB5_Ki_WT_ChEMBL.csv",
    "7kk3": "PARP1_Ki_WT_ChEMBL.csv",
}

# Legacy PDB codes that share a ligand table with the canonical panel.
LEGACY_PDB_TO_NODUBL_STEM: dict[str, str] = {
    **PDB_TO_NODUBL_STEM,
    "2z5x": "MAO-B_Ki_WT_ChEMBL_246",
    "3mjg": "PDGFRB_Ki_WT_ChEMBL_275",
    "7awe": "PSMB5_Ki_WT_ChEMBL_88",
}


def nodubl_csv_name(pdb_id: str) -> str:
    stem = LEGACY_PDB_TO_NODUBL_STEM[pdb_id.lower()]
    return f"{stem}_nodubl.csv"


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent.parent
