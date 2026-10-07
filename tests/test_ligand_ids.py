"""Stable ligand_id contract."""

from docking_benchmark2.ligand_ids import (
    dynamicbind_legacy_id,
    format_ligand_id,
    normalize_ligand_id,
    parse_dynamicbind_idx,
    vina_legacy_id,
)


def test_padding_and_legacy_aliases():
    assert format_ligand_id(1) == "ligand_0001"
    assert format_ligand_id(88) == "ligand_0088"
    assert normalize_ligand_id("ligand_1") == "ligand_0001"
    assert normalize_ligand_id("ligand_0001") == "ligand_0001"
    assert vina_legacy_id(1) == "ligand_1"
    assert dynamicbind_legacy_id(1) == "idx_0"
    assert parse_dynamicbind_idx("idx_0") == 1
    assert parse_dynamicbind_idx("idx_4") == 5
