"""Public pipeline config: 16-target pairing and dock preflight."""

from __future__ import annotations

from docking_benchmark2.pipeline import missing_dock_tools
from docking_benchmark2.utils.settings import get_protein_ligand_pairs, load_interaction_config


def test_interaction_json_has_sixteen_panel_pairs() -> None:
    pairs = get_protein_ligand_pairs(load_interaction_config())
    proteins = [p for p, _, _, _ in pairs]
    assert len(pairs) == 16
    assert "1g5m" in proteins
    assert "7kk3" in proteins


def test_missing_qvina_reported_when_binary_absent() -> None:
    missing = missing_dock_tools(
        ["qvina"],
        {"qvina": {"binary": "qvina02_not_installed_affitox"}},
    )
    assert missing
    assert "qvina" in missing[0]
