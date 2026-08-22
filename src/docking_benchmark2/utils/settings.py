"""Settings loading utilities."""

from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union
import json
import os
import yaml

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
_DEFAULT_INTERACTION_JSON = _REPO_ROOT / "config" / "interaction_protein_ligand_16target.json"


def _load_settings_file(path: Optional[Path]) -> Dict[str, Any]:
    """Load a YAML or JSON settings file if it exists."""
    if path is None:
        return {}
    
    resolved = Path(path).expanduser()
    if not resolved.exists():
        return {}
    
    if resolved.suffix.lower() in {".yaml", ".yml"}:
        with resolved.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
    elif resolved.suffix.lower() == ".json":
        with resolved.open("r", encoding="utf-8") as handle:
            data = json.load(handle) or {}
    else:
        raise ValueError(f"Unsupported settings file format: {resolved.suffix}")
    
    if not isinstance(data, dict):
        raise ValueError("Settings file must contain a mapping at the top level")
    return data


def load_protein_settings(path: Optional[Path] = None) -> Dict[str, Any]:
    """Return per-protein preparation overrides."""
    return _load_settings_file(path)


def load_box_settings(path: Optional[Path] = None) -> Dict[str, Any]:
    """Return per-protein box overrides."""
    return _load_settings_file(path)


def _resolve_interaction_path(path: Union[str, Path, None]) -> Path:
    if path is None or str(path).strip() in {"", "null", "None"}:
        return _DEFAULT_INTERACTION_JSON
    resolved = Path(path).expanduser()
    if resolved.is_absolute():
        return resolved
    candidates = []
    env_root = os.environ.get("TOXAFFINITY_ROOT")
    if env_root:
        candidates.append(Path(env_root) / resolved)
    candidates.append(Path.cwd() / resolved)
    candidates.append(_REPO_ROOT / resolved)
    for cand in candidates:
        if cand.exists():
            return cand.resolve()
    return candidates[0].resolve()


def load_interaction_config(path: Optional[Path] = None) -> Dict[str, List[str]]:
    """
    Load protein-ligand interaction configuration.

    Expected format:
    {
        "protein": ["1g5m", "2z5x"],
        "ligand": ["BCL2_Ki_WT_ChEMBL_252_nodubl", "MAO-B_Ki_WT_ChEMBL_246_nodubl"],
        "ref_ligand": ["", ""],
        "safe_chain": ["A", "A"]
    }
    """
    resolved = _resolve_interaction_path(path)
    if not resolved.exists():
        raise FileNotFoundError(
            f"Missing protein-ligand pairing file: {resolved}. "
            "Public checkout uses config/interaction_protein_ligand_16target.json."
        )

    with resolved.open("r", encoding="utf-8") as handle:
        data = json.load(handle) or {}

    if not isinstance(data, dict):
        raise ValueError(f"Interaction config must be a JSON object: {resolved}")
    return data


def get_protein_ligand_pairs(interaction_config: Optional[Dict[str, List[str]]] = None) -> List[Tuple[str, str, Optional[str], Optional[str]]]:
    """Get list of protein-ligand pairs from interaction config."""
    if interaction_config is None:
        interaction_config = load_interaction_config()
    
    if not interaction_config:
        return []
    
    proteins = interaction_config.get("protein", [])
    ligands = interaction_config.get("ligand", [])
    ref_ligands = interaction_config.get("ref_ligand", [])
    safe_chains = interaction_config.get("safe_chain", [])
    
    pairs = []
    for i in range(max(len(proteins), len(ligands))):
        protein = proteins[i] if i < len(proteins) else None
        ligand = ligands[i] if i < len(ligands) else None
        ref_ligand = ref_ligands[i] if i < len(ref_ligands) else None
        safe_chain = safe_chains[i] if i < len(safe_chains) else None
        
        if protein and ligand:
            pairs.append((protein, ligand, ref_ligand, safe_chain))
    
    return pairs


def get_proteins_for_ligand(ligand_dataset: str, interaction_config: Optional[Dict[str, List[str]]] = None) -> List[str]:
    """Get list of proteins associated with a ligand dataset."""
    pairs = get_protein_ligand_pairs(interaction_config)
    normalized_ligand = ligand_dataset.lower()
    return [protein for protein, ligand, _, _ in pairs if ligand.lower() == normalized_ligand]


def get_ligands_for_protein(protein: str, interaction_config: Optional[Dict[str, List[str]]] = None) -> List[str]:
    """Get list of ligand datasets associated with a protein."""
    pairs = get_protein_ligand_pairs(interaction_config)
    normalized_protein = protein.lower()
    return [ligand for prot, ligand, _, _ in pairs if prot.lower() == normalized_protein]

