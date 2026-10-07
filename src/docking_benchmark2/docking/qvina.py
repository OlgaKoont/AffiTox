"""QVina docking - simplified function-based approach."""

import json
import random
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from ..utils.env_utils import run_in_env
from ..utils.settings import load_interaction_config, get_protein_ligand_pairs


def dock_qvina(
    processed_dir: Path,
    output_dir: Path,
    config: Dict,
    interaction_config: Optional[Dict] = None,
) -> None:
    """
    Run QVina docking for all protein-ligand pairs.
    
    Results are saved in output_dir/{protein_name}/docking/qvina/
    
    Args:
        processed_dir: Directory with prepared proteins and ligands.
        output_dir: Directory to save docking results (organized by protein).
        config: QVina configuration (binary, exhaustiveness, etc.).
        interaction_config: Optional interaction config for protein-ligand pairs.
    """
    binary = config.get('binary', 'qvina02')
    exhaustiveness = config.get('exhaustiveness', 8)
    num_modes = config.get('num_modes', 9)
    docking_timeout = config.get('docking_timeout', None)
    random_seed = config.get('random_seed', 42)  # Fixed seed for reproducibility
    energy_range = config.get('energy_range', 3)  # Energy range for output
    verbose_energy = config.get('verbose_energy', False)  # Request detailed energy breakdown
    docking_env = config.get("conda_env") or None
    
    # Set random seed for reproducibility (QVina doesn't support seed directly,
    # but we set Python random and numpy random for any internal randomness)
    random.seed(random_seed)
    np.random.seed(random_seed)
    
    # Load interaction config
    if interaction_config is None:
        interaction_config = load_interaction_config()
    
    pairs = get_protein_ligand_pairs(interaction_config)
    if not pairs:
        print("  Warning: No protein-ligand pairs found. Skipping docking.")
        return
    
    print(f"  Running QVina docking for {len(pairs)} protein-ligand pairs...")
    
    # Directories
    protein_pdbqt_dir = processed_dir / "proteins"
    box_dir = processed_dir / "boxes"
    ligand_pdbqt_base = processed_dir / "ligands"
    
    docked_count = 0
    
    for protein_name, ligand_dataset, ref_ligand, safe_chain in pairs:
        normalized_protein_name = protein_name.lower()
        
        # Find protein PDBQT
        protein_pdbqt = protein_pdbqt_dir / f"{normalized_protein_name}.pdbqt"
        if not protein_pdbqt.exists():
            for pdbqt_file in protein_pdbqt_dir.glob("*.pdbqt"):
                if pdbqt_file.stem.lower() == normalized_protein_name:
                    protein_pdbqt = pdbqt_file
                    break
        
        if not protein_pdbqt.exists():
            print(f"    Warning: Protein {protein_name} not found, skipping")
            continue
        
        # Create results directory: output_dir/{protein_name}/docking/qvina/
        protein_results_dir = output_dir / protein_name / "docking" / "qvina"
        protein_results_dir.mkdir(parents=True, exist_ok=True)
        
        # Load box info
        box_file = box_dir / f"{normalized_protein_name}.json"
        box_info = None
        if box_file.exists():
            with open(box_file, 'r') as f:
                box_info = json.load(f)
        else:
            box_info = {"center": [0, 0, 0], "size": [20, 20, 20]}
        
        center = box_info['center']
        size = box_info['size']
        
        # Find ligands
        ligand_pdbqt_dir = ligand_pdbqt_base / normalized_protein_name / ligand_dataset
        if not ligand_pdbqt_dir.exists():
            print(f"    Warning: No ligands found for {protein_name}/{ligand_dataset}, skipping")
            continue
        
        # Dock each ligand
        for ligand_pdbqt in ligand_pdbqt_dir.glob("*.pdbqt"):
            ligand_name = ligand_pdbqt.stem
            output_pdbqt = protein_results_dir / f"{ligand_name}_out.pdbqt"
            log_file = protein_results_dir / f"{ligand_name}.log"
            
            if output_pdbqt.exists():
                continue
            
            cmd = [
                binary,
                '--receptor', str(protein_pdbqt),
                '--ligand', str(ligand_pdbqt),
                '--center_x', str(center[0]),
                '--center_y', str(center[1]),
                '--center_z', str(center[2]),
                '--size_x', str(size[0]),
                '--size_y', str(size[1]),
                '--size_z', str(size[2]),
                '--out', str(output_pdbqt),
                '--log', str(log_file),
                '--exhaustiveness', str(exhaustiveness),
                '--seed', str(random_seed),
                '--num_modes', str(num_modes),
                '--energy_range', str(energy_range),  # Show energy range for detailed analysis
            ]
            
            # If verbose_energy is enabled, increase energy_range to show more poses with energy details
            if verbose_energy:
                cmd[-1] = str(max(energy_range, 10))  # At least 10 kcal/mol range to see more detail
            
            try:
                start_time = time.time()
                run_in_env(
                    cmd,
                    env_name=docking_env,
                    check=True,
                    capture_output=True,
                    timeout=docking_timeout
                )
                elapsed_time = time.time() - start_time
                docked_count += 1
                print(f"    Docked: {protein_name}/{ligand_name} (took {elapsed_time:.1f}s)")
            except Exception as e:
                print(f"    Error docking {protein_name}/{ligand_name}: {e}")
                continue
    
    print(f"  QVina docking complete: {docked_count} dockings performed")


def _extract_affinity(log_file: Path) -> Optional[Dict]:
    """Extract affinity from QVina log file."""
    if not log_file.exists():
        return None
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        affinity_data = {}
        for i, line in enumerate(lines):
            if '-----+' in line and i + 1 < len(lines):
                for j in range(i + 1, min(i + 10, len(lines))):
                    parts = lines[j].split()
                    if len(parts) >= 2:
                        try:
                            mode = int(parts[0])
                            affinity = float(parts[1])
                            if mode == 1:
                                affinity_data['affinity'] = affinity
                                affinity_data['best_mode'] = mode
                                break
                        except (ValueError, IndexError):
                            continue
                break
        
        return affinity_data if affinity_data else None
    except Exception:
        return None


def extract_metrics_qvina(output_dir: Path) -> List[Dict]:
    """Extract metrics from QVina docking results."""
    metrics = []
    
    for protein_dir in output_dir.iterdir():
        if not protein_dir.is_dir() or protein_dir.name == 'global':
            continue
        
        protein_name = protein_dir.name
        qvina_dir = protein_dir / "docking" / "qvina"
        
        if not qvina_dir.exists():
            continue
        
        for log_file in qvina_dir.glob("*.log"):
            ligand_name = log_file.stem
            output_pdbqt = qvina_dir / f"{ligand_name}_out.pdbqt"
            
            affinity_data = _extract_affinity(log_file)
            
            metric = {
                'method': 'qvina',
                'protein': protein_name,
                'ligand': ligand_name,
                'output_file': str(output_pdbqt) if output_pdbqt.exists() else None,
            }
            
            if affinity_data:
                if 'affinity' in affinity_data:
                    metric['affinity'] = affinity_data['affinity']
                if 'best_mode' in affinity_data:
                    metric['best_mode'] = affinity_data['best_mode']
            
            metrics.append(metric)
    
    return metrics

