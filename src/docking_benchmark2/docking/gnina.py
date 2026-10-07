"""Gnina docking - simplified function-based approach."""

import json
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional

from ..utils.env_utils import run_in_env
from ..utils.settings import load_interaction_config, get_protein_ligand_pairs


def dock_gnina(
    processed_dir: Path,
    output_dir: Path,
    config: Dict,
    interaction_config: Optional[Dict] = None,
) -> None:
    """
    Run Gnina docking for all protein-ligand pairs.
    
    Results are saved in output_dir/{protein_name}/docking/gnina/
    
    Args:
        processed_dir: Directory with prepared proteins and ligands.
        output_dir: Directory to save docking results (organized by protein).
        config: Gnina configuration (binary, exhaustiveness, etc.).
        interaction_config: Optional interaction config for protein-ligand pairs.
    """
    binary = config.get('binary', 'gnina')
    exhaustiveness = config.get('exhaustiveness', 8)
    num_modes = config.get('num_modes', 9)
    use_cnn = config.get('use_cnn', False)
    # Default timeout: 10 minutes per ligand (600 seconds)
    # This prevents hanging on problematic ligands
    docking_timeout = config.get('docking_timeout', 600)
    random_seed = config.get('random_seed', 42)  # Fixed seed for reproducibility
    docking_env = config.get("conda_env") or None
    
    # Load interaction config
    if interaction_config is None:
        interaction_config = load_interaction_config()
    
    pairs = get_protein_ligand_pairs(interaction_config)
    if not pairs:
        print("  Warning: No protein-ligand pairs found. Skipping docking.")
        return
    
    print(f"  Running Gnina docking for {len(pairs)} protein-ligand pairs...")
    
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
            # Try case-insensitive
            for pdbqt_file in protein_pdbqt_dir.glob("*.pdbqt"):
                if pdbqt_file.stem.lower() == normalized_protein_name:
                    protein_pdbqt = pdbqt_file
                    break
        
        if not protein_pdbqt.exists():
            print(f"    Warning: Protein {protein_name} not found, skipping")
            continue
        
        # Create results directory: output_dir/{protein_name}/docking/gnina/
        protein_results_dir = output_dir / protein_name / "docking" / "gnina"
        protein_results_dir.mkdir(parents=True, exist_ok=True)
        
        # Load box info
        box_file = box_dir / f"{normalized_protein_name}.json"
        box_info = None
        if box_file.exists():
            with open(box_file, 'r') as f:
                box_info = json.load(f)
        else:
            print(f"    Warning: No box info for {protein_name}, using default")
            box_info = {"center": [0, 0, 0], "size": [20, 20, 20]}
        
        # Find ligands for this protein
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
            
            # Build command
            cmd = [
                binary,
                '--receptor', str(protein_pdbqt),
                '--ligand', str(ligand_pdbqt),
                '--out', str(output_pdbqt),
                '--log', str(log_file),
                '--exhaustiveness', str(exhaustiveness),
                '--num_modes', str(num_modes),
                '--seed', str(random_seed),  # Fixed seed for reproducibility
            ]
            
            # Add box parameters
            if box_info:
                center = box_info['center']
                size = box_info['size']
                cmd.extend([
                    '--center_x', str(center[0]),
                    '--center_y', str(center[1]),
                    '--center_z', str(center[2]),
                    '--size_x', str(size[0]),
                    '--size_y', str(size[1]),
                    '--size_z', str(size[2]),
                ])
            
            # AffiTox default (use_cnn=false): --no_gpu only. GNINA 1.3 still emits CNN
            # pose/affinity columns in the log; merge uses gnina_cnn_affinity_bestpose.
            # Optional use_cnn=true adds --cnn before --no_gpu (not used in published runs).
            use_cnn_actual = use_cnn
            if use_cnn:
                cmd.append('--cnn')
            cmd.append('--no_gpu')
            
            try:
                start_time = time.time()
                result = run_in_env(
                    cmd,
                    env_name=docking_env,
                    check=False,  # Don't raise on error, check return code manually
                    capture_output=True,
                    timeout=docking_timeout
                )
                # If failed with --cnn due to CUDA library issues, retry without --cnn
                if result.returncode != 0 and use_cnn_actual:
                    error_msg = result.stderr if isinstance(result.stderr, str) else (result.stderr.decode('utf-8', errors='ignore') if result.stderr else "")
                    if 'libcufft' in error_msg or 'libcublas' in error_msg or 'version' in error_msg.lower():
                        print(f"    Warning: CNN failed due to CUDA library issues, retrying without --cnn for {protein_name}/{ligand_name}")
                        # Remove --cnn flag and retry
                        cmd_no_cnn = [c for c in cmd if c != '--cnn']
                        result = run_in_env(
                            cmd_no_cnn,
                            env_name=docking_env,
                            check=False,
                            capture_output=True,
                            timeout=docking_timeout
                        )
                
                if result.returncode != 0:
                    error_msg = result.stderr if isinstance(result.stderr, str) else (result.stderr.decode('utf-8', errors='ignore') if result.stderr else "Unknown error")
                    print(f"    Error docking {protein_name}/{ligand_name}: exit code {result.returncode}")
                    if error_msg and error_msg != "Unknown error":
                        print(f"      Error details: {error_msg[:200]}")
                    continue
                elapsed_time = time.time() - start_time
                docked_count += 1
                print(f"    Docked: {protein_name}/{ligand_name} (took {elapsed_time:.1f}s)")
            except subprocess.TimeoutExpired:
                print(f"    Timeout docking {protein_name}/{ligand_name} (>{docking_timeout}s)")
                continue
            except Exception as e:
                print(f"    Error docking {protein_name}/{ligand_name}: {e}")
                continue
    
    print(f"  Gnina docking complete: {docked_count} dockings performed")


def _extract_affinity(log_file: Path) -> Optional[Dict]:
    """Extract affinity from Gnina log file."""
    if not log_file.exists():
        return None
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        affinity_data = {}
        for i, line in enumerate(lines):
            if '-----+' in line and i + 1 < len(lines):
                # Parse table
                for j in range(i + 1, min(i + 10, len(lines))):
                    parts = lines[j].split()
                    if len(parts) >= 2:
                        try:
                            mode = int(parts[0])
                            affinity = float(parts[1])
                            if mode == 1:  # Best mode
                                affinity_data['affinity'] = affinity
                                affinity_data['best_mode'] = mode
                                break
                        except (ValueError, IndexError):
                            continue
                break
        
        return affinity_data if affinity_data else None
    except Exception:
        return None


def extract_metrics_gnina(output_dir: Path) -> List[Dict]:
    """
    Extract metrics from Gnina docking results.
    
    Args:
        output_dir: Directory with results organized by protein.
    
    Returns:
        List of metric dictionaries.
    """
    metrics = []
    
    # Iterate over protein directories
    for protein_dir in output_dir.iterdir():
        if not protein_dir.is_dir() or protein_dir.name == 'global':
            continue
        
        protein_name = protein_dir.name
        gnina_dir = protein_dir / "docking" / "gnina"
        
        if not gnina_dir.exists():
            continue
        
        # Process log files
        for log_file in gnina_dir.glob("*.log"):
            ligand_name = log_file.stem
            output_pdbqt = gnina_dir / f"{ligand_name}_out.pdbqt"
            
            affinity_data = _extract_affinity(log_file)
            
            metric = {
                'method': 'gnina',
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

