"""DynamicBind docking - simplified function-based approach.

DynamicBind is a deep learning-based method for predicting protein-ligand binding
using diffusion models. It generates ligand poses and predicts binding affinity.

Reference: https://github.com/DeepGraphLearning/DynamicBind
"""

import json
import os
import subprocess
import sys
import time
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from ..utils.settings import load_interaction_config, get_protein_ligand_pairs
from ..utils.csv_utils import load_ligands_from_csv


def _check_dynamicbind_results(dataset_dir: Path, header: str, stdout: str, stderr: str) -> None:
    """
    Check and diagnose DynamicBind results after execution.
    
    Args:
        dataset_dir: Directory where results should be saved
        header: Header used for DynamicBind run
        stdout: Standard output from DynamicBind execution
        stderr: Standard error from DynamicBind execution
    """
    print(f"\n      === DynamicBind Results Diagnostic ===")
    
    # Check expected result directory structure
    header_dir = dataset_dir / header
    expected_result_patterns = [
        f"{dataset_dir}/{header}/index*_idx_*",
        f"{dataset_dir}/{header}/index*",
        f"{header_dir}/index*_idx_*",
        f"{header_dir}/index*",
    ]
    
    # Check if header directory exists
    if header_dir.exists():
        print(f"      ✓ Header directory exists: {header_dir}")
        
        # List all subdirectories
        subdirs = [d for d in header_dir.iterdir() if d.is_dir()]
        print(f"      Found {len(subdirs)} subdirectory(ies) in header dir:")
        for subdir in subdirs[:10]:  # Show first 10
            print(f"        - {subdir.name}")
            # Check what's inside
            files_in_subdir = list(subdir.glob("*"))
            print(f"          Contains {len(files_in_subdir)} item(s)")
            if files_in_subdir:
                for item in files_in_subdir[:5]:  # Show first 5 items
                    item_type = "DIR" if item.is_dir() else "FILE"
                    size = f" ({item.stat().st_size} bytes)" if item.is_file() else ""
                    print(f"            [{item_type}] {item.name}{size}")
    else:
        print(f"      ✗ Header directory NOT found: {header_dir}")
        print(f"      Checking parent directory: {dataset_dir}")
        if dataset_dir.exists():
            items = list(dataset_dir.iterdir())
            print(f"      Found {len(items)} item(s) in dataset_dir:")
            for item in items[:10]:
                item_type = "DIR" if item.is_dir() else "FILE"
                print(f"        [{item_type}] {item.name}")
    
    # Search for result directories (index*_idx_*)
    result_dirs = []
    if header_dir.exists():
        result_dirs = list(header_dir.glob("index*_idx_*"))
        if not result_dirs:
            result_dirs = list(header_dir.glob("index*"))
    
    print(f"\n      Result directories found: {len(result_dirs)}")
    if result_dirs:
        for result_dir in result_dirs[:5]:  # Show first 5
            print(f"        - {result_dir.name}")
            # Count files by type
            sdf_files = list(result_dir.glob("*.sdf"))
            csv_files = list(result_dir.glob("*.csv"))
            pdb_files = list(result_dir.glob("*.pdb"))
            other_files = [f for f in result_dir.iterdir() if f.is_file() and f.suffix not in ['.sdf', '.csv', '.pdb']]
            print(f"          SDF: {len(sdf_files)}, CSV: {len(csv_files)}, PDB: {len(pdb_files)}, Other: {len(other_files)}")
            
            # Check for affinity file
            affinity_file = result_dir / "affinity_prediction.csv"
            if affinity_file.exists():
                print(f"          ✓ Found affinity_prediction.csv")
                try:
                    import pandas as pd
                    df = pd.read_csv(affinity_file)
                    print(f"          Contains {len(df)} row(s)")
                    if len(df) > 0:
                        print(f"          Columns: {list(df.columns)}")
                except Exception as e:
                    print(f"          Warning: Could not read CSV: {e}")
            
            # List SDF files
            if sdf_files:
                print(f"          SDF files found:")
                for sdf_file in sdf_files[:5]:  # Show first 5
                    size = sdf_file.stat().st_size
                    print(f"            - {sdf_file.name} ({size} bytes)")
                if len(sdf_files) > 5:
                    print(f"            ... and {len(sdf_files) - 5} more")
    else:
        print(f"      ✗ No result directories found (expected: index*_idx_*)")
    
    # Check stdout/stderr for critical errors
    critical_errors = []
    if stdout:
        # Look for error messages or paths in stdout
        if "error" in stdout.lower() or "exception" in stdout.lower():
            print(f"\n      ⚠ Warnings/Errors found in stdout:")
            lines = stdout.split('\n')
            error_lines = [line for line in lines if any(word in line.lower() for word in ['error', 'exception', 'warning', 'failed', 'fail'])]
            for line in error_lines[:5]:
                print(f"        {line[:100]}")  # First 100 chars
                if 'FileNotFoundError' in line or 'model_parameters' in line:
                    critical_errors.append("Model file not found")
        
        # Look for output directory mentions
        if "out_dir" in stdout or "results" in stdout.lower():
            print(f"\n      Output directory mentioned in stdout:")
            lines = stdout.split('\n')
            relevant_lines = [line for line in lines if 'out_dir' in line or 'results' in line.lower()]
            for line in relevant_lines[:3]:
                print(f"        {line[:100]}")
    
    if stderr:
        print(f"\n      ⚠ STDERR output (showing critical errors):")
        error_lines = stderr.split('\n')
        shown_lines = 0
        for line in error_lines:
            if line.strip():
                # Show all tracebacks and critical errors
                if any(keyword in line for keyword in ['Traceback', 'FileNotFoundError', 'NameError', 'Error:', 'Exception']):
                    print(f"        {line[:150]}")
                    shown_lines += 1
                    if 'FileNotFoundError' in line and 'model_parameters' in line:
                        critical_errors.append("DynamicBind model not found - check model_parameters.yml")
                    if 'NameError' in line and 'PDBFixer' in line:
                        critical_errors.append("PDBFixer not available in DynamicBind environment")
                elif shown_lines < 10:  # Show first 10 non-traceback lines
                    print(f"        {line[:100]}")
                    shown_lines += 1
        
        if critical_errors:
            print(f"\n      ❌ CRITICAL ERRORS DETECTED:")
            for error in set(critical_errors):
                print(f"        - {error}")
            print(f"      These errors prevent DynamicBind from creating results.")
    
    # Search for SDF files in entire dataset_dir recursively
    all_sdf_files = list(dataset_dir.glob("**/*.sdf"))
    print(f"\n      Total SDF files found in dataset_dir (recursive): {len(all_sdf_files)}")
    if all_sdf_files:
        print(f"      Sample SDF files:")
        for sdf_file in all_sdf_files[:10]:
            rel_path = sdf_file.relative_to(dataset_dir)
            size = sdf_file.stat().st_size
            print(f"        - {rel_path} ({size} bytes)")
        if len(all_sdf_files) > 10:
            print(f"        ... and {len(all_sdf_files) - 10} more")
    
    print(f"      === End Diagnostic ===\n")


def _save_intermediate_poses(dataset_dir: Path, header: str, output_dataset_dir: Path, dynamicbind_path: Optional[Path] = None) -> None:
    """
    Save intermediate poses (all SDF files) to a structured directory.
    Also attempts to convert SDF to PDBQT format if obabel is available.
    
    Args:
        dataset_dir: Directory where DynamicBind saved results (may contain header subdirectory)
        header: Header used for DynamicBind run
        output_dataset_dir: Target directory for intermediate poses
    """
    intermediate_dir = output_dataset_dir / "intermediate_poses"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    
    # Look for results in dataset_dir - DynamicBind saves in {results}/{header}/index*_idx_*/
    search_dirs = [dataset_dir]
    header_subdir = dataset_dir / header
    if header_subdir.exists():
        search_dirs.append(header_subdir)
    
    # Also check DynamicBind results directory if provided
    # DynamicBind might save results in its own directory structure
    # But don't check data directory - it contains old test files
    if dynamicbind_path:
        dynamicbind_results = dynamicbind_path / "results" / header
        if dynamicbind_results.exists():
            search_dirs.append(dynamicbind_results)
    
    # Also check parent directories recursively - DynamicBind might save elsewhere
    # Check up to 3 levels up
    parent = dataset_dir.parent
    for _ in range(3):
        if parent.exists() and parent.is_dir():
            search_dirs.append(parent)
            parent = parent.parent
        else:
            break
    
    copied_count = 0
    converted_count = 0
    import shutil
    
    # First, try to find SDF files directly in search directories recursively
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue
        
        # Search recursively for SDF files
        sdf_files_found = list(search_dir.glob("**/*.sdf"))
        
        for sdf_file in sdf_files_found:
            # Skip files in intermediate_poses directory to avoid duplicates
            if 'intermediate_poses' in str(sdf_file):
                continue
            
            # Skip old test files from DynamicBind/data (they have specific names and old timestamps)
            filename = sdf_file.name
            if filename in ['1opj_STI_A.sdf', '1opl_P16_A.sdf']:
                continue
            
            # Only process files that look like DynamicBind results
            # Results are typically in index*_idx_* directories or have rank* prefix
            # Also check if file path contains header to ensure it's from current run
            path_str = str(sdf_file)
            
            # Skip files that are clearly old test files (check both filename and path)
            if any(old_file in path_str for old_file in ['1opj', '1opl', '1qg8', 'cmet', 'gpcr', 'SETD2', 'test_imatinib']):
                continue
            
            # Check file modification time - skip files older than 1 hour (likely old test files)
            try:
                file_age = time.time() - sdf_file.stat().st_mtime
                if file_age > 3600:  # Older than 1 hour
                    # Only skip if it's not in a result directory structure
                    if 'index' not in path_str and 'idx' not in path_str:
                        continue
            except Exception:
                pass  # If we can't check time, continue anyway
            
            # Priority: look for files in index*_idx_* directories (actual DynamicBind results)
            # or files that contain the header in their path
            if 'index' in path_str or 'idx' in path_str:
                # This is definitely a DynamicBind result - process it
                pass
            elif header in path_str:
                # File path contains header - likely from current run
                pass
            elif 'rank' in path_str.lower():
                # Rank files are also results
                pass
            else:
                # Skip files that don't look like results
                continue
            
            try:
                # Create unique filename to avoid collisions
                rel_path = sdf_file.relative_to(search_dir)
                safe_name = str(rel_path).replace('/', '_').replace('\\', '_')
                sdf_dest = intermediate_dir / safe_name
                
                # Skip if already copied
                if sdf_dest.exists():
                    continue
                
                # Copy SDF file to intermediate directory
                shutil.copy2(sdf_file, sdf_dest)
                copied_count += 1
                
                # Try to convert to PDBQT if obabel is available
                pdbqt_dest = intermediate_dir / f"{sdf_dest.stem}.pdbqt"
                if not pdbqt_dest.exists():  # Don't reconvert if already exists
                    try:
                        result = subprocess.run(
                            ['obabel', '-isdf', str(sdf_file), '-opdbqt', '-O', str(pdbqt_dest)],
                            capture_output=True,
                            check=True,
                            timeout=30
                        )
                        if pdbqt_dest.exists():
                            converted_count += 1
                    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
                        pass  # obabel not available or conversion failed - that's okay
                        
            except Exception as e:
                # Don't print warning for every failed file - too verbose
                continue
    
    if copied_count > 0:
        print(f"      Saved {copied_count} intermediate pose(s) (SDF format) to {intermediate_dir}")
        if converted_count > 0:
            print(f"      Converted {converted_count} pose(s) to PDBQT format")
    else:
        # Debug: print what directories were searched (only first few to avoid spam)
        print(f"      Debug: No SDF files found. Searched {len(search_dirs)} directory(ies)")
        if search_dirs:
            print(f"      Debug: First search dir: {search_dirs[0]} (exists: {search_dirs[0].exists()})")


def preprocess_dynamicbind(
    protein_dir: Path,
    processed_dir: Path,
    interaction_config: Optional[Dict] = None,
) -> None:
    """
    Preprocess proteins for DynamicBind docking.
    
    DynamicBind uses raw PDB files, so minimal preprocessing is needed.
    This function mainly validates that required files exist.
    
    Args:
        protein_dir: Directory with input protein PDB files.
        processed_dir: Directory for processed files (not heavily used for DynamicBind).
        interaction_config: Optional interaction config for protein-ligand pairs.
    """
    if interaction_config is None:
        interaction_config = load_interaction_config()
    
    pairs = get_protein_ligand_pairs(interaction_config)
    if not pairs:
        print("  Warning: No protein-ligand pairs found. Skipping preprocessing.")
        return
    
    # Validate that protein PDB files exist
    protein_names = set(pair[0] for pair in pairs)
    missing_proteins = []
    
    for protein_name in protein_names:
        normalized_name = protein_name.lower().replace('_', '')
        protein_pdb = None
        
        for pdb_file in protein_dir.glob("*.pdb"):
            if pdb_file.stem.lower().replace('_', '') == normalized_name:
                protein_pdb = pdb_file
                break
        
        if protein_pdb is None:
            missing_proteins.append(protein_name)
    
    if missing_proteins:
        print(f"  Warning: Could not find PDB files for proteins: {', '.join(missing_proteins)}")
    else:
        print(f"  Validated {len(protein_names)} protein PDB files for DynamicBind")


def dock_dynamicbind(
    processed_dir: Path,
    output_dir: Path,
    config: Dict,
    interaction_config: Optional[Dict] = None,
    ligand_dir: Optional[Path] = None,
    protein_dir: Optional[Path] = None,
) -> None:
    """
    Run DynamicBind docking for all protein-ligand pairs.
    
    Results are saved in output_dir/{protein_name}/docking/dynamicbind/
    
    Args:
        processed_dir: Directory with prepared proteins and ligands.
        output_dir: Directory to save docking results (organized by protein).
        config: DynamicBind configuration (dynamicbind_path, device, etc.).
        interaction_config: Optional interaction config for protein-ligand pairs.
        ligand_dir: Directory with ligand CSV files (required for DynamicBind).
    """
    dynamicbind_raw = os.environ.get('DYNAMICBIND_PATH') or config.get('dynamicbind_path') or ''
    if not str(dynamicbind_raw).strip():
        raise RuntimeError(
            "DynamicBind path is not set. Export DYNAMICBIND_PATH or set "
            "dynamicbind.dynamicbind_path in config/methods_config.yaml "
            "(see config/methods_config.hpc.example.yaml)."
        )
    dynamicbind_path = Path(dynamicbind_raw)
    device = config.get('device', 0)
    samples_per_complex = config.get('samples_per_complex', 10)
    savings_per_complex = config.get('savings_per_complex', 1)
    inference_steps = config.get('inference_steps', 20)
    num_workers = config.get('num_workers', 20)
    use_relax = config.get('use_relax', True)
    protein_dynamic = config.get('protein_dynamic', True)
    random_seed = config.get('random_seed', 42)
    python_env = os.environ.get('DYNAMICBIND_PYTHON') or config.get('python_env') or None
    relax_python_env = os.environ.get('RELAX_PYTHON') or config.get('relax_python_env') or None
    if python_env is not None and not str(python_env).strip():
        python_env = None
    if relax_python_env is not None and not str(relax_python_env).strip():
        relax_python_env = None
    
    if not dynamicbind_path.exists():
        raise RuntimeError(f"DynamicBind path does not exist: {dynamicbind_path}")
    
    # Load interaction config
    if interaction_config is None:
        interaction_config = load_interaction_config()
    
    pairs = get_protein_ligand_pairs(interaction_config)
    if not pairs:
        print("  Warning: No protein-ligand pairs found. Skipping docking.")
        return
    
    if ligand_dir is None:
        raise ValueError("ligand_dir is required for DynamicBind docking")
    
    ligand_dir = Path(ligand_dir)
    
    print(f"  Running DynamicBind for {len(set(pair[0] for pair in pairs))} proteins...")
    
    # Check if run script exists
    run_script = dynamicbind_path / "run_single_protein_inference.py"
    if not run_script.exists():
        raise RuntimeError(f"DynamicBind run script not found: {run_script}")
    
    # Check if model directory exists
    model_workdir = dynamicbind_path / "workdir" / "big_score_model_sanyueqi_with_time"
    model_params_file = model_workdir / "model_parameters.yml"
    if not model_params_file.exists():
        print(f"  WARNING: DynamicBind model not found: {model_params_file}")
        print(f"  This will cause DynamicBind to fail. Please ensure the model is properly set up.")
        print(f"  Expected model directory: {model_workdir}")
        # Don't raise error - let DynamicBind fail and we'll catch it
    
    # Determine Python executables
    # conda_base: root of miniconda3 (when CONDA_PREFIX is e.g. .../envs/dynamicbind, go up two levels)
    _conda_prefix = Path(os.environ.get("CONDA_PREFIX", "") or "")
    _conda_base_env = Path(os.environ.get("CONDA_BASE", "") or "")
    if _conda_prefix.exists() and "envs" in _conda_prefix.parts:
        conda_base = _conda_prefix.parent.parent  # .../miniconda3/envs/xyz -> miniconda3
    elif _conda_base_env.exists():
        conda_base = _conda_base_env
    elif _conda_prefix.exists():
        conda_base = _conda_prefix
    else:
        conda_base = Path()
    if python_env is None:
        dynamicbind_python = conda_base / "envs" / "dynamicbind" / "bin" / "python"
        if dynamicbind_python.exists():
            python_env = str(dynamicbind_python)
        else:
            python_env = "python"  # Will use current environment
    if relax_python_env is None:
        relax_python = conda_base / "envs" / "relax" / "bin" / "python"
        if relax_python.exists():
            relax_python_env = str(relax_python)
        else:
            relax_python_env = "python"  # Will use current environment
    
    # Group pairs by protein
    protein_ligand_map = {}
    for protein_name, ligand_dataset, ref_ligand, safe_chain in pairs:
        normalized_protein_name = protein_name.lower()
        if normalized_protein_name not in protein_ligand_map:
            protein_ligand_map[normalized_protein_name] = {
                'original_name': protein_name,
                'ligand_datasets': []
            }
        if ligand_dataset not in protein_ligand_map[normalized_protein_name]['ligand_datasets']:
            protein_ligand_map[normalized_protein_name]['ligand_datasets'].append(ligand_dataset)
    
    docked_count = 0
    
    # Process each protein
    for normalized_protein_name, protein_info in protein_ligand_map.items():
        protein_name = protein_info['original_name']
        print(f"  Processing protein: {protein_name}")
        
        # Find protein PDB file
        protein_pdb = None
        if protein_dir is None:
            # Try to infer from processed_dir structure
            protein_dir = processed_dir.parent / "input" / "proteins"
        else:
            protein_dir = Path(protein_dir)
        
        for pdb_file in protein_dir.glob("*.pdb"):
            if pdb_file.stem.lower() == normalized_protein_name:
                protein_pdb = pdb_file
                break
        
        if not protein_pdb or not protein_pdb.exists():
            print(f"    Warning: Protein PDB not found for {protein_name}, skipping")
            continue
        
        # Create results directory
        protein_results_dir = output_dir / protein_name / "docking" / "dynamicbind"
        protein_results_dir.mkdir(parents=True, exist_ok=True)
        
        # Process each ligand dataset
        for ligand_dataset in protein_info['ligand_datasets']:
            print(f"    Processing ligand dataset: {ligand_dataset}")
            
            # Load ligands from CSV
            # ligand_dataset may already include '.csv'
            if str(ligand_dataset).endswith(".csv"):
                ligand_csv = ligand_dir / str(ligand_dataset)
            else:
                ligand_csv = ligand_dir / f"{ligand_dataset}.csv"
            if not ligand_csv.exists():
                print(f"      Warning: Ligand CSV not found: {ligand_csv}")
                continue
            
            try:
                ligands = load_ligands_from_csv(ligand_csv)
            except Exception as e:
                print(f"      Error loading ligands: {e}")
                continue
            
            if not ligands:
                print(f"      Warning: No ligands found in {ligand_csv}")
                continue
            
            # Prepare CSV file for DynamicBind (needs 'ligand' column with SMILES)
            temp_csv = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
            temp_csv_path = Path(temp_csv.name)
            
            try:
                # Create DataFrame with 'ligand' column
                ligand_data = []
                for ligand_info in ligands:
                    ligand_data.append({
                        'ligand': ligand_info['smiles'],
                        'name': ligand_info.get('ligand_id', f"ligand_{len(ligand_data)+1}")
                    })
                
                df = pd.DataFrame(ligand_data)
                df.to_csv(temp_csv_path, index=False)
                temp_csv.close()
                
                # Create dataset subdirectory
                dataset_dir = protein_results_dir / ligand_dataset
                dataset_dir.mkdir(parents=True, exist_ok=True)
                
                # Build command
                header = f"{protein_name}_{ligand_dataset}"
                cmd = [
                    python_env,
                    str(run_script),
                    str(protein_pdb),
                    str(temp_csv_path),
                    '--header', header,
                    '--results', str(dataset_dir),
                    '--samples_per_complex', str(samples_per_complex),
                    '--savings_per_complex', str(savings_per_complex),
                    '--inference_steps', str(inference_steps),
                    '--seed', str(random_seed),
                    '--device', str(device),
                    '--python', python_env,
                    '--relax_python', relax_python_env,
                    '--num_workers', str(num_workers),
                ]
                
                if not use_relax:
                    cmd.append('--no_relax')
                
                if not protein_dynamic:
                    cmd.append('--rigid_protein')
                
                # Run DynamicBind
                try:
                    start_time = time.time()
                    
                    # Change to DynamicBind directory for execution
                    original_cwd = os.getcwd()
                    os.chdir(dynamicbind_path)
                    
                    # Save command to log
                    log_file = dataset_dir / f"{header}_dynamicbind.log"
                    with open(log_file, 'w') as log_f:
                        log_f.write(f"DynamicBind execution log for {header}\n")
                        log_f.write(f"Command: {' '.join(cmd)}\n")
                        log_f.write(f"Working directory: {os.getcwd()}\n")
                        log_f.write(f"Results directory: {dataset_dir}\n")
                        log_f.write("=" * 80 + "\n\n")
                    
                    result = subprocess.run(
                        cmd,
                        check=True,
                        capture_output=True,
                        text=True,
                        timeout=config.get('docking_timeout', None)
                    )
                    
                    os.chdir(original_cwd)
                    
                    # Save stdout and stderr to log file
                    with open(log_file, 'a') as log_f:
                        log_f.write("STDOUT:\n")
                        log_f.write("=" * 80 + "\n")
                        log_f.write(result.stdout)
                        log_f.write("\n\n")
                        log_f.write("STDERR:\n")
                        log_f.write("=" * 80 + "\n")
                        log_f.write(result.stderr)
                        log_f.write("\n")
                    
                    elapsed_time = time.time() - start_time
                    docked_count += len(ligands)
                    print(f"      Docked {len(ligands)} ligands (took {elapsed_time:.1f}s)")
                    print(f"      Full log saved to: {log_file}")
                    
                    # Check what files/directories were created
                    _check_dynamicbind_results(dataset_dir, header, result.stdout, result.stderr)
                    
                except subprocess.CalledProcessError as e:
                    os.chdir(original_cwd)
                    print(f"      Error running DynamicBind: {e}")
                    
                    # Save error to log file
                    log_file = dataset_dir / f"{header}_dynamicbind_error.log"
                    with open(log_file, 'w') as log_f:
                        log_f.write(f"DynamicBind Error Log for {header}\n")
                        log_f.write("=" * 80 + "\n")
                        log_f.write(f"Error: {e}\n")
                        log_f.write(f"Return code: {e.returncode}\n")
                        log_f.write("\nSTDOUT:\n")
                        log_f.write("=" * 80 + "\n")
                        log_f.write(e.stdout if e.stdout else "(empty)\n")
                        log_f.write("\nSTDERR:\n")
                        log_f.write("=" * 80 + "\n")
                        log_f.write(e.stderr if e.stderr else "(empty)\n")
                    
                    print(f"      Error log saved to: {log_file}")
                    if e.stdout:
                        print(f"      stdout: {e.stdout[:500]}")
                    if e.stderr:
                        print(f"      stderr: {e.stderr[:500]}")
                    
                    # Still check for any partial results
                    print(f"      Checking for partial results...")
                    _check_dynamicbind_results(dataset_dir, header, e.stdout or "", e.stderr or "")
                    continue
                except Exception as e:
                    os.chdir(original_cwd)
                    print(f"      Error: {e}")
                    continue
                finally:
                    # Save intermediate poses (all SDF files, not just rank1)
                    # This runs even if there was an error, to save any partial results
                    try:
                        # Wait a bit for files to be written
                        time.sleep(1)
                        _save_intermediate_poses(dataset_dir, header, dataset_dir, dynamicbind_path)
                    except Exception as e:
                        print(f"      Warning: Could not save intermediate poses: {e}")
                        import traceback
                        traceback.print_exc()
                    
                    # Clean up temp file
                    if temp_csv_path.exists():
                        temp_csv_path.unlink()
            
            except Exception as e:
                print(f"      Error preparing ligands: {e}")
                if temp_csv_path.exists():
                    temp_csv_path.unlink()
                continue
    
    print(f"  DynamicBind docking complete: {docked_count} dockings performed")


def extract_metrics_dynamicbind(output_dir: Path) -> List[Dict]:
    """Extract metrics from DynamicBind docking results."""
    metrics = []
    
    for protein_dir in output_dir.iterdir():
        if not protein_dir.is_dir() or protein_dir.name == 'global':
            continue
        
        protein_name = protein_dir.name
        dynamicbind_dir = protein_dir / "docking" / "dynamicbind"
        
        if not dynamicbind_dir.exists():
            continue
        
        # Look for ligand dataset subdirectories
        for dataset_dir in dynamicbind_dir.iterdir():
            if not dataset_dir.is_dir():
                continue
            
            # DynamicBind may create header subdirectories or save directly in dataset_dir
            search_dirs = [dataset_dir]
            
            # Check for header subdirectories (format: {protein}_{dataset})
            for header_subdir in dataset_dir.iterdir():
                if header_subdir.is_dir() and '_' in header_subdir.name:
                    search_dirs.append(header_subdir)
            
            # Search in all possible locations
            for search_dir in search_dirs:
                # Look for result subdirectories (index*_idx_* or similar patterns)
                for result_dir in search_dir.iterdir():
                    if not result_dir.is_dir():
                        continue
                    
                    # Check if this looks like a result directory
                    is_result_dir = (
                        result_dir.name.startswith('index') or
                        'idx' in result_dir.name or
                        result_dir.name.startswith('rank') or
                        any(result_dir.glob("*.sdf")) or
                        (result_dir / "affinity_prediction.csv").exists()
                    )
                    
                    if not is_result_dir:
                        continue
                    
                    # Check for affinity_prediction.csv
                    affinity_file = result_dir / "affinity_prediction.csv"
                    if affinity_file.exists():
                        try:
                            df = pd.read_csv(affinity_file)
                            
                            # Extract metrics for each ligand
                            for _, row in df.iterrows():
                                ligand_name = row.get('name', f"ligand_{len(metrics)+1}")
                                
                                metric = {
                                    'method': 'dynamicbind',
                                    'protein': protein_name,
                                    'ligand': ligand_name,
                                    'output_file': str(affinity_file),
                                }
                                
                                # Extract affinity if available
                                if 'affinity' in row:
                                    metric['affinity'] = float(row['affinity'])
                                elif 'predicted_affinity' in row:
                                    metric['affinity'] = float(row['predicted_affinity'])
                                
                                # Look for best pose SDF file
                                for sdf_file in result_dir.glob("rank1_*.sdf"):
                                    metric['output_file'] = str(sdf_file)
                                    break
                                
                                # If no rank1, look for any SDF file
                                if 'output_file' in metric and metric['output_file'] == str(affinity_file):
                                    for sdf_file in result_dir.glob("*.sdf"):
                                        metric['output_file'] = str(sdf_file)
                                        break
                                
                                metrics.append(metric)
                        except Exception as e:
                            print(f"    Warning: Could not read {affinity_file}: {e}")
                            continue
                    else:
                        # If no affinity file, look for SDF files directly
                        for sdf_file in result_dir.glob("*.sdf"):
                            # Try to extract ligand name from filename
                            ligand_name = sdf_file.stem
                            # Remove rank prefix if present
                            if ligand_name.startswith('rank'):
                                ligand_name = '_'.join(ligand_name.split('_')[1:])
                            # Extract first part as ligand name
                            ligand_name = ligand_name.split('_')[0] if '_' in ligand_name else ligand_name
                            
                            metric = {
                                'method': 'dynamicbind',
                                'protein': protein_name,
                                'ligand': ligand_name,
                                'output_file': str(sdf_file),
                            }
                            metrics.append(metric)
                
                # Also check intermediate_poses directory for saved poses
                intermediate_dir = search_dir / "intermediate_poses"
                if intermediate_dir.exists():
                    for sdf_file in intermediate_dir.glob("*.sdf"):
                        ligand_name = sdf_file.stem.split('_')[0] if '_' in sdf_file.stem else sdf_file.stem
                        metric = {
                            'method': 'dynamicbind',
                            'protein': protein_name,
                            'ligand': ligand_name,
                            'output_file': str(sdf_file),
                        }
                        metrics.append(metric)
    
    return metrics

