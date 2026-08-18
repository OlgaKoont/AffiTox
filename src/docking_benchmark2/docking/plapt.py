"""PLAPT docking - simplified function-based approach.

PLAPT (Protein-Ligand Affinity Prediction Using Pretrained Transformers) is a
deep learning-based method for predicting protein-ligand binding affinity.
It uses protein sequences and SMILES strings as input.

Reference: https://github.com/trrt-good/WELP-PLAPT
"""

import json
import os
import random
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# Suppress onnxruntime warnings
os.environ.setdefault('ORT_LOGGING_LEVEL', '3')
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*found in sys.modules.*')

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.simplefilter('ignore', PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    PDBParser = None

from ..utils.settings import load_interaction_config, get_protein_ligand_pairs
from ..utils.csv_utils import load_ligands_from_csv


def extract_protein_sequence(pdb_path: Path) -> Optional[str]:
    """
    Extract amino acid sequence from PDB file.
    
    Args:
        pdb_path: Path to PDB file.
    
    Returns:
        Protein sequence as string, or None if extraction fails.
    """
    if not BIOPYTHON_AVAILABLE:
        print("    Warning: BioPython not available, cannot extract sequence from PDB")
        return None
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', str(pdb_path))
        
        sequence = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname()
                    aa_map = {
                        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
                    }
                    if resname in aa_map:
                        sequence.append(aa_map[resname])
        
        return ''.join(sequence) if sequence else None
    except Exception as e:
        print(f"    Warning: Error extracting sequence from {pdb_path}: {e}")
        return None


def preprocess_plapt(
    protein_dir: Path,
    processed_dir: Path,
    interaction_config: Optional[Dict] = None,
) -> Dict[str, str]:
    """
    Preprocess PLAPT: extract protein sequences and save them.
    
    Args:
        protein_dir: Directory with protein PDB files.
        processed_dir: Directory to save processed sequences.
        interaction_config: Optional interaction config for protein-ligand pairs.
    
    Returns:
        Dictionary mapping protein names to sequences.
    """
    sequences_dir = processed_dir / "plapt_sequences"
    sequences_dir.mkdir(parents=True, exist_ok=True)
    
    sequences = {}
    
    # Load interaction config to get protein list
    if interaction_config is None:
        interaction_config = load_interaction_config()
    
    pairs = get_protein_ligand_pairs(interaction_config)
    if not pairs:
        print("  Warning: No protein-ligand pairs found. Skipping preprocessing.")
        return sequences
    
    # Get unique proteins from config
    config_proteins = {pair[0] for pair in pairs}
    normalized_config_proteins = {p.lower() for p in config_proteins}
    
    print(f"  Extracting sequences for {len(config_proteins)} proteins...")
    sequences_extracted = 0
    
    for pdb_file in protein_dir.glob("*.pdb"):
        protein_name = pdb_file.stem
        normalized_name = protein_name.lower()
        
        if normalized_name not in normalized_config_proteins:
            continue
        
        sequence_file = sequences_dir / f"{protein_name}.txt"
        
        # Skip if already extracted
        if sequence_file.exists():
            with open(sequence_file, 'r') as f:
                sequences[protein_name] = f.read().strip()
            continue
        
        print(f"    Extracting sequence for {protein_name}...")
        sequence = extract_protein_sequence(pdb_file)
        
        if sequence:
            sequences[protein_name] = sequence
            with open(sequence_file, 'w') as f:
                f.write(sequence)
            sequences_extracted += 1
            print(f"      Extracted sequence ({len(sequence)} residues)")
        else:
            print(f"      Warning: Failed to extract sequence for {protein_name}")
    
    print(f"  Extracted {sequences_extracted} protein sequences")
    return sequences


def dock_plapt(
    processed_dir: Path,
    output_dir: Path,
    config: Dict,
    interaction_config: Optional[Dict] = None,
    ligand_dir: Optional[Path] = None,
) -> None:
    """
    Run PLAPT affinity prediction for all protein-ligand pairs.
    
    Results are saved in output_dir/{protein_name}/docking/plapt/
    
    Args:
        processed_dir: Directory with prepared proteins and ligands.
        output_dir: Directory to save docking results (organized by protein).
        config: PLAPT configuration (plapt_path, device, batch_size, etc.).
        interaction_config: Optional interaction config for protein-ligand pairs.
        ligand_dir: Directory with ligand CSV files (required for PLAPT).
    """
    plapt_raw = os.environ.get('PLAPT_PATH') or config.get('plapt_path') or ''
    if not str(plapt_raw).strip():
        raise RuntimeError(
            "PLAPT path is not set. Export PLAPT_PATH or set plapt.plapt_path in "
            "config/methods_config.yaml (see config/methods_config.hpc.example.yaml)."
        )
    plapt_path = Path(plapt_raw)
    device = config.get('device', 'cuda')
    batch_size = config.get('batch_size', 16)
    affinity_batch_size = config.get('affinity_batch_size', 128)
    docking_timeout = config.get('docking_timeout', 300)
    random_seed = config.get('random_seed', 42)  # Fixed seed for reproducibility
    
    if not plapt_path.exists():
        raise RuntimeError(f"PLAPT path does not exist: {plapt_path}")
    
    # Set random seed for reproducibility
    random.seed(random_seed)
    np.random.seed(random_seed)
    # Also set torch random seed if available
    try:
        import torch
        torch.manual_seed(random_seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(random_seed)
    except ImportError:
        pass  # PyTorch not available
    
    # Load interaction config
    if interaction_config is None:
        interaction_config = load_interaction_config()
    
    pairs = get_protein_ligand_pairs(interaction_config)
    if not pairs:
        print("  Warning: No protein-ligand pairs found. Skipping docking.")
        return
    
    if ligand_dir is None:
        raise ValueError("ligand_dir is required for PLAPT docking")
    
    ligand_dir = Path(ligand_dir)
    
    print(f"  Running PLAPT for {len(set(pair[0] for pair in pairs))} proteins...")
    
    # Check if PLAPT module exists
    plapt_module_path = plapt_path / "plapt.py"
    if not plapt_module_path.exists():
        print(f"  Error: PLAPT module not found at {plapt_module_path}")
        print("  Please ensure PLAPT is installed at the specified path")
        return
    
    # Add PLAPT path to Python path
    plapt_parent = str(plapt_path)
    if plapt_parent not in sys.path:
        sys.path.insert(0, plapt_parent)
    
    # Import PLAPT
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from plapt import Plapt
    except ImportError as e:
        print(f"  Error: Could not import PLAPT: {e}")
        print(f"  Make sure PLAPT is installed and accessible from {plapt_path}")
        return
    
    # Initialize PLAPT
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plapt = Plapt(
                prediction_module_path=str(plapt_path / "models" / "affinity_predictor.onnx"),
                device=device,
                cache_dir=str(plapt_path / "embedding_cache"),
                use_tqdm=False
            )
    except Exception as e:
        print(f"  Error: Could not initialize PLAPT: {e}")
        return
    
    # Load sequences
    sequences_dir = processed_dir / "plapt_sequences"
    sequences = {}
    for seq_file in sequences_dir.glob("*.txt"):
        protein_name = seq_file.stem
        with open(seq_file, 'r') as f:
            sequences[protein_name.lower()] = f.read().strip()
    
    # Group pairs by protein for batch processing
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
        
        # Load protein sequence
        protein_sequence = sequences.get(normalized_protein_name)
        if not protein_sequence:
            print(f"    Warning: Sequence not found for {protein_name}, skipping")
            continue
        
        # Create results directory: output_dir/{protein_name}/docking/plapt/
        protein_results_dir = output_dir / protein_name / "docking" / "plapt"
        protein_results_dir.mkdir(parents=True, exist_ok=True)
        
        # Process each ligand dataset for this protein
        for ligand_dataset in protein_info['ligand_datasets']:
            print(f"    Processing ligand dataset: {ligand_dataset}")
            
            # Load ligands from CSV
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
            
            # Prepare batch
            ligand_ids = []
            smiles_list = []
            
            for ligand_data in ligands:
                ligand_id = ligand_data.get('ligand_id', f"ligand_{len(ligand_ids)+1}")
                smiles = ligand_data['smiles']
                ligand_ids.append(ligand_id)
                smiles_list.append(smiles)
            
            # Create dataset subdirectory
            dataset_dir = protein_results_dir / ligand_dataset
            dataset_dir.mkdir(parents=True, exist_ok=True)
            
            # Run PLAPT prediction using score_candidates
            try:
                start_time = time.time()
                
                # Use score_candidates for single protein, multiple ligands
                results = plapt.score_candidates(
                    target_protein=protein_sequence,
                    mol_smiles=smiles_list,
                    mol_batch_size=batch_size,
                    affinity_batch_size=affinity_batch_size
                )
                
                elapsed_time = time.time() - start_time
                
                # Save results
                for i, (ligand_id, result) in enumerate(zip(ligand_ids, results)):
                    result_file = dataset_dir / f"{ligand_id}.json"
                    result_data = {
                        'protein': protein_name,
                        'ligand_id': ligand_id,
                        'smiles': smiles_list[i],
                        'affinity': result.get('neg_log10_affinity_M'),
                        'affinity_uM': result.get('affinity_uM'),
                        'prediction_time': elapsed_time / len(ligands) if ligands else 0
                    }
                    
                    with open(result_file, 'w') as f:
                        json.dump(result_data, f, indent=2)
                
                docked_count += len(ligands)
                print(f"      Predicted affinity for {len(ligands)} ligands (took {elapsed_time:.1f}s)")
                
            except Exception as e:
                print(f"      Error running PLAPT for {protein_name}/{ligand_dataset}: {e}")
                continue
    
    print(f"  PLAPT prediction complete: {docked_count} predictions performed")


def extract_metrics_plapt(output_dir: Path) -> List[Dict]:
    """Extract metrics from PLAPT docking results."""
    metrics = []
    
    for protein_dir in output_dir.iterdir():
        if not protein_dir.is_dir() or protein_dir.name == 'global':
            continue
        
        protein_name = protein_dir.name
        plapt_dir = protein_dir / "docking" / "plapt"
        
        if not plapt_dir.exists():
            continue
        
        # Look for ligand dataset subdirectories
        for dataset_dir in plapt_dir.iterdir():
            if not dataset_dir.is_dir():
                continue
            
            for result_file in dataset_dir.glob("*.json"):
                ligand_id = result_file.stem
                
                try:
                    with open(result_file, 'r') as f:
                        result_data = json.load(f)
                    
                    affinity = result_data.get('affinity')  # neg_log10_affinity_M
                    
                    metric = {
                        'method': 'plapt',
                        'protein': protein_name,
                        'ligand': ligand_id,
                        'output_file': str(result_file),
                    }
                    
                    if affinity is not None:
                        metric['affinity'] = affinity
                    if 'affinity_uM' in result_data:
                        metric['affinity_uM'] = result_data['affinity_uM']
                    
                    metrics.append(metric)
                except Exception as e:
                    print(f"    Warning: Could not extract metrics from {result_file}: {e}")
                    continue
    
    return metrics
