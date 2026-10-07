"""Ligand preparation using Meeko Python API and RDKit - simplified function-based approach."""

import csv
import shutil
import tempfile
from pathlib import Path
from typing import Optional, Dict, List

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from meeko import pdbutils
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError as e:
    print(f"[ligand_prep] ERROR: Required packages not installed: {e}")
    print("[ligand_prep] Install with: pip install meeko rdkit scipy 'numpy<2.0'")
    raise

from ..utils.csv_utils import load_ligands_from_csv as load_ligands_from_csv_utils


def load_ligands_from_csv(csv_path: Path) -> List[Dict[str, str]]:
    """
    Load ligands from CSV file.
    
    Expected CSV format:
    - Must have a column with SMILES (can be named 'smiles', 'SMILES', 'ligand', etc.)
    - Optionally has 'ligand_id' or 'id' column for ligand identifiers
    
    Args:
        csv_path: Path to CSV file with SMILES.
    
    Returns:
        List of dictionaries with 'smiles' and 'ligand_id' keys.
    """
    # Use the improved version from csv_utils that handles None values and different delimiters
    return load_ligands_from_csv_utils(csv_path)


def _smiles_to_sdf(smiles: str, sdf_path: Path, mol_id: str = "mol") -> bool:
    """Convert SMILES to SDF file using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"    WARNING: failed to parse SMILES: {smiles[:50]}...")
            return False
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        writer = Chem.SDWriter(str(sdf_path))
        mol.SetProp("_Name", mol_id)
        writer.write(mol)
        writer.close()
        
        return True
    except Exception as e:
        print(f"    ERROR converting SMILES to SDF: {e}")
        return False


def _sdf_to_pdbqt(sdf_path: Path, pdbqt_path: Path) -> bool:
    """Convert SDF to PDBQT using Meeko Python API."""
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
        mol = supplier[0] if supplier else None
        
        if mol is None:
            print(f"    WARNING: failed to load molecule from {sdf_path}")
            return False
        
        # If molecule has multiple fragments, use only the largest
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) > 1:
            largest_frag = max(frags, key=lambda m: m.GetNumAtoms())
            mol = largest_frag
            print(f"    WARNING: molecule has {len(frags)} fragments, using largest fragment with {mol.GetNumAtoms()} atoms")
        
        # Prepare molecule through Meeko
        prep = MoleculePreparation.from_config({"charge_model": "gasteiger"})
        try:
            setups = prep.prepare(mol)
        except Exception as e:
            print(f"    WARNING: gasteiger charges failed for {sdf_path.name}, trying zero charges: {e}")
            prep = MoleculePreparation.from_config({"charge_model": "zero"})
            setups = prep.prepare(mol)
        
        # prep.prepare() can return a list or a single object
        if isinstance(setups, list):
            if len(setups) == 0:
                print(f"    WARNING: no setups generated for {sdf_path}")
                return False
            setup = setups[0]
        else:
            setup = setups
        
        # Get PDBQT string
        pdbqt_string, success, error = PDBQTWriterLegacy.write_string(setup, bad_charge_ok=True)
        
        if not success:
            print(f"    ERROR preparing ligand {sdf_path}: {error}")
            return False
        
        with open(pdbqt_path, 'w') as f:
            f.write(pdbqt_string)
        
        return True
        
    except Exception as e:
        print(f"    ERROR processing ligand {sdf_path}: {e}")
        import traceback
        traceback.print_exc()
        return False


def prepare_ligands(
    ligand_dir: Path,
    processed_dir: Path,
    interaction_config: Optional[Dict] = None,
    overwrite: bool = False,
) -> Dict[str, Dict[str, List[Path]]]:
    """
    Prepare ligands from CSV files, organized by protein.
    
    Args:
        ligand_dir: Directory with CSV files containing SMILES.
        processed_dir: Directory to save prepared ligands.
        interaction_config: Optional interaction config to determine protein-ligand pairs.
        overwrite: Whether to overwrite existing files.
    
    Returns:
        Dictionary: {protein_name: {ligand_dataset: [pdbqt_paths]}}
    """
    pdbqt_dir = processed_dir / "ligands"
    sdf_dir = processed_dir / "ligands_sdf"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    sdf_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # Get protein-ligand pairs from interaction config if provided
    if interaction_config:
        proteins = interaction_config.get("protein", [])
        ligands = interaction_config.get("ligand", [])
        pairs = list(zip(proteins, ligands))
    else:
        # If no config, process all CSV files for all proteins
        csv_files = list(ligand_dir.glob("*.csv"))
        pairs = [(None, csv_file.stem) for csv_file in csv_files]
    
    if not pairs:
        print(f"  Warning: No protein-ligand pairs found")
        return results
    
    print(f"  Preparing ligands for {len(pairs)} protein-ligand pairs...")
    
    temp_dir = Path(tempfile.gettempdir()) / "meeko_smiles_temp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    for protein_name, ligand_dataset in pairs:
        # ligand_dataset may already include '.csv' (e.g. FGFR1_Ki_WT_ChEMBL.csv)
        # Avoid constructing paths like '.../FGFR1_Ki_WT_ChEMBL.csv.csv'
        if str(ligand_dataset).endswith(".csv"):
            csv_path = ligand_dir / str(ligand_dataset)
        else:
            csv_path = ligand_dir / f"{ligand_dataset}.csv"
        if not csv_path.exists():
            print(f"    Warning: CSV file not found: {csv_path}")
            continue
        
        # Determine protein name (use from config or None)
        if protein_name:
            normalized_protein_name = protein_name.lower()
            ligand_sdf_dir = sdf_dir / normalized_protein_name / ligand_dataset
            ligand_pdbqt_dir = pdbqt_dir / normalized_protein_name / ligand_dataset
        else:
            # If no protein specified, use global directory
            ligand_sdf_dir = sdf_dir / ligand_dataset
            ligand_pdbqt_dir = pdbqt_dir / ligand_dataset
        
        ligand_sdf_dir.mkdir(parents=True, exist_ok=True)
        ligand_pdbqt_dir.mkdir(parents=True, exist_ok=True)
        
        # Load ligands from CSV
        try:
            ligands = load_ligands_from_csv(csv_path)
        except Exception as e:
            print(f"    Error loading ligands from {csv_path}: {e}")
            continue
        
        # Process each ligand
        pdbqt_paths = []
        for idx, ligand_data in enumerate(ligands):
            smiles = ligand_data['smiles']
            ligand_id = ligand_data.get("ligand_id") or f"ligand_{idx+1:04d}"
            
            sdf_path = ligand_sdf_dir / f"{ligand_id}.sdf"
            pdbqt_path = ligand_pdbqt_dir / f"{ligand_id}.pdbqt"
            
            # Skip if both files exist
            if sdf_path.exists() and pdbqt_path.exists() and not overwrite:
                pdbqt_paths.append(pdbqt_path)
                continue
            
            # Convert SMILES to SDF
            temp_sdf = temp_dir / f"{ligand_id}.sdf"
            if not _smiles_to_sdf(smiles, temp_sdf, ligand_id):
                print(f"    WARNING: failed to convert SMILES to 3D structure for {ligand_id}")
                continue
            
            # Copy SDF to output directory
            shutil.copy2(temp_sdf, sdf_path)
            
            # Convert SDF to PDBQT
            if _sdf_to_pdbqt(sdf_path, pdbqt_path):
                pdbqt_paths.append(pdbqt_path)
            
            # Clean up temp file
            if temp_sdf.exists():
                temp_sdf.unlink()
        
        # Store results
        if protein_name:
            if protein_name not in results:
                results[protein_name] = {}
            results[protein_name][ligand_dataset] = pdbqt_paths
        else:
            if 'global' not in results:
                results['global'] = {}
            results['global'][ligand_dataset] = pdbqt_paths
        
        print(f"    Prepared {len(pdbqt_paths)} ligands from {ligand_dataset}")
    
    print(f"  Successfully prepared ligands for {len(results)} proteins")
    return results

