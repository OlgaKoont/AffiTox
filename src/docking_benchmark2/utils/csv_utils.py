"""CSV utilities - no dependencies on meeko/gemmi."""

import csv
from pathlib import Path
from typing import Dict, List

from docking_benchmark2.ligand_ids import format_ligand_id, normalize_ligand_id


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
    ligands = []
    
    try:
        # Try to detect delimiter by reading first line
        with open(csv_path, 'r', encoding='utf-8') as f:
            first_line = f.readline()
            # Count delimiters
            comma_count = first_line.count(',')
            semicolon_count = first_line.count(';')
            tab_count = first_line.count('\t')
            
            # Choose delimiter with most occurrences
            if semicolon_count > comma_count and semicolon_count > tab_count:
                delimiter = ';'
            elif tab_count > comma_count:
                delimiter = '\t'
            else:
                delimiter = ','
        
        # Read with detected delimiter
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            rows = list(reader)
            
            # If rows are empty or first row has only one column, try alternative delimiters
            if not rows or (rows and len(rows[0]) == 1):
                for alt_delim in [';', '\t', ',']:
                    if alt_delim == delimiter:
                        continue
                    with open(csv_path, 'r', encoding='utf-8') as f2:
                        reader2 = csv.DictReader(f2, delimiter=alt_delim)
                        rows2 = list(reader2)
                        if rows2 and len(rows2[0]) > 1:
                            rows = rows2
                            break
    except Exception as e:
        raise ValueError(f"CSV file {csv_path} could not be parsed: {e}")
    
    if not rows:
        raise ValueError(f"CSV file {csv_path} is empty or has no data rows")
    
    # Find SMILES column
    smiles_col = None
    for col in rows[0].keys():
        if col is None:
            continue
        col_lower = col.lower() if col else ""
        if col_lower in ['smiles', 'smile', 'ligand', 'canonical_smiles', 
                       'canonical_smile', 'canonicalsmiles', 'canonicalsmile']:
            smiles_col = col
            break
    
    if smiles_col is None:
        raise ValueError(f"No SMILES column found in {csv_path}. Available columns: {list(rows[0].keys())}")
    
    # Find ID column
    id_col = None
    for col in rows[0].keys():
        if col is None:
            continue
        col_lower = col.lower() if col else ""
        if col_lower in ['ligand_id', 'id', 'name']:
            id_col = col
            break
    
    # Process rows
    for idx, row in enumerate(rows):
        smiles = row.get(smiles_col, "").strip() if smiles_col and smiles_col in row else ""
        if not smiles or (smiles and smiles.lower() == 'nan'):
            continue
        
        raw_id = row.get(id_col, "").strip() if id_col and id_col in row else ""
        ligand_id = normalize_ligand_id(raw_id) or format_ligand_id(idx + 1)
        ligands.append({
            'smiles': smiles,
            'ligand_id': ligand_id
        })
    
    print(f"    Loaded {len(ligands)} ligands from {csv_path.name}")
    return ligands



