"""Protein preparation using Meeko Python API - simplified function-based approach."""

import sys
import subprocess
import shutil
import re
from pathlib import Path
from typing import Dict, Optional, Any, List, Tuple

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy, Polymer, ResidueChemTemplates
    from meeko import pdbutils
    from rdkit import Chem
    
    # Patch for RDKit compatibility: add HasQuery method if missing
    if not hasattr(Chem.Mol, 'HasQuery'):
        def _has_query(self):
            """Compatibility method for older Meeko versions."""
            return False
        Chem.Mol.HasQuery = _has_query
except ImportError as e:
    print(f"[protein_prep] ERROR: Required packages not installed: {e}", file=sys.stderr)
    print("[protein_prep] Install with: pip install meeko scipy 'numpy<2.0'", file=sys.stderr)
    raise


# Mapping from REF_LIGAND names to PDB residue codes (3-letter codes)
# This maps common ligand names to their PDB HETATM residue names
# Format: "REF_LIGAND": "PDB_CODE" or None if not found
REF_LIGAND_TO_PDB_CODE = {
    "C8": "0S9",      # 4f65
    "BAK": None,      # 1g5m - no HETATM found, may need manual check
    "ONX": "S5K",     # 7awe (chain L)
    "ASC+NIL": "NIL", # 5mo4 - use NIL part (ASC is AY7 but we use NIL)
    "TLZ": "2YQ",     # 7kk3 (chain C)
    "RAG": "FAD",     # 2z5x - using FAD as best match
    "STI": "B49",     # 6jok - default, but see protein-specific mapping below
    "LEN": "1PE",     # 5jkv - using 1PE as best match (or could be ASD/HEM)
    "LZ0": "AV9",     # 4ase
    "TIV": "AV9",     # 4ase - same as LZ0
    "SCF": "F82",     # 6gqj
    "OSM": "YY3",     # 4zau
}

# Protein-specific mapping: (protein_id, ref_ligand) -> PDB code
# Use this when the same ref_ligand maps to different PDB codes for different proteins
PROTEIN_REF_LIGAND_TO_PDB_CODE = {
    ("4tz4", "STI"): "LVY",  # 4tz4 uses LVY in chain C (will be found even if safe_chain is B)
    ("3mjg", "STI"): None,   # 3mjg - no STI found, has NAG/NDG in chains X/Y
    ("1g5m", "BAK"): None,   # 1g5m - no BAK found, no HETATM in chain A
    ("6jok", "STI"): "B49",  # 6jok uses B49
}


def _normalize_protein_id(protein_id: str) -> str:
    """Normalize protein ID to lowercase for case-insensitive comparison."""
    return protein_id.lower()


def _get_setting(protein_id: str, settings: Dict[str, Any], key: str, default: Any = None) -> Any:
    """Get setting for protein with case-insensitive lookup."""
    normalized_id = _normalize_protein_id(protein_id)
    # Try exact match first
    if protein_id in settings:
        if isinstance(settings[protein_id], dict) and key in settings[protein_id]:
            return settings[protein_id].get(key, default)
    # Try case-insensitive match
    for setting_key, setting_value in settings.items():
        if _normalize_protein_id(setting_key) == normalized_id:
            if isinstance(setting_value, dict) and key in setting_value:
                return setting_value.get(key, default)
    # Try a global default section (optional)
    if "_default" in settings and isinstance(settings.get("_default"), dict):
        return settings["_default"].get(key, default)
    return default


def _detect_pdb2pqr_invocation() -> Optional[List[str]]:
    """
    Return an invocation prefix for pdb2pqr, e.g. ["/usr/bin/pdb2pqr"] or [sys.executable, "-m", "pdb2pqr"].
    Returns None if not available.
    """
    for exe in ("pdb2pqr30", "pdb2pqr"):
        path = shutil.which(exe)
        if path:
            return [path]
    # Try python module as a fallback
    try:
        import pdb2pqr  # type: ignore
        return [sys.executable, "-m", "pdb2pqr"]
    except Exception:
        return None


def _pdb2pqr_help(inv: List[str]) -> str:
    try:
        res = subprocess.run(inv + ["--help"], check=False, capture_output=True, text=True)
        return (res.stdout or "") + "\n" + (res.stderr or "")
    except Exception:
        return ""


def _pick_flag_style(help_text: str, flag: str) -> str:
    """
    Decide whether a CLI flag is shown as '--flag=VAL' or '--flag VAL'.
    Returns 'equals' or 'space'. Defaults to 'space'.
    """
    if f"{flag}=" in help_text:
        return "equals"
    return "space"


def _build_pdb2pqr_cmd(
    inv: List[str],
    inp_pdb: Path,
    out_pqr: Path,
    ph: float,
    ff: str = "AMBER",
    out_pdb: Optional[Path] = None,
) -> List[str]:
    """
    Build a pdb2pqr command line in a version-tolerant way (based on --help output).
    """
    help_text = _pdb2pqr_help(inv)

    # Force field flag
    ff_flag = "--ff"
    ff_style = _pick_flag_style(help_text, ff_flag)
    if ff_style == "equals":
        ff_arg = f"{ff_flag}={ff}"
        ff_parts = [ff_arg]
    else:
        ff_parts = [ff_flag, ff]

    # pH flag (varies across versions)
    # Common: --with-ph or --with-ph=7.0
    ph_flag = "--with-ph"
    if "with-ph" not in help_text and "with_ph" in help_text:
        ph_flag = "--with_ph"
    ph_style = _pick_flag_style(help_text, ph_flag)
    if ph_style == "equals":
        ph_parts = [f"{ph_flag}={ph}"]
    else:
        ph_parts = [ph_flag, str(ph)]

    cmd = inv + ff_parts + ph_parts

    # Optional PDB output (not always supported)
    if out_pdb is not None and ("--pdb-output" in help_text or "pdb-output" in help_text):
        pdb_flag = "--pdb-output"
        pdb_style = _pick_flag_style(help_text, pdb_flag)
        if pdb_style == "equals":
            cmd += [f"{pdb_flag}={out_pdb}"]
        else:
            cmd += [pdb_flag, str(out_pdb)]

    cmd += [str(inp_pdb), str(out_pqr)]
    return cmd


def _write_atom_only_pdb(source_pdb: Path, atom_only_pdb: Path) -> None:
    """
    Write only ATOM records from source PDB into a new file.
    This is used to run protonation tools that may not preserve HETATM cofactors/metals.
    """
    atom_only_pdb.parent.mkdir(parents=True, exist_ok=True)
    with source_pdb.open("r", encoding="utf-8") as fin, atom_only_pdb.open("w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("ATOM"):
                fout.write(line)


def _extract_hetatm_lines(source_pdb: Path) -> List[str]:
    het = []
    with source_pdb.open("r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith("HETATM"):
                het.append(line)
    return het


def _run_protonation_pdb2pqr(
    input_pdb: Path,
    output_pdb: Path,
    ph: float,
    ff: str = "AMBER",
) -> Tuple[bool, str]:
    """
    Run pdb2pqr to assign protonation at a given pH.
    Returns (success, message). On success, output_pdb is written.
    """
    inv = _detect_pdb2pqr_invocation()
    if inv is None:
        return False, "pdb2pqr not found (install with: pip install pdb2pqr)"

    out_pqr = output_pdb.with_suffix(".pqr")
    out_pdb_opt = output_pdb  # try to request PDB output if supported

    cmd = _build_pdb2pqr_cmd(inv, input_pdb, out_pqr, ph=ph, ff=ff, out_pdb=out_pdb_opt)
    try:
        res = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if res.returncode != 0:
            msg = (res.stdout or "") + "\n" + (res.stderr or "")
            return False, f"pdb2pqr failed (code={res.returncode}). Output:\n{msg}"

        # If pdb2pqr produced the requested PDB output, we are done.
        if out_pdb_opt.exists() and out_pdb_opt.stat().st_size > 0:
            return True, "ok"

        # Fallback: if only PQR exists, use it as PDB-like file (coordinates are in the same columns).
        if out_pqr.exists() and out_pqr.stat().st_size > 0:
            output_pdb.write_text(out_pqr.read_text(), encoding="utf-8")
            return True, "ok (used PQR as PDB-like input)"

        return False, "pdb2pqr produced no output files"
    except Exception as e:
        return False, f"pdb2pqr invocation error: {e}"


def _maybe_protonate_before_meeko(
    protein_id: str,
    cleaned_pdb: Path,
    include_cofactors: bool,
    protonate_cfg: Optional[Dict[str, Any]],
    overwrite: bool,
) -> Path:
    """
    Optionally run pH-based protonation on ATOM records prior to Meeko receptor prep,
    then merge back any HETATM lines from the cleaned PDB (cofactors/metals already filtered upstream).
    """
    if not protonate_cfg or not isinstance(protonate_cfg, dict) or not protonate_cfg.get("enabled", False):
        return cleaned_pdb

    ph = float(protonate_cfg.get("ph", 7.0))
    ff = str(protonate_cfg.get("ff", "AMBER"))
    tool = str(protonate_cfg.get("tool", "pdb2pqr")).strip().lower()
    if tool not in {"pdb2pqr"}:
        print(f"    [PREPARATION] Protonation tool '{tool}' not supported; skipping protonation", file=sys.stderr)
        return cleaned_pdb

    # We protonate only ATOM records and then re-attach HETATM lines from cleaned_pdb
    ph_tag = re.sub(r"[^0-9A-Za-z]+", "p", f"{ph:.2f}".rstrip("0").rstrip("."))
    atom_only = cleaned_pdb.with_name(f"{cleaned_pdb.stem}_atom_only.pdb")
    protonated_atom = cleaned_pdb.with_name(f"{cleaned_pdb.stem}_atom_only_ph{ph_tag}.pdb")
    protonated_merged = cleaned_pdb.with_name(f"{cleaned_pdb.stem}_ph{ph_tag}.pdb")

    if protonated_merged.exists() and not overwrite:
        return protonated_merged

    _write_atom_only_pdb(cleaned_pdb, atom_only)
    ok, msg = _run_protonation_pdb2pqr(atom_only, protonated_atom, ph=ph, ff=ff)
    if not ok:
        print(f"    [PREPARATION] Protonation skipped for {protein_id}: {msg}", file=sys.stderr)
        return cleaned_pdb

    hetatm_lines = _extract_hetatm_lines(cleaned_pdb) if include_cofactors else []
    # Merge: protonated ATOM + original HETATM from cleaned_pdb
    with protonated_atom.open("r", encoding="utf-8") as fin, protonated_merged.open("w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("ATOM"):
                fout.write(line)
        for line in hetatm_lines:
            fout.write(line)

    print(f"    [PREPARATION] Protonation (pH={ph}) applied via {tool} -> {protonated_merged.name} ({msg})")
    return protonated_merged


def _generate_clean_pdb(
    source: Path,
    destination: Path,
    chain_id: str,
    include_het: bool,
    include_cofactors: bool,
    include_waters: bool,
    keep_het_resnames: Optional[List[str]] = None,
    drop_het_resnames: Optional[List[str]] = None,
    ligand_destination: Optional[Path] = None,
    ref_ligand: Optional[str] = None,
    ligand_chain: Optional[str] = None,
    protein_id: Optional[str] = None,
) -> None:
    """
    Filter the original PDB to the requested chain and optional records.
    Also extract specific ligand (HETATM records) from specified chain to separate file.
    
    Args:
        source: Source PDB file
        destination: Output cleaned PDB file
        chain_id: Chain ID for protein
        include_het: Include HETATM in protein
        include_cofactors: Include cofactors in protein
        include_waters: Include waters in protein
        ligand_destination: Output file for extracted ligand
        ref_ligand: Residue name of the ligand to extract (e.g., "C8", "BAK")
        ligand_chain: Chain ID where the ligand is located
    """
    chain_id = chain_id.strip() or "A"
    ligand_chain = ligand_chain.strip() if ligand_chain else chain_id
    keep_het_resnames_set = set(r.strip().upper() for r in (keep_het_resnames or []) if str(r).strip())
    drop_het_resnames_set = set(r.strip().upper() for r in (drop_het_resnames or []) if str(r).strip())
    selected_lines = []
    ligand_lines_preferred = []  # Ligands from preferred chain
    ligand_lines_other = []       # Ligands from other chains
    chain_col = 21
    record_whitelist = {"ATOM"}

    if include_het or include_cofactors:
        record_whitelist.add("HETATM")

    with source.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            # Extract specific ligand from specified chain
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                resname_u = resname.upper()
                line_chain = line[chain_col].strip()
                
                # Extract ligand only if:
                # 1. Not water
                # 2. Matches ref_ligand name (if specified) - uses mapping for PDB codes
                # 3. Matches ligand_chain (if specified)
                if resname != "HOH":
                    match_ligand = False
                    if ref_ligand is None:
                        match_ligand = True
                    else:
                        # First, try to get PDB code from protein-specific mapping
                        pdb_codes = []
                        if protein_id and ref_ligand:
                            protein_key = (protein_id.lower(), ref_ligand)
                            if protein_key in PROTEIN_REF_LIGAND_TO_PDB_CODE:
                                mapped_code = PROTEIN_REF_LIGAND_TO_PDB_CODE[protein_key]
                                if mapped_code is not None:
                                    pdb_codes.append(mapped_code)
                        
                        # Then try general mapping
                        if ref_ligand in REF_LIGAND_TO_PDB_CODE:
                            mapped_code = REF_LIGAND_TO_PDB_CODE[ref_ligand]
                            if mapped_code is not None:
                                pdb_codes.append(mapped_code)
                        # Also check original ref_ligand (in case it's already a PDB code)
                        pdb_codes.append(ref_ligand)
                        
                        # Check if ref_ligand contains "+" (compound name like "ASC+NIL")
                        if "+" in ref_ligand:
                            # Split compound name and check each part
                            parts = [p.strip() for p in ref_ligand.split("+")]
                            for part in parts:
                                if part in REF_LIGAND_TO_PDB_CODE:
                                    mapped = REF_LIGAND_TO_PDB_CODE[part]
                                    if mapped is not None:
                                        pdb_codes.append(mapped)
                                pdb_codes.append(part)
                        
                        # Check if resname matches any of the possible PDB codes
                        if resname in pdb_codes:
                            match_ligand = True
                    
                    if match_ligand:
                        # Prefer ligands from specified chain, but also collect from other chains
                        if ligand_chain is None or line_chain == ligand_chain:
                            ligand_lines_preferred.append(line)
                        else:
                            # Store ligands from other chains as fallback
                            ligand_lines_other.append(line)
                
                # Continue processing for protein cleaning
                if "HETATM" not in record_whitelist:
                    continue
                
                if not include_waters and resname == "HOH":
                    continue

                # Optional: keep/drop specific HETATM residue names.
                #
                # This is useful to keep essential cofactors/metals (e.g., HEM/FAD/ZN/PTR)
                # while still excluding co-crystal ligands/solvents by default.
                if resname != "HOH":
                    if keep_het_resnames_set:
                        if resname_u not in keep_het_resnames_set:
                            continue
                    elif drop_het_resnames_set and resname_u in drop_het_resnames_set:
                        continue
                
                if line_chain != chain_id:
                    continue
                
                selected_lines.append(line)
            else:
                # ATOM records - only for the specified chain
                if line[chain_col].strip() != chain_id:
                    continue
                selected_lines.append(line)

    if not selected_lines:
        raise RuntimeError(
            f"No atoms selected for {source.name} with chain '{chain_id}'."
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        handle.writelines(selected_lines)
    
    # Save ligand if destination provided
    if ligand_destination is not None:
        # Use preferred chain ligands if available, otherwise use other chains
        ligand_lines = ligand_lines_preferred if ligand_lines_preferred else ligand_lines_other
        
        if ligand_lines:
            ligand_destination.parent.mkdir(parents=True, exist_ok=True)
            with ligand_destination.open("w", encoding="utf-8") as handle:
                handle.writelines(ligand_lines)
            
            ligand_info = f"resname={ref_ligand}" if ref_ligand else "all ligands"
            if ligand_lines_preferred:
                chain_info = f"chain={ligand_chain}"
                print(f"    Extracted {len(ligand_lines)} ligand atoms ({ligand_info}, {chain_info}) to {ligand_destination.name}")
            elif ligand_lines_other:
                # Found in different chain - extract chain info
                other_chains = set()
                for line in ligand_lines_other:
                    other_chains.add(line[chain_col].strip())
                chain_info = f"chain={ligand_chain} (not found), using chain(s)={','.join(sorted(other_chains))}"
                print(f"    WARNING: Extracted {len(ligand_lines)} ligand atoms ({ligand_info}, {chain_info}) to {ligand_destination.name}")
            else:
                chain_info = f"chain={ligand_chain}" if ligand_chain else "all chains"
                print(f"    Extracted {len(ligand_lines)} ligand atoms ({ligand_info}, {chain_info}) to {ligand_destination.name}")
        else:
            ligand_info = f"resname={ref_ligand}" if ref_ligand else "any ligand"
            chain_info = f"chain={ligand_chain}" if ligand_chain else "any chain"
            print(f"    Warning: No ligand atoms found in {source.name} ({ligand_info}, {chain_info})")


def _parse_meeko_padding_pairs(message: str) -> List[Tuple[str, str]]:
    """Residue pairs from Meeko 'Expected N paddings for (A:1, A:2)' errors."""
    return [
        (a.strip(), b.strip())
        for a, b in re.findall(r"paddings for \(([^,]+), ([^)]+)\)", message)
    ]


def _prepare_receptor_meeko(cleaned_pdb: Path, pdbqt_path: Path, pdb_path_out: Path, ligand_pdb_path: Optional[Path] = None) -> None:
    """
    Prepare receptor using Meeko Python API.
    Also process and save ligand if provided.
    """
    try:
        # Create templates and molecule preparation
        templates = ResidueChemTemplates.create_from_defaults()
        mk_prep = MoleculePreparation.from_config({})
        
        # Read PDB file and create Polymer
        with open(cleaned_pdb, 'r') as f:
            pdb_string = f.read()

        # pdb2pqr can place hydrogens in steric clash; Meeko then treats the
        # contact as a covalent bond and raises "Expected 2 paddings". Drop
        # those non-template bonds and retry.
        bonds_to_delete: List[Tuple[str, str]] = []
        polymer = None
        last_error: Optional[Exception] = None
        for _ in range(12):
            try:
                polymer = Polymer.from_pdb_string(
                    pdb_string,
                    templates,
                    mk_prep,
                    allow_bad_res=True,  # Automatically remove residues that don't match templates
                    default_altloc='A',   # Use first alternative location by default
                    bonds_to_delete=bonds_to_delete or None,
                )
                break
            except Exception as exc:
                last_error = exc
                pairs = _parse_meeko_padding_pairs(str(exc))
                new_pairs = [
                    pair for pair in pairs
                    if pair not in bonds_to_delete and (pair[1], pair[0]) not in bonds_to_delete
                ]
                if not new_pairs:
                    raise
                for pair in new_pairs:
                    print(
                        f"    [PREPARATION] Dropping spurious Meeko bond {pair[0]}-{pair[1]} "
                        "(pdb2pqr hydrogen clash, not a real covalent link)"
                    )
                    bonds_to_delete.append(pair)
        if polymer is None:
            raise last_error if last_error is not None else RuntimeError("Meeko Polymer.from_pdb_string failed")
        
        # Get PDBQT through PDBQTWriterLegacy
        pdbqt_tuple = PDBQTWriterLegacy.write_from_polymer(polymer)
        rigid_pdbqt, flex_pdbqt_dict = pdbqt_tuple
        
        # Save rigid part (main receptor) in PDBQT
        with open(pdbqt_path, 'w') as f:
            f.write(rigid_pdbqt)
        
        # Save cleaned and prepared PDB file
        pdb_string = polymer.to_pdb()
        with open(pdb_path_out, 'w') as f:
            f.write(pdb_string)
        
        # If there are flexible residues, save them separately
        if flex_pdbqt_dict:
            flex_out_path = pdbqt_path.parent / (pdbqt_path.stem + "_flex.pdbqt")
            all_flex_pdbqt = "".join(flex_pdbqt_dict.values())
            with open(flex_out_path, 'w') as f:
                f.write(all_flex_pdbqt)
        
        # Check if ligand was extracted (ligand file is already created by _generate_clean_pdb)
        if ligand_pdb_path is not None and ligand_pdb_path.exists():
            # Count ligand atoms to verify
            ligand_atom_count = 0
            with open(ligand_pdb_path, 'r') as f:
                for line in f:
                    if line.startswith("HETATM"):
                        ligand_atom_count += 1
            
            if ligand_atom_count > 0:
                print(f"    [PREPARATION] Ligand extracted: {ligand_pdb_path.name} ({ligand_atom_count} atoms)")
            else:
                print(f"    [PREPARATION] Warning: Ligand file exists but contains no atoms: {ligand_pdb_path.name}")
        
        print(f"    [PREPARATION] Protein {pdbqt_path.stem}: prepared using Meeko Python API")
        
    except Exception as e:
        raise RuntimeError(f"Meeko Python API receptor preparation failed: {e}") from e


def prepare_proteins(
    protein_dir: Path,
    processed_dir: Path,
    settings: Optional[Dict[str, Any]] = None,
    interaction_config: Optional[Dict] = None,
    overwrite: bool = False,
) -> Dict[str, Path]:
    """
    Prepare all proteins for docking.
    
    Args:
        protein_dir: Directory with input PDB files.
        processed_dir: Directory to save prepared proteins.
        settings: Optional per-protein settings (chain, include_het, etc.).
        interaction_config: Optional interaction config with ref_ligand and safe_chain info.
        overwrite: Whether to overwrite existing files.
    
    Returns:
        Dictionary mapping protein names to PDBQT paths.
    """
    settings = settings or {}
    protein_output_dir = processed_dir / "proteins"
    clean_dir = protein_output_dir / "cleaned"
    protein_output_dir.mkdir(parents=True, exist_ok=True)
    clean_dir.mkdir(parents=True, exist_ok=True)
    
    # Build mapping from interaction config: protein -> (ref_ligand, safe_chain)
    protein_ligand_map = {}
    if interaction_config:
        from ..utils.settings import get_protein_ligand_pairs
        pairs = get_protein_ligand_pairs(interaction_config)
        for protein_name, _, ref_ligand, safe_chain in pairs:
            normalized_protein = _normalize_protein_id(protein_name)
            protein_ligand_map[normalized_protein] = (ref_ligand, safe_chain)
    
    results = {}
    
    # Find all PDB files
    pdb_files = list(protein_dir.glob("*.pdb"))
    if not pdb_files:
        print(f"  Warning: No PDB files found in {protein_dir}")
        return results

    # If interaction_config specifies proteins, prepare only those (faster and avoids touching unrelated proteins).
    if interaction_config and isinstance(interaction_config, dict):
        proteins_filter = interaction_config.get("protein", [])
        if proteins_filter:
            wanted = set(str(p).lower() for p in proteins_filter if str(p).strip())
            pdb_files = [p for p in pdb_files if p.stem.lower() in wanted]
            if not pdb_files:
                print(f"  Warning: interaction_config provided proteins={sorted(wanted)} but no matching PDBs found in {protein_dir}")
                return results

    print(f"  Preparing {len(pdb_files)} proteins...")
    
    for pdb_path in pdb_files:
        protein_id = pdb_path.stem
        normalized_id = _normalize_protein_id(protein_id)
        
        # Get settings
        chain_id = _get_setting(protein_id, settings, "chain", "A") or "A"
        include_het = _get_setting(protein_id, settings, "include_ligands", False)
        include_cofactors = _get_setting(protein_id, settings, "include_cofactors", False)
        include_waters = _get_setting(protein_id, settings, "include_waters", False)
        keep_het_resnames = _get_setting(protein_id, settings, "keep_het_resnames", None)
        drop_het_resnames = _get_setting(protein_id, settings, "drop_het_resnames", None)
        protonate_cfg = _get_setting(protein_id, settings, "protonate", None)
        # Normalize single-string values to list
        if isinstance(keep_het_resnames, str):
            keep_het_resnames = [keep_het_resnames]
        if isinstance(drop_het_resnames, str):
            drop_het_resnames = [drop_het_resnames]
        
        # Get ref_ligand and safe_chain from interaction config
        ref_ligand = None
        ligand_chain = None
        if normalized_id in protein_ligand_map:
            ref_ligand, ligand_chain = protein_ligand_map[normalized_id]
            # Use safe_chain for ligand extraction, fallback to chain_id if not specified
            if ligand_chain:
                ligand_chain = ligand_chain.strip()
            else:
                ligand_chain = chain_id
        
        # Generate cleaned PDB
        cleaned_pdb = clean_dir / f"{normalized_id}_chain{chain_id}.pdb"
        ligand_pdb = clean_dir / f"{normalized_id}_ligand.pdb"
        pdbqt_path = protein_output_dir / f"{normalized_id}.pdbqt"
        pdb_path_out = protein_output_dir / f"{normalized_id}.pdb"
        
        # Check if already prepared
        ligand_exists = ligand_pdb.exists() if ligand_pdb else False
        if pdbqt_path.exists() and pdb_path_out.exists() and not overwrite:
            if not ligand_exists:
                print(f"    Protein {protein_id} prepared but ligand file missing, regenerating ligand extraction...")
                # Re-extract ligand from original PDB
                try:
                    _generate_clean_pdb(
                        pdb_path,
                        cleaned_pdb,
                        chain_id=chain_id,
                        include_het=include_het,
                        include_cofactors=include_cofactors,
                        include_waters=include_waters,
                        keep_het_resnames=keep_het_resnames,
                        drop_het_resnames=drop_het_resnames,
                        ligand_destination=ligand_pdb,
                        ref_ligand=ref_ligand,
                        ligand_chain=ligand_chain,
                        protein_id=protein_id,
                    )
                except Exception as e:
                    print(f"    Warning: Failed to extract ligand: {e}")
            else:
                print(f"    Protein {protein_id} already prepared, skipping")
                results[protein_id] = pdbqt_path
                continue
        
        try:
            # Generate cleaned PDB and extract ligand
            _generate_clean_pdb(
                pdb_path,
                cleaned_pdb,
                chain_id=chain_id,
                include_het=include_het,
                include_cofactors=include_cofactors,
                include_waters=include_waters,
                keep_het_resnames=keep_het_resnames,
                drop_het_resnames=drop_het_resnames,
                ligand_destination=ligand_pdb,
                ref_ligand=ref_ligand,
                ligand_chain=ligand_chain,
                protein_id=protein_id,
            )
            
            # Optional pH-based protonation before Meeko. This is best-effort and will
            # fall back to the original cleaned PDB if the tool is not available.
            meeko_input_pdb = _maybe_protonate_before_meeko(
                protein_id=protein_id,
                cleaned_pdb=cleaned_pdb,
                include_cofactors=include_cofactors,
                protonate_cfg=protonate_cfg,
                overwrite=overwrite,
            )
            
            # Prepare with Meeko (ligand will be processed if exists)
            _prepare_receptor_meeko(meeko_input_pdb, pdbqt_path, pdb_path_out, ligand_pdb_path=ligand_pdb)
            results[protein_id] = pdbqt_path
            
        except Exception as e:
            print(f"    Error preparing protein {protein_id}: {e}")
            continue
    
    print(f"  Successfully prepared {len(results)}/{len(pdb_files)} proteins")
    return results

