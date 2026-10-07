# AffiTox Zenodo record 04 — merged tables

One CSV per target: `merged_ligands_docking_<pdb>.csv`.

Join key: `(protein_pdb_id, ligand_id)` only. Missing method scores stay empty (never imputed).
Panel N after ligand_id merge (no SMILES collapse): **11180**.

Seed: 42 (documented for docking; QVina2 logs may print a different engine seed — see SOFTWARE_REGISTRY.md).
