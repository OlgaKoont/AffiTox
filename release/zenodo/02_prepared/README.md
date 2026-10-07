# AffiTox Zenodo record 02 — prepared inputs

Curated ligand tables with stable `ligand_id` (ligand_0001, …) in docked-table row order,
Meeko receptors, search boxes, PLAPT sequences, and engine-native ID maps.

Join key for every later record is `(pdb_id, ligand_id)`. Historical files named `ligand_1`
or `idx_0` are aliases listed in `id_maps/`.

Place Boltz-2, DynamicBind, and PLAPT weight files under `weights/` (see weights/README.md).
