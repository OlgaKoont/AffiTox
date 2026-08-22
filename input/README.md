# AffiTox input data

Canonical inputs for prepare/dock:

- `input/proteins/<pdb>.pdb` or `.cif`
- `input/ligands_nodubl/*_nodubl.csv`

Upstream BindingDB dumps in `input/bindingdb/` are provenance only and are gitignored. There is no curation CLI in `src/`. Rules: [docs/DATA_DICTIONARY.md](../docs/DATA_DICTIONARY.md).

Panel PDB IDs: `TARGETS` in `config/project.env.sh`. Extra structures may sit in `proteins/`; only selected targets are consumed when pairing is configured.

```bash
source config/project.env.sh
bash run_pipeline.sh prepare   # uses config/interaction_protein_ligand_16target.json
```

Paths in `config/toxdock_config.yaml`: `protein_dir: input/proteins`, `ligand_dir: input/ligands_nodubl`.

| Path / field | Meaning |
|--------------|---------|
| `proteins/<pdb>.pdb` or `.cif` | target structure |
| `ligands_nodubl/*_nodubl.csv` `canonical_smiles` | RDKit canonical SMILES |
| `standard_type` | endpoint (`Ki`) |
| `standard_value` | \(K_i\) nM |
| `pchembl_value` | p-scale potency when present |
| `*_nodubl_grouped.json` | grouped metadata |
| `duplicates_report.json` | removed replicate IDs |

The 16 panel CSVs currently contain 11,180 data rows. The manuscript reports 10,497 \(K_i\) records. Do not equate the two without an explicit filter.

Raw docking archives: [10.5281/zenodo.20825057](https://doi.org/10.5281/zenodo.20825057) (part 1 includes a BindingDB snapshot). License: MIT for code; BindingDB/ChEMBL/PDB terms for source records (`NOTICE`).
