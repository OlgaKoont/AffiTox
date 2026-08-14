# AffiTox Input Data

This directory contains the source input assets for the AffiTox benchmark pipeline.  
It is intended for users who want to run preparation and docking stages from raw structures and curated ligand tables.  
In this repository, the canonical input contract is: protein structures in `input/proteins`, per-target ligand/activity tables in `input/ligands_nodubl`, and upstream raw BindingDB dumps in `input/bindingdb`.  
The active benchmark panel used downstream is the curated 16-target set from `config/project.env.sh` (`TARGETS`).  
`input/bindingdb` is retained for provenance and regeneration workflows; routine benchmarking reads from curated `ligands_nodubl` files.

## Repository layout

```text
input/
├── proteins/                 raw target structures (.pdb/.cif)
├── ligands_nodubl/           curated per-target ligand/activity tables
│   ├── *_nodubl.csv
│   ├── *_nodubl_grouped.json
│   └── duplicates_report.json
└── bindingdb/                upstream raw BindingDB TSV exports
```

## Minimal Runnable Example

From repository root:

```bash
source config/project.env.sh
bash run_pipeline.sh prepare
```

Expected inputs used by this stage:
- proteins: `input/proteins/<pdb>.pdb` or `input/proteins/<pdb>.cif`
- ligands: `input/ligands_nodubl/*_nodubl.csv`

Main output location:
- prepared assets: `processed/`

Check config resolving these paths:

```bash
python -c "import yaml;print(yaml.safe_load(open('config/toxdock_config.yaml'))['protein_dir']);print(yaml.safe_load(open('config/toxdock_config.yaml'))['ligand_dir'])"
```

## Input Contract

| Path / field | Meaning | Units / scale | Used by |
|---|---|---|---|
| `proteins/<pdb>.pdb` or `proteins/<pdb>.cif` | target 3D structure | structural coordinates | `prepare` stage |
| `ligands_nodubl/*_nodubl.csv` + `canonical_smiles` | ligand chemical structure | SMILES string | ligand prep and docking |
| `ligands_nodubl/*_nodubl.csv` + `standard_type` | assay endpoint label (`Ki`, `IC50`, etc.) | categorical | filtering/traceability |
| `ligands_nodubl/*_nodubl.csv` + `standard_value` | reported activity value | typically nM in source | activity normalization |
| `ligands_nodubl/*_nodubl.csv` + `pchembl_value` | transformed potency (`-log10(M)`) | p-scale | ranking/scoring reference |
| `ligands_nodubl/*_nodubl_grouped.json` | grouped ligands by activity bins | grouped records | reproducible grouping metadata |
| `ligands_nodubl/duplicates_report.json` | removed duplicate summary | counts + IDs | curation QC |
| `bindingdb/*.tsv` | raw upstream BindingDB exports | source table | provenance only |

Interpretation scope:
- `ligands_nodubl` is the benchmark-ready curation layer and should be used for production runs.
- `bindingdb` is not benchmark-ready without additional filtering/curation.
- Filename prefix (target name) is metadata; canonical target selection is controlled by `TARGETS`.

## Reproducibility

- Canonical path mapping is defined in `config/toxdock_config.yaml`:
  - `protein_dir: input/proteins`
  - `ligand_dir: input/ligands_nodubl`
- Canonical target panel is defined in `config/project.env.sh` (`TARGETS`).
- `random_state: 42` is set in `config/toxdock_config.yaml` for deterministic sampling where applicable.
- Mixed file presence is allowed in `proteins/` (panel targets + extra structures); only selected `TARGETS` are consumed in standard runs.

## Evaluation Status

- Implemented now:
  - benchmark-ready protein structures and curated ligand tables are present;
  - duplicate-removed ligand datasets and grouped JSON descriptors are present;
  - raw BindingDB dumps are retained for provenance.
- Planned / evolving:
  - optional extension of curated target set beyond current 16-target panel;
  - additional curation reports may be added if preprocessing rules evolve.

## Zenodo Artifacts

| Zenodo link | What it contains | Notes |
|---|---|---|
| [10.5281/zenodo.20825057](https://doi.org/10.5281/zenodo.20825057) | BindingDB snapshot + part of raw docking | prereserved draft DOI |
| [docs/DATA_DICTIONARY.md](../docs/DATA_DICTIONARY.md) | curated `ligands_nodubl` schema | use this for column definitions |

## Citation, License, Contribution

### Citation

If you reuse these input assets or curation conventions, cite:
- BindingDB data source,
- ChEMBL source where applicable in curated files,
- benchmark manuscript (when DOI is finalized).

### License

Follow the repository-level license and the original terms of third-party datasets (BindingDB/ChEMBL).  
When redistributing subsets, preserve attribution and source identifiers.

### Contribution

- Keep `input/` changes traceable (what changed, why, and source).
- Do not overwrite curated files without updating curation metadata/report artifacts.
- When adding new targets, update both `input/` contents and path/target config in `config/`.
