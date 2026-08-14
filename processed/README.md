# AffiTox Prepared Assets (`processed/`)

This directory stores deterministic intermediate artifacts produced by the `prepare` stage and consumed by docking methods.  
It is intended for users who already have raw inputs in `input/` and need method-ready protein, ligand, and box files.  
What is prepared here: (i) receptor structures (`.pdb`/`.pdbqt`), (ii) per-target ligand ensembles (`.sdf`/`.pdbqt`), (iii) docking boxes (`center`/`size` JSON), and (iv) optional method-specific helper artifacts (for example, `plapt_sequences`).  
In the pipeline contract, this directory is the bridge between `input/` and `results/`.

## Repository layout

```text
processed/
├── proteins/                 prepared receptors (.pdb/.pdbqt)
│   └── cleaned/              chain/ligand extracted intermediates
├── ligands/                  prepared ligand .pdbqt libraries
├── ligands_sdf/              intermediate ligand .sdf libraries
├── boxes/                    docking box JSON files per target
└── plapt_sequences/          extracted sequences for PLAPT
```

## Minimal Runnable Example

From repository root:

```bash
source config/project.env.sh
bash run_pipeline.sh prepare
```

Expected input:
- proteins: `input/proteins/<pdb>.pdb` (or `.cif` as source)
- ligands: `input/ligands_nodubl/*.csv`

Primary output:
- prepared assets: `processed/`

CLI help (stage runner):

```bash
python -m docking_benchmark2.cli.run_benchmark --help
```

## Output Contract

| File / field | Meaning | Units / scale | Use case |
|---|---|---|---|
| `proteins/<target>.pdbqt` | receptor prepared for docking engines | PDBQT coordinates + atom typing/charges | GNINA/QVina docking input |
| `proteins/<target>.pdb` | cleaned/prepared receptor in PDB | structural coordinates | inspection/debugging/reproducibility |
| `proteins/cleaned/<target>_chain<id>.pdb` | chain-filtered receptor before final conversion | structural coordinates | provenance of chain/cofactor filtering |
| `proteins/cleaned/<target>_ligand.pdb` | co-crystal ligand extracted from complex (if found) | ligand atom coordinates | reference pocket context, QA |
| `ligands/<target>/<dataset>/ligand_<n>.pdbqt` | prepared ligand conformers for docking | PDBQT coordinates + charges | method input libraries |
| `ligands_sdf/<target>/<dataset>/ligand_<n>.sdf` | intermediate 3D ligand structures | SDF coordinates | traceability and ligand QA |
| `boxes/<target>.json` + `center`, `size` | docking search region | Angstrom coordinates (`[x,y,z]`) | docking grid definition |
| `plapt_sequences/<target>.txt` | amino acid sequence used by PLAPT pre-step | one-letter sequence string | PLAPT-specific preprocessing |

Interpretation scope:
- `processed/` artifacts are **intermediate operational outputs**, not final benchmarking results.
- File presence does not imply method success; method execution outcomes are tracked under `results/` and `analysis/`.
- Docking box values are geometric setup parameters and should not be interpreted as quality metrics.

## Reproducibility

- Path contract is fixed by `config/toxdock_config.yaml` (`processed_dir: "processed"`).
- Preparation logic runs through `pipeline/prepare/run_prepare.sh` and `docking_benchmark2` preprocessing modules.
- Deterministic components:
  - directory/file naming contract,
  - target/labeled dataset mapping,
  - ligand embedding seed (`randomSeed=42`) in ligand preparation.
- Potentially environment-sensitive components:
  - protonation tool availability/options,
  - external binary/library versions affecting atom typing or charge assignment.
- Reproducible out-of-the-box artifacts:
  - `processed/proteins/`, `processed/ligands/`, `processed/ligands_sdf/`, `processed/boxes/`, and method helper subfolders when enabled.

## Evaluation Status

- Implemented now:
  - receptor preparation to `.pdb`/`.pdbqt`;
  - ligand preparation from curated CSVs to `.sdf` and `.pdbqt`;
  - per-target box generation via configured LaBOX-style workflow;
  - optional helper artifacts for enabled methods (for example, PLAPT sequences).
- Planned / evolving:
  - additional method-specific preprocessors may add new subfolders under `processed/`;
  - preparation heuristics (protonation/cofactor handling) can be refined as benchmark scope expands.

## Citation, License, Contribution

### Citation

If you use the prepared-asset workflow, cite the benchmark repository/manuscript and the upstream data sources:

```bibtex
@article{liu2007bindingdb,
  title   = {BindingDB: a web-accessible database of experimentally determined protein-ligand binding affinities},
  author  = {Liu, Tiqing and Lin, Yuhmei and Wen, Xian and Jorissen, Ryan N. and Gilson, Michael K.},
  journal = {Nucleic Acids Research},
  year    = {2007},
  volume  = {35},
  pages   = {D198--D201},
  doi     = {10.1093/nar/gkl999}
}

@article{burley2019rcsb,
  title   = {RCSB Protein Data Bank: biological macromolecular structures enabling research and education in fundamental biology, biomedicine, biotechnology and energy},
  author  = {Burley, Stephen K. and others},
  journal = {Nucleic Acids Research},
  year    = {2019},
  volume  = {47},
  number  = {D1},
  pages   = {D464--D474},
  doi     = {10.1093/nar/gky1004}
}
```

### License

Follow the repository-level license and the terms of original upstream datasets/tools used to generate these artifacts.

### Contribution

- Keep `processed/` schema stable (subfolder names and file naming conventions).
- If adding a new preprocessing artifact family, document it here in `Output Contract`.
- When changing preparation logic, include reproducibility notes (version, seed, and expected artifact diffs).
