# Prepare stage

```bash
source config/project.env.sh
bash pipeline/prepare/run_prepare.sh
# or
bash run_pipeline.sh prepare
```

Calls `python -m docking_benchmark2.cli.run_benchmark --stage preparation` with `TOXDOCK_CONFIG` and `METHODS_CONFIG`.

Expected outputs when protein-ligand pairing is configured: `processed/proteins/*.pdbqt`, `processed/ligands/<target>/<dataset>/*.pdbqt`, `processed/boxes/*.json`, `processed/plapt_sequences/*.txt`.

Pairing file (tracked): [`config/interaction_protein_ligand_16target.json`](../../config/interaction_protein_ligand_16target.json), set in [`config/toxdock_config.yaml`](../../config/toxdock_config.yaml) as `interaction_config_file`. Override with `--interaction-config`.

Settings that change downstream docking: `config/protein_settings_keep_cofactors_v2.yaml`, `labox` in `toxdock_config.yaml`, `random_state: 42`.
