# Dock stage

```bash
source config/project.env.sh
bash pipeline/dock/run_docking.sh
# or
bash run_pipeline.sh dock
```

Calls `python -m docking_benchmark2.cli.run_benchmark --stage docking`. Methods in tracked `config/toxdock_config.yaml`: `qvina`, `gnina`, `plapt`, `dynamicbind`. Boltz-2 is **not** in this list.

`METHODS` from `project.env.sh` is for analysis. It does not subset this wrapper. To subset engines:

```bash
bash pipeline/dock/run_docking.sh --methods gnina qvina
```

Outputs: `results/<target>/docking/<method>/` and per-method metrics CSVs used by merge.

Hyperparameters: `config/methods_config.yaml`. Log-recovered pins (which can disagree with YAML): [docs/SOFTWARE_REGISTRY.md](../../docs/SOFTWARE_REGISTRY.md). Empty `plapt_path` / `dynamicbind_path` fail unless `PLAPT_PATH` / `DYNAMICBIND_PATH` are set.
