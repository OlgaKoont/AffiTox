# AffiTox example (mini pipeline)

Intended layout: 16 targets \(\times\) 2 ligands, writing under `example/` rather than the production tree.

## Incomplete on GitHub

Git currently tracks:

- `example/README.md` (this file)
- `example/scripts/run_boltz2_example.sh`

It does **not** track `example/run_example_pipeline.sh`, `example/run_example_pipeline.sbatch`, `config/project.env.example.sh`, `config/toxdock_config.example.yaml`, or the mini `example/input/` tree. Commands below describe the **local** driver when those files are present. They are not a verified public-clone recipe.

For manuscript statistics, use [docs/INSTALL.md](../docs/INSTALL.md) (analysis from `analysis/tables/`).

## Local driver stages (when `run_example_pipeline.sh` exists)

```text
prepare-inputs   2 ligands/target + protein links
install          pip install -e . (skipped when SKIP_INSTALL=1)
prepare          protein/ligand/box preparation
dock             qvina, gnina, plapt, dynamicbind
boltz-prepare    stage mini ligands + CIF + MSA
boltz-run        boltz predict (GPU; skipped when BOLTZ_SKIP_RUN=1)
boltz-sync       import into example/results/
merge            merged tables + pValue
posebusters      PoseBusters
analysis         scoring, ranking, screening, pose aggregations, figures
all              the above (install still skipped unless SKIP_INSTALL=0)
```

```bash
# only if the untracked driver is on disk
SKIP_INSTALL=1 bash example/run_example_pipeline.sh all
```

Boltz GPU (tracked script, but it sources untracked `config/project.env.example.sh` and `config/boltz_example.env.sh`):

```bash
bash example/scripts/run_boltz2_example.sh
```

## Outputs (local tree)

- `example/analysis/tables/`
- `example/analysis/tables/posebuster/`
- `example/analysis/figures/`
