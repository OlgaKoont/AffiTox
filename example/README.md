# AffiTox example (smoke-test pipeline)

Mini-dataset: **16 targets × 2 ligands** → full pipeline into `example/`.

## Layout

```text
example/
├── input/           mini CSV + protein PDB symlinks
├── boltz/           Boltz-2 staging (CIF, MSA, mini ligands)
├── processed/       prepared structures
├── results/         docking outputs (all methods incl. Boltz after sync)
├── analysis/        merged tables + figures
├── logs/            SLURM logs
├── run_example_pipeline.sh
└── run_example_pipeline.sbatch
```

Code and configs stay in repo root (`src/`, `config/`, `pipeline/`).

## Quick start (local)

```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
export PYTHON="$(command -v python)"

# 1) Inputs only
bash example/run_example_pipeline.sh prepare-inputs boltz-prepare

# 2) Full pipeline (skip pip install)
SKIP_INSTALL=1 bash example/run_example_pipeline.sh all
```

## SLURM (recommended for dock + analysis)

```bash
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)"

# Dock + postprocess (prepare already done)
sbatch --export=ALL,EXAMPLE_STAGES="dock boltz-sync merge posebusters analysis" \
  example/run_example_pipeline.sbatch

# Full run from scratch + Boltz GPU predict
sbatch --export=ALL,BOLTZ_SKIP_RUN=0 \
  example/run_example_pipeline.sbatch
```

## Boltz-2

| Variable | Default | Meaning |
|----------|---------|---------|
| `BOLTZ_ROOT` | `.../boltz/data` | External Boltz workspace |
| `BOLTZ_WRITE_DIR` | `${BOLTZ_ROOT}/results` | Where `boltz predict` writes |
| `BOLTZ_SKIP_RUN` | `1` | Skip GPU predict in `all` / sbatch |
| `BOLTZ_SYNC_MODE` | `symlink` | Import into `example/results/` |
| `BOLTZ_USE_EXTERNAL=1` | off | Merge reads `boltz/data/results` directly |

Manual Boltz on GPU node:

```bash
bash example/scripts/run_boltz2_example.sh
bash example/run_example_pipeline.sh boltz-sync merge posebusters analysis
```

## Outputs

- Merged tables: `example/analysis/tables/`
- PoseBusters: `example/analysis/tables/posebuster/`
- Figures: `example/analysis/figures/`
