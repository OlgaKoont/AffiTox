# AffiTox pipeline

The public switch is [`run_pipeline.sh`](../run_pipeline.sh). It sources [`config/project.env.sh`](../config/project.env.sh).

```text
Usage: bash run_pipeline.sh [stage ...]

Stages (default: analysis):
  install      pip install -e .
  prepare      protein/ligand/box preparation
  dock         run docking methods
  merge        merge tables + pValue
  posebusters  PoseBusters pass-rate tables
  analysis     article statistics + figures
  all          prepare -> dock -> merge -> posebusters -> analysis
```

With no arguments the driver runs **analysis**, not `all`.

Manuscript tables and figures come from `src/analysis/run_pipeline.py` via [`pipeline/postprocess/run_article_analysis.sh`](../pipeline/postprocess/run_article_analysis.sh). The older CLI stages `aggregation` / `analysis` on `python -m docking_benchmark2.cli.run_benchmark` are a different, narrower code path. Do not use them for the paper figures.

## Flow

```text
input/proteins + input/ligands_nodubl
        |  prepare
        v
    processed/
        |  dock: qvina, gnina, plapt, dynamicbind
        v
    results/                    Boltz-2 (external `boltz predict`)
        |                                |
        +---------------+----------------+
                        |  merge
                        v
        analysis/tables/merged_ligands_docking_<pdb>.csv
                        |  posebusters (`bust`)
                        v
        analysis/tables/posebuster/
                        |  analysis
                        v
        analysis/tables/*  and  analysis/figures/*
```

## Stage notes

**install.** Editable pip install only. [`pipeline/install/README.md`](../pipeline/install/README.md).

**prepare.** Shared structures for all methods. Needs an interaction JSON that lists PDB codes and ligand CSV stems. Tracked `config/toxdock_config.yaml` has `interaction_config_file: null`. The loader defaults to `config/interaction_protein_ligand.json`. If that file is missing, pairing is empty. [`pipeline/prepare/README.md`](../pipeline/prepare/README.md).

**dock.** Methods in `config/toxdock_config.yaml`: `qvina`, `gnina`, `plapt`, `dynamicbind`. Boltz-2 is not in that list. `METHODS` in `project.env.sh` is consumed by **analysis**, not by `run_docking.sh`. Subset engines with `--methods`. [`pipeline/dock/README.md`](../pipeline/dock/README.md).

**merge.** `src/analysis/merge_ligands_docking_from_dir.py` then `add_pvalue_column.py`. Optional `--boltz-results-dir` / `BOLTZ_RESULTS_DIR`. DynamicBind subdirectory default: `dynamicbind_new`. [`pipeline/postprocess/README.md`](../pipeline/postprocess/README.md).

**posebusters.** `prepare_and_run_posebusters.py`; Boltz helper only if `BOLTZ_RESULTS_DIR` is set. PLAPT excluded. Needs `bust`.

**analysis.** [`analysis/config/defaults.sh`](../analysis/config/defaults.sh). Skip pieces with `RUN_FIGURES=0` and similar. SI: `bash pipeline/postprocess/run_si_full.sh`. [`analysis/README.md`](../analysis/README.md).

**all.** prepare through analysis. No install. No Boltz inference.

## Config files

| File | Role |
|------|------|
| [`config/project.env.sh`](../config/project.env.sh) | paths, `TARGETS`, `METHODS` (analysis), `BOLTZ_RESULTS_DIR` |
| [`config/toxdock_config.yaml`](../config/toxdock_config.yaml) | input dirs, labox, dock method list |
| [`config/methods_config.yaml`](../config/methods_config.yaml) | per-engine flags and empty tool paths |
| [`config/protein_settings_keep_cofactors_v2.yaml`](../config/protein_settings_keep_cofactors_v2.yaml) | receptor preparation |
| [`analysis/config/defaults.sh`](../analysis/config/defaults.sh) | bootstrap, colours, stage switches |
