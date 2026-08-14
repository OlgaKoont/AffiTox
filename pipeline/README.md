# AffiTox Pipeline

This directory contains executable stage wrappers for reproducible end-to-end benchmark runs.  
It is intended for users who need a stable operational entrypoint from environment setup to manuscript-grade analysis artifacts.  
The stage axis covered here is: install -> prepare -> dock -> postprocess (merge, PoseBusters, inferential analysis and figures).  
Methodological anchors used downstream include [Boltz-2](https://doi.org/10.1101/2025.06.14.659707), [PoseBusters](https://doi.org/10.1038/s41586-024-07487-1), [EF/BEDROC](https://doi.org/10.1021/ci600426e), [ROC](https://doi.org/10.1016/j.patrec.2005.10.010), [BH-FDR](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x), and [Holm correction](https://www.jstor.org/stable/4615733).

## Repository layout

```text
pipeline/
├── install/                 environment/bootstrap wrappers
├── prepare/                 structure/ligand/box preparation wrappers
├── dock/                    method execution wrappers
└── postprocess/             merge, PoseBusters, analysis wrappers
```

## Minimal run

```bash
source config/project.env.sh
bash run_pipeline.sh all
```

Common stage-by-stage execution:

```bash
bash run_pipeline.sh install prepare dock merge posebusters analysis
```

Show top-level help:

```bash
bash run_pipeline.sh --help
```

## Output contract

| Stage | Main input | Main output | Use case |
|---|---|---|---|
| `install` | `environment.yml`, project package | active environment + editable install | reproducible runtime |
| `prepare` | `input/` | `processed/` | method-ready structures/ligands/boxes |
| `dock` | `processed/` | `results/` | per-method predictions/poses and metrics |
| `merge` | `results/`, optional Boltz-2 external results | `analysis/tables/merged_ligands_docking_<pdb>.csv` | unified panel table |
| `posebusters` | `results/` and `processed/proteins` | `analysis/tables/posebuster/*.csv` | pose plausibility assessment |
| `analysis` | merged tables + PoseBusters tables | `analysis/tables/*`, `analysis/figures/*` | manuscript-grade statistics/plots |

## Interpretation scope

- `pipeline/` orchestrates execution; scientific interpretation belongs to `analysis/` outputs.
- Stage success means artifacts were generated, not that methods are scientifically superior.
- Compare methods using unified downstream statistics, not raw stage logs alone.

## Reproducibility

Core environment/path variables (from `config/project.env.sh`):
- `TOXAFFINITY_ROOT`, `PYTHON`, `TOXDOCK_CONFIG`, `METHODS_CONFIG`
- `INPUT_DIR`, `PROCESSED_DIR`, `RESULTS_DIR`, `ANALYSIS_ROOT`
- `TARGETS`, `METHODS`, `BOLTZ_RESULTS_DIR`

Parameter variants and expected impact:
- `TARGETS` subset (e.g., `TARGETS="1g5m 3eyg"`): reduces panel size; faster but changes all aggregate statistics.
- `METHODS` subset (e.g., `METHODS="gnina qvina"`): removes method comparisons for omitted methods.
- `PYTHON` executable override: can change dependency stack and numeric behavior if environments differ.
- Alternative `TOXDOCK_CONFIG` / `METHODS_CONFIG`: changes preprocessing, method hyperparameters, and raw score distributions.
- `BOLTZ_RESULTS_DIR` unset/set: controls whether Boltz-2 results participate in merged/analysis tables.

## Evaluation status

- Implemented now:
  - complete stage wrappers for install/prepare/dock/postprocess;
  - end-to-end run via `run_pipeline.sh all`.
- Planned / evolving:
  - additional docking methods and postprocessing adapters can extend stage outputs.

## Citation, license, contribution

### Citation

Cite method-specific papers for methods you run and statistics you report (see `analysis/README.md` for full metric/statistics references).

### License

Follow the MIT license in `LICENSE` and third-party terms in `NOTICE`.

### Contribution

- Keep wrappers shell-only and stage-focused.
- Document any new environment variable or CLI pass-through in the corresponding stage README.
