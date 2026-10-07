# Postprocess stage

```bash
source config/project.env.sh
bash pipeline/postprocess/run_merge.sh
bash pipeline/postprocess/run_posebusters.sh
bash pipeline/postprocess/run_article_analysis.sh
```

SI resampling profile:

```bash
bash pipeline/postprocess/run_si_full.sh   # N_BOOTSTRAP=10000 N_PERMUTATION=10000, logs under analysis/logs/
```

| Script | Input | Output |
|--------|-------|--------|
| `run_merge.sh` | `results/`, `processed/id_maps/`, optional `BOLTZ_RESULTS_DIR` | `analysis/excluding_2z5x_3mjg/tables/merged_ligands_docking_<pdb>.csv` joined on `ligand_id` |
| `run_posebusters.sh` | poses in `results/` plus `bust`, or deposited CSVs if those are missing | `analysis/tables/posebuster/*.csv` |
| `run_article_analysis.sh` | merged CSVs + PoseBusters CSVs | correlations, enrichment, inferential tables, figures |
| `run_si_full.sh` | same | same with forced 10k/10k |

Analysis overrides (`analysis/config/defaults.sh`):

```bash
N_BOOTSTRAP=300 N_PERMUTATION=300 bash pipeline/postprocess/run_article_analysis.sh
RUN_FIGURES=0 bash pipeline/postprocess/run_article_analysis.sh
TARGETS="1g5m 3eyg 6gqj" METHODS="gnina qvina" \
  bash pipeline/postprocess/run_article_analysis.sh
```

`DYNAMICBIND_POSEBUSTERS_LABEL` default is `dynamicbind_new`. Mismatch with the on-disk subdirectory drops DynamicBind pose rows.

Tracked Python helpers in this folder (`rename_zenodo_affitox.py`, `upload_zenodo.py`) are for archive naming/upload, not for manuscript metrics.

Details: [docs/PIPELINE.md](../../docs/PIPELINE.md), [analysis/README.md](../../analysis/README.md).
