# Files to stage for the AffiTox v1.0.0 reproducibility commit

Do **not** `git add .`. The working tree also contains unrelated robustness WIP,
regenerated figures, and cluster-only scripts.

## Must include (analysis blockers)

- `analysis/config/defaults.sh`
- `src/analysis/constants.py`
- `src/analysis/plots/style.py`
- `src/analysis/run_pipeline.py`
- `src/analysis/property_baselines.py`
- `src/analysis/pki_range_audit.py`
- `src/analysis/config.py`
- `src/analysis/inferential.py`
- `src/analysis/plots/correlations.py`
- `src/analysis/plots/posebusters.py`
- `tests/test_analysis_imports.py`
- `.github/workflows/ci.yml`
- `setup.py`

## Must include (portable paths)

- `config/methods_config.yaml`
- `config/methods_config.hpc.example.yaml`
- `config/example_toolchain.env.sh`
- `config/boltz_example.env.sh`
- `example/scripts/run_boltz2_example.sh`
- `src/docking_benchmark2/docking/plapt.py`
- `src/docking_benchmark2/docking/dynamicbind.py`
- `src/analysis/cognate_rmsd_redocking.py`
- `src/analysis/prepare_and_run_posebusters_boltz2.py`

## Must include (docs / citation)

- `README.md`
- `NOTICE`
- `CITATION.cff`
- `docs/INSTALL.md`
- `docs/METHODS.md`
- `docs/SOFTWARE_REGISTRY.md`
- `docs/DATA_AVAILABILITY.md`
- `docs/RELEASE_FILELIST.md`
- `docs/ZENODO_RELEASE.md`
- `docs/zenodo/README.md`
- `docs/zenodo/affitox_data_manifest.tsv`
- `analysis/README.md`
- `input/README.md`

Do not create `.zenodo.json`. Do not publish GitHub `v1.0.0` until CI is green.
