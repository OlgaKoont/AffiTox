# AffiTox installation

Two reproduction modes exist. Do not expect `conda env create -f environment.yml` to install GNINA, DynamicBind, Boltz-2, PLAPT, or PoseBusters.

| Mode | What you can reproduce | Hardware |
|------|------------------------|----------|
| **Analysis-only** | Pearson / Spearman / Kendall, nEF, Wilcoxon/Holm, PoseBusters aggregations, manuscript figures from deposited tables | Linux x86_64 CPU, several GB RAM (10k bootstrap is slow) |
| **Full docking** | prepare → dock/infer → merge → PoseBusters → analysis | Linux HPC; GPU for Boltz-2, DynamicBind, PLAPT |

Windows is unsupported.

The environment variable `TOXAFFINITY_ROOT` is the AffiTox repository root (historical name; do not rename it or scripts will break).

## 1. Analysis-only (recommended first)

```bash
git clone https://github.com/OlgaKoont/AffiTox.git
cd AffiTox
conda env create -f environment.yml   # env name: docking
conda activate docking
export PYTHON="$(command -v python)"
bash run_pipeline.sh install
```

`environment.yml` pins Python 3.10, RDKit, OpenBabel, NumPy, pandas, matplotlib, Meeko, and CUDA toolkit headers. It does **not** install docking method binaries.

If CUDA packages are unwanted on a CPU machine, install the Python stack from `setup.py` instead:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[analysis]"
export PYTHON=python
```

Then, with `analysis/tables/` present (Git clone or Zenodo):

```bash
export TOXAFFINITY_ROOT="$(pwd)"
bash pipeline/postprocess/run_article_analysis.sh
```

Full SI profile (`N_BOOTSTRAP=10000`, `N_PERMUTATION=10000`):

```bash
bash pipeline/postprocess/run_si_full.sh
```

Smoke tests (no GPU, no docking):

```bash
PYTHONPATH=src pytest tests -q
```

## 2. Method-specific docking software

Install each tool into its own conda env or prefix, then point `config/methods_config.yaml`
or environment variables (`PLAPT_PATH`, `DYNAMICBIND_PATH`, `GNINA_BIN`, `BOLTZ_ROOT`) at
**your** paths. Cluster-local examples live in `config/methods_config.hpc.example.yaml`.
Exact third-party pins: [docs/SOFTWARE_REGISTRY.md](SOFTWARE_REGISTRY.md).

### QVina2 (classical)

- Binary expected: `qvina02` (`config/methods_config.yaml` → `qvina.binary`)
- Related conda pins: `vina=1.2.6`, `mgltools=1.5.7`, `openbabel=3.1.1` in `environment.yml`
- Published settings: `exhaustiveness=8`, `num_modes=10`, `random_seed=42`
- Primary score: `qvina_affinity_bestpose` (lower-better; sign-flipped for correlation)

### GNINA 1.3 (hybrid)

- Binary expected: `gnina`
- Published AffiTox runs: `--no_gpu` only (`use_cnn: false`). GNINA 1.3 still writes CNN columns.
- Primary score: `gnina_cnn_affinity_bestpose` (higher-better; **not** sign-flipped)
- Timeout: 600 s/ligand
- PoseBusters: included

### DynamicBind (learned poses)

- Separate env and checkpoint tree; set `dynamicbind.dynamicbind_path` and Python interpreters in `methods_config.yaml`
- Primary score: `dynamicbind_affinity_bestpose` (higher-better)
- GPU recommended
- PoseBusters: included (`dynamicbind_new` CSV label)

### Boltz-2 (learned poses + affinity)

- External `boltz predict` workspace. Set `BOLTZ_RESULTS_DIR` before merge.
- Primary score: `boltz2_affinity_pred_value` (lower-better; sign-flipped)
- GPU required for inference
- PoseBusters: included (mmCIF export)

### PLAPT (affinity only)

- Sequence + canonical SMILES; **no 3D pose**
- Set `plapt.plapt_path` to a WELP-PLAPT checkout
- Primary score: `plapt_affinity` (higher-better)
- **Excluded from PoseBusters** (`POSEBUSTERS_EXCLUDE`)

### PoseBusters

- CLI must be on `PATH` for `pipeline/postprocess/run_posebusters.sh`
- Deposited CSVs use the 20-check PoseBusters **dock** schema (≥ 0.3). The exact `bust --version` is not stored in those CSVs; see [SOFTWARE_REGISTRY.md](SOFTWARE_REGISTRY.md).

## 3. Full pipeline

```bash
source config/project.env.sh
export BOLTZ_RESULTS_DIR=/path/to/boltz/data/results   # if using Boltz-2

bash run_pipeline.sh install
bash run_pipeline.sh prepare
bash run_pipeline.sh dock
bash run_pipeline.sh merge
bash run_pipeline.sh posebusters
bash run_pipeline.sh analysis
```

Mini smoke dataset (16 targets × 2 ligands): `example/README.md`.

## 4. Data archives

Raw docking outputs exceed Git. The six Zenodo records below are **published and open** (not drafts):

| Part | DOI | Contents |
|------|-----|----------|
| 1/6 | [10.5281/zenodo.20825057](https://doi.org/10.5281/zenodo.20825057) | BindingDB snapshot + 7awe, 3mjg, 1g5m, 4f65, 4zau, 4tz4, 2z5x |
| 2/6 | [10.5281/zenodo.20825059](https://doi.org/10.5281/zenodo.20825059) | 6jok, 3lxk, 5jkv |
| 3/6 | [10.5281/zenodo.20825061](https://doi.org/10.5281/zenodo.20825061) | 4ase, 6gqj |
| 4/6 | [10.5281/zenodo.20825063](https://doi.org/10.5281/zenodo.20825063) | 7kk3, 5mo4 |
| 5/6 | [10.5281/zenodo.20825065](https://doi.org/10.5281/zenodo.20825065) | 3jy9 |
| 6/6 | [10.5281/zenodo.20825067](https://doi.org/10.5281/zenodo.20825067) | 3eyg |

Community: [zenodo.org/communities/affitox](https://zenodo.org/communities/affitox/).
File-level MD5 checksums: [docs/zenodo/affitox_data_manifest.tsv](zenodo/affitox_data_manifest.tsv).
Release procedure: [docs/ZENODO_RELEASE.md](ZENODO_RELEASE.md). Data availability wording: [docs/DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).

Analysis tables and figures are intended to ship with GitHub (`analysis/tables/`, `analysis/figures/`).
