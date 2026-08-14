# AffiTox

AffiTox is a toxicity-oriented docking and binding-affinity benchmark: **16 BindingDB targets**, experimental \(K_i\), and **five** methods (QVina2, GNINA 1.3, DynamicBind, Boltz-2, PLAPT).

Code: [github.com/OlgaKoont/AffiTox](https://github.com/OlgaKoont/AffiTox). Historical identifiers kept so existing scripts still run: Python package `docking_benchmark2`, env `TOXAFFINITY_ROOT`, CLI `toxdock-pipeline` (alias `affitox-pipeline`).

## Paper

Konovalova et al., *Limited and Target-Dependent Transferability of Docking Methods to Toxicity-Linked Targets*, Journal of Chemical Information and Modeling (in review). See `CITATION.cff`.

## What this repository can reproduce

- **Analysis-only (supported without authors):** scoring, ranking, screening, and PoseBusters *aggregations* from `analysis/tables/` via `pipeline/postprocess/run_article_analysis.sh`.
- **Full docking (not from `environment.yml` alone):** protein/ligand prep and method binaries/weights. See [docs/INSTALL.md](docs/INSTALL.md).

## Methods

| Method | Role | Primary score | pKi correlation | Poses |
|--------|------|---------------|-----------------|-------|
| QVina2 | classical docking | `qvina_affinity_bestpose` | sign-flipped (lower-better) | yes |
| GNINA 1.3 | hybrid CNN affinity | `gnina_cnn_affinity_bestpose` | as reported (higher-better) | yes |
| DynamicBind | learned poses | `dynamicbind_affinity_bestpose` | as reported | yes |
| Boltz-2 | learned poses + affinity | `boltz2_affinity_pred_value` | sign-flipped | yes |
| PLAPT | seq+SMILES affinity | `plapt_affinity` | as reported | **no** (excluded from PoseBusters) |

Details: [docs/METHODS.md](docs/METHODS.md).

## Panel

Targets (`config/project.env.sh`): `1g5m 2z5x 3eyg 3jy9 3lxk 3mjg 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 7awe 7kk3`.

Activity: \(K_i\) nM → \(\mathrm{p}K_i=9-\log_{10}(K_i)\). Screening: active \(K_i<1000\) nM, inactive \(\ge 1000\) nM. Curation: [docs/DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md).

## Quick start (analysis)

Linux x86_64. Requires Python 3.10+ and the deposited tables.

```bash
git clone https://github.com/OlgaKoont/AffiTox.git
cd AffiTox
python3 -m pip install -e ".[analysis]"
export TOXAFFINITY_ROOT="$(pwd)"
bash pipeline/postprocess/run_article_analysis.sh
```

Outputs: `analysis/tables/`, `analysis/figures/`. SI bootstrap: `bash pipeline/postprocess/run_si_full.sh`.

Tests: `PYTHONPATH=src pytest tests -q`.

## Full pipeline

```bash
source config/project.env.sh
bash run_pipeline.sh install    # pip install -e .  — not GNINA/Boltz
bash run_pipeline.sh prepare
bash run_pipeline.sh dock
bash run_pipeline.sh merge
bash run_pipeline.sh posebusters
bash run_pipeline.sh analysis
```

Method install, GPU, and Zenodo raw results: [docs/INSTALL.md](docs/INSTALL.md). Mini run: [example/README.md](example/README.md).

## Metrics (confirmatory)

- Scoring: Pearson \(r\) vs experimental pKi (sign-harmonized).
- Screening: active \(\mathrm{nEF}_{10}\). Cross-target Wilcoxon + Holm.
- Exploratory: Spearman, Kendall, \(\mathrm{nEF}_{10,\mathrm{low}}\), PoseBusters pass-all vs mean compliance (mean ≠ valid-pose rate).

## Layout

```text
AffiTox/
├── config/                 path contract (TOXAFFINITY_ROOT) + YAML
├── input/                  proteins + ligands_nodubl
├── processed/              prepared structures
├── results/                raw docking (gitignored; Zenodo)
├── src/
│   ├── docking_benchmark2/ docking adapters (import path unchanged)
│   └── analysis/           metrics and plots
├── pipeline/               install / prepare / dock / postprocess
├── analysis/tables|figures manuscript bundle
├── docs/                   INSTALL, METHODS, DATA_DICTIONARY, FIGURES
├── tests/                  scientific invariant tests
└── run_pipeline.sh
```

## Data

| Location | Role |
|----------|------|
| `input/ligands_nodubl/` | curated Ki panel |
| `analysis/tables/` | merged scores + statistics |
| [Zenodo AffiTox community](https://zenodo.org/communities/affitox/) | raw docking parts 1–6 (prereserved DOIs `10.5281/zenodo.20825057` … `20825067`) |

## Requirements

- OS: Linux x86_64 (HPC). Windows unsupported.
- Analysis: CPU; 10k bootstrap/permutation is the expensive step.
- Docking: QVina2/GNINA CPU; Boltz-2, DynamicBind, PLAPT typically GPU; tens of GB for `results/`.

## Limitations

MW ≤ 500 Da; Ki-only nM records; no cognate RMSD (PoseBusters instead); PLAPT is not docking; complete-case \(N\) differs by method.

## Licence and citation

This repository’s source code is released under the [MIT License](LICENSE). BindingDB/ChEMBL/PDB records and third-party method weights are **not** covered by MIT; see [NOTICE](NOTICE).

Please cite the AffiTox manuscript (JCIM, in review) using [CITATION.cff](CITATION.cff). Bug reports and questions: [github.com/OlgaKoont/AffiTox/issues](https://github.com/OlgaKoont/AffiTox/issues).
