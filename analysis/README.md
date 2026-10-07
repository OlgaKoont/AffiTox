# AffiTox analysis

Post-docking statistics and figures for the 16-target panel. Entry point:

```bash
bash pipeline/postprocess/run_article_analysis.sh   # canonical root: analysis/excluding_2z5x_3mjg/
```

CLI (flags required if you call Python directly):

```bash
PYTHONPATH=src python src/analysis/run_pipeline.py --help
```

Inputs:

- `analysis/tables/merged_ligands_docking_<pdb>.csv`
- `analysis/tables/posebuster/posebusters_results_<pdb>_<method>.csv`

Outputs: `analysis/tables/{correlations,enrichment,inferential,posebusters}/` and `analysis/figures/`.

Do not mix axes. Scoring is Pearson \(r\). Ranking is Spearman/Kendall. Screening confirmatory estimand is active \(\mathrm{nEF}_{10}\). Pose validity is PoseBusters. ROC-AUC and BEDROC20 are computed in `src/analysis/enrichment.py` but are exploratory.

## Layout

```text
analysis/
├── config/defaults.sh        N_BOOTSTRAP, N_PERMUTATION, colours, stage switches
├── tables/                   canonical CSVs (tracked)
└── figures/                  PNG + SVG (tracked)
```

`analysis/logs/` is gitignored.

## Output contract

| File / field | Axis | Meaning |
|--------------|------|---------|
| `tables/correlations/summary_all_proteins.csv` (`pearson_r`) | scoring | within-target Pearson vs pKi |
| same (`spearman_rho`, `kendall_tau`) | ranking | rank correlations |
| `tables/correlations/correlations_with_ci_perm.csv` | scoring/ranking | bootstrap CI, permutation p |
| `tables/enrichment/ef_summary_all_proteins.csv` (`nEF10_active`) | screening | confirmatory early enrichment of actives |
| same (`nEF10_low`) | screening | exploratory weak-binding-end enrichment |
| `tables/posebusters/pass_rates_by_target_method.csv` (`pass_rate_all`) | pose validity | fraction passing all checks |
| `tables/inferential/pairwise_method_tests.csv` | all families | Wilcoxon + Holm across 16 targets |
| `tables/inferential/per_target_pairwise_tests.csv` | all families | per-target permutation + BH |

Cognate co-crystal RMSD (`src/analysis/cognate_rmsd_redocking.py`) is optional and is **not** run by `run_article_analysis.sh`. Figure generators: [docs/FIGURES.md](../docs/FIGURES.md). Metric definitions: [docs/METHODS.md](../docs/METHODS.md).

## Reproducibility

- Seed: `RANDOM_SEED=42`
- SI: `bash pipeline/postprocess/run_si_full.sh` (`N_BOOTSTRAP=10000`, `N_PERMUTATION=10000`)
- Fast check: `N_BOOTSTRAP=10 N_PERMUTATION=10 TARGETS=1g5m bash pipeline/postprocess/run_article_analysis.sh`

GitHub Actions uses that fast check with `ANALYSIS_ROOT` pointed at a scratch directory and `MERGED_DATA_DIR` at the deposited tables.

## Data

Raw docking (not required for this folder): [10.5281/zenodo.20825057](https://doi.org/10.5281/zenodo.20825057) through [10.5281/zenodo.20825067](https://doi.org/10.5281/zenodo.20825067). Community: [zenodo.org/communities/affitox](https://zenodo.org/communities/affitox/).

## Citation

Cite the AffiTox manuscript ([CITATION.cff](../CITATION.cff)), the methods you compare, and the metrics you report.

```bibtex
@article{passaro2025boltz2,
  title   = {Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction},
  year    = {2025},
  journal = {bioRxiv},
  doi     = {10.1101/2025.06.14.659707}
}

@article{buttenschoen2024posebusters,
  title   = {PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences},
  year    = {2024},
  journal = {Nature},
  doi     = {10.1038/s41586-024-07487-1}
}

@article{truchon2007evaluation,
  author  = {Truchon, Jean-Francois and Bayly, Chris I.},
  title   = {Evaluating Virtual Screening Methods: Good and Bad Metrics for the ``Early Recognition'' Problem},
  journal = {Journal of Chemical Information and Modeling},
  year    = {2007},
  volume  = {47},
  number  = {2},
  pages   = {488--508},
  doi     = {10.1021/ci600426e}
}

@article{fawcett2006roc,
  author  = {Fawcett, Tom},
  title   = {An Introduction to ROC Analysis},
  journal = {Pattern Recognition Letters},
  year    = {2006},
  volume  = {27},
  number  = {8},
  pages   = {861--874},
  doi     = {10.1016/j.patrec.2005.10.010}
}

@article{holm1979simple,
  author  = {Holm, Sture},
  title   = {A Simple Sequentially Rejective Multiple Test Procedure},
  journal = {Scandinavian Journal of Statistics},
  year    = {1979},
  volume  = {6},
  number  = {2},
  pages   = {65--70}
}

@article{benjamini1995controlling,
  author  = {Benjamini, Yoav and Hochberg, Yosef},
  title   = {Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing},
  journal = {Journal of the Royal Statistical Society: Series B},
  year    = {1995},
  volume  = {57},
  number  = {1},
  pages   = {289--300},
  doi     = {10.1111/j.2517-6161.1995.tb02031.x}
}

@article{pearson1895note,
  author  = {Pearson, Karl},
  title   = {Note on Regression and Inheritance in the Case of Two Parents},
  journal = {Proceedings of the Royal Society of London},
  year    = {1895},
  volume  = {58},
  pages   = {240--242},
  doi     = {10.1098/rspl.1895.0041}
}

@article{spearman1904proof,
  author  = {Spearman, Charles},
  title   = {The Proof and Measurement of Association between Two Things},
  journal = {The American Journal of Psychology},
  year    = {1904},
  volume  = {15},
  number  = {1},
  pages   = {72--101},
  doi     = {10.2307/1412159}
}

@article{kendall1938new,
  author  = {Kendall, Maurice G.},
  title   = {A New Measure of Rank Correlation},
  journal = {Biometrika},
  year    = {1938},
  volume  = {30},
  number  = {1/2},
  pages   = {81--93},
  doi     = {10.1093/biomet/30.1-2.81}
}

@article{wilcoxon1945individual,
  author  = {Wilcoxon, Frank},
  title   = {Individual Comparisons by Ranking Methods},
  journal = {Biometrics Bulletin},
  year    = {1945},
  volume  = {1},
  number  = {6},
  pages   = {80--83},
  doi     = {10.2307/3001968}
}

@book{efron1993bootstrap,
  author    = {Efron, Bradley and Tibshirani, Robert J.},
  title     = {An Introduction to the Bootstrap},
  publisher = {Chapman and Hall/CRC},
  year      = {1993},
  doi       = {10.1007/978-1-4899-4541-9}
}

@book{good2005permutation,
  author    = {Good, Phillip},
  title     = {Permutation, Parametric and Bootstrap Tests of Hypotheses},
  publisher = {Springer},
  year      = {2005},
  edition   = {3},
  doi       = {10.1007/b138696}
}
```

License: MIT (`LICENSE`). Third-party data and weights: `NOTICE`.
