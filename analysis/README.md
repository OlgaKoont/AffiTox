# AffiTox Analysis

This directory contains the post-docking analysis layer for AffiTox: manuscript-grade tables and figures for a 16-target toxicity panel across 5 methods.  
It is designed for users who already have merged per-target prediction tables and PoseBusters outputs and need reproducible, publication-ready scoring/ranking/screening/pose-quality analysis.  
What is analyzed here: (i) affinity agreement and rank preservation (Pearson/Spearman/Kendall), (ii) early recognition metrics for actives/inactives (EF/nEF, ROC-AUC, BEDROC), (iii) pose plausibility profiles, and (iv) cross-method inferential statistics with uncertainty and multiplicity correction (bootstrap/permutation, Wilcoxon, BH/Holm).  
Methodological anchors and model context: [Boltz-2](https://doi.org/10.1101/2025.06.14.659707), [PoseBusters](https://doi.org/10.1038/s41586-024-07487-1), [Early-recognition metrics (EF/BEDROC)](https://doi.org/10.1021/ci600426e), [ROC analysis](https://doi.org/10.1016/j.patrec.2005.10.010), [BH-FDR](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x), [Holm correction](https://www.jstor.org/stable/4615733).


## Repository layout

```text
analysis/
├── config/                   analysis-specific defaults/env
├── tables/                   canonical CSV outputs
│   ├── correlations/
│   ├── enrichment/
│   ├── inferential/
│   └── posebuster/
├── figures/                  generated figures for manuscript/SI
└── logs/                     run logs (including SI/full runs)
```

## Minimal Runnable Example

```bash
bash ../pipeline/postprocess/run_article_analysis.sh
```

- Input directories:
  - merged tables: `analysis/tables/merged_ligands_docking_<pdb>.csv`
  - pose checks: `analysis/tables/posebuster/posebusters_results_<pdb>_<method>.csv`
- Output directories:
  - tables: `analysis/tables/*`
  - figures: `analysis/figures/*`
- Help for full CLI:

```bash
python ../src/analysis/run_pipeline.py --help
```

## Output Contract

| File / field | Meaning | Units / scale | Use case |
|---|---|---|---|
| `tables/correlations/summary_all_proteins.csv` + `pearson_r`, `spearman_rho`, `kendall_tau` | score-vs-pKi association per target/method | correlation coefficients | scoring/ranking quality |
| `tables/correlations/correlations_with_ci_perm.csv` + `*_ci_low/high`, `*_perm_p` | uncertainty and nonparametric significance | CI bounds, permutation p | robustness of associations |
| `tables/enrichment/ef_summary_all_proteins.csv` + `nEF10_active`, `nEF10_low` | early enrichment for active/inactive compounds | normalized EF in [0,1] | screening triage behavior |
| `tables/posebusters/pass_rates_by_target_method.csv` + `pass_rate_all` | fraction of poses passing all checks | proportion in [0,1] | docking plausibility |
| `tables/inferential/pairwise_method_tests.csv` + `p_holm`, `median_diff_a_minus_b` | cross-target method-vs-method inference | adjusted p + effect | confirmatory global comparisons |
| `tables/inferential/per_target_pairwise_tests.csv` + `perm_p_bh_target`, `effect_obs` | per-target pairwise inference | BH-adjusted p + effect | target-specific method behavior |

Interpretation scope:
- Use correlation outputs for continuous affinity fidelity.
- Use nEF outputs for early retrieval/deprioritization scenarios.
- Use PoseBusters to assess geometric/chemical plausibility of predicted poses.
- Do not use a single metric as a stand-alone quality verdict.

## Reproducibility

- Default seed: `RANDOM_SEED=42`.
- SI profile: `N_BOOTSTRAP=10000`, `N_PERMUTATION=10000` (`pipeline/postprocess/run_si_full.sh`).
- Deterministic components: table aggregation and plotting given fixed inputs and seed.
- Partially stochastic components: bootstrap/permutation statistics (seed-controlled).
- Reproducible out-of-the-box artifacts:
  - all CSV tables in `analysis/tables/`
  - all figure assets in `analysis/figures/` (PNG/SVG).

## Evaluation Status

- Implemented now:
  - scoring/ranking correlations + CI/permutation;
  - screening EF/nEF summaries;
  - PoseBusters pass-rate summaries and per-check heatmaps;
  - cross-target and per-target pairwise inferential tests.
- Planned / evolving:
  - additional endpoint families can be added in `src/analysis/inferential.py`;
  - figure layout/style refinements continue in `src/analysis/plots/`.
- Relevant scripts:
  - pipeline wrappers: `pipeline/postprocess/run_article_analysis.sh`, `pipeline/postprocess/run_si_full.sh`
  - core code: `src/analysis/`
  - postprocess bridge: `pipeline/postprocess/`.

## GitHub release bundle

The repository tracks **`analysis/tables/`** and **`analysis/figures/`** (~111 MB total) as the canonical manuscript/SI analysis bundle. Logs under `analysis/logs/` stay local.

| Path | Contents | Size (approx.) |
|---|---|---|
| `analysis/tables/` | CSV summaries, merged ligand scores, inferential tests, overlap audit | ~20 MB |
| `analysis/figures/` | PNG + SVG for heatmaps, scatters, nEF, PoseBusters, Wilcoxon | ~92 MB |

SI LaTeX PNGs are a **flat copy** of figures from `analysis/figures/`:

```bash
bash manuscript/supplementary/build_si.sh   # copies *.png → manuscript/supplementary/images/
```

To publish on GitHub after updating analysis:

```bash
git add analysis/tables analysis/figures analysis/README.md analysis/config
git status analysis/
git commit -m "Add AffiTox analysis tables and figure bundle for manuscript reproduction."
git push origin main
```

Optional Zenodo zip of the same tree:

```bash
bash pipeline/postprocess/package_zenodo.sh   # creates toxdock-analysis-artifacts.zip
```


| Zenodo link | What it contains | Notes |
|---|---|---|
| [10.5281/zenodo.20825057](https://doi.org/10.5281/zenodo.20825057) (part 1/6, includes BindingDB snapshot) | raw inputs + selected `results/` | draft/prereserved; publish before citing as a stable record |
| [AffiTox Zenodo community](https://zenodo.org/communities/affitox/) | parts 2–6 docking zips | DOIs `10.5281/zenodo.20825059` … `20825067` |

## Citation, License, Contribution

### Citation

If you use this analysis pipeline, cite:
- the benchmark manuscript (update DOI once published),
- the model papers you benchmark,
- and the metric/statistics references corresponding to the endpoints you report (EF/nEF, ROC-AUC/BEDROC, correlation, Wilcoxon, BH/Holm, bootstrap/permutation).

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


### License

MIT (`LICENSE`). Third-party data and method weights: `NOTICE`.

### Contribution

- Open an issue/PR in this repository.
- Keep metric names and method labels aligned with manuscript/SI terminology.
- Add tests or reproducibility notes for any new endpoint or plotting logic.
