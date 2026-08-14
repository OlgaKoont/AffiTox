# AffiTox methods and metrics

Python package import remains `docking_benchmark2` (historical). Display name is AffiTox.

## Methods (not equivalent docking engines)

| Method | Class | Primary column (`src/analysis/constants.py`) | Score direction for correlation | PoseBusters |
|--------|-------|-----------------------------------------------|---------------------------------|-------------|
| QVina2 | classical docking | `qvina_affinity_bestpose` | lower-better → flipped | yes |
| GNINA 1.3 | hybrid docking + CNN | `gnina_cnn_affinity_bestpose` | higher-better | yes |
| DynamicBind | learned pose generator | `dynamicbind_affinity_bestpose` | higher-better | yes |
| Boltz-2 | learned pose + affinity | `boltz2_affinity_pred_value` | lower-better → flipped | yes |
| PLAPT | sequence–SMILES affinity | `plapt_affinity` | higher-better | **no** |

Sign flips are **fixed per metric**, not chosen per target (`SIGN_FLIP_METRICS` in `src/analysis/constants.py`). Screening uses `detect_score_direction` in `src/analysis/enrichment.py` (inactive set reverses ranking).

## Confirmatory vs exploratory

Confirmatory (manuscript): target-wise Pearson \(r\) and active \(\mathrm{nEF}_{10}\), paired Wilcoxon across 16 targets, Holm within each family.

Exploratory: Spearman, Kendall, inactive \(\mathrm{nEF}_{10,\mathrm{low}}\), PoseBusters profiles, RMSE/MAE.

## nEF

\[
\mathrm{EF}_f = \frac{H/T}{A/N},\quad
\mathrm{nEF}_f = \frac{\mathrm{EF}_f}{\mathrm{EF}_f^{\max}} = \frac{H}{\min(A,T)}
\]

with \(T=\lceil f N\rceil\). If \(A\ge T\), nEF is precision at fraction \(f\); if \(A<T\), nEF is recall. Active: \(K_i<1000\) nM, stronger-first. Inactive: \(K_i\ge 1000\) nM, weaker-first.

## PoseBusters summaries (do not mix)

1. Pass-all: fraction of poses with all 20 checks True (`pass_rate_all`).
2. At least \(N\) checks: cumulative `% ligands` with `passed >= N`.
3. Mean per-check compliance: average fraction of checks passed (`mean_frac_pass`). **Not** “valid pose rate”.
4. Per-check pass rates: one heatmap row per check.

## RMSE/MAE

OLS in-sample calibration of score → experimental pKi, then residual RMSE/MAE. Not raw-score error.
