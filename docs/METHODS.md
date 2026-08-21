# AffiTox methods and metrics

Python package import remains `docking_benchmark2`. Display name is AffiTox.

## Four axes (keep separate)

| Axis | Estimand | Confirmatory in manuscript | Code |
|------|----------|----------------------------|------|
| Scoring | Pearson \(r\) of method score vs experimental \(\mathrm{p}K_i\) | yes (target-wise \(r\); panel summary = **median** of within-target \(r\)) | `src/analysis/correlations.py`, `plots/correlations.py` |
| Ranking | Spearman \(\rho\), Kendall \(\tau\) | no (exploratory) | same |
| Screening | \(\mathrm{nEF}_{f}\) at fractions 1/5/10% | yes for active \(\mathrm{nEF}_{10}\); inactive \(\mathrm{nEF}_{10,\mathrm{low}}\) is exploratory | `src/analysis/enrichment.py` |
| Pose validity | PoseBusters dock-schema checks | exploratory profiles; not a substitute for scoring | `src/analysis/posebusters.py` |

ROC-AUC and BEDROC20 are **computed** in `enrichment.py` and appear in inferential families. They are not the confirmatory screening estimand. Confirmatory screening is active \(\mathrm{nEF}_{10}\).

Sign flips for Pearson/Spearman/Kendall are **fixed per metric** (`SIGN_FLIP_METRICS` in `src/analysis/constants.py`). Screening ranking direction uses `detect_score_direction` in `enrichment.py` (inactive set reverses ranking).

## Methods (not equivalent engines)

| Method | Class | Primary column | Score direction for correlation | PoseBusters |
|--------|-------|----------------|---------------------------------|-------------|
| QVina2 | classical docking | `qvina_affinity_bestpose` | lower-better \(\rightarrow\) flipped | yes |
| GNINA 1.3 | hybrid docking + CNN | `gnina_cnn_affinity_bestpose` | higher-better | yes |
| DynamicBind | learned pose generator | `dynamicbind_affinity_bestpose` | higher-better | yes |
| Boltz-2 | learned pose + affinity | `boltz2_affinity_pred_value` | lower-better \(\rightarrow\) flipped | yes |
| PLAPT | sequence-SMILES affinity | `plapt_affinity` | higher-better | **no** |

Third-party pins: [SOFTWARE_REGISTRY.md](SOFTWARE_REGISTRY.md).

## Inferential tests

Confirmatory: target-wise Pearson \(r\) and active \(\mathrm{nEF}_{10}\), paired Wilcoxon across 16 targets, Holm within each family (`src/analysis/inferential.py`).

Exploratory: Spearman, Kendall, inactive \(\mathrm{nEF}_{10,\mathrm{low}}\), PoseBusters profiles, cognate RMSD redocking (QVina2/GNINA), RMSE/MAE, ROC-AUC, BEDROC20.

## nEF

\[
\mathrm{EF}_f = \frac{H/T}{A/N},\quad
\mathrm{nEF}_f = \frac{\mathrm{EF}_f}{\mathrm{EF}_f^{\max}} = \frac{H}{\min(A,T)}
\]

with \(T=\lceil f N\rceil\). If \(A\ge T\), nEF is precision at fraction \(f\); if \(A<T\), nEF is recall. Active: \(K_i<1000\) nM, stronger-first. Inactive: \(K_i\ge 1000\) nM, weaker-first.

## PoseBusters summaries (do not mix)

1. Pass-all: fraction of poses with all 20 checks True (`pass_rate_all`).
2. At least \(N\) checks: cumulative % ligands with `passed >= N`.
3. Mean per-check compliance: average fraction of checks passed (`mean_frac_pass`). This is **not** a valid-pose rate.
4. Per-check pass rates: one heatmap row per check.

## Cognate RMSD redocking anchor

CASF-style geometry check on drug-like co-crystal ligands (`src/analysis/cognate_rmsd_redocking.py`):

- Tables (when generated): `analysis/tables/cognate_rmsd_redocking.csv`, `cognate_rmsd_summary.csv`
- **Not** invoked by `run_article_analysis.sh`
- Primary AffiTox pose-validity endpoint remains PoseBusters on the BindingDB panel
- DynamicBind / Boltz-2 cognate redock deferred; PLAPT has no poses

## RMSE/MAE

OLS in-sample calibration of score \(\rightarrow\) experimental pKi, then residual RMSE/MAE. Not raw-score error.
