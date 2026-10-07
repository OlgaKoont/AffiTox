"""Pearson / Spearman / Kendall correlations with multiple-testing correction."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr

from .config import AnalysisConfig
from .data import load_merged, prepare_xy
from .inferential import benjamini_hochberg, bootstrap_ci, bootstrap_error_ci, calibrated_error, permutation_p


def compute_correlation_row(
    target: str,
    method_id: str,
    metric_col: str,
    x: np.ndarray,
    y: np.ndarray,
) -> dict:
    n = len(x)
    row = {
        "target": target.lower(),
        "method_id": method_id,
        "metric": metric_col,
        "n_points": n,
        "pearson_r": np.nan,
        "pearson_p": np.nan,
        "spearman_rho": np.nan,
        "spearman_p": np.nan,
        "kendall_tau": np.nan,
        "kendall_p": np.nan,
    }
    if n < 3:
        return row

    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)
    kt, kp = kendalltau(x, y)
    row.update(
        pearson_r=float(pr),
        pearson_p=float(pp),
        spearman_rho=float(sr),
        spearman_p=float(sp),
        kendall_tau=float(kt),
        kendall_p=float(kp),
    )
    return row


def run_correlations(cfg: AnalysisConfig) -> pd.DataFrame:
    rows: list[dict] = []
    for target in cfg.targets:
        df = load_merged(cfg.merged_dir, target)
        per_target_rows: list[dict] = []
        for method_id in cfg.methods:
            metric = cfg.primary_metric(method_id)
            if metric not in df.columns:
                continue
            x, y = prepare_xy(df, metric, cfg.exp_col)
            row = compute_correlation_row(target, method_id, metric, x, y)
            row["method_label"] = cfg.method_label(method_id)
            per_target_rows.append(row)
            rows.append(row)

        if per_target_rows:
            pd.DataFrame(per_target_rows).to_csv(
                cfg.tables_dir / "correlations" / "per_target" / f"correlations_{target}.csv",
                index=False,
            )

    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary

    for pcol, qcol in (
        ("pearson_p", "pearson_p_bh"),
        ("spearman_p", "spearman_p_bh"),
        ("kendall_p", "kendall_p_bh"),
    ):
        mask = summary[pcol].notna()
        summary[qcol] = np.nan
        if mask.any():
            summary.loc[mask, qcol] = benjamini_hochberg(summary.loc[mask, pcol].values)

    out = cfg.tables_dir / "correlations" / "summary_all_proteins.csv"
    summary.to_csv(out, index=False)
    return summary


def add_bootstrap_permutation(cfg: AnalysisConfig, summary: pd.DataFrame) -> pd.DataFrame:
    """Extend correlation summary with bootstrap CI and permutation p-values."""
    ext_rows: list[dict] = []
    for _, row in summary.iterrows():
        target = row["target"]
        metric = row["metric"]
        df = load_merged(cfg.merged_dir, target)
        xy = prepare_xy(df, metric, cfg.exp_col)
        if len(xy[0]) < 3:
            ext_rows.append(dict(row))
            continue

        x, y = xy
        rec = dict(row)
        p_lo, p_hi, _ = bootstrap_ci(
            x, y, method="pearson", n_boot=cfg.n_bootstrap, seed=cfg.random_seed
        )
        s_lo, s_hi, _ = bootstrap_ci(
            x, y, method="spearman", n_boot=cfg.n_bootstrap, seed=cfg.random_seed + 1
        )
        k_lo, k_hi, _ = bootstrap_ci(
            x, y, method="kendall", n_boot=cfg.n_bootstrap, seed=cfg.random_seed + 4
        )
        rmse_lo, rmse_hi, rmse_mean = bootstrap_error_ci(
            x, y, metric="rmse", n_boot=cfg.n_bootstrap, seed=cfg.random_seed + 5
        )
        mae_lo, mae_hi, mae_mean = bootstrap_error_ci(
            x, y, metric="mae", n_boot=cfg.n_bootstrap, seed=cfg.random_seed + 6
        )
        rec["pearson_ci_low"] = p_lo
        rec["pearson_ci_high"] = p_hi
        rec["pearson_perm_p"] = permutation_p(
            x, y, method="pearson", n_perm=cfg.n_permutation, seed=cfg.random_seed + 2
        )
        rec["spearman_ci_low"] = s_lo
        rec["spearman_ci_high"] = s_hi
        rec["spearman_perm_p"] = permutation_p(
            x, y, method="spearman", n_perm=cfg.n_permutation, seed=cfg.random_seed + 3
        )
        rec["kendall_ci_low"] = k_lo
        rec["kendall_ci_high"] = k_hi
        rec["kendall_perm_p"] = permutation_p(
            x, y, method="kendall", n_perm=cfg.n_permutation, seed=cfg.random_seed + 7
        )
        rec["calibrated_rmse"] = calibrated_error(x, y, "rmse")
        rec["calibrated_rmse_ci_low"] = rmse_lo
        rec["calibrated_rmse_ci_high"] = rmse_hi
        rec["calibrated_rmse_boot_mean"] = rmse_mean
        rec["calibrated_mae"] = calibrated_error(x, y, "mae")
        rec["calibrated_mae_ci_low"] = mae_lo
        rec["calibrated_mae_ci_high"] = mae_hi
        rec["calibrated_mae_boot_mean"] = mae_mean
        ext_rows.append(rec)

    extended = pd.DataFrame(ext_rows)
    out = cfg.tables_dir / "correlations" / "correlations_with_ci_perm.csv"
    extended.to_csv(out, index=False)
    return extended
