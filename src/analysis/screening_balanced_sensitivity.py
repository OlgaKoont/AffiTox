#!/usr/bin/env python3
"""Part D: balanced-class resampling sensitivity for nEF10 (changes the estimand)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.constants import PRIMARY_METRICS
from analysis.data import activity_mask, load_merged
from analysis.enrichment import (
    compute_auc_roc,
    compute_bedroc,
    compute_ef,
    detect_score_direction,
)
from analysis.robustness_common import (
    ALL_METHODS,
    MERGED_DIR,
    ROBUST_TABLES,
    TARGETS,
    RobustnessConfig,
    method_label,
    target_label,
)


def _active_inactive_masks(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    act = activity_mask(df, "active")
    low = activity_mask(df, "low")
    return np.asarray(act, dtype=bool), np.asarray(low, dtype=bool)


def _metrics(is_active: np.ndarray, scores: np.ndarray, set_name: str, metric: str) -> dict:
    direction = detect_score_direction(metric, set_name)
    ef, nef, n, a, h = compute_ef(is_active, scores, 0.10, direction)
    top_n = max(1, int(np.ceil(0.10 * n))) if n else 0
    precision = (h / top_n) if top_n else np.nan
    recall = (h / a) if a else np.nan
    return {
        "nEF10": nef,
        "precision_at_10": precision,
        "recall_at_10": recall,
        "roc_auc": compute_auc_roc(is_active, scores, direction),
        "bedroc20": compute_bedroc(is_active, scores, direction, alpha=20.0),
        "N": n,
        "A": a,
        "H_top": h,
        "direction": direction,
    }


def run_balanced_screening(cfg: RobustnessConfig) -> Path:
    cfg.ensure_dirs()
    rng = np.random.default_rng(cfg.seed)
    rows = []
    for target in TARGETS:
        df = load_merged(MERGED_DIR, target)
        for method_id in ALL_METHODS:
            metric = PRIMARY_METRICS[method_id]
            if metric not in df.columns:
                continue
            scores = pd.to_numeric(df[metric], errors="coerce").to_numpy(dtype=float)
            act, low = _active_inactive_masks(df)
            ok = np.isfinite(scores) & (act | low)
            scores = scores[ok]
            act = act[ok]
            low = low[ok]
            n_act = int(act.sum())
            n_inact = int(low.sum())
            unb_act = _metrics(act, scores, "active", metric)
            unb_low = _metrics(low, scores, "low", metric)
            n_bal = min(n_act, n_inact)
            rec = {
                "target": target,
                "target_label": target_label(target),
                "method_id": method_id,
                "method": method_label(method_id),
                "n_active_unbalanced": n_act,
                "n_inactive_unbalanced": n_inact,
                "n_balanced_per_class": n_bal,
                "unbalanced_nEF10_active": unb_act["nEF10"],
                "unbalanced_nEF10_low": unb_low["nEF10"],
                "unbalanced_precision_at_10_active": unb_act["precision_at_10"],
                "unbalanced_recall_at_10_active": unb_act["recall_at_10"],
                "unbalanced_roc_auc_active": unb_act["roc_auc"],
                "unbalanced_bedroc20_active": unb_act["bedroc20"],
                "active_direction": unb_act["direction"],
                "inactive_direction": unb_low["direction"],
                "analysis_label": "class_balance_sensitivity_changes_estimand",
            }
            if n_bal < 5:
                rec.update(
                    {
                        "balanced_nEF10_active_mean": np.nan,
                        "balanced_nEF10_active_ci_low": np.nan,
                        "balanced_nEF10_active_ci_high": np.nan,
                        "balanced_nEF10_low_mean": np.nan,
                        "balanced_nEF10_low_ci_low": np.nan,
                        "balanced_nEF10_low_ci_high": np.nan,
                        "balanced_precision_at_10_active_mean": np.nan,
                        "balanced_recall_at_10_active_mean": np.nan,
                        "balanced_roc_auc_active_mean": np.nan,
                        "balanced_bedroc20_active_mean": np.nan,
                        "n_successful_replicates": 0,
                        "skip_reason": "min(n_active,n_inactive)<5",
                    }
                )
                rows.append(rec)
                continue
            act_idx = np.flatnonzero(act)
            inact_idx = np.flatnonzero(low)
            store = {k: [] for k in (
                "nEF10_active",
                "nEF10_low",
                "precision_at_10_active",
                "recall_at_10_active",
                "roc_auc_active",
                "bedroc20_active",
            )}
            for _ in range(cfg.n_balanced):
                ia = rng.choice(act_idx, size=n_bal, replace=False)
                ii = rng.choice(inact_idx, size=n_bal, replace=False)
                idx = np.concatenate([ia, ii])
                s = scores[idx]
                a = act[idx]
                low_s = low[idx]
                m_act = _metrics(a, s, "active", metric)
                m_low = _metrics(low_s, s, "low", metric)
                store["nEF10_active"].append(m_act["nEF10"])
                store["nEF10_low"].append(m_low["nEF10"])
                store["precision_at_10_active"].append(m_act["precision_at_10"])
                store["recall_at_10_active"].append(m_act["recall_at_10"])
                store["roc_auc_active"].append(m_act["roc_auc"])
                store["bedroc20_active"].append(m_act["bedroc20"])
            def _summ(key: str) -> tuple[float, float, float]:
                arr = np.asarray(store[key], dtype=float)
                arr = arr[np.isfinite(arr)]
                if len(arr) == 0:
                    return np.nan, np.nan, np.nan
                return float(np.mean(arr)), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))

            mu, lo, hi = _summ("nEF10_active")
            mu_l, lo_l, hi_l = _summ("nEF10_low")
            rec.update(
                {
                    "balanced_nEF10_active_mean": mu,
                    "balanced_nEF10_active_ci_low": lo,
                    "balanced_nEF10_active_ci_high": hi,
                    "balanced_nEF10_low_mean": mu_l,
                    "balanced_nEF10_low_ci_low": lo_l,
                    "balanced_nEF10_low_ci_high": hi_l,
                    "balanced_precision_at_10_active_mean": _summ("precision_at_10_active")[0],
                    "balanced_recall_at_10_active_mean": _summ("recall_at_10_active")[0],
                    "balanced_roc_auc_active_mean": _summ("roc_auc_active")[0],
                    "balanced_bedroc20_active_mean": _summ("bedroc20_active")[0],
                    "n_successful_replicates": cfg.n_balanced,
                    "skip_reason": "",
                    "n_balanced": cfg.n_balanced,
                    "seed": cfg.seed,
                }
            )
            rows.append(rec)
    out = pd.DataFrame(rows)
    path = ROBUST_TABLES / "balanced_class_nef_sensitivity.csv"
    out.to_csv(path, index=False)
    note = """# Balanced-class nEF10 sensitivity

This is a **sensitivity analysis**. Sampling equal numbers of actives and inactives
without replacement **changes the target population and the estimand**.
Balanced nEF10 is therefore not a drop-in replacement for the confirmatory
unbalanced AffiTox nEF10.

Active enrichment uses stronger-predicted-binding-first ranking; inactive
enrichment uses the reversed score direction (`detect_score_direction`).
"""
    (ROBUST_TABLES / "balanced_class_nef_sensitivity.md").write_text(note)
    return path


if __name__ == "__main__":
    run_balanced_screening(RobustnessConfig())
