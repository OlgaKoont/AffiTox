"""PoseBusters pass rates: 100%, 90%, 50% and per-check breakdown."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .constants import POSEBUSTERS_EXCLUDE
from .data import resolve_posebusters_csv


def _boolean_check_columns(df: pd.DataFrame) -> list[str]:
    skip = {"file", "molecule", "position"}
    cols: list[str] = []
    for c in df.columns:
        if c in skip:
            continue
        if pd.api.types.is_bool_dtype(df[c]):
            cols.append(c)
        else:
            s = df[c].dropna().astype(str).str.lower()
            if len(s) > 0 and s.isin(["true", "false"]).all():
                cols.append(c)
    return cols


def _to_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series
    return series.astype(str).str.lower().map({"true": True, "false": False})


def run_posebusters(cfg: AnalysisConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict] = []
    check_rows: list[dict] = []

    for target in cfg.targets:
        for method_id in cfg.methods:
            if method_id in POSEBUSTERS_EXCLUDE:
                continue
            csv_path = resolve_posebusters_csv(
                cfg.posebusters_dir,
                cfg.posebusters_boltz2_dir,
                target,
                method_id,
                cfg.dynamicbind_posebusters_label,
            )
            if csv_path is None:
                continue

            df = pd.read_csv(csv_path)
            bool_cols = _boolean_check_columns(df)
            if not bool_cols:
                continue

            bool_df = pd.DataFrame({c: _to_bool(df[c]) for c in bool_cols})
            frac_pass = bool_df.mean(axis=1)
            n_poses = len(bool_df)
            n_pass_all = int((frac_pass >= 1.0).sum())
            n_pass_90pct = int((frac_pass >= 0.90).sum())
            n_pass_50pct = int((frac_pass >= 0.50).sum())

            summary_rows.append(
                {
                    "target": target.lower(),
                    "method_id": method_id,
                    "method_label": cfg.method_label(method_id),
                    "n_poses": n_poses,
                    "n_checks": len(bool_cols),
                    "n_pass_all": n_pass_all,
                    "n_pass_90pct": n_pass_90pct,
                    "n_pass_50pct": n_pass_50pct,
                    "pass_rate_all": float(n_pass_all / n_poses) if n_poses else float("nan"),
                    "pass_rate_90pct": float(n_pass_90pct / n_poses) if n_poses else float("nan"),
                    "pass_rate_50pct": float(n_pass_50pct / n_poses) if n_poses else float("nan"),
                    "mean_frac_pass": float(frac_pass.mean()),
                }
            )

            for check in bool_cols:
                check_rows.append(
                    {
                        "target": target.lower(),
                        "method_id": method_id,
                        "method_label": cfg.method_label(method_id),
                        "check": check,
                        "pass_rate": float(_to_bool(df[check]).mean()),
                    }
                )

    summary = pd.DataFrame(summary_rows)
    checks = pd.DataFrame(check_rows)

    summary.to_csv(cfg.tables_dir / "posebusters" / "pass_rates_by_target_method.csv", index=False)
    checks.to_csv(cfg.tables_dir / "posebusters" / "pass_rates_by_check_target_method.csv", index=False)
    return summary, checks


def canonical_posebusters_checks(reference_csv: Path | None = None) -> list[str]:
    """Boolean PoseBusters check names from a results CSV (21 in the current deposit)."""
    if reference_csv is None:
        root = Path(__file__).resolve().parents[2]
        for candidate in (
            root / "analysis/excluding_2z5x_3mjg/tables/posebuster/posebusters_results_1g5m_boltz2.csv",
            root / "release/zenodo/05_posebusters/tables/posebusters_results_1g5m_boltz2.csv",
            root / "example/05_posebusters/tables/posebusters_results_1g5m_boltz2.csv",
        ):
            if candidate.exists():
                reference_csv = candidate
                break
    if reference_csv is None or not reference_csv.exists():
        raise FileNotFoundError("Missing PoseBusters reference CSV")
    df = pd.read_csv(reference_csv, nrows=1)
    return _boolean_check_columns(df)


def collect_pooled_pass_counts(cfg: AnalysisConfig) -> tuple[int, dict[str, list[int]]]:
    """Pool pass counts per pose across all targets (All targets).

    Each pose is scored on the boolean columns present in its own CSV, matching
    ``run_posebusters`` / ``n_pass_all``. Current panel CSVs have 21 checks.
    """
    pooled: dict[str, list[int]] = {
        m: [] for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE
    }
    n_checks_seen: list[int] = []

    for target in cfg.targets:
        for method_id in pooled:
            csv_path = resolve_posebusters_csv(
                cfg.posebusters_dir,
                cfg.posebusters_boltz2_dir,
                target,
                method_id,
                cfg.dynamicbind_posebusters_label,
            )
            if csv_path is None:
                continue
            df = pd.read_csv(csv_path)
            bool_cols = _boolean_check_columns(df)
            if not bool_cols:
                continue
            n_checks_seen.append(len(bool_cols))
            bool_df = pd.DataFrame({c: _to_bool(df[c]) for c in bool_cols})
            pooled[method_id].extend(bool_df.sum(axis=1).astype(int).tolist())

    unique_n = sorted(set(n_checks_seen))
    if not unique_n:
        return 0, pooled
    if len(unique_n) > 1:
        raise ValueError(
            f"PoseBusters CSVs mix n_checks={unique_n}; refuse to pool pass-count curves"
        )
    return unique_n[0], pooled


def pass_count_distribution(counts: list[int], n_tests: int, *, mode: str) -> pd.DataFrame:
    """Build % ligands vs number of tests passed (exactly or at least)."""
    if not counts:
        return pd.DataFrame(columns=["n_tests_passed", "pct_ligands"])
    arr = np.asarray(counts, dtype=int)
    n = len(arr)
    rows: list[dict] = []
    for k in range(1, n_tests + 1):
        if mode == "exactly":
            pct = 100.0 * float((arr == k).sum()) / n
        elif mode == "at_least":
            pct = 100.0 * float((arr >= k).sum()) / n
        else:
            raise ValueError(mode)
        rows.append({"n_tests_passed": k, "pct_ligands": pct})
    return pd.DataFrame(rows)
