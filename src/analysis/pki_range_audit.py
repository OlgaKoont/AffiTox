"""Within-target pKi-range audit for scoring-power interpretation (reviewer W5).

Reports per-target pKi SD/IQR, effective N, Bemis–Murcko scaffold count, and
joins method Pearson r so range attenuation can be inspected against the
target-dependence claim.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import pearsonr, spearmanr

from .config import AnalysisConfig
from .data import ensure_pki, load_merged


def _murcko_count(smiles: pd.Series) -> int:
    scaffolds: set[str] = set()
    for smi in smiles.dropna().astype(str):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        except Exception:
            continue
        if scaf:
            scaffolds.add(scaf)
    return len(scaffolds)


def _pki_stats(pki: np.ndarray) -> dict[str, float]:
    x = np.asarray(pki, dtype=float)
    x = x[np.isfinite(x)]
    n = int(len(x))
    if n == 0:
        return {
            "n_pki": 0,
            "pki_min": np.nan,
            "pki_max": np.nan,
            "pki_range": np.nan,
            "pki_sd": np.nan,
            "pki_iqr": np.nan,
            "pki_mean": np.nan,
            "pki_median": np.nan,
        }
    q25, q75 = np.percentile(x, [25, 75])
    return {
        "n_pki": n,
        "pki_min": float(np.min(x)),
        "pki_max": float(np.max(x)),
        "pki_range": float(np.max(x) - np.min(x)),
        "pki_sd": float(np.std(x, ddof=1)) if n > 1 else 0.0,
        "pki_iqr": float(q75 - q25),
        "pki_mean": float(np.mean(x)),
        "pki_median": float(np.median(x)),
    }


def run_pki_range_audit(cfg: AnalysisConfig, corr: pd.DataFrame | None = None) -> pd.DataFrame:
    """Write per-target pKi spread + scaffold counts; optionally join method Pearson r."""
    rows: list[dict] = []
    for target in cfg.targets:
        df = load_merged(cfg.merged_dir, target)
        pki = ensure_pki(df, cfg.exp_col).to_numpy(dtype=float)
        stats = _pki_stats(pki)
        n_scaffolds = (
            _murcko_count(df["canonical_smiles"])
            if "canonical_smiles" in df.columns
            else np.nan
        )
        rows.append(
            {
                "target": target.lower(),
                "n_scaffolds_murcko": n_scaffolds,
                **stats,
            }
        )

    panel = pd.DataFrame(rows)

    if corr is None:
        corr_path = cfg.tables_dir / "correlations" / "summary_all_proteins.csv"
        if corr_path.exists():
            corr = pd.read_csv(corr_path)

    if corr is not None and not corr.empty:
        # Wide join: one row per target with method r and n_points
        for method_id in cfg.methods:
            sub = corr.loc[corr["method_id"] == method_id, ["target", "n_points", "pearson_r"]].copy()
            sub["target"] = sub["target"].str.lower()
            label = cfg.method_label(method_id)
            sub = sub.rename(
                columns={
                    "n_points": f"n_points_{method_id}",
                    "pearson_r": f"pearson_r_{method_id}",
                }
            )
            # also store under display label keys for plotting convenience
            sub[f"pearson_r_{label}"] = sub[f"pearson_r_{method_id}"]
            panel = panel.merge(sub, on="target", how="left")

    out = cfg.tables_dir / "correlations" / "pki_range_by_target.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(out, index=False)

    # Compact association of Boltz-2 r with spread (null that target-dependence == range)
    assoc_rows = []
    if "pearson_r_boltz2" in panel.columns:
        for spread_col in ("pki_sd", "pki_iqr", "pki_range", "n_scaffolds_murcko", "n_pki"):
            x = panel[spread_col].to_numpy(dtype=float)
            y = panel["pearson_r_boltz2"].to_numpy(dtype=float)
            mask = np.isfinite(x) & np.isfinite(y)
            if mask.sum() < 4:
                continue
            pr, pp = pearsonr(x[mask], y[mask])
            sr, sp = spearmanr(x[mask], y[mask])
            assoc_rows.append(
                {
                    "method_id": "boltz2",
                    "spread_metric": spread_col,
                    "n_targets": int(mask.sum()),
                    "pearson_r_vs_spread": float(pr),
                    "pearson_p": float(pp),
                    "spearman_rho_vs_spread": float(sr),
                    "spearman_p": float(sp),
                }
            )
    if assoc_rows:
        assoc = pd.DataFrame(assoc_rows)
        assoc.to_csv(
            cfg.tables_dir / "correlations" / "pki_range_vs_pearson_association.csv",
            index=False,
        )

    return panel
