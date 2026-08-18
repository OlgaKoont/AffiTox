"""Within-target molecular-property null baselines for scoring-power heatmaps.

Computes Pearson r(pKi, property) for heavy-atom count, MW, and cLogP, plus a
permutation-floor |r| under shuffled pairings. Within-target estimands avoid the
cross-target size confound that inflates pooled MW/HAC–affinity correlations.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
from scipy.stats import pearsonr

from .config import AnalysisConfig
from .data import ensure_pki, load_merged

BASELINE_LABELS: dict[str, str] = {
    "heavy_atom_count": "Heavy-atom count",
    "mol_wt": "MW",
    "clogp": "cLogP",
    "permutation_floor": "Permutation floor",
}


def _descriptor_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Return pKi + HAC/MW/cLogP columns aligned to ligand rows."""
    pki = ensure_pki(df)
    smiles = df["canonical_smiles"].astype(str)
    if "clogp" in df.columns:
        clogp = pd.to_numeric(df["clogp"], errors="coerce")
    else:
        clogp = pd.Series(np.nan, index=df.index, dtype=float)

    hac = np.full(len(df), np.nan, dtype=float)
    mw = np.full(len(df), np.nan, dtype=float)
    clogp_rdkit = np.full(len(df), np.nan, dtype=float)
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        hac[i] = float(Descriptors.HeavyAtomCount(mol))
        mw[i] = float(Descriptors.MolWt(mol))
        clogp_rdkit[i] = float(Crippen.MolLogP(mol))

    out = pd.DataFrame(
        {
            "pki": pki.to_numpy(dtype=float),
            "heavy_atom_count": hac,
            "mol_wt": mw,
            "clogp": clogp.to_numpy(dtype=float),
            "clogp_rdkit": clogp_rdkit,
        }
    )
    # Prefer curated clogp; fall back to RDKit Crippen when missing.
    miss = ~np.isfinite(out["clogp"])
    out.loc[miss, "clogp"] = out.loc[miss, "clogp_rdkit"]
    return out.drop(columns=["clogp_rdkit"])


def _pearson_safe(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan")
    xv = x[mask]
    yv = y[mask]
    if np.std(xv) == 0.0 or np.std(yv) == 0.0:
        return float("nan")
    r, _ = pearsonr(xv, yv)
    return float(r)


def permutation_floor_abs_r(
    pki: np.ndarray,
    *,
    n_perm: int,
    seed: int,
) -> float:
    """Mean |Pearson r| under random re-pairings of pKi (noise floor for this N)."""
    x = np.asarray(pki, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 3:
        return float("nan")
    rng = np.random.default_rng(seed)
    abs_rs = np.empty(n_perm, dtype=float)
    for b in range(n_perm):
        yp = rng.permutation(x)
        r, _ = pearsonr(x, yp)
        abs_rs[b] = abs(float(r))
    return float(np.mean(abs_rs))


def run_property_baselines(cfg: AnalysisConfig, *, n_perm: int | None = None) -> pd.DataFrame:
    """Per-target property and permutation baselines; write CSV under correlations/."""
    n_perm = int(n_perm if n_perm is not None else min(max(cfg.n_permutation, 200), 2000))
    rows: list[dict] = []
    for t_i, target in enumerate(cfg.targets):
        df = load_merged(cfg.merged_dir, target)
        if "canonical_smiles" not in df.columns:
            continue
        props = _descriptor_frame(df)
        pki = props["pki"].to_numpy(dtype=float)
        n_ok = int(np.isfinite(pki).sum())
        floor = permutation_floor_abs_r(
            pki, n_perm=n_perm, seed=cfg.random_seed + 17 + t_i
        )
        for key, label in BASELINE_LABELS.items():
            if key == "permutation_floor":
                r = floor
            else:
                r = _pearson_safe(pki, props[key].to_numpy(dtype=float))
            rows.append(
                {
                    "target": target.lower(),
                    "baseline_id": key,
                    "baseline_label": label,
                    "n_points": n_ok,
                    "pearson_r": r,
                    "n_permutation": n_perm if key == "permutation_floor" else np.nan,
                    "note": (
                        "mean_|r|_under_pKi_shuffle"
                        if key == "permutation_floor"
                        else "within_target_r_property_vs_pKi"
                    ),
                }
            )

    summary = pd.DataFrame(rows)
    out = cfg.tables_dir / "correlations" / "property_null_baselines.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out, index=False)
    return summary


def baseline_pearson_matrix(baselines: pd.DataFrame, cfg: AnalysisConfig) -> pd.DataFrame:
    """Rows = baseline labels, columns = target PDB ids (uppercase)."""
    idx = list(BASELINE_LABELS.values())
    cols = [t.upper() for t in cfg.targets]
    mat = pd.DataFrame(index=idx, columns=cols, dtype=float)
    for _, row in baselines.iterrows():
        label = row["baseline_label"]
        tgt = str(row["target"]).upper()
        if label in mat.index and tgt in mat.columns:
            mat.loc[label, tgt] = float(row["pearson_r"])
    return mat


def baseline_summary_medians(baselines: pd.DataFrame) -> pd.Series:
    """Median across targets of within-target baseline r (not pooled)."""
    out = {}
    for key, label in BASELINE_LABELS.items():
        sub = baselines.loc[baselines["baseline_id"] == key, "pearson_r"]
        out[label] = float(sub.median()) if len(sub) else float("nan")
    return pd.Series(out, dtype=float)
