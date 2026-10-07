"""Shared paths, inventory helpers, and resampling utilities for robustness analyses."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, rankdata, spearmanr

from .constants import METHOD_LABELS, POSEBUSTERS_EXCLUDE, PRIMARY_METRICS
from .data import resolve_posebusters_csv
from .inferential import benjamini_hochberg
from .posebusters import _boolean_check_columns, _to_bool
from .training_overlap_audit import LIGANDS, PDB_TO_CSV, TARGET_LABEL

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TABLES = PROJECT_ROOT / "analysis" / "tables"
FIGURES = PROJECT_ROOT / "analysis" / "figures"
ROBUST_TABLES = TABLES / "robustness"
ROBUST_FIGURES = FIGURES / "robustness"
POSE_DIR = TABLES / "posebuster"
MERGED_DIR = TABLES

POSE_METHODS = [m for m in PRIMARY_METRICS if m not in POSEBUSTERS_EXCLUDE]
ALL_METHODS = list(PRIMARY_METRICS.keys())
TARGETS = [p.lower() for p in PDB_TO_CSV]
DEFAULT_SEED = 42
DEFAULT_N_BOOT = 10_000
DEFAULT_N_PERM = 10_000
DEFAULT_N_BALANCED = 2_000
MIN_N_ELIGIBLE = 10
MIN_PKI_RANGE = 1.0

CHEMBL_RE = re.compile(r"(CHEMBL\d+)", re.IGNORECASE)
LIGAND_RE = re.compile(r"ligand_(\d+)", re.IGNORECASE)

EVIDENCE_LEVELS = (
    "exact_manifest_intersection",
    "processed_source_table_intersection",
    "temporal_eligibility_proxy",
    "source_domain_exposure",
    "similarity_based_exposure_proxy",
    "indeterminate",
    "not_applicable",
)


@dataclass(frozen=True)
class RobustnessConfig:
    seed: int = DEFAULT_SEED
    n_boot: int = DEFAULT_N_BOOT
    n_perm: int = DEFAULT_N_PERM
    n_balanced: int = DEFAULT_N_BALANCED
    min_n: int = MIN_N_ELIGIBLE
    min_pki_range: float = MIN_PKI_RANGE

    def ensure_dirs(self) -> None:
        ROBUST_TABLES.mkdir(parents=True, exist_ok=True)
        ROBUST_FIGURES.mkdir(parents=True, exist_ok=True)


def file_md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def inventory_row(path: Path, *, unit: str, key_cols: list[str] | None = None) -> dict:
    try:
        file_path = str(path.resolve().relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        file_path = str(path)
    rec = {
        "file_path": file_path,
        "exists": path.exists(),
        "row_count": np.nan,
        "n_columns": np.nan,
        "columns": "",
        "unique_targets": "",
        "n_unique_targets": np.nan,
        "unique_methods": "",
        "n_unique_methods": np.nan,
        "n_missing_cells": np.nan,
        "unit_of_analysis": unit,
        "duplicate_key_rows": np.nan,
        "key_columns": ",".join(key_cols or []),
        "notes": "",
    }
    if not path.exists():
        rec["notes"] = "missing"
        return rec
    if path.suffix.lower() not in {".csv", ".tsv", ".txt"}:
        rec["notes"] = f"not a table ({path.suffix})"
        rec["row_count"] = np.nan
        return rec
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} and path.name.endswith(".tsv") else ","
    if path.suffix.lower() == ".txt":
        rec["row_count"] = sum(1 for _ in path.open() if _.strip())
        rec["unit_of_analysis"] = unit
        rec["notes"] = "line-oriented text"
        return rec
    try:
        df = pd.read_csv(path, sep=sep if path.suffix.lower() == ".tsv" else ",")
    except Exception as exc:
        rec["notes"] = f"unreadable: {exc}"
        return rec
    rec["row_count"] = int(len(df))
    rec["n_columns"] = int(df.shape[1])
    rec["columns"] = "|".join(map(str, df.columns))
    rec["n_missing_cells"] = int(df.isna().sum().sum())
    target_col = next((c for c in ("target", "pdb", "protein") if c in df.columns), None)
    method_col = next((c for c in ("method_id", "method") if c in df.columns), None)
    if target_col:
        vals = sorted({str(v).lower() for v in df[target_col].dropna().unique()})
        rec["unique_targets"] = "|".join(vals)
        rec["n_unique_targets"] = len(vals)
    if method_col:
        vals = sorted({str(v) for v in df[method_col].dropna().unique()})
        rec["unique_methods"] = "|".join(vals)
        rec["n_unique_methods"] = len(vals)
    if key_cols and all(c in df.columns for c in key_cols):
        rec["duplicate_key_rows"] = int(df.duplicated(key_cols).sum())
    else:
        rec["duplicate_key_rows"] = np.nan
        if key_cols:
            rec["notes"] = "key columns absent"
    return rec


def load_nodubl_index(target: str) -> pd.DataFrame:
    stem = PDB_TO_CSV[target.lower()]
    path = LIGANDS / f"{stem}_nodubl.csv"
    df = pd.read_csv(path, sep=";")
    df = df.copy()
    df["csv_index"] = np.arange(len(df))
    df["ligand_token"] = [f"ligand_{i + 1}" for i in df["csv_index"]]
    if "molecule_chembl_id" in df.columns:
        df["molecule_chembl_id"] = df["molecule_chembl_id"].astype(str)
    return df


def pb_ligand_token(file_col: object) -> str | None:
    if not isinstance(file_col, str):
        return None
    m = CHEMBL_RE.search(file_col)
    if m:
        return m.group(1).upper()
    m = LIGAND_RE.search(file_col)
    if m:
        return f"ligand_{int(m.group(1))}"
    return None


def collapse_posebusters(
    df: pd.DataFrame,
    nodubl: pd.DataFrame,
    merged: pd.DataFrame,
) -> pd.DataFrame:
    """One row per AffiTox merged ligand that has an evaluable PoseBusters pose.

    Duplicate pose rows (e.g. two Boltz-2 SDF copies) are collapsed by taking
    ``position==0`` when present, otherwise the first row. Unmatched PoseBusters
    ligands are dropped from the merged-panel denominators and counted separately.
    """
    work = df.copy()
    work["_token"] = work["file"].map(pb_ligand_token) if "file" in work.columns else None
    bool_cols = _boolean_check_columns(work)
    if not bool_cols:
        return pd.DataFrame()
    for c in bool_cols:
        work[c] = _to_bool(work[c])
    work["_n_checks_defined"] = work[bool_cols].notna().sum(axis=1)
    work["_n_pass"] = work[bool_cols].sum(axis=1, min_count=len(bool_cols))
    work["_pass_all"] = work["_n_pass"] == len(bool_cols)
    work["_frac_pass"] = work[bool_cols].mean(axis=1)
    if "position" in work.columns:
        work["_pos"] = pd.to_numeric(work["position"], errors="coerce").fillna(10**9)
        work = work.sort_values(["_token", "_pos"])
    else:
        work = work.sort_values(["_token"])
    work = work.drop_duplicates("_token", keep="first")

    chembl_from_token = {}
    if "molecule_chembl_id" in nodubl.columns:
        for tok, chembl in zip(nodubl["ligand_token"], nodubl["molecule_chembl_id"]):
            chembl_from_token[tok] = str(chembl).upper() if pd.notna(chembl) else None
    merged = merged.copy()
    merged["_chembl"] = (
        merged["molecule_chembl_id"].astype(str).str.upper()
        if "molecule_chembl_id" in merged.columns
        else ""
    )
    merged["_smiles"] = merged["canonical_smiles"].astype(str)
    chembl_to_smiles = dict(zip(merged["_chembl"], merged["_smiles"]))
    token_to_smiles = {}
    if "canonical_smiles" in nodubl.columns:
        token_to_smiles = dict(zip(nodubl["ligand_token"], nodubl["canonical_smiles"].astype(str)))
    merged_smiles = set(merged["_smiles"])

    rows = []
    for _, row in work.iterrows():
        token = row["_token"]
        if not isinstance(token, str):
            continue
        chembl = token if token.startswith("CHEMBL") else chembl_from_token.get(token)
        smiles = chembl_to_smiles.get(chembl) if chembl else None
        if smiles is None:
            cand = token_to_smiles.get(token)
            if cand in merged_smiles:
                smiles = cand
        rows.append(
            {
                "token": token,
                "molecule_chembl_id": chembl,
                "canonical_smiles": smiles,
                "matched_to_merged": bool(smiles in merged_smiles) if smiles else False,
                "n_checks": int(len(bool_cols)),
                "n_checks_defined": int(row["_n_checks_defined"]),
                "evaluable": bool(row["_n_checks_defined"] == len(bool_cols)),
                "pass_all": bool(row["_pass_all"]) if pd.notna(row["_pass_all"]) else False,
                "frac_checks_passed": float(row["_frac_pass"]) if pd.notna(row["_frac_pass"]) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def two_way_dummies(targets: np.ndarray, methods: np.ndarray) -> np.ndarray:
    targets = np.asarray(targets)
    methods = np.asarray(methods)
    t_levels = pd.unique(targets)
    m_levels = pd.unique(methods)
    t_mat = (
        np.column_stack([(targets == lev).astype(float) for lev in t_levels[1:]])
        if len(t_levels) > 1
        else np.zeros((len(targets), 0))
    )
    m_mat = (
        np.column_stack([(methods == lev).astype(float) for lev in m_levels[1:]])
        if len(m_levels) > 1
        else np.zeros((len(methods), 0))
    )
    intercept = np.ones((len(targets), 1))
    return np.column_stack([intercept, t_mat, m_mat])


def residualize_with_design(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    resid = np.full_like(y, np.nan, dtype=float)
    ok = np.isfinite(y)
    if x.ndim != 2 or len(x) != len(y):
        raise ValueError("design matrix must align with y")
    ok = ok & np.isfinite(x).all(axis=1)
    if ok.sum() < 3:
        return resid
    beta, *_ = np.linalg.lstsq(x[ok], y[ok], rcond=None)
    resid[ok] = y[ok] - x[ok] @ beta
    return resid


def residualize_two_way_means(
    y: np.ndarray,
    t_codes: np.ndarray,
    m_codes: np.ndarray,
    n_t: int,
    n_m: int,
) -> np.ndarray:
    """Additive two-way residual y_ij − ȳ_i. − ȳ_.j + ȳ_.. (balanced or unbalanced)."""
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(y)
    resid = np.full_like(y, np.nan, dtype=float)
    if ok.sum() < 3:
        return resid
    w = ok.astype(float)
    yw = np.where(ok, y, 0.0)
    grand = yw.sum() / w.sum()
    t_sum = np.bincount(t_codes, weights=yw, minlength=n_t)
    t_n = np.bincount(t_codes, weights=w, minlength=n_t)
    m_sum = np.bincount(m_codes, weights=yw, minlength=n_m)
    m_n = np.bincount(m_codes, weights=w, minlength=n_m)
    t_mean = t_sum / np.maximum(t_n, 1.0)
    m_mean = m_sum / np.maximum(m_n, 1.0)
    resid[ok] = y[ok] - t_mean[t_codes[ok]] - m_mean[m_codes[ok]] + grand
    return resid


def residualize_two_way(y: np.ndarray, targets: np.ndarray, methods: np.ndarray) -> np.ndarray:
    t_codes, t_unq = pd.factorize(np.asarray(targets))
    m_codes, m_unq = pd.factorize(np.asarray(methods))
    return residualize_two_way_means(y, t_codes, m_codes, len(t_unq), len(m_unq))


def spearman_safe(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 3:
        return np.nan
    xs, ys = x[mask], y[mask]
    if np.unique(xs).size < 2 or np.unique(ys).size < 2:
        return np.nan
    rx = rankdata(xs)
    ry = rankdata(ys)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if denom == 0:
        return np.nan
    return float((rx * ry).sum() / denom)


def kendall_safe(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    tau, _ = kendalltau(x[mask], y[mask])
    return float(tau) if np.isfinite(tau) else np.nan


def cluster_bootstrap_spearman(
    x: np.ndarray,
    y: np.ndarray,
    clusters: np.ndarray,
    n_boot: int,
    seed: int,
) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    uniq = np.unique(clusters)
    vals: list[float] = []
    index_map = {c: np.flatnonzero(clusters == c) for c in uniq}
    for _ in range(n_boot):
        draw = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([index_map[c] for c in draw])
        rho = spearman_safe(x[idx], y[idx])
        if np.isfinite(rho):
            vals.append(rho)
    if not vals:
        return np.nan, np.nan, np.nan
    arr = np.asarray(vals)
    return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5)), float(np.mean(arr))


def restricted_permutation_p(
    pose: np.ndarray,
    affinity: np.ndarray,
    targets: np.ndarray,
    methods: np.ndarray,
    n_perm: int,
    seed: int,
    *,
    adjusted: bool,
) -> tuple[float, float]:
    """Shuffle pose among methods within each target; residualize inside each permutation."""
    pose = np.asarray(pose, dtype=float)
    affinity = np.asarray(affinity, dtype=float)
    t_codes, t_unq = pd.factorize(np.asarray(targets))
    m_codes, m_unq = pd.factorize(np.asarray(methods))
    n_t, n_m = len(t_unq), len(m_unq)
    groups = [np.flatnonzero(t_codes == i) for i in range(n_t)]
    if adjusted:
        y_aff = residualize_two_way_means(affinity, t_codes, m_codes, n_t, n_m)
        rho_obs = spearman_safe(residualize_two_way_means(pose, t_codes, m_codes, n_t, n_m), y_aff)
    else:
        y_aff = affinity
        rho_obs = spearman_safe(pose, affinity)
    if not np.isfinite(rho_obs):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    count = 0
    n_valid = 0
    pose_perm = np.empty_like(pose)
    for _ in range(n_perm):
        pose_perm[:] = pose
        for idx in groups:
            pose_perm[idx] = rng.permutation(pose[idx])
        xp = residualize_two_way_means(pose_perm, t_codes, m_codes, n_t, n_m) if adjusted else pose_perm
        rho = spearman_safe(xp, y_aff)
        if not np.isfinite(rho):
            continue
        n_valid += 1
        if abs(rho) >= abs(rho_obs):
            count += 1
    p = (count + 1) / (n_valid + 1) if n_valid else np.nan
    return float(rho_obs), float(p)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n")


_MORGAN_GEN = None


def morgan_bitvect(mol, radius: int = 2, n_bits: int = 2048):
    """Deterministic ECFP4-style Morgan fingerprint (RDKit bit vector)."""
    global _MORGAN_GEN
    try:
        from rdkit.Chem import rdFingerprintGenerator

        if _MORGAN_GEN is None or getattr(_MORGAN_GEN, "_fp_size", None) != n_bits:
            _MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
            _MORGAN_GEN._fp_size = n_bits
        return _MORGAN_GEN.GetFingerprint(mol)
    except Exception:
        from rdkit.Chem import AllChem

        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def canonicalize_evidence_level(raw: object, quantity: object = "") -> str:
    """Map source-table labels onto the seven prespecified evidence levels."""
    text = "" if raw is None or (isinstance(raw, float) and np.isnan(raw)) else str(raw).strip().lower()
    qty = "" if quantity is None or (isinstance(quantity, float) and np.isnan(quantity)) else str(quantity).lower()
    if text in {"not_applicable", "n/a", "na"}:
        return "not_applicable"
    if text in {"indeterminate", "unknown"}:
        return "indeterminate"
    if text in {"exact_types", "exact_manifest_intersection", "exact_manifest"}:
        return "exact_manifest_intersection"
    if text in {"processed_source_table", "processed_source_table_intersection"}:
        return "processed_source_table_intersection"
    if text in {"similarity_based_exposure_proxy", "similarity"}:
        return "similarity_based_exposure_proxy"
    if text in {"source_domain_exposure"}:
        return "source_domain_exposure"
    if text in {"temporal_eligibility_proxy"}:
        return "temporal_eligibility_proxy"
    if text in {"proxy"}:
        if "source-domain" in qty or "source domain" in qty:
            return "source_domain_exposure"
        if "year" in qty or "deposition" in qty or "temporal" in qty or "pubyear" in qty:
            return "temporal_eligibility_proxy"
        return "temporal_eligibility_proxy"
    if text:
        return text
    return "indeterminate"


def method_label(method_id: str) -> str:
    return METHOD_LABELS.get(method_id, method_id)


def target_label(pdb: str) -> str:
    return TARGET_LABEL.get(pdb.upper(), pdb.upper())
