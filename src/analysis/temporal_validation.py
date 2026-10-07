#!/usr/bin/env python3
"""Part C: temporal sensitivity and external-validation readiness (not prospective validation)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.constants import PRIMARY_METRICS, SIGN_FLIP_METRICS
from analysis.data import ensure_pki, load_merged
from analysis.leakage_stratification import (
    STRATUM_GE,
    STRATUM_LE,
    STRATUM_UNK,
    parent_mol,
)
from analysis.robustness_common import (
    ALL_METHODS,
    MERGED_DIR,
    ROBUST_TABLES,
    TABLES,
    TARGETS,
    RobustnessConfig,
    method_label,
    morgan_bitvect,
    target_label,
)
from analysis.training_overlap_audit import PDB_DEPOSITION_YEAR, PDB_TO_UNIPROT
from rdkit import Chem, RDLogger
from rdkit.Chem import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

OVERLAP = TABLES / "overlap"

# Method-specific documented cutoffs. Boltz-2 dates are NOT copied onto other engines.
CUTOFFS = [
    {
        "method_id": "boltz2",
        "documented_training_cutoff": "affinity-head snapshot reconstructed independently as BindingDB/PubChem ~June 2023; official pair manifest unpublished",
        "receptor_structure_cutoff": "PDB <= 2023-06-01 (structure module; not automatically the affinity cutoff)",
        "bioactivity_data_cutoff": "unknown (no official affinity-training dump)",
        "checkpoint_release_date": "public Boltz-2 affinity release (see SI software registry)",
        "unknown_cutoff": False,
        "temporal_proxy_year_used_in_affitox": 2022,
        "temporal_proxy_note": "publication year <=2022 is a conservative temporal proxy, not membership",
    },
    {
        "method_id": "dynamicbind",
        "documented_training_cutoff": "PDBbind v2020 / timesplit depositions before 2019",
        "receptor_structure_cutoff": "deposition year < 2019 for the timesplit training list",
        "bioactivity_data_cutoff": "not a bioactivity-supervised affinity head; crystal complexes",
        "checkpoint_release_date": "2024 public DynamicBind release",
        "unknown_cutoff": False,
        "temporal_proxy_year_used_in_affitox": 2018,
        "temporal_proxy_note": "deposition-year proxy is not exact PDB occurrence in the processed table",
    },
    {
        "method_id": "gnina",
        "documented_training_cutoff": "CrossDocked2020 v1.3 default ensemble (full v1.3 training types)",
        "receptor_structure_cutoff": "PDB-derived CrossDocked2020 complexes (compilation ~2020)",
        "bioactivity_data_cutoff": "not applicable (pose/score CNN on docked complexes)",
        "checkpoint_release_date": "GNINA 1.3",
        "unknown_cutoff": False,
        "temporal_proxy_year_used_in_affitox": 2019,
        "temporal_proxy_note": "Boltz-2 2023 cutoff is not applied to GNINA",
    },
    {
        "method_id": "plapt",
        "documented_training_cutoff": "HuggingFace binding_affinity first 100k rows (snapshot-dependent)",
        "receptor_structure_cutoff": "not applicable (sequence + SMILES; no 3D receptor cutoff)",
        "bioactivity_data_cutoff": "unknown snapshot / row order",
        "checkpoint_release_date": "preprint checkpoint",
        "unknown_cutoff": True,
        "temporal_proxy_year_used_in_affitox": np.nan,
        "temporal_proxy_note": "Boltz-2 cutoff is not applied to PLAPT; exact membership indeterminate",
    },
    {
        "method_id": "qvina",
        "documented_training_cutoff": "not applicable",
        "receptor_structure_cutoff": "not applicable",
        "bioactivity_data_cutoff": "not applicable",
        "checkpoint_release_date": "qvina02 / Vina 1.2 lineage",
        "unknown_cutoff": False,
        "temporal_proxy_year_used_in_affitox": np.nan,
        "temporal_proxy_note": "no supervised bioactivity training set",
    },
]


def _pairs() -> pd.DataFrame:
    return pd.read_csv(OVERLAP / "affinity_pairs_uniprot_inchikey.csv")


def _score_xy(df: pd.DataFrame, metric: str) -> tuple[np.ndarray, np.ndarray]:
    pki = ensure_pki(df)
    score = pd.to_numeric(df.get(metric), errors="coerce")
    mask = pki.notna() & score.notna() & np.isfinite(pki) & np.isfinite(score)
    x = pki[mask].to_numpy(dtype=float)
    y = score[mask].to_numpy(dtype=float)
    if metric in SIGN_FLIP_METRICS:
        y = -y
    return x, y


def _pearson_safe(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x, ddof=1) == 0 or np.std(y, ddof=1) == 0:
        return np.nan
    r, _ = pearsonr(x, y)
    return float(r) if np.isfinite(r) else np.nan


def _scaffold(smiles: str) -> str | None:
    mol = parent_mol(smiles)
    if mol is None:
        return None
    try:
        sc = MurckoScaffold.GetScaffoldForMol(mol)
        if sc is None or sc.GetNumHeavyAtoms() == 0:
            return None
        return Chem.MolToSmiles(sc)
    except Exception:
        return None


def eligibility_audit(cfg: RobustnessConfig) -> pd.DataFrame:
    flags = pd.read_csv(OVERLAP / "pair_membership_flags.csv")
    rows = []
    for method_id in ALL_METHODS:
        for target in TARGETS:
            merged = load_merged(MERGED_DIR, target)
            metric = PRIMARY_METRICS[method_id]
            merged = merged.copy()
            if "canonical_smiles" not in flags.columns:
                raise AssertionError("pair_membership_flags.csv missing canonical_smiles")
            subf = flags.loc[flags["pdb"].astype(str).str.upper() == target.upper()]
            year_map = {}
            stratum_col = "temporal_stratum" if "temporal_stratum" in subf.columns else None
            if stratum_col:
                year_map = dict(zip(subf["canonical_smiles"], subf[stratum_col]))
            merged["_stratum"] = merged["canonical_smiles"].map(year_map).fillna(STRATUM_UNK)
            for stratum, lab in (
                (STRATUM_LE, "pre_cutoff_proxy_le2022"),
                (STRATUM_GE, "post_cutoff_proxy_ge2023"),
            ):
                part = merged.loc[merged["_stratum"] == stratum]
                if method_id != "boltz2":
                    # Still report the AffiTox 2022/2023 split as a shared ligand-era table,
                    # but label it as not that method's cutoff.
                    lab_use = stratum
                else:
                    lab_use = lab
                x, y = _score_xy(part, metric)
                pki_part = ensure_pki(part).to_numpy(dtype=float)
                pki_part = pki_part[np.isfinite(pki_part)]
                pki_range = float(np.nanmax(pki_part) - np.nanmin(pki_part)) if len(pki_part) else np.nan
                n = int(len(x))
                eligible = bool(n >= cfg.min_n and np.isfinite(pki_range) and pki_range >= cfg.min_pki_range)
                n_act = int((pd.to_numeric(part["standard_value"], errors="coerce") < 1000).sum())
                n_inact = int((pd.to_numeric(part["standard_value"], errors="coerce") >= 1000).sum())
                scafs = {_scaffold(s) for s in part["canonical_smiles"].astype(str)}
                scafs.discard(None)
                rows.append(
                    {
                        "method_id": method_id,
                        "method": method_label(method_id),
                        "target": target,
                        "target_label": target_label(target),
                        "stratum": stratum,
                        "stratum_label": lab_use,
                        "cutoff_applied": (
                            "AffiTox publication-year split 2022/2023 (shared ligand-era table; "
                            "not this method's training cutoff unless method_id=boltz2 proxy)"
                        ),
                        "n_ligands_stratum": int(len(part)),
                        "n_scored": n,
                        "pki_range": pki_range,
                        "pki_sd": float(np.nanstd(pki_part, ddof=1)) if len(pki_part) > 1 else np.nan,
                        "n_active": n_act,
                        "n_inactive": n_inact,
                        "n_scaffolds": int(len(scafs)),
                        "eligible_min_n_and_pki_range": eligible,
                        "min_n_rule": cfg.min_n,
                        "min_pki_range_rule": cfg.min_pki_range,
                        "pearson_r": _pearson_safe(x, y),
                        "grain": "target_method_pubyear_stratum",
                        "cutoff_policy": (
                            "AffiTox publication-year split is a ligand-era table. "
                            "It is a Boltz-2 temporal proxy only; it is not the training cutoff "
                            "for PLAPT, DynamicBind, or GNINA 1.3."
                        ),
                    }
                )
    return pd.DataFrame(rows)


def max_tanimoto_to_pre(post_smiles: list[str], pre_smiles: list[str]) -> list[float]:
    fps_pre = []
    for smi in pre_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fps_pre.append(morgan_bitvect(mol))
    out = []
    for smi in post_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None or not fps_pre:
            out.append(np.nan)
            continue
        fp = morgan_bitvect(mol)
        out.append(float(max(DataStructs.BulkTanimotoSimilarity(fp, fps_pre))))
    return out


def within_target_delta_r(cfg: RobustnessConfig, audit: pd.DataFrame) -> pd.DataFrame:
    flags = pd.read_csv(OVERLAP / "pair_membership_flags.csv")
    rows = []
    rng = np.random.default_rng(cfg.seed)
    dual_targets = []
    for target in TARGETS:
        pre_ok = audit.loc[
            (audit.target == target)
            & (audit.stratum == STRATUM_LE)
            & (audit.eligible_min_n_and_pki_range)
        ]
        post_ok = audit.loc[
            (audit.target == target)
            & (audit.stratum == STRATUM_GE)
            & (audit.eligible_min_n_and_pki_range)
        ]
        if pre_ok.empty or post_ok.empty:
            continue
        methods_both = set(pre_ok.method_id) & set(post_ok.method_id)
        if methods_both:
            dual_targets.append(target)
        merged = load_merged(MERGED_DIR, target)
        subf = flags.loc[flags["pdb"].astype(str).str.upper() == target.upper()]
        smap = dict(zip(subf["canonical_smiles"], subf["temporal_stratum"]))
        merged = merged.copy()
        merged["_stratum"] = merged["canonical_smiles"].map(smap).fillna(STRATUM_UNK)
        pre = merged.loc[merged["_stratum"] == STRATUM_LE]
        post = merged.loc[merged["_stratum"] == STRATUM_GE]
        pre_sc = {_scaffold(s) for s in pre["canonical_smiles"].astype(str)}
        pre_sc.discard(None)
        sims = max_tanimoto_to_pre(
            post["canonical_smiles"].astype(str).tolist(),
            pre["canonical_smiles"].astype(str).tolist(),
        )
        max_sim_median = float(np.nanmedian(np.asarray(sims, dtype=float))) if sims else np.nan
        for method_id in ALL_METHODS:
            metric = PRIMARY_METRICS[method_id]
            x_pre, y_pre = _score_xy(pre, metric)
            x_post, y_post = _score_xy(post, metric)
            elig_pre = bool(
                len(x_pre) >= cfg.min_n
                and (np.nanmax(x_pre) - np.nanmin(x_pre) >= cfg.min_pki_range if len(x_pre) else False)
            )
            elig_post = bool(
                len(x_post) >= cfg.min_n
                and (np.nanmax(x_post) - np.nanmin(x_post) >= cfg.min_pki_range if len(x_post) else False)
            )
            r_pre = _pearson_safe(x_pre, y_pre)
            r_post = _pearson_safe(x_post, y_post)
            delta = r_post - r_pre if np.isfinite(r_pre) and np.isfinite(r_post) else np.nan
            lo = hi = np.nan
            if elig_pre and elig_post:
                deltas = []
                for _ in range(cfg.n_boot):
                    i = rng.integers(0, len(x_pre), len(x_pre))
                    j = rng.integers(0, len(x_post), len(x_post))
                    r1 = _pearson_safe(x_pre[i], y_pre[i])
                    r2 = _pearson_safe(x_post[j], y_post[j])
                    if np.isfinite(r1) and np.isfinite(r2):
                        deltas.append(r2 - r1)
                if deltas:
                    arr = np.asarray(deltas)
                    lo, hi = float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))
            n_novel_sc = 0
            for smi in post["canonical_smiles"].astype(str):
                sc = _scaffold(smi)
                if sc and sc not in pre_sc:
                    n_novel_sc += 1
            rows.append(
                {
                    "target": target,
                    "target_label": target_label(target),
                    "method_id": method_id,
                    "method": method_label(method_id),
                    "n_pre": int(len(x_pre)),
                    "n_post": int(len(x_post)),
                    "eligible_pre": elig_pre,
                    "eligible_post": elig_post,
                    "both_eligible": elig_pre and elig_post,
                    "r_pre": r_pre,
                    "r_post": r_post,
                    "delta_r": delta,
                    "delta_r_ci_low": lo,
                    "delta_r_ci_high": hi,
                    "n_boot": cfg.n_boot,
                    "median_max_tanimoto_post_to_pre": max_sim_median,
                    "n_post_novel_scaffolds": n_novel_sc,
                    "analysis_label": "temporal_sensitivity",
                    "not_prospective_validation": True,
                }
            )
    return pd.DataFrame(rows)


def method_postcutoff_counts() -> pd.DataFrame:
    """Method-specific post-cutoff pair counts. Does not copy the Boltz-2 cutoff onto other engines."""
    flags = pd.read_csv(OVERLAP / "pair_membership_flags.csv")
    gnina_pdbs = set()
    gnina_path = OVERLAP / "external" / "gnina" / "crossdocked_v1.3_train_pdbs.txt"
    if gnina_path.exists():
        gnina_pdbs = {line.strip().upper() for line in gnina_path.read_text().splitlines() if line.strip()}
    rows = []

    def _append(method_id, quantity, n_items, n_post, note, computable=True, evidence="proxy"):
        rows.append(
            {
                "method_id": method_id,
                "method": method_label(method_id),
                "quantity": quantity,
                "n_ligand_target_pairs": n_items,
                "n_post_cutoff": n_post,
                "n_eligible_targets": np.nan,
                "computable": computable,
                "evidence": evidence,
                "note": note,
                "grain": "method_panel",
            }
        )

    n_pairs = int(len(flags))
    n_ge = int((flags["temporal_stratum"] == STRATUM_GE).sum())
    n_unk = int((flags["temporal_stratum"] == STRATUM_UNK).sum())
    _append(
        "boltz2",
        "publication year >= 2023 ligand–target pairs (temporal eligibility proxy)",
        n_pairs,
        n_ge,
        f"Unknown-year pairs excluded from the post-cutoff count (n_unknown={n_unk}). Not exact Boltz-2 membership.",
        True,
        "temporal_eligibility_proxy",
    )
    flags["deposition_year"] = pd.to_numeric(flags["deposition_year"], errors="coerce")
    n_dyn_post = int((flags["deposition_year"] >= 2019).sum())
    _append(
        "dynamicbind",
        "receptor deposition year >= 2019 ligand–target pairs",
        n_pairs,
        n_dyn_post,
        "Timesplit-style receptor-era proxy. Not exact processed-table occurrence.",
        True,
        "temporal_eligibility_proxy",
    )
    if "dynamicbind_processed_pair" in flags.columns:
        n_not_pair = int((flags["dynamicbind_processed_pair"] == False).sum())  # noqa: E712
        _append(
            "dynamicbind",
            "pairs absent from processed-table (UniProt, parent InChIKey)",
            n_pairs,
            n_not_pair,
            "Complement of processed-source-table pair intersection; not a calendar cutoff.",
            True,
            "processed_source_table_intersection",
        )
    n_gnina_notypes = int((~flags["pdb"].astype(str).str.upper().isin(gnina_pdbs)).sum()) if gnina_pdbs else np.nan
    _append(
        "gnina",
        "ligand–target pairs whose AffiTox receptor PDB is absent from CrossDocked v1.3 training types",
        n_pairs,
        n_gnina_notypes,
        "Receptor PDB non-membership in types files. Ligand/pair types membership remains indeterminate.",
        bool(gnina_pdbs),
        "exact_manifest_intersection",
    )
    _append(
        "plapt",
        "post-cutoff ligand–target pairs vs PLAPT training snapshot",
        n_pairs,
        np.nan,
        "PLAPT cutoff is unknown; Boltz-2 year split is not applied.",
        False,
        "indeterminate",
    )
    _append(
        "qvina",
        "supervised bioactivity post-cutoff pairs",
        n_pairs,
        np.nan,
        "Not applicable: QVina2 has no supervised bioactivity training set.",
        False,
        "not_applicable",
    )
    return pd.DataFrame(rows)


def candidate_manifest() -> pd.DataFrame:
    flags = pd.read_csv(OVERLAP / "pair_membership_flags.csv")
    if "temporal_stratum" not in flags.columns:
        raise AssertionError("Need temporal_stratum in pair_membership_flags.csv")
    post = flags.loc[flags["temporal_stratum"] == STRATUM_GE].copy()
    merged_parts = []
    for target in TARGETS:
        m = load_merged(MERGED_DIR, target).copy()
        m["pdb"] = target.upper()
        cols = [
            c
            for c in (
                "pdb",
                "canonical_smiles",
                "molecule_chembl_id",
                "assay_chembl_id",
                "standard_value",
                "type",
                "units",
                "uo_units",
                "document_chembl_id",
            )
            if c in m.columns
        ]
        merged_parts.append(m[cols])
    merged_all = pd.concat(merged_parts, ignore_index=True)
    if merged_all.duplicated(["pdb", "canonical_smiles"]).any():
        raise AssertionError("Duplicate pdb×smiles in concatenated merged tables")
    out = post.merge(merged_all, on=["pdb", "canonical_smiles"], how="left", suffixes=("", "_merged"))
    out["publication_year"] = pd.to_numeric(out.get("min_document_year"), errors="coerce")
    out["assay_identifier"] = out["assay_chembl_id"] if "assay_chembl_id" in out.columns else np.nan
    out["parent_inchikey"] = out["parent_inchikey"] if "parent_inchikey" in out.columns else np.nan
    out["murcko_scaffold"] = out["scaffold"] if "scaffold" in out.columns else np.nan
    out["uniprot"] = out["uniprot"] if "uniprot" in out.columns else out["pdb"].map(PDB_TO_UNIPROT)
    ki_type = out["type"].astype(str).str.upper().eq("KI") if "type" in out.columns else True
    out["exact_ki"] = ki_type
    out["single_human_protein_target"] = True
    out["accepted_relation"] = True
    out["standardized_nm_units"] = (
        out["uo_units"].astype(str).eq("UO_0000065") if "uo_units" in out.columns else True
    )
    out["first_publication_strictly_after_2022"] = True
    out["relevant_model_cutoff"] = (
        "AffiTox publication year >= 2023 (Boltz-2 temporal proxy only; "
        "not the PLAPT, DynamicBind, or GNINA training cutoff)"
    )
    out["provenance"] = "AffiTox BindingDB/ChEMBL curated panel (already in-sample)"
    out["novelty_relative_to_current_affitox"] = "post-2022 within AffiTox; not an external holdout"
    out["n_missing_assay_id"] = int(out["assay_identifier"].isna().sum()) if "assay_identifier" in out.columns else np.nan
    out["note"] = (
        "These records are already inside AffiTox. A stronger external temporal test "
        "requires new Ki labels collected after model freeze, not a re-split of the same panel."
    )
    return out


def feasibility(delta: pd.DataFrame, cfg: RobustnessConfig) -> pd.DataFrame:
    """Power/feasibility for a cross-target median Δr under the observed effect-size spread."""
    obs = delta.loc[delta["both_eligible"] & delta["method_id"].eq("boltz2"), "delta_r"].to_numpy(dtype=float)
    obs = obs[np.isfinite(obs)]
    rng = np.random.default_rng(cfg.seed)
    rows = []
    if len(obs) == 0:
        return pd.DataFrame(
            [{"note": "no dual-eligible Boltz-2 targets; feasibility not estimated from data"}]
        )
    sigma = float(np.std(obs, ddof=1)) if len(obs) > 1 else float(np.abs(obs).mean())
    mu = float(np.mean(obs))
    for n_targets in (4, 8, 16, 24, 32):
        for n_lig in (20, 50, 100, 200):
            se_r = 1.0 / np.sqrt(max(n_lig - 3, 1))
            se_delta = float(np.sqrt(2.0) * se_r)
            widths = []
            for _ in range(400):
                true = rng.normal(mu, sigma if sigma > 0 else 0.2, size=n_targets)
                noisy = true + rng.normal(0.0, se_delta, size=n_targets)
                meds = [np.median(rng.choice(noisy, size=n_targets, replace=True)) for __ in range(200)]
                widths.append(float(np.percentile(meds, 97.5) - np.percentile(meds, 2.5)))
            rows.append(
                {
                    "n_targets": n_targets,
                    "n_ligands_per_stratum": n_lig,
                    "observed_boltz2_dual_n": int(len(obs)),
                    "observed_mean_delta_r": mu,
                    "observed_sd_delta_r": sigma,
                    "approx_se_delta_r_from_n_lig": se_delta,
                    "median_ci_width_of_cross_target_median_delta_r": float(np.median(widths)),
                    "label": "feasibility_simulation_not_sample_size_rule",
                }
            )
    return pd.DataFrame(rows)


def write_limitations() -> None:
    text = """# Temporal / external-validation limitations

A retrospective post-cutoff analysis cannot become a true prospective
validation merely through statistical reanalysis. True prospective validation
requires frozen models and analysis rules before new experimental labels are
obtained.

Results in this directory are labelled **temporal sensitivity** (and, where
the post-cutoff slice is treated as an external-like check,
**external temporal validation**). They are **not** prospective validation.

## Cutoff policy

The Boltz-2 structure cutoff (PDB ≤ 2023-06-01) and the AffiTox publication-year
split (≤2022 vs ≥2023) are **not** applied as training cutoffs for PLAPT,
DynamicBind, or GNINA 1.3. Method-specific documented cutoffs are in
`temporal_cutoffs.csv`.

## Eligibility

A target–method–stratum cell is inferentially eligible only if
`n_scored ≥ 10` and `pKi range ≥ 1.0`. A 16-target pre-cutoff median is not
compared with a four-target post-cutoff median.

## Candidate manifest

`external_validation_candidate_manifest.csv` lists AffiTox records whose
document year is ≥2023. They are already inside the present panel. A stronger
test needs new Ki measurements after a model freeze.
"""
    (ROBUST_TABLES / "temporal_validation_limitations.md").write_text(text)


def run_temporal(cfg: RobustnessConfig) -> dict[str, Path]:
    cfg.ensure_dirs()
    cut = pd.DataFrame(CUTOFFS)
    cut["method"] = cut["method_id"].map(method_label)
    audit = eligibility_audit(cfg)
    method_counts = method_postcutoff_counts()
    audit = pd.concat([audit, method_counts], ignore_index=True, sort=False)
    delta = within_target_delta_r(cfg, audit)
    manifest = candidate_manifest()
    feas = feasibility(delta, cfg)
    paths = {
        "cutoffs": ROBUST_TABLES / "temporal_cutoffs.csv",
        "audit": ROBUST_TABLES / "temporal_eligibility_audit.csv",
        "delta": ROBUST_TABLES / "temporal_within_target_delta_r.csv",
        "manifest": ROBUST_TABLES / "external_validation_candidate_manifest.csv",
        "feasibility": ROBUST_TABLES / "external_validation_feasibility.csv",
    }
    cut.to_csv(paths["cutoffs"], index=False)
    audit.to_csv(paths["audit"], index=False)
    delta.to_csv(paths["delta"], index=False)
    manifest.to_csv(paths["manifest"], index=False)
    feas.to_csv(paths["feasibility"], index=False)
    write_limitations()
    return paths


if __name__ == "__main__":
    run_temporal(RobustnessConfig())
