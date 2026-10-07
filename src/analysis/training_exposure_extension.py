#!/usr/bin/env python3
"""Part B: method-specific training-exposure evidence (no proxy→membership relabeling)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.constants import METHOD_LABELS, PRIMARY_METRICS
from analysis.data import load_merged
from analysis.robustness_common import (
    MERGED_DIR,
    ROBUST_TABLES,
    TABLES,
    TARGETS,
    RobustnessConfig,
    canonicalize_evidence_level,
    method_label,
    morgan_bitvect,
    target_label,
)

OVERLAP = TABLES / "overlap"
LP_CACHE = OVERLAP / "external" / "LP_PDBBind_ligand_features.csv"
DYN_CSV = OVERLAP / "external" / "dynamicbind" / "d3_with_clash_info.csv"
GNINA_PDBS = OVERLAP / "external" / "gnina" / "crossdocked_v1.3_train_pdbs.txt"
FASTA = OVERLAP / "external" / "uniprot_panel.fasta"

NOT_COMPUTABLE = [
    {
        "quantity": "Boltz-2 exact affinity-training pair membership (InChIKey × UniProt against official filters)",
        "method": "Boltz-2",
        "reason": "No official training/evaluation manifest is bundled with the public checkpoint.",
        "required_files": "Official Boltz-2 affinity-training export (assay filters, snapshot date, pair keys).",
        "evidence_level": "indeterminate",
    },
    {
        "quantity": "PLAPT exact first-100k pair membership",
        "method": "PLAPT",
        "reason": "HuggingFace first-100k composition depends on dataset snapshot and parquet row order; no immutable snapshot is in this repository.",
        "required_files": "Frozen parquet + documented row-order hash from the 2024-01-27 access cited in the preprint.",
        "evidence_level": "indeterminate",
    },
    {
        "quantity": "GNINA 1.3 ligand- and pair-level CrossDocked types membership",
        "method": "GNINA 1.3",
        "reason": "Training types name docked poses, not parent InChIKeys. Molecular structures / molcache are absent.",
        "required_files": "CrossDocked2020 v1.3 molcache or SDF for every types-file pose.",
        "evidence_level": "indeterminate",
    },
    {
        "quantity": "Maximum ECFP4 Tanimoto vs Boltz-2 / PLAPT / GNINA training ligands",
        "method": "Boltz-2, PLAPT, GNINA 1.3",
        "reason": "No training SMILES manifests for these engines are locally available.",
        "required_files": "Training SMILES or InChIKeys for each engine.",
        "evidence_level": "indeterminate",
    },
    {
        "quantity": "Maximum pairwise protein sequence identity vs each engine's training receptors (Needleman–Wunsch / BLAST)",
        "method": "all ML methods",
        "reason": "No complete training-receptor FASTA is available except LP-PDBBind sequences usable as a DynamicBind crystal-ligand companion. A proper identity matrix was not substituted from unaligned prefixes.",
        "required_files": "Training-receptor FASTA per engine, or a documented alignment protocol.",
        "evidence_level": "indeterminate",
    },
    {
        "quantity": "New post-cutoff Ki records outside the current AffiTox panel",
        "method": "all",
        "reason": "No frozen external ChEMBL/BindingDB dump dated after model freeze is bundled.",
        "required_files": "Post-freeze assay export with exact Ki, human single-protein targets, nM units, parent InChIKey, UniProt, document year, assay ID.",
        "evidence_level": "indeterminate",
    },
]


def _read_gnina_pdbs() -> set[str]:
    if not GNINA_PDBS.exists():
        return set()
    return {line.strip().upper() for line in GNINA_PDBS.read_text().splitlines() if line.strip()}


def _dynamicbind_train_pdbs() -> set[str]:
    if not DYN_CSV.exists():
        return set()
    d3 = pd.read_csv(DYN_CSV, usecols=["pdb", "group"])
    return set(d3.loc[d3["group"].astype(str) == "train", "pdb"].astype(str).str.upper())


def evidence_catalog() -> pd.DataFrame:
    frac = pd.read_csv(OVERLAP / "overlap_fractions_by_method.csv")
    rows = []
    for _, r in frac.iterrows():
        src_level = r["evidence_level"]
        canonical = canonicalize_evidence_level(src_level, r["overlap_dimension"])
        exact_or_proxy = {
            "exact_manifest_intersection": "exact",
            "processed_source_table_intersection": "exact",
            "not_applicable": "not_applicable",
            "indeterminate": "indeterminate",
        }.get(canonical, "proxy")
        n_items = r["n_items"]
        n_missing = np.nan if pd.isna(n_items) else 0
        rows.append(
            {
                "method": r["method"],
                "method_id": r["method_id"],
                "quantity": r["overlap_dimension"],
                "n_items": n_items,
                "n_overlap": r["n_overlap"],
                "frac": r["frac_overlap"],
                "status": r["status"],
                "evidence_level_source": src_level,
                "evidence_level": canonical,
                "exact_or_proxy": exact_or_proxy,
                "denominator": n_items,
                "n_missing": n_missing,
                "evidence_source": "analysis/tables/overlap/overlap_fractions_by_method.csv",
                "definition": r["definition"],
            }
        )
    rows.append(
        {
            "method": "Boltz-2",
            "method_id": "boltz2",
            "quantity": "official affinity-training manifest intersection",
            "n_items": np.nan,
            "n_overlap": np.nan,
            "frac": np.nan,
            "status": "not_computable",
            "evidence_level_source": "indeterminate",
            "evidence_level": "indeterminate",
            "exact_or_proxy": "indeterminate",
            "denominator": np.nan,
            "n_missing": np.nan,
            "evidence_source": "none",
            "definition": "Exact pair membership is indeterminate unless an official manifest is present.",
        }
    )
    return pd.DataFrame(rows)


def similarity_dynamicbind(cfg: RobustnessConfig) -> pd.DataFrame:
    """ECFP4 max Tanimoto vs processed-table crystal ligands (similarity proxy, not membership)."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import DataStructs

    RDLogger.DisableLog("rdApp.*")

    if not LP_CACHE.exists() or not DYN_CSV.exists():
        return pd.DataFrame()
    train_pdbs = _dynamicbind_train_pdbs()
    lp = pd.read_csv(LP_CACHE)
    lp["pdb"] = lp["pdb"].astype(str).str.upper()
    train_smiles = lp.loc[lp["pdb"].isin(train_pdbs), "smiles"].dropna().astype(str).unique().tolist()
    fps_train = []
    n_fail_train = 0
    for smi in train_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_fail_train += 1
            continue
        fps_train.append(morgan_bitvect(mol))
    rows = []
    for target in TARGETS:
        merged = load_merged(MERGED_DIR, target)
        max_sims = []
        n_fail = 0
        for smi in merged["canonical_smiles"].astype(str):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                n_fail += 1
                continue
            fp = morgan_bitvect(mol)
            if not fps_train:
                max_sims.append(np.nan)
                continue
            sims = DataStructs.BulkTanimotoSimilarity(fp, fps_train)
            max_sims.append(float(max(sims)))
        arr = np.asarray(max_sims, dtype=float)
        rows.append(
            {
                "method": "DynamicBind",
                "method_id": "dynamicbind",
                "target": target,
                "target_label": target_label(target),
                "quantity": "max ECFP4 Tanimoto vs processed-table crystal ligands",
                "evidence_level": "similarity_based_exposure_proxy",
                "exact_or_proxy": "proxy",
                "n_affitox_ligands": int(len(merged)),
                "n_missing_fingerprints": int(n_fail),
                "n_train_fingerprints": int(len(fps_train)),
                "n_train_smiles_unparsed": int(n_fail_train),
                "median_max_tanimoto": float(np.nanmedian(arr)) if np.isfinite(arr).any() else np.nan,
                "mean_max_tanimoto": float(np.nanmean(arr)) if np.isfinite(arr).any() else np.nan,
                "frac_max_tanimoto_eq_1": float(np.nanmean(arr == 1.0)) if np.isfinite(arr).any() else np.nan,
                "evidence_source": "LP_PDBBind_ligand_features.csv ∩ d3 group==train",
            }
        )
    return pd.DataFrame(rows)


def by_method_table(evidence: pd.DataFrame, sim: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method_id, label in METHOD_LABELS.items():
        sub = evidence.loc[evidence["method_id"] == method_id]
        levels = sorted(set(sub["evidence_level"].dropna().astype(str)))
        rows.append(
            {
                "method": label,
                "method_id": method_id,
                "evidence_levels_present": "|".join(levels),
                "n_quantities": int(len(sub)),
                "has_exact_manifest": int((sub["evidence_level"] == "exact_manifest_intersection").any()),
                "has_processed_table": int(
                    (sub["evidence_level"] == "processed_source_table_intersection").any()
                ),
                "has_only_proxy_or_indeterminate": int(
                    method_id in {"boltz2", "plapt"}
                ),
                "supervised_bioactivity_overlap": (
                    "not_applicable" if method_id == "qvina" else "see evidence table"
                ),
                "similarity_proxy_computed": int(len(sim) > 0 and method_id == "dynamicbind"),
            }
        )
    return pd.DataFrame(rows)


def exposure_performance(cfg: RobustnessConfig) -> pd.DataFrame:
    """Descriptive target-level performance vs PDB-level exposure flags. Not causal inference."""
    corr = pd.read_csv(TABLES / "correlations" / "correlations_with_ci_perm.csv")
    gnina_pdbs = _read_gnina_pdbs()
    dyn_pdbs = _dynamicbind_train_pdbs()
    rng = np.random.default_rng(cfg.seed)
    rows = []

    def _median_ci(vals: np.ndarray) -> tuple[float, float, float]:
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            return np.nan, np.nan, np.nan
        boots = []
        for _ in range(cfg.n_boot):
            boots.append(float(np.median(rng.choice(vals, size=len(vals), replace=True))))
        arr = np.asarray(boots)
        return float(np.median(vals)), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))

    for method_id, flag_fn, level, name in (
        (
            "gnina",
            lambda t: t.upper() in gnina_pdbs,
            "exact_manifest_intersection",
            "GNINA types PDB intersection",
        ),
        (
            "dynamicbind",
            lambda t: t.upper() in dyn_pdbs,
            "processed_source_table_intersection",
            "DynamicBind processed-table PDB intersection",
        ),
    ):
        sub = corr.loc[corr["method_id"] == method_id].copy()
        sub["exposed"] = sub["target"].map(flag_fn)
        for exposed, lab in ((True, "exposed"), (False, "not_exposed")):
            vals = sub.loc[sub["exposed"] == exposed, "pearson_r"].to_numpy(dtype=float)
            med, lo, hi = _median_ci(vals)
            rows.append(
                {
                    "method": method_label(method_id),
                    "method_id": method_id,
                    "exposure_flag": name,
                    "evidence_level": level,
                    "stratum": lab,
                    "n_targets": int(np.isfinite(vals).sum()),
                    "median_pearson_r": med,
                    "median_pearson_r_ci_low": lo,
                    "median_pearson_r_ci_high": hi,
                    "inference": "descriptive_target_level_only",
                    "note": "PDB-level flag; ligand/pair membership may still be absent. Not a leakage test.",
                }
            )
    return pd.DataFrame(rows)


def write_limitations(evidence: pd.DataFrame, sim: pd.DataFrame) -> None:
    lines = [
        "# Training-exposure limitations",
        "",
        "Proxy-based exposure estimates are **not** exact training-set membership.",
        "Processed-source-table and types-file intersections are exact only for the",
        "released record that was actually intersected (DynamicBind processed CSV;",
        "GNINA CrossDocked v1.3 training types PDB IDs).",
        "",
        "## Evidence-level vocabulary used here",
        "",
        "1. exact manifest intersection",
        "2. processed-source-table intersection",
        "3. temporal eligibility proxy",
        "4. source-domain exposure",
        "5. similarity-based exposure proxy",
        "6. indeterminate",
        "7. not applicable",
        "",
        "QVina2 supervised bioactivity-training overlap is **not applicable**, not zero.",
        "",
        "## Computable in this repository",
        "",
        f"- {len(evidence)} evidence rows written to training_exposure_evidence.csv",
        f"- DynamicBind ECFP4 similarity proxy rows: {len(sim)}",
        "",
        "## Not computable",
        "",
        "See NOT_COMPUTABLE.md in this directory. No guessed numeric substitutes were inserted.",
        "",
        "## Performance versus exposure",
        "",
        "Target-level median Pearson r by PDB exposure flag is descriptive.",
        "Within-target ligand-level DynamicBind pair hits are too sparse for a powered",
        "exposed-vs-unseen Pearson comparison (see overlap_fractions_by_method.csv).",
        "",
    ]
    (ROBUST_TABLES / "training_exposure_limitations.md").write_text("\n".join(lines))
    nc = "# NOT_COMPUTABLE\n\n"
    for item in NOT_COMPUTABLE:
        nc += (
            f"## {item['quantity']}\n\n"
            f"- Method: {item['method']}\n"
            f"- Evidence level: {item['evidence_level']}\n"
            f"- Why: {item['reason']}\n"
            f"- Required: {item['required_files']}\n\n"
        )
    (ROBUST_TABLES / "NOT_COMPUTABLE.md").write_text(nc)


def run_training_exposure(cfg: RobustnessConfig) -> dict[str, Path]:
    cfg.ensure_dirs()
    evidence = evidence_catalog()
    sim = similarity_dynamicbind(cfg)
    if len(sim):
        evidence = pd.concat(
            [
                evidence,
                pd.DataFrame(
                    [
                        {
                            "method": "DynamicBind",
                            "method_id": "dynamicbind",
                            "quantity": "max ECFP4 Tanimoto vs processed-table crystal ligands (panel median of per-target medians)",
                            "n_items": float(sim["n_affitox_ligands"].sum()),
                            "n_overlap": np.nan,
                            "frac": float(sim["median_max_tanimoto"].median()),
                            "status": "similarity_proxy",
                            "evidence_level_source": "similarity_based_exposure_proxy",
                            "evidence_level": "similarity_based_exposure_proxy",
                            "exact_or_proxy": "proxy",
                            "denominator": float(sim["n_affitox_ligands"].sum()),
                            "n_missing": int(sim["n_missing_fingerprints"].sum()),
                            "evidence_source": "LP_PDBBind_ligand_features.csv ∩ d3 group==train",
                            "definition": "Maximum ECFP4 Tanimoto to a processed-table crystal ligand. Similarity proxy, not training-set membership.",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
    by_method = by_method_table(evidence, sim)
    perf = exposure_performance(cfg)
    p1 = ROBUST_TABLES / "training_exposure_evidence.csv"
    p2 = ROBUST_TABLES / "training_exposure_by_method.csv"
    p3 = ROBUST_TABLES / "exposure_performance_sensitivity.csv"
    evidence.to_csv(p1, index=False)
    by_method.to_csv(p2, index=False)
    perf.to_csv(p3, index=False)
    if len(sim):
        sim.to_csv(ROBUST_TABLES / "dynamicbind_ecfp4_similarity_proxy.csv", index=False)
    write_limitations(evidence, sim)
    return {"evidence": p1, "by_method": p2, "performance": p3}


if __name__ == "__main__":
    run_training_exposure(RobustnessConfig())
