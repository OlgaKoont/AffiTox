#!/usr/bin/env python3
"""W3: method-specific training-exposure sensitivity analysis.

Evidence is not uniform across methods. We report the most specific
publicly verifiable record for each engine:

  processed_source_table  DynamicBind Zenodo clash-annotated table ∩ timesplit
  exact_types             PDB IDs in CrossDocked2020 v1.3 training types (GNINA)
  proxy                   source-domain / temporal eligibility, not membership
  indeterminate           no immutable pair-level manifest (Boltz-2, PLAPT;
                          GNINA ligand/pair without molecule cache)
  not_applicable          no supervised bioactivity training set (QVina2)

These quantities are not equivalent measures of training-set membership.
Exact intersections show that a benchmark entity occurs in a released
training record. Proxy-based estimates show only potential availability.

Affinity-pair identity is (UniProt, parent InChIKey), aggregated *before*
temporal labelling. Publication year is a conservative temporal proxy:
document_year != database ingestion date != training membership.

Strata:
  pubyear_le_2022   min ChEMBL document_year <= 2022
  pubyear_ge_2023   min ChEMBL document_year >= 2023
  unknown_year      missing document_year (excluded from both strata)
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tarfile
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import pearsonr

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.constants import METHOD_LABELS, PRIMARY_METRICS, SIGN_FLIP_METRICS  # noqa: E402
from analysis.data import ensure_pki, load_merged  # noqa: E402
from analysis.inferential import bootstrap_ci, holm, wilcoxon_paired  # noqa: E402
from analysis.training_overlap_audit import (  # noqa: E402
    LIGANDS,
    PDB_DEPOSITION_YEAR,
    PDB_TO_CSV,
    PDB_TO_UNIPROT,
    TARGET_LABEL,
)

TABLES = PROJECT_ROOT / "analysis" / "tables"
OUT = TABLES / "overlap"
EXTERNAL = OUT / "external"
LP_PATH = EXTERNAL / "LP_PDBBind.csv"
LP_CACHE = EXTERNAL / "LP_PDBBind_ligand_features.csv"
SIFTS_PATH = EXTERNAL / "pdb_chain_uniprot.tsv.gz"
DYN_CSV = EXTERNAL / "dynamicbind" / "d3_with_clash_info.csv"
SPLIT_TRAIN = EXTERNAL / "splits" / "timesplit_no_lig_overlap_train"
SPLIT_VAL = EXTERNAL / "splits" / "timesplit_no_lig_overlap_val"
GNINA_TYPES_TGZ = EXTERNAL / "gnina" / "CrossDocked2020_v1.3_types.tgz"
GNINA_PDB_CACHE = EXTERNAL / "gnina" / "crossdocked_v1.3_train_pdbs.txt"
PDB_ID_RE = re.compile(r"/([0-9][A-Za-z0-9]{3})_", re.IGNORECASE)

PUBYEAR_LE_MAX = 2022
PUBYEAR_GE_MIN = 2023
MIN_N = 10
N_BOOT = 2000
BOOT_SEED = 42
METHODS = list(PRIMARY_METRICS.keys())
BOLTZ2_COMPARATORS = [m for m in METHODS if m != "boltz2"]
STRATUM_LE = "pubyear_le_2022"
STRATUM_GE = "pubyear_ge_2023"
STRATUM_UNK = "unknown_year"


def _largest_fragment(mol: Chem.Mol) -> Chem.Mol | None:
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def parent_mol(smiles: str) -> Chem.Mol | None:
    """Salt/fragment-stripped parent molecule; stereo retained; no tautomer enum."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = _largest_fragment(mol)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    except Exception:
        pass
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def parent_inchikey(smiles: str) -> str | None:
    mol = parent_mol(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol) or None
    except Exception:
        return None


def murcko_smiles(smiles: str) -> tuple[str | None, bool]:
    """Return (Murcko SMILES, empty_flag). Empty scaffolds are not overlap keys."""
    mol = parent_mol(smiles)
    if mol is None:
        return None, False
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:
        return None, False
    if scaf is None or scaf.GetNumHeavyAtoms() == 0:
        return None, True
    return Chem.MolToSmiles(scaf, canonical=True), False


def temporal_stratum(year: float) -> str:
    if pd.isna(year):
        return STRATUM_UNK
    if year <= PUBYEAR_LE_MAX:
        return STRATUM_LE
    return STRATUM_GE


def load_affitox_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return smiles-level rows, properly aggregated affinity pairs, exclusions.

    Affinity pairs are grouped by (UniProt, parent InChIKey) *before* the
    temporal label is assigned. Smiles-level rows inherit that pair stratum
    so Pearson uses the same membership definition as overlap counts.
    """
    frames = []
    for pdb, stem in PDB_TO_CSV.items():
        path = LIGANDS / f"{stem}_nodubl.csv"
        df = pd.read_csv(path, sep=";")
        keep = df[
            [
                "canonical_smiles",
                "molecule_chembl_id",
                "parent_molecule_chembl_id",
                "target_chembl_id",
                "document_year",
            ]
        ].copy()
        keep["pdb"] = pdb.upper()
        keep["document_year"] = pd.to_numeric(keep["document_year"], errors="coerce")
        frames.append(keep)
    raw = pd.concat(frames, ignore_index=True)
    raw["uniprot"] = raw["pdb"].map(PDB_TO_UNIPROT)
    raw["parent_inchikey"] = raw["canonical_smiles"].map(parent_inchikey)
    scaf = raw["canonical_smiles"].map(murcko_smiles)
    raw["scaffold"] = [s for s, _ in scaf]
    raw["empty_murcko"] = [e for _, e in scaf]

    smiles_rows = (
        raw.groupby(["pdb", "canonical_smiles"], as_index=False)
        .agg(
            molecule_chembl_id=("molecule_chembl_id", "first"),
            parent_molecule_chembl_id=("parent_molecule_chembl_id", "first"),
            target_chembl_id=("target_chembl_id", "first"),
            uniprot=("uniprot", "first"),
            parent_inchikey=("parent_inchikey", "first"),
            scaffold=("scaffold", "first"),
            empty_murcko=("empty_murcko", "first"),
            min_document_year_smiles=("document_year", "min"),
            n_records=("document_year", "size"),
        )
    )
    smiles_rows["deposition_year"] = smiles_rows["pdb"].map(PDB_DEPOSITION_YEAR)

    affinity_pairs = (
        raw.dropna(subset=["uniprot", "parent_inchikey"])
        .groupby(["uniprot", "parent_inchikey"], as_index=False)
        .agg(
            min_document_year=("document_year", "min"),
            n_records=("document_year", "size"),
            n_pdb=("pdb", "nunique"),
            n_smiles=("canonical_smiles", "nunique"),
            scaffold=("scaffold", "first"),
            empty_murcko=("empty_murcko", "max"),
            pdb_list=("pdb", lambda s: ";".join(sorted(set(s)))),
        )
    )
    affinity_pairs["temporal_stratum"] = affinity_pairs["min_document_year"].map(
        temporal_stratum
    )
    affinity_pairs["affinity_pair_id"] = (
        affinity_pairs["uniprot"] + "|" + affinity_pairs["parent_inchikey"]
    )

    pair_stratum = affinity_pairs.set_index(["uniprot", "parent_inchikey"])[
        ["temporal_stratum", "min_document_year", "affinity_pair_id"]
    ]
    smiles_rows = smiles_rows.merge(
        pair_stratum, on=["uniprot", "parent_inchikey"], how="left"
    )
    smiles_rows["temporal_stratum"] = smiles_rows["temporal_stratum"].fillna(STRATUM_UNK)
    smiles_rows["affinity_pair_id"] = smiles_rows["uniprot"].fillna("") + "|" + smiles_rows[
        "parent_inchikey"
    ].fillna("")

    exclusions = _exclusion_counts(raw, smiles_rows, affinity_pairs)
    return smiles_rows, affinity_pairs, exclusions


def _exclusion_counts(
    raw: pd.DataFrame, smiles_rows: pd.DataFrame, affinity_pairs: pd.DataFrame
) -> pd.DataFrame:
    n_raw = len(raw)
    n_smiles = len(smiles_rows)
    n_fail_key = int(smiles_rows["parent_inchikey"].isna().sum())
    n_aff = len(affinity_pairs)
    n_le = int((affinity_pairs["temporal_stratum"] == STRATUM_LE).sum())
    n_ge = int((affinity_pairs["temporal_stratum"] == STRATUM_GE).sum())
    n_unk = int((affinity_pairs["temporal_stratum"] == STRATUM_UNK).sum())
    n_empty = int(smiles_rows["empty_murcko"].sum())
    n_collapse = n_smiles - n_fail_key - n_aff
    rows = [
        {"step": "raw ChEMBL/BindingDB records after per-target nodubl CSVs", "n": n_raw, "reason": "input"},
        {"step": "unique (PDB, canonical_smiles) evaluation rows", "n": n_smiles, "reason": "score tables are stored at this grain"},
        {"step": "evaluation rows with failed parent InChIKey", "n": n_fail_key, "reason": "RDKit parse / InChI failure; excluded from affinity-pair key"},
        {"step": "unique (UniProt, parent InChIKey) affinity pairs", "n": n_aff, "reason": "primary overlap unit; min document_year across all contributing records"},
        {"step": "affinity pairs collapsed from multiple SMILES or PDBs", "n": max(n_collapse, 0), "reason": "same parent InChIKey after salt/fragment stripping"},
        {"step": "affinity pairs with unknown publication year", "n": n_unk, "reason": "excluded from both temporal strata"},
        {"step": "affinity pairs with publication year <= 2022", "n": n_le, "reason": "temporal proxy stratum pubyear_le_2022"},
        {"step": "affinity pairs with publication year >= 2023", "n": n_ge, "reason": "temporal proxy stratum pubyear_ge_2023"},
        {"step": "evaluation rows with empty Murcko scaffold", "n": n_empty, "reason": "acyclic parents; excluded from scaffold keys"},
    ]
    return pd.DataFrame(rows)


def _map_lp_row(smiles: str) -> tuple[str | None, str | None, bool]:
    key = parent_inchikey(smiles)
    scaf, empty = murcko_smiles(smiles)
    return key, scaf, empty


def load_lp_features() -> pd.DataFrame:
    if LP_CACHE.exists():
        feat = pd.read_csv(LP_CACHE)
        feat["date"] = pd.to_datetime(feat["date"], errors="coerce")
        if "year" not in feat.columns:
            feat["year"] = feat["date"].dt.year
        return feat
    if not LP_PATH.exists():
        raise FileNotFoundError(f"Missing LP-PDBBind table: {LP_PATH}")
    lp = pd.read_csv(LP_PATH)
    pdb_col = lp.columns[0]
    feat = pd.DataFrame(
        {
            "pdb": lp[pdb_col].astype(str).str.upper(),
            "smiles": lp["smiles"],
            "seq": lp["seq"].astype(str),
            "date": pd.to_datetime(lp["date"], errors="coerce"),
            "new_split": lp.get("new_split"),
        }
    )
    mapped = feat["smiles"].map(_map_lp_row)
    feat["parent_inchikey"] = [m[0] for m in mapped]
    feat["scaffold"] = [m[1] for m in mapped]
    feat["empty_murcko"] = [m[2] for m in mapped]
    feat["year"] = feat["date"].dt.year
    LP_CACHE.parent.mkdir(parents=True, exist_ok=True)
    feat.to_csv(LP_CACHE, index=False)
    return feat


def _read_split(path: Path) -> set[str]:
    return {line.strip().upper() for line in path.read_text().splitlines() if line.strip()}


def load_dynamicbind_processed(lp: pd.DataFrame) -> dict:
    """Processed DynamicBind training pool: Zenodo clash-annotated table ∩ timesplit.

    The published loader does not drop complexes by a clashScore cutoff. The
    CSV ``group`` labels already equal the EquiBind/DiffDock timesplit after
    preprocessing: ``train`` / ``valid`` / ``test``, plus ``remove``.
    Ligand InChIKeys are taken from LP-PDBBind crystal SMILES of the same PDB
    because the Zenodo table stores CCD codes, not SMILES.
    """
    if not DYN_CSV.exists():
        raise FileNotFoundError(f"Missing DynamicBind processed table: {DYN_CSV}")
    if not SPLIT_TRAIN.exists():
        raise FileNotFoundError(f"Missing timesplit: {SPLIT_TRAIN}")
    d3 = pd.read_csv(DYN_CSV)
    d3["pdb_u"] = d3["pdb"].astype(str).str.upper()
    d3["uid"] = d3["uid"].astype(str)
    timesplit_train = _read_split(SPLIT_TRAIN)
    timesplit_val = _read_split(SPLIT_VAL) if SPLIT_VAL.exists() else set()
    train = d3.loc[d3["group"].astype(str) == "train"].copy()
    train_pdbs = set(train["pdb_u"])
    expected = timesplit_train & set(d3["pdb_u"])
    if train_pdbs != expected:
        raise AssertionError(
            f"DynamicBind group==train ({len(train_pdbs)}) does not equal "
            f"timesplit_train ∩ CSV ({len(expected)})"
        )
    lp_idx = lp.drop_duplicates("pdb").set_index("pdb")
    inchikeys = []
    scaffolds = []
    for pdb in train["pdb_u"]:
        if pdb in lp_idx.index:
            inchikeys.append(lp_idx.at[pdb, "parent_inchikey"] if "parent_inchikey" in lp_idx.columns else None)
            scaffolds.append(lp_idx.at[pdb, "scaffold"] if "scaffold" in lp_idx.columns else None)
        else:
            inchikeys.append(None)
            scaffolds.append(None)
    train = train.copy()
    train["parent_inchikey"] = inchikeys
    train["scaffold"] = scaffolds
    pairs = {
        f"{uid}|{key}"
        for uid, key in zip(train["uid"], train["parent_inchikey"])
        if uid and isinstance(key, str) and key
    }
    affitox_pdbs = {p.upper() for p in PDB_TO_CSV}
    hit = sorted(affitox_pdbs & train_pdbs)
    return {
        "pdb": train_pdbs,
        "inchikey": set(k for k in train["parent_inchikey"] if isinstance(k, str) and k),
        "scaffold": set(s for s in train["scaffold"] if isinstance(s, str) and s),
        "pair": pairs,
        "uniprot": set(train["uid"].dropna()),
        "n_pool": int(len(train)),
        "n_csv": int(len(d3)),
        "n_timesplit_train": int(len(timesplit_train)),
        "n_timesplit_val": int(len(timesplit_val)),
        "n_remove": int((d3["group"].astype(str) == "remove").sum()),
        "n_absent_from_csv": int(len(timesplit_train - set(d3["pdb_u"]))),
        "affitox_pdb_hits": hit,
        "clashscore_note": (
            "clashScore is annotation; training pool is group==train "
            "(timesplit training list after preprocessing). No clashScore cutoff in the loader."
        ),
    }


def extract_crossdocked_train_pdbs(tgz: Path, cache: Path) -> set[str]:
    """Unique receptor PDB IDs from CrossDocked2020 v1.3 *train*.types files."""
    pdbs: set[str] = set()
    n_files = 0
    with tarfile.open(tgz, "r:gz") as tar:
        for member in tar:
            name = member.name
            base = Path(name).name.lower()
            if not base.endswith(".types"):
                continue
            if "train" not in base:
                continue
            n_files += 1
            handle = tar.extractfile(member)
            if handle is None:
                continue
            for raw in handle:
                try:
                    line = raw.decode("utf-8", errors="ignore")
                except Exception:
                    continue
                parts = line.split()
                rec = parts[3] if len(parts) >= 4 else line
                m = PDB_ID_RE.search(rec)
                if m:
                    pdbs.add(m.group(1).upper())
    if not pdbs:
        raise RuntimeError(f"No PDB IDs extracted from {tgz} ({n_files} train types files)")
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text("\n".join(sorted(pdbs)) + "\n")
    (cache.parent / "crossdocked_v1.3_train_types_file_count.txt").write_text(
        f"n_train_types_files={n_files}\nn_unique_pdb={len(pdbs)}\n"
    )
    return pdbs


def load_gnina_types_pdbs() -> set[str]:
    if GNINA_PDB_CACHE.exists():
        return _read_split(GNINA_PDB_CACHE)
    if not GNINA_TYPES_TGZ.exists():
        raise FileNotFoundError(
            f"Missing CrossDocked v1.3 types archive {GNINA_TYPES_TGZ}. "
            "Download from http://bits.csb.pitt.edu/files/crossdock2020/"
            "CrossDocked2020_v1.3_types.tgz"
        )
    print(f"Extracting train PDB IDs from {GNINA_TYPES_TGZ.name} (once)...")
    return extract_crossdocked_train_pdbs(GNINA_TYPES_TGZ, GNINA_PDB_CACHE)


def load_sifts_panel_map() -> dict[str, set[str]]:
    """PDB ID -> AffiTox-panel UniProt accessions present in any chain (SIFTS)."""
    if not SIFTS_PATH.exists():
        raise FileNotFoundError(
            f"Missing SIFTS mapping {SIFTS_PATH}; download pdb_chain_uniprot.tsv.gz from PDBe"
        )
    panel = set(PDB_TO_UNIPROT.values())
    sifts = pd.read_csv(SIFTS_PATH, sep="\t", comment="#")
    sifts["PDB"] = sifts["PDB"].astype(str).str.upper()
    sifts = sifts[sifts["SP_PRIMARY"].isin(panel)]
    out: dict[str, set[str]] = {}
    for pdb, sub in sifts.groupby("PDB"):
        out[str(pdb)] = set(sub["SP_PRIMARY"].astype(str))
    return out


def assign_lp_uniprot_sifts(feat: pd.DataFrame, sifts_map: dict[str, set[str]]) -> pd.DataFrame:
    """Unique SIFTS mapping onto the AffiTox UniProt panel.

    Ambiguous PDBs (two or more panel UniProts) are not distributed across
    proteins and do not contribute to pair-level source exposure.
    """
    mapping = []
    matched = []
    n_hits = []
    panel_list = []
    for pdb in feat["pdb"].astype(str):
        hits = sifts_map.get(pdb, set())
        n_hits.append(len(hits))
        panel_list.append(";".join(sorted(hits)))
        if len(hits) == 1:
            acc = next(iter(hits))
            mapping.append("unique")
            matched.append(acc)
        elif len(hits) > 1:
            mapping.append("ambiguous")
            matched.append("")
        else:
            mapping.append("unmapped")
            matched.append("")
    out = feat.copy()
    out["sifts_mapping"] = mapping
    out["matched_uniprot"] = matched
    out["n_panel_uniprot_sifts"] = n_hits
    out["panel_uniprots"] = panel_list
    return out


def corpus_sets(feat: pd.DataFrame, train_mask: pd.Series) -> dict[str, set[str]]:
    sub = feat.loc[train_mask]
    unique = sub.loc[sub["sifts_mapping"] == "unique"]
    pair_ids = {
        f"{row.matched_uniprot}|{row.parent_inchikey}"
        for row in unique.itertuples(index=False)
        if row.matched_uniprot and row.parent_inchikey
    }
    uniprot_present = set()
    for val in sub["panel_uniprots"].fillna(""):
        for acc in str(val).split(";"):
            if acc:
                uniprot_present.add(acc)
    return {
        "pdb": set(sub["pdb"].dropna().astype(str)),
        "inchikey": set(sub["parent_inchikey"].dropna()),
        "scaffold": set(sub["scaffold"].dropna()),
        "pair": pair_ids,
        "uniprot": uniprot_present,
        "n_ambiguous_pdb_rows": int((sub["sifts_mapping"] == "ambiguous").sum()),
        "n_unmapped_pdb_rows": int((sub["sifts_mapping"] == "unmapped").sum()),
        "n_unique_pdb_rows": int((sub["sifts_mapping"] == "unique").sum()),
    }


def overlap_row(
    method: str,
    method_id: str,
    dimension: str,
    n_items: int | None,
    n_overlap: int | None,
    definition: str,
    status: str,
    evidence_level: str,
) -> dict:
    frac = np.nan
    if n_items and n_overlap is not None and n_items > 0:
        frac = float(n_overlap) / float(n_items)
    return {
        "method": method,
        "method_id": method_id,
        "overlap_dimension": dimension,
        "n_items": n_items if n_items is not None else np.nan,
        "n_overlap": n_overlap if n_overlap is not None else np.nan,
        "frac_overlap": frac,
        "status": status,
        "evidence_level": evidence_level,
        "definition": definition,
    }


def overlap_tables(
    smiles_rows: pd.DataFrame,
    affinity_pairs: pd.DataFrame,
    dyn_sets: dict,
    gnina_pdbs: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    n_smiles_pairs = len(smiles_rows)
    n_aff_pairs = len(affinity_pairs)
    pre_smiles = smiles_rows["temporal_stratum"] == STRATUM_LE
    pre_aff = affinity_pairs["temporal_stratum"] == STRATUM_LE

    scafs = set(smiles_rows["scaffold"].dropna())
    pre_scafs = set(smiles_rows.loc[pre_smiles, "scaffold"].dropna())

    n_dyn_pdb = sum(pdb.upper() in dyn_sets["pdb"] for pdb in PDB_TO_CSV)
    n_gnina_pdb = sum(pdb.upper() in gnina_pdbs for pdb in PDB_TO_CSV)
    dyn_pdb_ids = dyn_sets.get("affitox_pdb_hits") or sorted(
        p.upper() for p in PDB_TO_CSV if p.upper() in dyn_sets["pdb"]
    )
    gnina_pdb_ids = sorted(p.upper() for p in PDB_TO_CSV if p.upper() in gnina_pdbs)
    plapt_note = (
        "PLAPT was trained on the first 100,000 rows of HuggingFace "
        "jglaser/binding_affinity (accessed 2024-01-27 per the preprint). "
        "That subset is not an immutable pair-level manifest: composition "
        "depends on the dataset snapshot and parquet row order. Exact "
        "membership is therefore indeterminate; we do not reconstruct an HF-100k proxy."
    )

    method_rows = [
        overlap_row(
            "Boltz-2",
            "boltz2",
            "source-domain exposure (BindingDB/ChEMBL families)",
            n_smiles_pairs,
            n_smiles_pairs,
            "All AffiTox records originate from database families reported as Boltz-2 affinity-training sources. Complete source-domain overlap, not pair membership in the filtered training set.",
            "complete_source_exposure",
            "proxy",
        ),
        overlap_row(
            "Boltz-2",
            "boltz2",
            "publication year <= 2022 affinity pairs (UniProt, parent InChIKey)",
            n_aff_pairs,
            int(pre_aff.sum()),
            f"min document_year <= {PUBYEAR_LE_MAX} after aggregating on (UniProt, parent InChIKey). Publication year is a temporal proxy, not ingestion date or training membership.",
            STRATUM_LE,
            "proxy",
        ),
        overlap_row(
            "Boltz-2",
            "boltz2",
            "internal AffiTox scaffold coverage among publication year <= 2022 pairs",
            len(scafs),
            len(scafs & pre_scafs),
            "Bemis-Murcko scaffold appears in any AffiTox pair whose aggregated document_year is <= 2022. Internal temporal coverage, not AffiTox ∩ Boltz-2 training scaffolds.",
            "internal_temporal_coverage",
            "proxy",
        ),
        overlap_row(
            "Boltz-2",
            "boltz2",
            "exact filtered training-set pair membership",
            None,
            None,
            "Boltz-2 training/evaluation manifests are unpublished (repository: coming soon). Exact pair-level membership cannot be recovered.",
            "unknown",
            "indeterminate",
        ),
        overlap_row(
            "PLAPT",
            "plapt",
            "exact pair/scaffold exposure vs PLAPT training set",
            None,
            None,
            plapt_note,
            "indeterminate",
            "indeterminate",
        ),
        overlap_row(
            "DynamicBind",
            "dynamicbind",
            "exact PDB intersection with processed training table",
            16,
            n_dyn_pdb,
            (
                f"AffiTox PDB ID in Zenodo d3_with_clash_info.csv group==train "
                f"(n_pool={dyn_sets.get('n_pool')}; timesplit train n="
                f"{dyn_sets.get('n_timesplit_train')}; {dyn_sets.get('n_absent_from_csv')} "
                f"timesplit complexes absent from the processed CSV; "
                f"{dyn_sets.get('n_remove')} labelled remove). Hits: "
                f"{', '.join(dyn_pdb_ids) or 'none'}. "
                f"{dyn_sets.get('clashscore_note')}"
            ),
            "processed_source_table_pdb",
            "processed_source_table",
        ),
        overlap_row(
            "DynamicBind",
            "dynamicbind",
            "receptor deposition-year proxy (<2019)",
            16,
            int(sum(y < 2019 for y in PDB_DEPOSITION_YEAR.values())),
            "Calendar deposition year < 2019. Proxy only: this PDB ID need not appear in the processed table, and a different PDB of the same protein may.",
            "deposition_year_proxy",
            "proxy",
        ),
        overlap_row(
            "DynamicBind",
            "dynamicbind",
            "ligand InChIKey intersection with processed training table",
            n_aff_pairs,
            int(affinity_pairs["parent_inchikey"].isin(dyn_sets["inchikey"]).sum()),
            "AffiTox parent InChIKey present among crystal ligands of processed-table PDBs (InChIKey from LP-PDBBind SMILES of the same PDB; Zenodo table stores CCD codes).",
            "processed_source_table_ligand",
            "processed_source_table",
        ),
        overlap_row(
            "DynamicBind",
            "dynamicbind",
            "pair intersection with processed training table (UniProt, parent InChIKey)",
            n_aff_pairs,
            int(affinity_pairs["affinity_pair_id"].isin(dyn_sets["pair"]).sum()),
            "UniProt from the Zenodo table uid column; parent InChIKey of the co-crystal ligand. Not a SMILES dump from DynamicBind and not a guarantee that the AffiTox assay record was used for fitting.",
            "processed_source_table_pair",
            "processed_source_table",
        ),
        overlap_row(
            "DynamicBind",
            "dynamicbind",
            "scaffold intersection with processed training table",
            len(scafs),
            len(scafs & dyn_sets["scaffold"]),
            "AffiTox non-empty Bemis-Murcko scaffolds intersecting crystal-ligand scaffolds of processed-table PDBs.",
            "processed_source_table_scaffold",
            "processed_source_table",
        ),
        overlap_row(
            "DynamicBind",
            "dynamicbind",
            "UniProt present in processed training table",
            16,
            sum(PDB_TO_UNIPROT[p.upper()] in dyn_sets["uniprot"] for p in PDB_TO_CSV),
            "AffiTox UniProt appears as uid in at least one processed-table training complex (possibly a different PDB of the same protein).",
            "processed_source_table_uniprot",
            "processed_source_table",
        ),
        overlap_row(
            "GNINA 1.3",
            "gnina",
            "exact PDB intersection with CrossDocked2020 v1.3 training types",
            16,
            n_gnina_pdb,
            (
                f"Receptor PDB identifiers extracted from CrossDocked2020 v1.3 "
                f"*train*.types files (n_types_pdb={len(gnina_pdbs)}). Hits: "
                f"{', '.join(gnina_pdb_ids) or 'none'}. Ligand- and pair-level "
                "membership require molecular structures / molcache and remain indeterminate. "
                "Primary AffiTox score: CNNaffinity of the max-CNNscore pose."
            ),
            "exact_types_pdb",
            "exact_types",
        ),
        overlap_row(
            "GNINA 1.3",
            "gnina",
            "receptor deposition-year proxy (<2020)",
            16,
            int(sum(y < 2020 for y in PDB_DEPOSITION_YEAR.values())),
            "Calendar deposition year < 2020 as a CrossDocked2020-era proxy. Not equivalent to types-file membership.",
            "deposition_year_proxy",
            "proxy",
        ),
        overlap_row(
            "GNINA 1.3",
            "gnina",
            "ligand / target-ligand pair membership vs CrossDocked types",
            None,
            None,
            "Types files name docked poses, not standardized parent InChIKeys. Exact ligand- and pair-level membership cannot be established without the corresponding molecular structures and molecule cache.",
            "indeterminate",
            "indeterminate",
        ),
        overlap_row(
            "QVina2",
            "qvina",
            "supervised bioactivity-training overlap",
            None,
            None,
            "QVina2 does not use a supervised ligand-target affinity training set and is a non-ML anchor, not a numeric zero from pair matching.",
            "not_applicable",
            "not_applicable",
        ),
    ]
    method_overlap = pd.DataFrame(method_rows)

    by_target_rows = []
    for pdb, sub in smiles_rows.groupby("pdb"):
        n = len(sub)
        n_pre = int((sub["temporal_stratum"] == STRATUM_LE).sum())
        n_ge = int((sub["temporal_stratum"] == STRATUM_GE).sum())
        n_unk = int((sub["temporal_stratum"] == STRATUM_UNK).sum())
        sc = set(sub["scaffold"].dropna())
        pre_sc = set(sub.loc[sub["temporal_stratum"] == STRATUM_LE, "scaffold"].dropna())
        pair_ids = sub.loc[sub["parent_inchikey"].notna(), "affinity_pair_id"]
        aff_sub = affinity_pairs.loc[affinity_pairs["uniprot"] == PDB_TO_UNIPROT[pdb]]
        by_target_rows.append(
            {
                "pdb": pdb,
                "target": TARGET_LABEL[pdb],
                "uniprot": PDB_TO_UNIPROT[pdb],
                "n_smiles_pairs": n,
                "n_affinity_pairs": int(len(aff_sub)),
                "n_pubyear_le_2022": n_pre,
                "n_pubyear_ge_2023": n_ge,
                "n_unknown_year": n_unk,
                "frac_pubyear_le_2022": n_pre / n if n else np.nan,
                "frac_pubyear_ge_2023": n_ge / n if n else np.nan,
                "frac_unknown_year": n_unk / n if n else np.nan,
                "n_scaffolds": len(sc),
                "n_internal_le2022_scaffolds": len(pre_sc),
                "frac_internal_le2022_scaffold_coverage": (len(sc & pre_sc) / len(sc)) if sc else np.nan,
                "n_empty_murcko": int(sub["empty_murcko"].sum()),
                "dynamicbind_processed_pdb": pdb in dyn_sets["pdb"],
                "gnina_types_pdb": pdb in gnina_pdbs,
                "frac_dynamicbind_ligand_inchikey": float(sub["parent_inchikey"].isin(dyn_sets["inchikey"]).mean()) if n else np.nan,
                "frac_dynamicbind_processed_pair": float(pair_ids.isin(dyn_sets["pair"]).mean()) if len(pair_ids) else np.nan,
            }
        )
    by_target = pd.DataFrame(by_target_rows).sort_values("pdb")

    pair_flags = smiles_rows.copy()
    pair_flags["boltz2_source_domain_exposed"] = True
    pair_flags["boltz2_pubyear_le_2022"] = pre_smiles
    pair_flags["plapt_exposure"] = "indeterminate"
    pair_flags["qvina_supervised_overlap"] = "not_applicable"
    pair_flags["dynamicbind_processed_pdb"] = pair_flags["pdb"].isin(dyn_sets["pdb"])
    pair_flags["gnina_types_pdb"] = pair_flags["pdb"].isin(gnina_pdbs)
    pair_flags["dynamicbind_ligand_inchikey"] = pair_flags["parent_inchikey"].isin(dyn_sets["inchikey"])
    pair_flags["dynamicbind_processed_pair"] = pair_flags["affinity_pair_id"].isin(dyn_sets["pair"])
    pair_flags["gnina_ligand_pair_membership"] = "indeterminate"
    pair_flags["internal_le2022_scaffold"] = False
    for pdb, sub in smiles_rows.groupby("pdb"):
        seen_sc = set(sub.loc[sub["temporal_stratum"] == STRATUM_LE, "scaffold"].dropna())
        idx = pair_flags["pdb"] == pdb
        pair_flags.loc[idx, "internal_le2022_scaffold"] = pair_flags.loc[idx, "scaffold"].isin(seen_sc)

    return method_overlap, by_target, pair_flags


def common_mask(df: pd.DataFrame) -> pd.Series:
    pki = ensure_pki(df)
    ok = pki.notna() & np.isfinite(pki)
    for metric in PRIMARY_METRICS.values():
        if metric not in df.columns:
            return pd.Series(False, index=df.index)
        s = pd.to_numeric(df[metric], errors="coerce")
        ok &= s.notna() & np.isfinite(s)
    return ok


def _score_xy(df: pd.DataFrame, metric: str) -> tuple[np.ndarray, np.ndarray]:
    pki = ensure_pki(df)
    score = pd.to_numeric(df.get(metric), errors="coerce")
    mask = pki.notna() & score.notna() & np.isfinite(pki) & np.isfinite(score)
    x = pki[mask].to_numpy(dtype=float)
    y = score[mask].to_numpy(dtype=float)
    if metric in SIGN_FLIP_METRICS:
        y = -y
    return x, y


def _pearson_safe(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 3 or np.std(x, ddof=1) == 0 or np.std(y, ddof=1) == 0:
        return np.nan, np.nan
    r, p = pearsonr(x, y)
    return float(r), float(p)


def stratified_pearson(smiles_rows: pd.DataFrame) -> pd.DataFrame:
    flags = smiles_rows.set_index(["pdb", "canonical_smiles"])["temporal_stratum"]
    rows = []
    for pdb in sorted(PDB_TO_CSV):
        pdb_u = pdb.upper()
        merged = load_merged(TABLES, pdb)
        if "canonical_smiles" not in merged.columns:
            continue
        merged = merged.copy()
        merged["_stratum"] = merged["canonical_smiles"].map(
            lambda s, t=pdb_u: flags.get((t, s), STRATUM_UNK)
        )
        for stratum in (STRATUM_LE, STRATUM_GE):
            sub = merged.loc[merged["_stratum"] == stratum]
            mask = common_mask(sub)
            common = sub.loc[mask]
            n_common = int(len(common))
            for method_id, metric in PRIMARY_METRICS.items():
                x, y = _score_xy(common, metric)
                n = len(x)
                if n_common > 0 and n != n_common:
                    raise AssertionError(
                        f"Common-mask length mismatch {pdb_u} {method_id}: "
                        f"len(x)={n} n_common={n_common}"
                    )
                r, p = _pearson_safe(x, y)
                lo = hi = np.nan
                if n >= MIN_N and np.isfinite(r):
                    lo, hi, _ = bootstrap_ci(
                        x, y, method="pearson", n_boot=N_BOOT, seed=BOOT_SEED
                    )
                rows.append(
                    {
                        "pdb": pdb_u,
                        "target": TARGET_LABEL[pdb_u],
                        "uniprot": PDB_TO_UNIPROT[pdb_u],
                        "stratum": stratum,
                        "method_id": method_id,
                        "method": METHOD_LABELS[method_id],
                        "metric": metric,
                        "n_stratum_rows": int(len(sub)),
                        "n_common_mask": n_common,
                        "n_points": n,
                        "pearson_r": r,
                        "pearson_p": p,
                        "pearson_ci_low": lo,
                        "pearson_ci_high": hi,
                        "eligible": n >= MIN_N and np.isfinite(r),
                        "mask": "all_five_methods_finite",
                    }
                )
    return pd.DataFrame(rows)


def dual_stratum_targets(corr: pd.DataFrame) -> list[str]:
    """Targets with n_common >= MIN_N in both publication-year strata."""
    ok = {}
    for stratum in (STRATUM_LE, STRATUM_GE):
        sub = corr.loc[corr["stratum"] == stratum]
        ok[stratum] = set(sub.loc[sub["n_common_mask"] >= MIN_N, "pdb"])
    return sorted(ok[STRATUM_LE] & ok[STRATUM_GE])


def stratum_delta_r(corr: pd.DataFrame, targets: list[str]) -> pd.DataFrame:
    """Within-target difference in Pearson r between publication-year strata.

    Not a paired ligand analysis: the two strata are independent compound sets.
    """
    rows = []
    pre = corr.loc[corr["stratum"] == STRATUM_LE]
    post = corr.loc[corr["stratum"] == STRATUM_GE]
    for pdb in targets:
        for method_id in METHODS:
            a = pre.loc[(pre.pdb == pdb) & (pre.method_id == method_id)]
            b = post.loc[(post.pdb == pdb) & (post.method_id == method_id)]
            if a.empty or b.empty:
                continue
            r_pre = float(a["pearson_r"].iloc[0]) if pd.notna(a["pearson_r"].iloc[0]) else np.nan
            r_post = float(b["pearson_r"].iloc[0]) if pd.notna(b["pearson_r"].iloc[0]) else np.nan
            rows.append(
                {
                    "pdb": pdb,
                    "target": TARGET_LABEL[pdb],
                    "method_id": method_id,
                    "method": METHOD_LABELS[method_id],
                    "n_pre": int(a["n_points"].iloc[0]),
                    "n_post": int(b["n_points"].iloc[0]),
                    "eligible_pre": bool(a["eligible"].iloc[0]),
                    "eligible_post": bool(b["eligible"].iloc[0]),
                    "r_pre": r_pre,
                    "r_post": r_post,
                    "r_pre_ci_low": float(a["pearson_ci_low"].iloc[0]) if pd.notna(a["pearson_ci_low"].iloc[0]) else np.nan,
                    "r_pre_ci_high": float(a["pearson_ci_high"].iloc[0]) if pd.notna(a["pearson_ci_high"].iloc[0]) else np.nan,
                    "r_post_ci_low": float(b["pearson_ci_low"].iloc[0]) if pd.notna(b["pearson_ci_low"].iloc[0]) else np.nan,
                    "r_post_ci_high": float(b["pearson_ci_high"].iloc[0]) if pd.notna(b["pearson_ci_high"].iloc[0]) else np.nan,
                    "delta_r_ge2023_minus_le2022": (
                        r_post - r_pre if np.isfinite(r_pre) and np.isfinite(r_post) else np.nan
                    ),
                    "both_eligible": bool(a["eligible"].iloc[0] and b["eligible"].iloc[0]),
                    "ci_excludes_zero_descriptive": np.nan,
                }
            )
    return pd.DataFrame(rows)


def bootstrap_delta_ci(
    smiles_rows: pd.DataFrame, corr_delta: pd.DataFrame, targets: list[str]
) -> pd.DataFrame:
    """Independent-stratum percentile bootstrap CI for within-target Δr."""
    flags = smiles_rows.set_index(["pdb", "canonical_smiles"])["temporal_stratum"]
    out_rows = []
    for pdb in targets:
        merged = load_merged(TABLES, pdb.lower())
        merged = merged.copy()
        merged["_stratum"] = merged["canonical_smiles"].map(
            lambda s, t=pdb: flags.get((t, s), STRATUM_UNK)
        )
        pre_sub = merged.loc[merged["_stratum"] == STRATUM_LE]
        post_sub = merged.loc[merged["_stratum"] == STRATUM_GE]
        pre = pre_sub.loc[common_mask(pre_sub)]
        post = post_sub.loc[common_mask(post_sub)]
        for method_id, metric in PRIMARY_METRICS.items():
            x_pre, y_pre = _score_xy(pre, metric)
            x_post, y_post = _score_xy(post, metric)
            if len(x_pre) != len(pre) or len(x_post) != len(post):
                raise AssertionError(
                    f"Delta-r mask mismatch {pdb} {method_id}: "
                    f"pre {len(x_pre)}/{len(pre)} post {len(x_post)}/{len(post)}"
                )
            if len(x_pre) < MIN_N or len(x_post) < MIN_N:
                lo = hi = np.nan
            else:
                rng = np.random.default_rng(BOOT_SEED)
                deltas = []
                for _ in range(N_BOOT):
                    i = rng.integers(0, len(x_pre), len(x_pre))
                    j = rng.integers(0, len(x_post), len(x_post))
                    r1, _ = _pearson_safe(x_pre[i], y_pre[i])
                    r2, _ = _pearson_safe(x_post[j], y_post[j])
                    if np.isfinite(r1) and np.isfinite(r2):
                        deltas.append(r2 - r1)
                if deltas:
                    arr = np.asarray(deltas)
                    lo, hi = float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))
                else:
                    lo = hi = np.nan
            match = corr_delta.loc[
                (corr_delta.pdb == pdb) & (corr_delta.method_id == method_id)
            ]
            rec = dict(match.iloc[0]) if len(match) else {
                "pdb": pdb,
                "method_id": method_id,
                "method": METHOD_LABELS[method_id],
            }
            rec["delta_r_ci_low"] = lo
            rec["delta_r_ci_high"] = hi
            rec["n_boot"] = N_BOOT
            rec["ci_excludes_zero_descriptive"] = bool(
                np.isfinite(lo) and np.isfinite(hi) and (hi < 0 or lo > 0)
            )
            out_rows.append(rec)
    return pd.DataFrame(out_rows)


def wilcoxon_exploratory(corr: pd.DataFrame) -> pd.DataFrame:
    """Post-hoc exploratory Wilcoxon on publication-year<=2022 targets.

    Primary Holm family: Boltz-2 versus the four comparators (four tests),
    matching the scientific question reported in Methods. All ten pairwise
    raw p-values are retained; Holm over those ten is not used for claims.
    """
    rows = []
    stratum = STRATUM_LE
    sub = corr.loc[(corr["stratum"] == stratum) & (corr["eligible"])].copy()
    family = (
        "Exploratory post-hoc Pearson r (publication year <= 2022, common mask, n>=10); "
        "Holm family = Boltz-2 vs four comparators; not confirmatory"
    )
    for m1, m2 in combinations(METHODS, 2):
        a, b = [], []
        for pdb in sorted(sub["pdb"].unique()):
            v1 = sub.loc[(sub.pdb == pdb) & (sub.method_id == m1), "pearson_r"]
            v2 = sub.loc[(sub.pdb == pdb) & (sub.method_id == m2), "pearson_r"]
            if len(v1) and len(v2) and pd.notna(v1.iloc[0]) and pd.notna(v2.iloc[0]):
                a.append(float(v1.iloc[0]))
                b.append(float(v2.iloc[0]))
        p, n_nz, n_tgt, med, mean = wilcoxon_paired(a, b)
        rows.append(
            {
                "family": family,
                "stratum": stratum,
                "comparison": f"{METHOD_LABELS[m1]} vs {METHOD_LABELS[m2]}",
                "method_a": m1,
                "method_b": m2,
                "n_targets": n_tgt,
                "n_nonzero_diff": n_nz,
                "median_diff_a_minus_b": med,
                "mean_diff_a_minus_b": mean,
                "p_raw": p,
                "in_boltz2_holm_family": m1 == "boltz2" or m2 == "boltz2",
            }
        )
    result = pd.DataFrame(rows)
    result["p_holm"] = np.nan
    fam = result["in_boltz2_holm_family"]
    result.loc[fam, "p_holm"] = holm(result.loc[fam, "p_raw"].fillna(1.0).to_numpy(dtype=float))
    return result


def median_summary(corr: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (stratum, method_id), sub in corr.groupby(["stratum", "method_id"]):
        elig = sub.loc[sub["eligible"]]
        rows.append(
            {
                "stratum": stratum,
                "method_id": method_id,
                "method": METHOD_LABELS[method_id],
                "n_targets_all": int(sub["pdb"].nunique()),
                "n_targets_n_ge_min": int(elig["pdb"].nunique()),
                "median_r_all": float(sub["pearson_r"].median()) if len(sub) else np.nan,
                "median_r_eligible": float(elig["pearson_r"].median()) if len(elig) else np.nan,
                "mean_n_eligible": float(elig["n_points"].mean()) if len(elig) else np.nan,
                "note": (
                    "Medians across strata are not a drop: eligible target sets differ."
                    if stratum == STRATUM_GE
                    else "Common five-method mask; exploratory sensitivity analysis."
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(["stratum", "method"])


def provenance_metadata(
    n_empty: int,
    dual_targets: list[str],
    dyn_sets: dict,
    gnina_pdbs: set[str],
) -> dict:
    return {
        "rdkit_version": rdBase.rdkitVersion,
        "standardization": (
            "largest fragment; Uncharger; parent InChIKey via RDKit MolToInchiKey; "
            "stereo retained; no explicit tautomer enumeration (InChI canonicalization only)"
        ),
        "scaffold": "RDKit Bemis-Murcko on parent mol; empty (acyclic) scaffolds excluded from overlap keys",
        "n_empty_murcko_rows": n_empty,
        "pair_key": "(canonical UniProt accession, parent-ligand InChIKey); aggregated before temporal labelling",
        "temporal_proxy": (
            f"{STRATUM_LE}: min document_year <= {PUBYEAR_LE_MAX}; "
            f"{STRATUM_GE}: >= {PUBYEAR_GE_MIN}; {STRATUM_UNK} excluded from strata. "
            "Publication year is not database ingestion date or training membership."
        ),
        "evidence_policy": (
            "Method-specific: processed-source-table intersections (DynamicBind), "
            "exact types-file PDB intersections (GNINA 1.3), and proxy-based exposure "
            "estimates (Boltz-2 source-domain/temporal; PLAPT indeterminate). These are "
            "not equivalent membership measures across methods."
        ),
        "protein_mapping_dynamicbind": (
            "UniProt from Zenodo d3_with_clash_info.csv uid; ligand InChIKey from "
            "LP-PDBBind crystal SMILES of the same PDB (Zenodo table stores CCD codes)."
        ),
        "boltz2_affinity_cutoff_note": (
            "Boltz-2 structure cutoff 2023-06-01 is not automatically the affinity-head cutoff. "
            "Affinity-training manifests are unpublished. Independent reconstructions use "
            "ChEMBL v34 and BindingDB/PubChem as of June 2023 with InChIKey and UniProt "
            "(Shimizu et al., Chem-Bio Informatics Journal 2026)."
        ),
        "dynamicbind_source": (
            f"Zenodo 10.5281/zenodo.10429051 d3_with_clash_info.csv group==train "
            f"(n={dyn_sets.get('n_pool')}) equals timesplit_no_lig_overlap_train ∩ CSV. "
            f"Timesplit train n={dyn_sets.get('n_timesplit_train')}; "
            f"{dyn_sets.get('n_absent_from_csv')} timesplit complexes never entered the "
            f"processed CSV; {dyn_sets.get('n_remove')} labelled remove. "
            f"AffiTox PDB hits: {', '.join(dyn_sets.get('affitox_pdb_hits') or [])}. "
            f"{dyn_sets.get('clashscore_note')}"
        ),
        "gnina_source": (
            f"CrossDocked2020 v1.3 training types PDB identifiers "
            f"(n_unique_pdb={len(gnina_pdbs)}). Ligand/pair membership indeterminate "
            "without molcache/SDF. Primary AffiTox score: CNNaffinity of the max-CNNscore pose."
        ),
        "plapt": (
            "First-100k HuggingFace jglaser/binding_affinity subset is snapshot- and "
            "row-order-dependent; exact pair membership left indeterminate (no HF-100k reconstruction)."
        ),
        "qvina": "Supervised bioactivity-training overlap not applicable (non-ML anchor).",
        "pearson": (
            f"Identical five-method finite-score masks with assert len(x)==n_common; "
            f"ensure_pki used for both mask and Pearson; bootstrap B={N_BOOT}, seed={BOOT_SEED}; "
            f"eligible n>={MIN_N}. Dual-stratum targets auto-selected by n_common>={MIN_N} "
            f"in both strata: {dual_targets}."
        ),
        "wilcoxon_holm_family": "Boltz-2 versus four comparators (four tests), exploratory post-hoc",
        "delta_r": (
            "Within-target difference in Pearson r between independent publication-year "
            "compound sets; not a paired-ligand analysis. Percentile bootstrap intervals "
            "that exclude zero are descriptive and are not multiplicity-adjusted tests."
        ),
        "analysis_status": (
            "method-specific training-exposure analysis; not a complete leakage audit "
            "and not a prospective seen-unseen evaluation"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-pearson",
        action="store_true",
        help="Write overlap tables only; reuse existing Pearson CSVs if present.",
    )
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    EXTERNAL.mkdir(parents=True, exist_ok=True)

    print("Loading AffiTox pairs; aggregating on (UniProt, parent InChIKey)...")
    smiles_rows, affinity_pairs, exclusions = load_affitox_tables()
    n_empty = int(smiles_rows["empty_murcko"].sum())
    affinity_pairs.to_csv(OUT / "affinity_pairs_uniprot_inchikey.csv", index=False)
    exclusions.to_csv(OUT / "exclusion_audit.csv", index=False)

    print("Loading DynamicBind processed table ∩ timesplit...")
    lp = load_lp_features()
    dyn_sets = load_dynamicbind_processed(lp)
    print("Loading GNINA CrossDocked v1.3 training-types PDB IDs...")
    gnina_pdbs = load_gnina_types_pdbs()

    method_overlap, by_target, pair_flags = overlap_tables(
        smiles_rows, affinity_pairs, dyn_sets, gnina_pdbs
    )
    method_overlap.to_csv(OUT / "overlap_fractions_by_method.csv", index=False)
    by_target.to_csv(OUT / "overlap_by_target.csv", index=False)
    pair_flags.to_csv(OUT / "pair_membership_flags.csv", index=False)

    corr_path = OUT / "pearson_by_bioactivity_stratum.csv"
    if args.skip_pearson and corr_path.exists():
        print("Skipping Pearson recompute; loading existing CSVs...")
        corr = pd.read_csv(corr_path)
        med = pd.read_csv(OUT / "pearson_stratum_median_summary.csv")
        delta = pd.read_csv(OUT / "pearson_paired_delta_dual_stratum.csv")
        tests = pd.read_csv(OUT / "wilcoxon_pearson_by_stratum.csv")
        dual = dual_stratum_targets(corr)
    else:
        print("Computing common-mask stratified Pearson r...")
        corr = stratified_pearson(smiles_rows)
        corr.to_csv(corr_path, index=False)
        med = median_summary(corr)
        med.to_csv(OUT / "pearson_stratum_median_summary.csv", index=False)

        dual = dual_stratum_targets(corr)
        print(f"Auto dual-stratum targets (n_common>={MIN_N} in both strata): {dual}")
        print("Computing within-target stratum Δr and exploratory Wilcoxon...")
        delta = stratum_delta_r(corr, dual)
        delta = bootstrap_delta_ci(smiles_rows, delta, dual)
        delta.to_csv(OUT / "pearson_paired_delta_dual_stratum.csv", index=False)
        tests = wilcoxon_exploratory(corr)
        tests.to_csv(OUT / "wilcoxon_pearson_by_stratum.csv", index=False)

    meta = provenance_metadata(n_empty, dual, dyn_sets, gnina_pdbs)
    (OUT / "exposure_analysis_metadata.json").write_text(json.dumps(meta, indent=2))

    print("\n=== Exclusion audit ===")
    print(exclusions.to_string(index=False))
    print("\n=== Exposure / overlap (method-specific) ===")
    show = method_overlap[
        ["method", "overlap_dimension", "n_items", "n_overlap", "frac_overlap", "evidence_level"]
    ]
    print(show.to_string(index=False))
    print("\n=== Median Pearson r (common mask; do not compare strata as a drop) ===")
    print(med.to_string(index=False))
    print("\n=== Within-target Δr (publication year >=2023 minus <=2022) ===")
    cols = [
        "pdb",
        "method",
        "n_pre",
        "n_post",
        "r_pre",
        "r_post",
        "delta_r_ge2023_minus_le2022",
        "delta_r_ci_low",
        "delta_r_ci_high",
        "ci_excludes_zero_descriptive",
        "both_eligible",
    ]
    print(delta[cols].to_string(index=False) if len(delta) else "(none)")
    print("\n=== Exploratory Wilcoxon (Boltz-2 Holm family of four) ===")
    print(
        tests.loc[
            tests.in_boltz2_holm_family,
            ["comparison", "n_targets", "p_raw", "p_holm"],
        ].to_string(index=False)
    )
    print(f"\nRDKit {meta['rdkit_version']}; wrote {OUT}")


if __name__ == "__main__":
    main()
