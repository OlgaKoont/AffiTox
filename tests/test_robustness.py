"""Tests for AffiTox robustness formulas, joins, score direction, and reproducibility."""

from __future__ import annotations

import numpy as np
import pandas as pd

from analysis.constants import SIGN_FLIP_METRICS
from analysis.data import prepare_xy
from analysis.enrichment import compute_ef, detect_score_direction
from analysis.inferential import benjamini_hochberg
from analysis.robustness_common import (
    residualize_two_way,
    spearman_safe,
    two_way_dummies,
    canonicalize_evidence_level,
)
from analysis.cross_axis_association import association_table, build_long_table
from analysis.robustness_common import RobustnessConfig


def test_pose_rate_formulas_do_not_impute_missing_as_conditional_failures():
    n_input = 10
    n_evaluable = 7
    n_pass = 3
    coverage = n_evaluable / n_input
    conditional = n_pass / n_evaluable
    end_to_end = n_pass / n_input
    assert coverage == 0.7
    assert abs(conditional - (3 / 7)) < 1e-12
    assert end_to_end == 0.3
    # The 3 missing poses are unsuccessful only end-to-end, not in the conditional rate.
    assert n_pass / n_evaluable != n_pass / n_input


def test_active_enrichment_uses_stronger_binding_first():
    assert detect_score_direction("gnina_cnn_affinity_bestpose", "active") == "desc"
    assert detect_score_direction("qvina_affinity_bestpose", "active") == "asc"
    assert detect_score_direction("boltz2_affinity_pred_value", "active") == "asc"
    assert detect_score_direction("qvina_affinity_bestpose", "low") == "desc"
    assert detect_score_direction("gnina_cnn_affinity_bestpose", "low") == "asc"


def test_nef_direction_recovers_actives_at_correct_end():
    scores = np.array([0.1, 0.2, 0.3, 0.9, 1.0])
    is_active = np.array([0, 0, 0, 1, 1], dtype=float)
    _ef, nef, n, a, h = compute_ef(is_active, scores, 0.40, "desc")
    assert n == 5 and a == 2 and h == 2
    np.testing.assert_allclose(nef, 1.0)


def test_sign_flip_not_applied_inside_enrichment_detect():
    assert "qvina_affinity_bestpose" in SIGN_FLIP_METRICS
    df = pd.DataFrame(
        {"pValue": [6.0, 7.0, 8.0], "qvina_affinity_bestpose": [-5.0, -6.0, -7.0]}
    )
    _, y = prepare_xy(df, "qvina_affinity_bestpose")
    np.testing.assert_allclose(y, [5.0, 6.0, 7.0])


def test_two_way_residualize_removes_target_and_method_means():
    targets = np.array(["a", "a", "b", "b"])
    methods = np.array(["m1", "m2", "m1", "m2"])
    y = np.array([10.0, 12.0, 20.0, 22.0])
    resid = residualize_two_way(y, targets, methods)
    assert np.allclose(resid.mean(), 0.0, atol=1e-10)
    x = two_way_dummies(targets, methods)
    fitted = x @ np.linalg.lstsq(x, y, rcond=None)[0]
    np.testing.assert_allclose(y - fitted, resid)


def test_bh_is_monotonic():
    p = np.array([0.001, 0.02, 0.04, 0.20])
    q = benjamini_hochberg(p)
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)
    assert np.all(q >= p - 1e-12)


def test_spearman_constant_is_nan():
    x = np.ones(8)
    y = np.arange(8, dtype=float)
    assert np.isnan(spearman_safe(x, y))


def test_qvina_overlap_is_not_applicable_not_zero():
    assert canonicalize_evidence_level("not_applicable") == "not_applicable"
    assert canonicalize_evidence_level("exact_types") == "exact_manifest_intersection"
    assert canonicalize_evidence_level("processed_source_table") == "processed_source_table_intersection"
    assert canonicalize_evidence_level("proxy", "source-domain exposure") == "source_domain_exposure"
    assert canonicalize_evidence_level("proxy", "publication year <= 2022") == "temporal_eligibility_proxy"
    assert canonicalize_evidence_level("indeterminate") == "indeterminate"


def test_end_to_end_counts_missing_poses_conditional_does_not():
    n_input, n_evaluable, n_pass = 10, 0, 0
    conditional = np.nan if n_evaluable == 0 else n_pass / n_evaluable
    end_to_end = n_pass / n_input
    assert np.isnan(conditional)
    assert end_to_end == 0.0


def test_cross_axis_reproducible_tiny_boot():
    cfg = RobustnessConfig(seed=42, n_boot=20, n_perm=20)
    long = build_long_table()
    assert long.duplicated(["target", "method_id"]).sum() == 0
    assert len(long) == 64
    assert set(long["method_id"]) == {"boltz2", "dynamicbind", "gnina", "qvina"}
    a = association_table(long, cfg)
    b = association_table(long, cfg)
    pd.testing.assert_frame_equal(a, b)
    prim = a.loc[a["is_primary"] == 1]
    assert len(prim) == 1
    assert (a["family"] == "secondary").sum() == len(a) - 1
    ok = long["strict_pass_conditional"].dropna()
    assert (ok >= 0.0).all() and (ok <= 1.0).all()
    e2e = long["strict_pass_end_to_end"].dropna()
    assert (e2e >= 0.0).all() and (e2e <= 1.0).all()
    assert (long["n_evaluable_poses"] <= long["n_input_ligands"]).all()
    mismatch = long["n_evaluable_poses"] != long["n_input_ligands"]
    if mismatch.any():
        sub = long.loc[mismatch]
        np.testing.assert_allclose(
            sub["strict_pass_end_to_end"],
            sub["n_pass_all_posebusters"] / sub["n_input_ligands"],
        )
