"""Scientific invariants for AffiTox metrics (does not rerun docking)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from analysis.constants import POSEBUSTERS_EXCLUDE, SIGN_FLIP_METRICS
from analysis.data import ensure_pki, prepare_xy
from analysis.enrichment import compute_ef
from analysis.inferential import holm


def test_pki_from_ki_nm():
    df = pd.DataFrame({"standard_value": [10.0, 1.0, 1000.0]})
    pki = ensure_pki(df, exp_col="missing")
    np.testing.assert_allclose(pki.to_numpy(), [8.0, 9.0, 6.0])


def test_nef_precision_when_a_ge_t():
    # N=10, f=0.10 → T=1; A=5 ≥ T → nEF = H/T (precision)
    scores = np.linspace(0.0, 1.0, 10)
    is_active = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1], dtype=float)
    _ef, n_ef, n, a, h_top = compute_ef(is_active, scores, 0.10, "desc")
    assert n == 10 and a == 5 and h_top == 1
    np.testing.assert_allclose(n_ef, 1.0)


def test_nef_recall_when_a_lt_t():
    # N=10, f=0.50 → T=5; A=2 < T → nEF = H/A (recall)
    scores = np.array([0, 1, 2, 3, 4, 5, 6, 7, 10, 11], dtype=float)
    is_active = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1], dtype=float)
    _ef, n_ef, n, a, h_top = compute_ef(is_active, scores, 0.50, "desc")
    assert n == 10 and a == 2 and h_top == 2
    np.testing.assert_allclose(n_ef, 1.0)


def test_nef_partial_recall():
    scores = np.array([10, 9, 8, 7, 6, 5, 4, 3, 2, 1], dtype=float)
    is_active = np.array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype=float)
    _ef, n_ef, n, a, h_top = compute_ef(is_active, scores, 0.50, "desc")
    assert a == 2 and h_top == 1
    np.testing.assert_allclose(n_ef, 0.5)


def test_holm_monotonic_and_capped():
    p = np.array([0.01, 0.04, 0.03])
    adj = holm(p)
    np.testing.assert_allclose(adj, [0.03, 0.06, 0.06])
    assert np.all(adj >= p)
    assert np.all(adj <= 1.0)


def test_sign_flip_qvina_not_gnina_cnn():
    assert "qvina_affinity_bestpose" in SIGN_FLIP_METRICS
    assert "boltz2_affinity_pred_value" in SIGN_FLIP_METRICS
    assert "gnina_cnn_affinity_bestpose" not in SIGN_FLIP_METRICS
    df = pd.DataFrame(
        {
            "pValue": [6.0, 7.0],
            "qvina_affinity_bestpose": [1.0, 2.0],
            "gnina_cnn_affinity_bestpose": [1.0, 2.0],
        }
    )
    _, y_q = prepare_xy(df, "qvina_affinity_bestpose")
    _, y_g = prepare_xy(df, "gnina_cnn_affinity_bestpose")
    np.testing.assert_allclose(y_q, [-1.0, -2.0])
    np.testing.assert_allclose(y_g, [1.0, 2.0])


def test_plapt_excluded_from_posebusters():
    assert POSEBUSTERS_EXCLUDE == {"plapt"}


def test_posebusters_pairwise_without_count_columns():
    from analysis.inferential import run_posebusters_pairwise_tests

    summary = pd.DataFrame(
        {
            "target": ["1g5m", "1g5m"],
            "method_id": ["gnina", "qvina"],
            "method_label": ["GNINA 1.3", "QVina2"],
            "n_poses": [10, 10],
            "pass_rate_all": [0.5, 0.4],
            "pass_rate_90pct": [0.8, 0.7],
            "pass_rate_50pct": [0.9, 0.85],
        }
    )
    result = run_posebusters_pairwise_tests(summary)
    assert not result.empty
    assert result["p_raw"].notna().any()
