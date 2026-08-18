"""Shared constants: method labels, primary score columns, sign flips."""

from __future__ import annotations

# Internal method id -> article display label
METHOD_LABELS: dict[str, str] = {
    "boltz2": "Boltz-2",
    "dynamicbind": "DynamicBind",
    "gnina": "GNINA 1.3",
    "plapt": "PLAPT",
    "qvina": "QVina2",
}

# Primary affinity column per method (article convention)
PRIMARY_METRICS: dict[str, str] = {
    "boltz2": "boltz2_affinity_pred_value",
    "dynamicbind": "dynamicbind_affinity_bestpose",
    "gnina": "gnina_cnn_affinity_bestpose",
    "plapt": "plapt_affinity",
    "qvina": "qvina_affinity_bestpose",
}

# Methods excluded from PoseBusters (no 3D poses)
POSEBUSTERS_EXCLUDE: set[str] = {"plapt"}

# PoseBusters / scatter line colors: Pearson-r anchors on correlation heatmap scale.
METHOD_COLOR_CORRELATION_ANCHORS: dict[str, float] = {
    "boltz2": 1.0,
    "dynamicbind": 0.5,
    "gnina": -0.35,
    "qvina": -1.0,
    "plapt": 0.0,
}

# Fallback hex (r = 1, 0.5, -0.35, -1, 0); overridden dynamically in method_color().
DEFAULT_METHOD_COLORS: dict[str, str] = {
    "boltz2": "#E5885F",
    "dynamicbind": "#FDBBAA",
    "gnina": "#AFC6F2",
    "qvina": "#6584E1",
    "plapt": "#DEDAD7",
}

# PoseBusters CSV method keys (may differ from internal id)
POSEBUSTERS_METHOD_KEYS: dict[str, str] = {
    "boltz2": "boltz2",
    "dynamicbind": "dynamicbind",
    "gnina": "gnina",
    "qvina": "qvina",
}

SIGN_FLIP_METRICS: set[str] = {
    "qvina_affinity_min",
    "qvina_affinity_bestpose",
    "qvina_max_affinity",
    "qvina_lb_affinity",
    "qvina_ub_affinity",
    "gnina_affinity_min",
    "gnina_affinity_bestpose",
    "boltz2_affinity_pred_value",
    "boltz2_affinity_pred_value1",
    "boltz2_affinity_pred_value2",
}

TOP_FRACS: dict[str, float] = {"1": 0.01, "5": 0.05, "10": 0.10}

ACTIVITY_SETS: tuple[str, ...] = ("high_active", "active", "low")
