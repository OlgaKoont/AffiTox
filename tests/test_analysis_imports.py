"""Import smoke tests so CI fails if the analysis pipeline cannot start."""

from __future__ import annotations


def test_plot_style_imports() -> None:
    from analysis.constants import DEFAULT_METHOD_COLORS, METHOD_COLOR_CORRELATION_ANCHORS
    from analysis.plots.style import apply_style

    assert callable(apply_style)
    assert "gnina" in DEFAULT_METHOD_COLORS
    assert "qvina" in METHOD_COLOR_CORRELATION_ANCHORS
