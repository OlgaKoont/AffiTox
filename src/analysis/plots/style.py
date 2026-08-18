"""Figure styling: palette, DPI, save helpers, heatmap layout with summary row/column."""

from __future__ import annotations

from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

from ..config import AnalysisConfig
from ..constants import DEFAULT_METHOD_COLORS, METHOD_COLOR_CORRELATION_ANCHORS

ALL_TARGETS_LABEL = "All targets"
# Vertical/horizontal gap between summary (All targets) and per-target blocks in GridSpec.
GAP_RATIO = 0.025  # was 0.10; 4× tighter summary-to-body spacing
SUMMARY_BODY_HSPACE = GAP_RATIO * 5  # combined multi-panel heatmaps (e.g. correlations_combined)

# Article diverging palette: blue (negative) → gray/white (0) → salmon (positive).
# Positions are in colormap space [0, 1] matching data range [-1, +1] with center=0.
CORRELATION_DIVERGING_STOPS: list[tuple[float, str]] = [
    (0.0, "#6584E1"),    # -1.0  saturated blue
    (0.25, "#A2BFFF"),   # -0.5  periwinkle
    (0.375, "#B7CAEA"),  # -0.25 light blue
    (0.5, "#DEDAD7"),    #  0.0  light gray / near-white
    (0.625, "#F0D1BF"),  # +0.25 peach
    (0.75, "#FDBBAB"),   # +0.5  light salmon
    (0.875, "#F18C6E"),  # +0.75 salmon-orange
    (1.0, "#E5885F"),    # +1.0  saturated salmon
]

# Sequential variant for bounded positive metrics (nEF 0–1, pass rates 0–100).
POSITIVE_SEQUENTIAL_STOPS: list[tuple[float, str]] = [
    (0.0, "#6584E1"),
    (0.35, "#B7CAEA"),
    (0.55, "#DEDAD7"),
    (0.75, "#FDBBAB"),
    (1.0, "#E5885F"),
]


def round_half_up(value: float, digits: int = 2) -> float:
    quant = Decimal("1").scaleb(-digits)
    return float(Decimal(str(value)).quantize(quant, rounding=ROUND_HALF_UP))


def correlation_diverging_cmap(
    cfg: AnalysisConfig | None = None, name: str = "article_correlation"
) -> LinearSegmentedColormap:
    """Salmon (positive) — gray (0) — blue (negative), article heatmap palette."""
    _ = cfg  # reserved for future per-run overrides
    return LinearSegmentedColormap.from_list(name, CORRELATION_DIVERGING_STOPS, N=256)


def correlation_value_color(value: float) -> str:
    """Hex color from the article correlation heatmap scale at Pearson r = value."""
    cmap = correlation_diverging_cmap()
    t = (float(value) + 1.0) / 2.0
    t = min(max(t, 0.0), 1.0)
    return mcolors.to_hex(cmap(t))


def positive_sequential_cmap(
    cfg: AnalysisConfig | None = None, name: str = "article_positive"
) -> LinearSegmentedColormap:
    """Low values blue → high values salmon (nEF, PoseBusters pass rates)."""
    _ = cfg
    return LinearSegmentedColormap.from_list(name, POSITIVE_SEQUENTIAL_STOPS, N=256)


def diverging_cmap(cfg: AnalysisConfig, name: str = "article_diverging") -> LinearSegmentedColormap:
    return correlation_diverging_cmap(cfg, name=name)


def sequential_cmap(cfg: AnalysisConfig, name: str = "article_sequential") -> LinearSegmentedColormap:
    return positive_sequential_cmap(cfg, name=name)


def apply_style(cfg: AnalysisConfig) -> None:
    plt.rcParams.update(
        {
            "font.family": cfg.font_family,
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "figure.dpi": cfg.figure_dpi,
            "savefig.dpi": cfg.figure_dpi,
            "svg.fonttype": "none",
        }
    )


def save_figure(fig: plt.Figure, path_base: Path) -> None:
    path_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_base.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(path_base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def method_color(cfg: AnalysisConfig, method_id: str, default: str = "#555555") -> str:
    if method_id in cfg.method_colors:
        return cfg.method_colors[method_id]
    if method_id in METHOD_COLOR_CORRELATION_ANCHORS:
        return correlation_value_color(METHOD_COLOR_CORRELATION_ANCHORS[method_id])
    return DEFAULT_METHOD_COLORS.get(method_id, default)


def _format_annot(mat: pd.DataFrame, fmt: str) -> np.ndarray:
    annot = np.empty(mat.shape, dtype=object)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if pd.isna(val):
                annot[i, j] = ""
            elif fmt == ".2f":
                annot[i, j] = f"{round_half_up(float(val), 2):.2f}"
            elif fmt == ".1f":
                annot[i, j] = f"{round_half_up(float(val), 1):.1f}"
            else:
                annot[i, j] = f"{float(val):{fmt}}"
    return annot


def plot_heatmap_with_summary_row(
    summary_row: pd.Series,
    body: pd.DataFrame,
    *,
    title: str,
    cbar_label: str,
    cmap: LinearSegmentedColormap,
    vmin: float | None,
    vmax: float | None,
    center: float | None,
    fmt: str,
    path_base: Path,
    figsize_scale: tuple[float, float] = (0.55, 0.35),
    summary_label: str | None = None,
) -> None:
    """Heatmap: one summary row + gap + per-target rows."""
    summary_df = summary_row.to_frame().T
    summary_df.index = [summary_label or ALL_TARGETS_LABEL]
    n_cols = len(body.columns)
    n_body = len(body)

    fig_w = max(6.5, figsize_scale[0] * n_cols)
    fig_h = max(4.5, figsize_scale[1] * (1 + n_body) + 0.8)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(
        2, 2,
        figure=fig,
        height_ratios=[1.0, max(n_body, 1)],
        width_ratios=[1.0, 0.045],
        hspace=GAP_RATIO,
        wspace=0.05,
    )
    ax_top = fig.add_subplot(gs[0, 0])
    ax_bot = fig.add_subplot(gs[1, 0])
    cbar_ax = fig.add_subplot(gs[:, 1])

    heatmap_kw = dict(
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=center,
        linewidths=0.5,
        linecolor="white",
        annot_kws={"size": 9, "color": "black"},
        cbar_ax=cbar_ax,
        cbar_kws={"label": cbar_label},
    )

    sns.heatmap(
        summary_df,
        ax=ax_top,
        annot=_format_annot(summary_df, fmt),
        fmt="",
        cbar=False,
        **heatmap_kw,
    )
    sns.heatmap(
        body,
        ax=ax_bot,
        annot=_format_annot(body, fmt) if fmt else None,
        fmt="" if fmt else "",
        **heatmap_kw,
    )

    ax_top.set_ylabel("")
    ax_top.set_xlabel("")
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_bot.set_xlabel("Docking method")
    ax_bot.set_ylabel("Target")
    ax_top.tick_params(axis="y", rotation=0)
    ax_bot.tick_params(axis="x", rotation=0)
    ax_bot.tick_params(axis="y", rotation=0)
    fig.suptitle(title, y=1.02, fontsize=11)
    save_figure(fig, path_base)


def plot_heatmap_with_summary_col(
    summary_col: pd.Series,
    body: pd.DataFrame,
    *,
    title: str,
    cbar_label: str,
    cmap: LinearSegmentedColormap,
    vmin: float,
    vmax: float,
    fmt: str,
    path_base: Path,
    figsize_scale: tuple[float, float] = (0.45, 0.22),
    figsize: tuple[float, float] | None = None,
    summary_panel_ratio: float = 1.35,
    label_panel_ratio: float = 2.2,
    label_summary_wspace: float = 0.02,
    cbar_width_ratio: float = 0.045,
    summary_body_gap_scale: float = 1.0,
    summary_body_wspace: float | None = None,
    equal_summary_width: bool = False,
    body_annot: bool = False,
    body_fmt: str | None = None,
    y_label: str = "Check",
    y_tick_rotation: int = 0,
    center: float | None = None,
    summary_label: str | None = None,
) -> None:
    """Heatmap: one summary column + gap + per-target columns."""
    summary_df = summary_col.to_frame()
    summary_df.columns = [summary_label or ALL_TARGETS_LABEL]
    n_rows = len(body.index)
    n_cols = len(body.columns)

    if figsize is not None:
        fig_w, fig_h = figsize
    else:
        fig_w = max(8.0, figsize_scale[0] * (1 + n_cols) + 1.2)
        fig_h = max(6.0, figsize_scale[1] * n_rows)
    fig = plt.figure(figsize=(fig_w, fig_h))

    if equal_summary_width:
        # [method labels | 1-col summary | gap | n-col body | colorbar]
        gap_wspace = 0.23 if summary_body_wspace is None else summary_body_wspace
        gs = GridSpec(
            1, 2,
            figure=fig,
            width_ratios=[label_panel_ratio + 1.0, n_cols + 0.045],
            wspace=gap_wspace,
        )
        gs_left = gs[0].subgridspec(
            1, 2, width_ratios=[label_panel_ratio, 1.0], wspace=label_summary_wspace
        )
        gs_right = gs[1].subgridspec(
            1, 2, width_ratios=[n_cols, cbar_width_ratio], wspace=0.05
        )
        ax_left = fig.add_subplot(gs_left[0, 1])
        ax_label = fig.add_subplot(gs_left[0, 0])
        ax_mid = fig.add_subplot(gs_right[0, 0])
        cbar_ax = fig.add_subplot(gs_right[0, 1])
    else:
        gs = GridSpec(
            1, 2,
            figure=fig,
            width_ratios=[summary_panel_ratio, max(n_cols, 1)],
            wspace=GAP_RATIO if summary_body_wspace is None else summary_body_wspace,
        )
        gs_inner = gs[1].subgridspec(1, 2, width_ratios=[1.0, 0.045], wspace=0.05)
        ax_label = None
        ax_left = fig.add_subplot(gs[0, 0])
        ax_mid = fig.add_subplot(gs_inner[0, 0])
        cbar_ax = fig.add_subplot(gs_inner[0, 1])

    heatmap_kw = dict(
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=center,
        linewidths=0.3,
        linecolor="white",
        cbar_ax=cbar_ax,
        cbar_kws={"label": cbar_label},
        annot_kws={"size": 8, "color": "black"},
    )

    sns.heatmap(
        summary_df,
        ax=ax_left,
        annot=_format_annot(summary_df, fmt) if fmt else None,
        fmt="" if fmt else "",
        cbar=False,
        yticklabels=(ax_label is None),
        **heatmap_kw,
    )
    sns.heatmap(
        body,
        ax=ax_mid,
        annot=_format_annot(body, body_fmt or fmt) if body_annot and (body_fmt or fmt) else False,
        fmt="" if body_annot else "",
        yticklabels=False,
        **heatmap_kw,
    )

    if ax_label is not None:
        tick_pos = np.arange(len(body.index)) + 0.5
        ax_label.set_ylim(len(body.index), 0)
        ax_label.set_yticks(tick_pos)
        ax_label.set_yticklabels(list(body.index))
        ax_label.yaxis.tick_right()
        ax_label.tick_params(
            axis="y",
            rotation=y_tick_rotation,
            labelsize=9,
            pad=2,
            length=0,
            labelleft=False,
            labelright=True,
            right=True,
            left=False,
        )
        plt.setp(
            ax_label.get_yticklabels(),
            ha="right" if y_tick_rotation == 0 else "center",
            rotation=y_tick_rotation,
            visible=True,
        )
        ax_label.set_xticks([])
        ax_label.set_xlabel("")
        for spine in ax_label.spines.values():
            spine.set_visible(False)
        ax_label.set_ylabel(y_label, rotation=90, labelpad=28)
        ax_left.set_yticklabels([])
        ax_left.tick_params(axis="y", left=False, labelleft=False)
    else:
        ax_left.set_ylabel(y_label, rotation=90, labelpad=28)
        ax_left.set_xlabel("")
        ax_left.tick_params(axis="y", rotation=y_tick_rotation, labelsize=9, pad=4)
        plt.setp(
            ax_left.get_yticklabels(),
            ha="right" if y_tick_rotation == 0 else "center",
            rotation=y_tick_rotation,
        )

    ax_mid.set_ylabel("")
    ax_mid.set_xlabel("Target")
    ax_mid.tick_params(axis="x", rotation=0)
    plt.setp(ax_mid.get_xticklabels(), ha="center")
    ax_mid.tick_params(axis="y", left=False, labelleft=False)

    if equal_summary_width and ax_label is not None:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        sum_pos = ax_left.get_position()
        body_pos = ax_mid.get_position()
        cb_pos = cbar_ax.get_position()
        gap = body_pos.x0 - sum_pos.x1
        if gap > 0 and summary_body_gap_scale != 1.0:
            shift = gap * (1.0 - summary_body_gap_scale)
            ax_left.set_position(
                [sum_pos.x0 + shift, sum_pos.y0, sum_pos.width, sum_pos.height]
            )
            sum_pos = ax_left.get_position()
            body_pos = ax_mid.get_position()
            cb_pos = cbar_ax.get_position()
        label_bboxes = [
            t.get_window_extent(renderer).transformed(fig.transFigure.inverted())
            for t in ax_label.get_yticklabels()
        ]
        max_label_w = max((bb.width for bb in label_bboxes), default=0.08)
        label_pad = 0.008
        lab_pos = ax_label.get_position()
        ax_label.set_position(
            [
                sum_pos.x0 - label_pad - max_label_w,
                lab_pos.y0,
                max_label_w,
                lab_pos.height,
            ]
        )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        sum_pos = ax_left.get_position()
        body_pos = ax_mid.get_position()
        cb_pos = cbar_ax.get_position()
        col_w = body_pos.width / max(n_cols, 1)
        if ax_mid.texts:
            text_bbox = ax_mid.texts[0].get_window_extent(renderer).transformed(
                fig.transFigure.inverted()
            )
            val_len = max(len(ax_mid.texts[0].get_text().strip("*")), 1)
            two_digit_w = 2.0 * text_bbox.width / val_len
        else:
            two_digit_w = col_w * 0.45
        two_digit_w = max(two_digit_w, col_w * 0.35)
        cbar_ax.set_position(
            [cb_pos.x1 - two_digit_w, cb_pos.y0, two_digit_w, cb_pos.height]
        )

    fig.suptitle(title, y=1.02, fontsize=11)
    save_figure(fig, path_base)


def plot_heatmap_with_summary_row_and_col(
    summary_row: pd.Series,
    summary_col: pd.Series,
    body: pd.DataFrame,
    *,
    title: str,
    cbar_label: str,
    cmap: LinearSegmentedColormap,
    vmin: float,
    vmax: float,
    fmt: str,
    path_base: Path,
    figsize_scale: tuple[float, float] = (0.90, 0.44),
    x_label: str = "Target",
    y_label: str = "Tests",
) -> None:
    """Heatmap with both All tests row and All targets column summaries."""
    top_left = pd.DataFrame(
        [[float(body.stack().mean()) if not body.stack().empty else np.nan]],
        index=["All"],
        columns=["All"],
    )
    top_row = summary_row.to_frame().T
    top_row.index = ["All tests"]
    left_col = summary_col.to_frame()
    left_col.columns = [ALL_TARGETS_LABEL]

    n_rows = len(body.index)
    n_cols = len(body.columns)
    fig_w = max(10.0, figsize_scale[0] * (1 + n_cols) + 2.0)
    fig_h = max(6.8, figsize_scale[1] * (1 + n_rows) + 1.0)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(
        2,
        3,
        figure=fig,
        height_ratios=[1.0, max(n_rows, 1)],
        width_ratios=[1.35, max(n_cols, 1), 0.045],
        hspace=GAP_RATIO,
        wspace=GAP_RATIO,
    )
    ax_corner = fig.add_subplot(gs[0, 0])
    ax_top = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax_body = fig.add_subplot(gs[1, 1])
    cbar_ax = fig.add_subplot(gs[:, 2])

    heatmap_kw = dict(
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.3,
        linecolor="white",
        annot_kws={"size": 8, "color": "black"},
    )

    sns.heatmap(
        top_left,
        ax=ax_corner,
        annot=_format_annot(top_left, fmt),
        fmt="",
        cbar=False,
        **heatmap_kw,
    )
    sns.heatmap(
        top_row,
        ax=ax_top,
        annot=_format_annot(top_row, fmt),
        fmt="",
        cbar=False,
        **heatmap_kw,
    )
    sns.heatmap(
        left_col,
        ax=ax_left,
        annot=_format_annot(left_col, fmt),
        fmt="",
        cbar=False,
        **heatmap_kw,
    )
    sns.heatmap(
        body,
        ax=ax_body,
        annot=_format_annot(body, fmt),
        fmt="",
        cbar_ax=cbar_ax,
        cbar_kws={"label": cbar_label},
        **heatmap_kw,
    )

    ax_corner.set_xlabel("")
    ax_corner.set_ylabel("")
    ax_corner.tick_params(axis="x", rotation=0)
    ax_corner.tick_params(axis="y", rotation=0)

    ax_top.set_xlabel("")
    ax_top.set_ylabel("")
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_top.tick_params(axis="y", rotation=0)

    ax_left.set_xlabel("")
    ax_left.set_ylabel(y_label)
    ax_left.tick_params(axis="x", rotation=0)
    ax_left.tick_params(axis="y", labelleft=True, pad=2)
    for tick in ax_left.get_yticklabels():
        tick.set_fontsize(8)

    ax_body.set_xlabel(x_label)
    ax_body.set_ylabel("")
    ax_body.tick_params(axis="x", rotation=0)
    ax_body.tick_params(axis="y", left=False, labelleft=False)

    fig.suptitle(title, y=1.02, fontsize=11)
    save_figure(fig, path_base)


def plot_combined_heatmaps_with_summary_row(
    panels: list[tuple[str, pd.Series, pd.DataFrame]],
    *,
    cbar_label: str,
    cmap: LinearSegmentedColormap,
    vmin: float | None,
    vmax: float | None,
    center: float | None,
    fmt: str,
    path_base: Path,
    x_label: str = "Docking method",
    y_label: str = "Target",
    col_width: float = 0.85,
    row_height: float = 0.36,
    summary_label: str | None = None,
) -> None:
    """Side-by-side heatmaps sharing one colorbar and one bottom x-axis label."""
    if not panels:
        return

    n_panels = len(panels)
    n_cols = len(panels[0][2].columns)
    n_body = len(panels[0][2].index)
    row_label = summary_label or ALL_TARGETS_LABEL

    panel_w = max(3.8, col_width * n_cols + 0.6)
    fig_w = panel_w * n_panels + 1.1
    fig_h = max(5.8, row_height * (1 + n_body) + 1.6)
    width_ratios = [panel_w] * n_panels + [0.045]

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(
        2,
        n_panels + 1,
        figure=fig,
        height_ratios=[1.0, max(n_body, 1)],
        width_ratios=width_ratios,
        hspace=SUMMARY_BODY_HSPACE,
        wspace=0.10,
        left=0.06,
        right=0.94,
        top=0.90,
        bottom=0.14,
    )
    cbar_ax = fig.add_subplot(gs[:, n_panels])

    heatmap_kw = dict(
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=center,
        linewidths=0.5,
        linecolor="white",
        annot_kws={"size": 8, "color": "black"},
    )

    for idx, (title, summary_row, body) in enumerate(panels):
        summary_df = summary_row.to_frame().T
        summary_df.index = [row_label]

        ax_top = fig.add_subplot(gs[0, idx])
        ax_bot = fig.add_subplot(gs[1, idx], sharex=ax_top)

        sns.heatmap(
            summary_df,
            ax=ax_top,
            annot=_format_annot(summary_df, fmt),
            fmt="",
            cbar=False,
            **heatmap_kw,
        )
        sns.heatmap(
            body,
            ax=ax_bot,
            annot=_format_annot(body, fmt),
            fmt="",
            cbar=(idx == n_panels - 1),
            cbar_ax=cbar_ax if idx == n_panels - 1 else None,
            cbar_kws={"label": cbar_label},
            **heatmap_kw,
        )

        ax_top.set_title(title, fontsize=11, pad=8)
        ax_top.set_xlabel("")
        ax_top.set_ylabel("")
        ax_top.tick_params(axis="x", labelbottom=False)
        ax_top.tick_params(axis="y", rotation=0)

        ax_bot.set_xlabel("")
        ax_bot.tick_params(axis="x", rotation=0)
        ax_bot.tick_params(axis="y", rotation=0)
        if idx == 0:
            ax_top.set_ylabel("")
            ax_bot.set_ylabel(y_label)
        else:
            ax_top.set_ylabel("")
            ax_bot.set_ylabel("")
            ax_top.tick_params(axis="y", left=False, labelleft=False)
            ax_bot.tick_params(axis="y", left=False, labelleft=False)

    fig.supxlabel(x_label, fontsize=10, y=0.01)
    save_figure(fig, path_base)
