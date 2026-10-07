"""PoseBusters heatmaps with All targets summary row/column."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

from ..config import AnalysisConfig
from ..constants import POSEBUSTERS_EXCLUDE
from ..posebusters import collect_pooled_pass_counts, pass_count_distribution
from .style import (
    ALL_TARGETS_LABEL,
    SUMMARY_BODY_HSPACE,
    apply_style,
    method_color,
    plot_combined_heatmaps_with_summary_row,
    plot_heatmap_with_summary_row,
    plot_heatmap_with_summary_row_and_col,
    save_figure,
    sequential_cmap,
)
from .style import _format_annot  # shared heatmap cell labels

ALL_TESTS_HEATMAP_VMIN = 76.0  # 40–100 base; color at 75 → now at 90 (vmax=100)
ALL_TESTS_HEATMAP_VMAX = 100.0


def _style_heatmap_text_contrast(
    ax: plt.Axes, vmin: float, vmax: float, *, fontsize: float = 8
) -> None:
    """White/dark annotation text tuned for the 76–100% pass-rate palette."""
    for text in ax.texts:
        raw = text.get_text().strip()
        if not raw:
            continue
        try:
            val = float(raw)
        except ValueError:
            continue
        t = (val - vmin) / (vmax - vmin) if vmax > vmin else 0.5
        text.set_color("#ffffff" if t >= 0.50 else "#222222")
        text.set_fontsize(fontsize)


def _add_panel_letter(ax: plt.Axes, letter: str, fontsize: float) -> None:
    """Place A/B the same number of points above the axes top (not axes-fraction)."""
    ax.annotate(
        letter,
        xy=(0.0, 1.0),
        xycoords=ax.transAxes,
        xytext=(-2, 6),
        textcoords="offset points",
        fontsize=fontsize,
        fontweight="bold",
        ha="right",
        va="bottom",
        clip_on=False,
        annotation_clip=False,
    )


def _draw_all_tests_heatmap_panel(
    fig: plt.Figure,
    gs_slot,
    *,
    summary_df: pd.DataFrame,
    body: pd.DataFrame,
    n_body: int,
    cmap,
    title: str,
    adaptive_text: bool,
    vmin: float = ALL_TESTS_HEATMAP_VMIN,
    vmax: float = ALL_TESTS_HEATMAP_VMAX,
    cbar_ticks: list[float] | None = None,
    panel_letter: str | None = None,
    title_fontsize: float | None = None,
    annot_size: float = 8,
) -> plt.Axes:
    inner = gs_slot.subgridspec(
        2,
        2,
        height_ratios=[1.0, max(n_body, 1)],
        width_ratios=[1.0, 0.052],
        hspace=SUMMARY_BODY_HSPACE / 4,
        wspace=0.05,
    )
    ax_top = fig.add_subplot(inner[0, 0])
    ax_bot = fig.add_subplot(inner[1, 0], sharex=ax_top)
    cbar_ax = fig.add_subplot(inner[:, 1])

    if cbar_ticks is None:
        cbar_ticks = [76, 80, 85, 90, 95, 100]
    heatmap_kw = dict(
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        annot_kws={"size": annot_size, "color": "black"},
    )
    sns.heatmap(
        summary_df,
        ax=ax_top,
        annot=_format_annot(summary_df, ".1f"),
        fmt="",
        cbar=False,
        **heatmap_kw,
    )
    sns.heatmap(
        body,
        ax=ax_bot,
        annot=_format_annot(body, ".1f"),
        fmt="",
        cbar_ax=cbar_ax,
        cbar_kws={"label": "Pass rate (%)", "ticks": cbar_ticks},
        **heatmap_kw,
    )
    if adaptive_text:
        _style_heatmap_text_contrast(ax_top, vmin, vmax, fontsize=annot_size)
        _style_heatmap_text_contrast(ax_bot, vmin, vmax, fontsize=annot_size)

    title_fs = (
        plt.rcParams["axes.titlesize"] if title_fontsize is None else title_fontsize
    )
    if title:
        ax_top.set_title(title, fontsize=title_fs, pad=6)
    ax_top.set_xlabel("")
    ax_top.set_ylabel("")
    ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_top.tick_params(axis="y", rotation=0)

    ax_bot.set_xlabel("")
    ax_bot.set_ylabel("Target")
    ax_bot.tick_params(axis="x", rotation=0)
    ax_bot.tick_params(axis="y", rotation=0, pad=4)
    ax_top.tick_params(axis="y", pad=4)
    cbar_ax.tick_params(labelsize=8)
    if panel_letter:
        _add_panel_letter(ax_top, panel_letter, title_fs)
    return ax_top


def _prepare_at_least_all_tests_panels(
    summary: pd.DataFrame,
    checks: pd.DataFrame,
    cfg: AnalysisConfig,
) -> tuple[int, dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame, int] | None:
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]
    n_tests, pooled = collect_pooled_pass_counts(cfg)
    if n_tests == 0 or not any(pooled.values()):
        return None

    at_least: dict[str, pd.DataFrame] = {}
    for method_id, counts in pooled.items():
        if counts:
            at_least[method_id] = pass_count_distribution(counts, n_tests, mode="at_least")

    summary_row = _all_tests_summary_row(checks, summary, cfg, methods)
    body = _all_tests_body(checks, cfg, methods)
    if body.isna().all().all():
        return None

    summary_df = summary_row.to_frame().T
    summary_df.index = [ALL_TARGETS_LABEL]
    return n_tests, at_least, summary_df, body.astype(float), len(body)

def _pretty_check(name: str) -> str:
    return (
        name.replace("minimum_distance_to_", "min_dist_")
        .replace("volume_overlap_with_", "overlap_")
        .replace("protein-ligand_", "prot_lig_")
    )


def _pass_rate_body(summary: pd.DataFrame, rate_col: str, cfg: AnalysisConfig, methods: list[str]) -> pd.DataFrame:
    labels = [cfg.method_label(m) for m in methods]
    idx = [cfg.target_label(t) for t in cfg.targets]
    body = pd.DataFrame(index=idx, columns=labels, dtype=float)
    for _, row in summary.iterrows():
        if row["method_id"] not in methods:
            continue
        body.loc[
            cfg.target_label(str(row["target"])),
            cfg.method_label(row["method_id"]),
        ] = 100.0 * row[rate_col]
    return body


def _weighted_summary_row(summary: pd.DataFrame, rate_col: str, cfg: AnalysisConfig, methods: list[str]) -> pd.Series:
    sub = summary[summary["method_id"].isin(methods)].copy()
    sub["weighted"] = sub[rate_col] * sub["n_poses"]
    summary_row = pd.Series(dtype=float)
    for method_id in methods:
        label = cfg.method_label(method_id)
        part = sub[sub["method_id"] == method_id]
        if part.empty or part["n_poses"].sum() == 0:
            summary_row[label] = np.nan
        else:
            summary_row[label] = 100.0 * part["weighted"].sum() / part["n_poses"].sum()
    return summary_row


def _all_tests_body(checks: pd.DataFrame, cfg: AnalysisConfig, methods: list[str]) -> pd.DataFrame:
    """Mean per-check pass rate by target and method (All tests row from per-check heatmaps)."""
    labels = [cfg.method_label(m) for m in methods]
    idx = [cfg.target_label(t) for t in cfg.targets]
    body = pd.DataFrame(index=idx, columns=labels, dtype=float)
    for method_id in methods:
        sub = checks[checks.method_id == method_id]
        if sub.empty:
            continue
        label = cfg.method_label(method_id)
        for target in cfg.targets:
            part = sub[sub.target == target.lower()]
            if part.empty:
                continue
            body.loc[cfg.target_label(target), label] = 100.0 * part["pass_rate"].mean()
    return body


def _all_tests_summary_row(checks: pd.DataFrame, summary: pd.DataFrame, cfg: AnalysisConfig, methods: list[str]) -> pd.Series:
    if "mean_frac_pass" in summary.columns:
        sub = summary[summary["method_id"].isin(methods)].copy()
        summary_row = pd.Series(dtype=float)
        for method_id in methods:
            label = cfg.method_label(method_id)
            part = sub[sub["method_id"] == method_id]
            if part.empty or part["n_poses"].sum() == 0:
                summary_row[label] = np.nan
            else:
                summary_row[label] = 100.0 * (
                    part["mean_frac_pass"] * part["n_poses"]
                ).sum() / part["n_poses"].sum()
        return summary_row

    body = _all_tests_body(checks, cfg, methods)
    return body.mean(axis=0, skipna=True)


def plot_posebusters_summary(summary: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters" / "heatmaps"
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]

    for rate_col, title in (
        ("pass_rate_all", "Pass all checks (100%)"),
        ("pass_rate_90pct", "Pass >=90% checks"),
        ("pass_rate_50pct", "Pass >=50% checks"),
    ):
        body = _pass_rate_body(summary, rate_col, cfg, methods)
        if body.isna().all().all():
            continue

        summary_row = _weighted_summary_row(summary, rate_col, cfg, methods)
        plot_heatmap_with_summary_row(
            summary_row,
            body.astype(float),
            title=title,
            cbar_label="Pass rate (%)",
            cmap=cmap,
            vmin=0.0,
            vmax=100.0,
            center=None,
            fmt=".1f",
            path_base=out_dir / f"heatmap_{rate_col}",
        )


def plot_posebusters_combined(summary: pd.DataFrame, checks: pd.DataFrame, cfg: AnalysisConfig) -> None:
    """Three-panel figure: pass all, pass >=90%, and mean All-tests pass rate."""
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters" / "heatmaps"
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]

    panels: list[tuple[str, pd.Series, pd.DataFrame]] = []
    for rate_col, title in (
        ("pass_rate_all", "Pass all checks (100%)"),
        ("pass_rate_90pct", "Pass >=90% checks"),
    ):
        body = _pass_rate_body(summary, rate_col, cfg, methods)
        if body.isna().all().all():
            continue
        panels.append((title, _weighted_summary_row(summary, rate_col, cfg, methods), body.astype(float)))

    all_tests_body = _all_tests_body(checks, cfg, methods)
    if not all_tests_body.isna().all().all():
        panels.append(
            (
                "All tests (mean check pass rate)",
                _all_tests_summary_row(checks, summary, cfg, methods),
                all_tests_body.astype(float),
            )
        )

    if len(panels) == 3:
        plot_combined_heatmaps_with_summary_row(
            panels,
            cbar_label="Pass rate (%)",
            cmap=cmap,
            vmin=0.0,
            vmax=100.0,
            center=None,
            fmt=".1f",
            path_base=out_dir / "heatmap_posebusters_combined",
            col_width=0.95,
            row_height=0.36,
        )


def plot_posebusters_per_check(checks: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters" / "heatmaps"
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]

    for method_id in methods:
        sub = checks[checks.method_id == method_id].copy()
        if sub.empty:
            continue
        label = cfg.method_label(method_id)
        body = sub.pivot(index="check", columns="target", values="pass_rate")
        body.index = [_pretty_check(c) for c in body.index]
        body.columns = [cfg.target_label(str(c)) for c in body.columns]
        body = 100.0 * body

        summary_col = body.mean(axis=1, skipna=True)
        summary_row = body.mean(axis=0, skipna=True)

        plot_heatmap_with_summary_row_and_col(
            summary_row,
            summary_col,
            body.astype(float),
            title=f"PoseBusters per-check pass rate — {label}",
            cbar_label="Pass rate (%)",
            cmap=cmap,
            vmin=0.0,
            vmax=100.0,
            fmt=".1f",
            path_base=out_dir / f"heatmap_checks_{method_id}",
            figsize_scale=(0.90, 0.44),
            y_label="Tests",
        )


def _draw_pass_count_curve_on_ax(
    ax: plt.Axes,
    curves: dict[str, pd.DataFrame],
    cfg: AnalysisConfig,
    n_tests: int,
    *,
    x_label: str,
    title: str,
    panel_letter: str | None = None,
    title_fontsize: float | None = None,
) -> None:
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]
    for method_id in methods:
        df = curves.get(method_id)
        if df is None or df.empty:
            continue
        ax.plot(
            df["n_tests_passed"],
            df["pct_ligands"],
            label=cfg.method_label(method_id),
            color=method_color(cfg, method_id),
            linestyle="-",
            linewidth=2.0,
            marker="o",
            markersize=3.5,
        )

    ax.set_xlim(0.5, n_tests + 0.5)
    ax.set_ylim(0, 105)
    ax.set_xticks(range(1, n_tests + 1))
    ax.set_yticks(range(0, 101, 10))
    title_fs = (
        plt.rcParams["axes.titlesize"] if title_fontsize is None else title_fontsize
    )
    ax.set_xlabel(x_label)
    ax.set_ylabel("% ligands")
    if title:
        ax.set_title(title, fontsize=title_fs, pad=6)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(loc="upper left", frameon=True, framealpha=0.9, fontsize=9)
    if panel_letter:
        _add_panel_letter(ax, panel_letter, title_fs)


def _plot_pass_count_curve(
    curves: dict[str, pd.DataFrame],
    *,
    cfg: AnalysisConfig,
    n_tests: int,
    title: str,
    x_label: str,
    path_base,
) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    _draw_pass_count_curve_on_ax(ax, curves, cfg, n_tests, x_label=x_label, title=title)
    save_figure(fig, path_base)


def plot_posebusters_at_least_all_tests_combined(
    summary: pd.DataFrame,
    checks: pd.DataFrame,
    cfg: AnalysisConfig,
) -> None:
    """Left: at-least-N curves; right: pass-all heatmap (same data as heatmap_pass_rate_all)."""
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters"

    prepared = _prepare_at_least_all_tests_panels(summary, checks, cfg)
    if prepared is None:
        return
    n_tests, at_least, _mean_summary, _mean_body, _n_mean = prepared

    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]
    body = _pass_rate_body(summary, "pass_rate_all", cfg, methods)
    if body.isna().all().all():
        return
    summary_row = _weighted_summary_row(summary, "pass_rate_all", cfg, methods)
    summary_df = summary_row.to_frame().T
    summary_df.index = [ALL_TARGETS_LABEL]
    n_body = len(body)

    title_fs = float(plt.rcParams["axes.titlesize"]) * 1.4
    annot_fs = 8 * 1.3

    fig = plt.figure(figsize=(16.2, 6.8))
    outer = GridSpec(
        1,
        2,
        figure=fig,
        width_ratios=[1.0, 1.12],
        # 0.28 → 0.16 ≈ half a right-panel heatmap cell closer.
        wspace=0.16,
        left=0.06,
        right=0.97,
        top=0.86,
        bottom=0.08,
    )

    ax_line = fig.add_subplot(outer[0, 0])
    _draw_pass_count_curve_on_ax(
        ax_line,
        at_least,
        cfg,
        n_tests,
        x_label=f"Minimum N tests passed (of {n_tests})",
        title="",
        title_fontsize=title_fs,
    )

    ax_heat_top = _draw_all_tests_heatmap_panel(
        fig,
        outer[0, 1],
        summary_df=summary_df,
        body=body.astype(float),
        n_body=n_body,
        cmap=cmap,
        title="",
        adaptive_text=False,
        vmin=0.0,
        vmax=100.0,
        cbar_ticks=[0, 20, 40, 60, 80, 100],
        title_fontsize=title_fs,
        annot_size=annot_fs,
    )

    # A, B and both titles on one figure y (axes-fraction offset made B sit lower).
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    y_shared = max(
        inv.transform(ax.transAxes.transform((0.0, 1.0)))[1]
        for ax in (ax_line, ax_heat_top)
    )
    dy = 8.0 / 72.0 / fig.get_figheight()
    y_label = y_shared + dy
    for ax, letter, heading in (
        (ax_line, "A", "All targets: ligands passing at least N tests"),
        (ax_heat_top, "B", "Pass all checks"),
    ):
        x0, _ = inv.transform(ax.transAxes.transform((0.0, 1.0)))
        x1, _ = inv.transform(ax.transAxes.transform((1.0, 1.0)))
        fig.text(
            x0,
            y_label,
            letter,
            fontsize=title_fs,
            fontweight="bold",
            ha="right",
            va="bottom",
        )
        fig.text(
            0.5 * (x0 + x1),
            y_label,
            heading,
            fontsize=title_fs,
            ha="center",
            va="bottom",
        )

    save_figure(fig, out_dir / "pass_count_at_least_with_all_tests")


def plot_posebusters_at_least_all_tests_combined_var2(
    summary: pd.DataFrame,
    checks: pd.DataFrame,
    cfg: AnalysisConfig,
) -> None:
    """Refined layout: same structure as v1, adaptive heatmap labels, explicit colorbar ticks."""
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters"

    prepared = _prepare_at_least_all_tests_panels(summary, checks, cfg)
    if prepared is None:
        return
    n_tests, at_least, summary_df, body, n_body = prepared

    fig = plt.figure(figsize=(15.5, 6.8))
    outer = GridSpec(
        1,
        2,
        figure=fig,
        width_ratios=[1.08, 1.0],
        wspace=0.30,
        left=0.07,
        right=0.97,
        top=0.92,
        bottom=0.12,
    )

    ax_line = fig.add_subplot(outer[0, 0])
    _draw_pass_count_curve_on_ax(
        ax_line,
        at_least,
        cfg,
        n_tests,
        x_label=f"Minimum N tests passed (of {n_tests})",
        title="All targets: ligands passing at least N tests",
    )

    _draw_all_tests_heatmap_panel(
        fig,
        outer[0, 1],
        summary_df=summary_df,
        body=body,
        n_body=n_body,
        cmap=cmap,
        title="All tests (mean check pass rate)",
        adaptive_text=True,
    )

    save_figure(fig, out_dir / "pass_count_at_least_with_all_tests_var2")


def plot_posebusters_pass_count_curves(cfg: AnalysisConfig) -> None:
    """Line plots: exactly / at least N PoseBusters tests passed (All targets pooled)."""
    apply_style(cfg)
    out_dir = cfg.figures_dir / "posebusters"
    out_dir.mkdir(parents=True, exist_ok=True)

    n_tests, pooled = collect_pooled_pass_counts(cfg)
    if n_tests == 0 or not any(pooled.values()):
        return

    exactly: dict[str, pd.DataFrame] = {}
    at_least: dict[str, pd.DataFrame] = {}
    for method_id, counts in pooled.items():
        if not counts:
            continue
        exactly[method_id] = pass_count_distribution(counts, n_tests, mode="exactly")
        at_least[method_id] = pass_count_distribution(counts, n_tests, mode="at_least")

    tables_dir = cfg.tables_dir / "posebusters"
    tables_dir.mkdir(parents=True, exist_ok=True)
    for method_id in exactly:
        label = cfg.method_label(method_id)
        exactly[method_id].assign(method_id=method_id, method_label=label).to_csv(
            tables_dir / f"pass_count_exactly_{method_id}.csv", index=False
        )
        at_least[method_id].assign(method_id=method_id, method_label=label).to_csv(
            tables_dir / f"pass_count_at_least_{method_id}.csv", index=False
        )

    x_label = f"Number of PoseBusters tests passed (of {n_tests})"
    _plot_pass_count_curve(
        exactly,
        cfg=cfg,
        n_tests=n_tests,
        title="All targets: ligands passing exactly N tests",
        x_label=x_label,
        path_base=out_dir / "pass_count_exactly_n_tests",
    )
    _plot_pass_count_curve(
        at_least,
        cfg=cfg,
        n_tests=n_tests,
        title="All targets: ligands passing at least N tests",
        x_label=f"Minimum N tests passed (of {n_tests})",
        path_base=out_dir / "pass_count_at_least_n_tests",
    )
