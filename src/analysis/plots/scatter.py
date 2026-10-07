"""Per-target scatter plots: pKi vs docking score, all methods in one row."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress, pearsonr, spearmanr

from ..config import AnalysisConfig
from ..data import load_merged, prepare_xy
from .style import apply_style, method_color, save_figure


def plot_target_scatter(target: str, cfg: AnalysisConfig) -> None:
    df = load_merged(cfg.merged_dir, target)
    n_methods = len(cfg.methods)
    apply_style(cfg)

    fig, axes = plt.subplots(1, n_methods, figsize=(3.2 * n_methods, 3.4), squeeze=False)
    axes = axes[0]

    for ax, method_id in zip(axes, cfg.methods):
        label = cfg.method_label(method_id)
        metric = cfg.primary_metric(method_id)
        color = method_color(cfg, method_id)

        if metric not in df.columns:
            ax.text(0.5, 0.5, f"{label}\nno data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(label)
            continue

        x, y = prepare_xy(df, metric, cfg.exp_col)
        if len(x) < 3:
            ax.text(0.5, 0.5, f"{label}\ninsufficient data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(label)
            continue

        ax.scatter(x, y, s=14, alpha=0.55, color=color, edgecolors="white", linewidths=0.2)
        slope, intercept, _, _, _ = linregress(x, y)
        xx = np.linspace(x.min(), x.max(), 100)
        ax.plot(xx, slope * xx + intercept, color="black", lw=1.2, alpha=0.8)

        pr, _ = pearsonr(x, y)
        sr, _ = spearmanr(x, y)
        stats = f"r={pr:.3f}\n$\\rho$={sr:.3f}\nn={len(x)}"
        ax.text(
            0.03, 0.97, stats, transform=ax.transAxes, va="top", ha="left", fontsize=6.5,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.92, edgecolor="#cccccc"),
        )
        ax.set_title(label)
        ax.set_xlabel("Experimental pKi")
        ax.set_ylabel("")
        ax.grid(True, alpha=0.25)

    ylabel_fs = 10
    char_w_pt = ylabel_fs * 0.55
    axes[0].set_ylabel("Predicted score", fontsize=ylabel_fs, labelpad=2 * char_w_pt)
    fig.suptitle(
        f"{cfg.target_label(target)}: experimental pKi vs docking scores",
        y=1.02,
        fontsize=11,
    )
    out = cfg.figures_dir / "correlations" / "scatter" / f"scatter_{target.lower()}"
    save_figure(fig, out)


def plot_all_scatters(cfg: AnalysisConfig) -> None:
    for target in cfg.targets:
        plot_target_scatter(target, cfg)
