"""nEF violin figures (article style, active + inactive panels)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.gridspec import GridSpec

from ..config import AnalysisConfig
from .style import apply_style, save_figure

# Article blue → cyan at nEF = 0.5 (same for every target) → ivory → vermilion.
# Ivory after 0.5 keeps the 0.5 band from mixing into purple.
NEF_POINT_CMAP = LinearSegmentedColormap.from_list(
    "nef_points",
    [
        (0.00, "#3558C5"),
        (0.50, "#7BA3E8"),
        (0.62, "#FFF9EE"),
        (1.00, "#EF4938"),
    ],
    N=256,
)


def _panel_plot(
    ax: plt.Axes,
    df: pd.DataFrame,
    cfg: AnalysisConfig,
    col: str,
    title: str,
    panel_letter: str,
    order: list[str],
    cmap,
    seed: int,
    frac: int,
) -> None:
    rng = np.random.default_rng(seed)
    x_positions = np.arange(len(order))
    values_by_method = [
        df.loc[df["method_label"] == method, col].astype(float).dropna().to_numpy()
        for method in order
    ]

    violin = ax.violinplot(
        values_by_method,
        positions=x_positions,
        widths=0.82,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body in violin["bodies"]:
        body.set_facecolor("#DEDAD7")
        body.set_edgecolor("#2A2A2A")
        body.set_linewidth(0.6)
        body.set_alpha(0.22)

    label_thr = 0.45
    for i, vals in enumerate(values_by_method):
        if len(vals) == 0:
            continue
        method = order[i]
        df_method = (
            df.loc[df["method_label"] == method, ["target", col]]
            .dropna()
            .copy()
            .reset_index(drop=True)
        )
        proteins = (
            df_method["target"].astype(str).map(cfg.target_label).to_numpy()
        )
        jitter_x = x_positions[i] + rng.uniform(-0.18, 0.18, size=len(vals))
        jitter_y = vals + rng.uniform(-0.015, 0.015, size=len(vals))
        jitter_y = np.clip(jitter_y, 0.0, 1.0)
        low_mask = vals < label_thr
        jitter_x[low_mask] = x_positions[i]
        jitter_y[low_mask] = vals[low_mask]
        ax.scatter(
            jitter_x,
            jitter_y,
            s=34,
            c=[cmap(float(np.clip(v, 0.0, 1.0))) for v in vals],
            alpha=0.82,
            edgecolors="white",
            linewidths=0.55,
            zorder=3,
        )
        label_idx = sorted(
            [j for j, yi in enumerate(vals) if yi < label_thr],
            key=lambda j: vals[j],
        )[:4]
        last_y_by_side = {"right": -10.0, "left": -10.0}
        min_gap = 0.08
        for rank, j in enumerate(label_idx):
            side = "right" if rank % 2 == 0 else "left"
            ha = "left" if side == "right" else "right"
            x_text = (jitter_x[j] + 0.10) if side == "right" else (jitter_x[j] - 0.10)
            y_text = float(jitter_y[j])
            if (y_text - last_y_by_side[side]) < min_gap:
                y_text = last_y_by_side[side] + min_gap
            y_text = float(np.clip(y_text, 0.0, 1.0))
            last_y_by_side[side] = y_text
            ax.annotate(
                proteins[j],
                (jitter_x[j], jitter_y[j]),
                textcoords="data",
                xytext=(x_text, y_text),
                ha=ha,
                va="center",
                fontsize=6,
                alpha=0.9,
                zorder=5,
                arrowprops=dict(arrowstyle="-", color="#666666", lw=0.4, shrinkA=0, shrinkB=2),
            )
        ax.scatter(
            [x_positions[i]],
            [float(np.median(vals))],
            s=50,
            marker="D",
            c=["#222222"],
            edgecolors="white",
            linewidths=0.7,
            zorder=4,
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(order, fontsize=9)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel(f"nEF{frac}", fontsize=11)
    ax.grid(axis="y", linestyle="--", alpha=0.35, linewidth=0.7)
    ax.axhline(0.5, color="#666666", linestyle=":", linewidth=1.0, alpha=0.7)
    ax.set_title(title, fontsize=12, pad=10)
    ax.text(
        0.03, 0.98, panel_letter,
        transform=ax.transAxes,
        fontsize=18, fontweight="bold", va="top", ha="left",
    )


def plot_nef_violins(ef: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    out_dir = cfg.figures_dir / "enrichment" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    for old in out_dir.glob("heatmap_nEF*"):
        old.unlink()

    order = [cfg.method_label(m) for m in cfg.methods if m in set(ef["method_id"])]

    for frac in ("1", "5", "10"):
        active_col = f"nEF{frac}_active"
        low_col = f"nEF{frac}_low"
        if active_col not in ef.columns or low_col not in ef.columns:
            continue

        fig = plt.figure(figsize=(13.4, 5.6))
        gs = GridSpec(
            1, 3, figure=fig,
            width_ratios=[1.0, 1.0, 0.045],
            wspace=0.22,
            left=0.07, right=0.94, top=0.88, bottom=0.12,
        )
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)
        cax = fig.add_subplot(gs[0, 2])
        _panel_plot(
            ax0, ef, cfg, active_col,
            f"nEF{frac} (actives: Ki < 1000 nM)", "A", order, NEF_POINT_CMAP,
            seed=cfg.random_seed, frac=int(frac),
        )
        _panel_plot(
            ax1, ef, cfg, low_col,
            f"nEF{frac},low (inactives: Ki ≥ 1000 nM)", "B", order, NEF_POINT_CMAP,
            seed=cfg.random_seed + 1, frac=int(frac),
        )
        sm = ScalarMappable(norm=Normalize(0.0, 1.0), cmap=NEF_POINT_CMAP)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_ticks([0.0, 0.5, 1.0])
        cbar.set_label("nEF")
        fig.canvas.draw()
        pos = ax0.get_position()
        cpos = cax.get_position()
        cax.set_position([cpos.x0, pos.y0, cpos.width, pos.height])
        save_figure(fig, out_dir / f"summary_nEF{frac}_combined")
