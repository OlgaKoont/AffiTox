#!/usr/bin/env python3
"""Render the corrected Pearson heatmap in the ultramarine/ivory/vermilion style."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.collections import QuadMesh
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec


ROOT = Path("/mnt/tank/scratch/okonovalova/AffiTox")
sys.path.insert(0, str(ROOT))

from pipeline.postprocess.export_bmc_vector_pdfs import _convert, _prepare_svg


ANALYSIS_ROOT = ROOT / "analysis" / "excluding_2z5x_3mjg"
TARGETS = [
    "1g5m",
    "2v5z",
    "3eyg",
    "3jy9",
    "3lxk",
    "11ue",
    "4ase",
    "4f65",
    "4tz4",
    "4zau",
    "5jkv",
    "5mo4",
    "6gqj",
    "6jok",
    "5lf3",
    "7kk3",
]
METHODS = ["boltz2", "dynamicbind", "gnina", "plapt", "qvina"]
METHOD_LABELS = {
    "boltz2": "Boltz-2",
    "dynamicbind": "DynamicBind",
    "gnina": "GNINA 1.3",
    "plapt": "PLAPT",
    "qvina": "QVina2",
}


def _annotations(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.map(lambda value: "" if pd.isna(value) else f"{value:.2f}")


def main() -> None:
    table_path = ANALYSIS_ROOT / "tables" / "correlations" / "summary_all_proteins.csv"
    correlations = pd.read_csv(table_path)
    pearson = correlations[
        correlations["method_id"].isin(METHODS)
    ].pivot(index="method_id", columns="target", values="pearson_r")
    pearson = pearson.reindex(index=METHODS, columns=TARGETS).astype(float)
    pearson.index = [METHOD_LABELS[method] for method in pearson.index]
    pearson.columns = [target.upper() for target in pearson.columns]
    median = pearson.median(axis=1, skipna=True).to_frame("Median")

    cmap = LinearSegmentedColormap.from_list(
        "ultramarine_ivory_vermilion",
        ["#3558C5", "#FFF9EE", "#EF4938"],
        N=256,
    )
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(10.54, 3.48))
    grid = GridSpec(
        1,
        4,
        figure=fig,
        width_ratios=[1.0, 0.38, 16.0, 0.55],
        wspace=0.08,
        left=0.13,
        right=0.96,
        bottom=0.22,
        top=0.95,
    )
    median_ax = fig.add_subplot(grid[0, 0])
    body_ax = fig.add_subplot(grid[0, 2])
    colorbar_ax = fig.add_subplot(grid[0, 3])

    common = {
        "cmap": cmap,
        "vmin": -1.0,
        "vmax": 1.0,
        "center": 0.0,
        "linewidths": 0.3,
        "linecolor": "white",
        "square": True,
        "annot_kws": {"size": 8, "color": "black"},
    }
    sns.heatmap(
        median,
        ax=median_ax,
        annot=_annotations(median),
        fmt="",
        cbar=False,
        **common,
    )
    sns.heatmap(
        pearson,
        ax=body_ax,
        annot=_annotations(pearson),
        fmt="",
        yticklabels=False,
        cbar=True,
        cbar_ax=colorbar_ax,
        cbar_kws={"label": "Correlation coefficient"},
        **common,
    )

    median_ax.set_xlabel("")
    median_ax.set_ylabel("Docking method")
    median_ax.tick_params(axis="x", rotation=0)
    median_ax.tick_params(axis="y", rotation=0)
    body_ax.set_xlabel("Target")
    body_ax.set_ylabel("")
    body_ax.tick_params(axis="x", rotation=0)
    body_ax.tick_params(axis="y", left=False, labelleft=False)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    mesh = next(c for c in body_ax.collections if isinstance(c, QuadMesh))
    table_bbox = mesh.get_window_extent(renderer)
    inv = fig.transFigure.inverted()
    _, y0 = inv.transform((table_bbox.x0, table_bbox.y0))
    _, y1 = inv.transform((table_bbox.x1, table_bbox.y1))
    cpos = colorbar_ax.get_position()
    colorbar_ax.set_position([cpos.x0, y0, cpos.width, y1 - y0])

    output_base = (
        ANALYSIS_ROOT
        / "figures"
        / "correlations"
        / "heatmaps"
        / "pearson_heatmap_ultramarine_ivory_vermilion_square_cells"
    )
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    pdf_path, width_mm, height_mm = _convert(output_base.with_suffix(".svg"))
    vector_svg, _, _ = _prepare_svg(output_base.with_suffix(".svg"))
    output_base.with_suffix(".svg").write_bytes(vector_svg)
    requested_svg = (
        ROOT
        / "analysis"
        / "figures"
        / "correlations"
        / "heatmaps"
        / "pearson_heatmap_ultramarine_ivory_vermilion_square_cells (1).svg"
    )
    requested_svg.write_bytes(vector_svg)
    requested_pdf = (
        ROOT
        / "analysis"
        / "figures"
        / "correlations"
        / "heatmaps"
        / "pearson_heatmap_ultramarine_ivory_vermilion_square_cells (1).pdf"
    )
    shutil.copy2(pdf_path, requested_pdf)
    for variant in ("by_protein", "by_pdb"):
        dest = (
            ANALYSIS_ROOT
            / "figures"
            / variant
            / "correlations"
            / "heatmaps"
            / output_base.name
        )
        dest.parent.mkdir(parents=True, exist_ok=True)
        for suffix in (".png", ".svg", ".pdf"):
            src = output_base.with_suffix(suffix)
            if src.exists():
                shutil.copy2(src, dest.with_suffix(suffix))
    print(
        f"{pdf_path} ({width_mm:.1f} x {height_mm:.1f} mm)\n"
        f"{requested_svg}\n"
        f"{requested_pdf}"
    )


if __name__ == "__main__":
    main()
