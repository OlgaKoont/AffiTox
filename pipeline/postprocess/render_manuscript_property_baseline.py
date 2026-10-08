#!/usr/bin/env python3
"""Render the manuscript Pearson/property-baseline heatmap from corrected data."""

from __future__ import annotations

import sys
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/tank/scratch/okonovalova/AffiTox")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from analysis.config import AnalysisConfig
from analysis.constants import DEFAULT_METHOD_COLORS
from analysis.correlations import run_correlations
from analysis.plots.correlations import _pivot
from analysis.plots.style import (
    apply_style,
    diverging_cmap,
    plot_heatmap_with_summary_col,
)
from analysis.property_baselines import (
    baseline_pearson_matrix,
    baseline_summary_medians,
    run_property_baselines,
)
from pipeline.postprocess.export_bmc_vector_pdfs import _convert


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


def main() -> None:
    analysis_root = ROOT / "analysis" / "excluding_2z5x_3mjg"
    cfg = AnalysisConfig(
        analysis_root=analysis_root,
        merged_dir=analysis_root / "tables",
        posebusters_dir=analysis_root / "tables" / "posebuster",
        posebusters_boltz2_dir=analysis_root / "tables" / "posebuster",
        dynamicbind_posebusters_label="dynamicbind_new",
        targets=TARGETS,
        methods=METHODS,
        exp_col="pValue",
        color_low="#6584E1",
        color_mid="#DEDAD7",
        color_high="#E5885F",
        method_colors=DEFAULT_METHOD_COLORS.copy(),
        figure_dpi=300,
        font_family="DejaVu Sans",
        n_bootstrap=10_000,
        n_permutation=10_000,
        random_seed=42,
    )
    cfg.ensure_dirs()

    correlation_path = analysis_root / "tables" / "correlations" / "summary_all_proteins.csv"
    if correlation_path.exists():
        correlations = pd.read_csv(correlation_path)
    else:
        correlations = run_correlations(cfg)
    method_body = _pivot(correlations, "pearson_r", cfg).T.astype(float)
    method_summary = method_body.median(axis=1, skipna=True)

    baseline_path = analysis_root / "tables" / "correlations" / "property_null_baselines.csv"
    if baseline_path.exists():
        baselines = pd.read_csv(baseline_path)
    else:
        baselines = run_property_baselines(cfg)
    baseline_body = baseline_pearson_matrix(baselines, cfg)
    baseline_summary = baseline_summary_medians(baselines)

    combined_body = pd.concat([method_body, baseline_body], axis=0)
    combined_summary = pd.concat([method_summary, baseline_summary]).reindex(
        combined_body.index
    )

    apply_style(cfg)
    output_base = (
        analysis_root
        / "figures"
        / "correlations"
        / "heatmaps"
        / "pearson_heatmap_with_property_baselines"
    )
    plot_heatmap_with_summary_col(
        combined_summary,
        combined_body,
        title="Pearson r by method/target with property null baselines",
        cbar_label="Correlation coefficient",
        cmap=diverging_cmap(cfg),
        vmin=-1.0,
        vmax=1.0,
        center=0.0,
        fmt=".2f",
        body_annot=True,
        y_label="Method / null baseline",
        figsize=(12.4, max(5.2, 0.42 * len(combined_body.index) + 1.2)),
        equal_summary_width=True,
        label_panel_ratio=1.55,
        label_summary_wspace=0.005,
        cbar_width_ratio=0.35,
        summary_body_gap_scale=0.5,
        summary_body_wspace=0.14,
        y_tick_rotation=0,
        path_base=output_base,
        summary_label="Median",
    )
    pdf_path, width_mm, height_mm = _convert(output_base.with_suffix(".svg"))
    manuscript_pdf = ROOT / "manuscript" / "pearson_heatmap_with_property_baselines.pdf"
    shutil.copy2(pdf_path, manuscript_pdf)
    print(
        f"{pdf_path} ({width_mm:.1f} x {height_mm:.1f} mm)\n"
        f"{manuscript_pdf}"
    )


if __name__ == "__main__":
    main()
