#!/usr/bin/env python3
"""Orchestrate ToxAffinity article analysis pipeline."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# Allow running without pip install
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.config import build_parser, config_from_args
from analysis.correlations import add_bootstrap_permutation, run_correlations
from analysis.enrichment import run_enrichment
from analysis.property_baselines import run_property_baselines
from analysis.pki_range_audit import run_pki_range_audit
from analysis.inferential import (
    run_friedman_tests,
    run_pairwise_method_tests,
    run_per_target_pairwise_tests,
    run_posebusters_pairwise_tests,
)
from analysis.posebusters import run_posebusters
from analysis.constants import METHOD_LABELS, PRIMARY_METRICS
from analysis.plots.correlations import plot_correlation_heatmaps
from analysis.plots.scatter import plot_all_scatters
from analysis.plots.enrichment import plot_nef_violins
from analysis.plots.posebusters import (
    plot_posebusters_at_least_all_tests_combined,
    plot_posebusters_at_least_all_tests_combined_var2,
    plot_posebusters_combined,
    plot_posebusters_pass_count_curves,
    plot_posebusters_per_check,
    plot_posebusters_summary,
)
from analysis.plots.inferential import (
    plot_inferential_heatmaps,
    plot_per_target_inferential_heatmaps,
)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    cfg = config_from_args(args)
    cfg.ensure_dirs()

    print("=" * 72)
    print("ToxAffinity analysis pipeline")
    print(f"  analysis root : {cfg.analysis_root}")
    print(f"  merged data   : {cfg.merged_dir}")
    print(f"  targets       : {', '.join(cfg.targets)}")
    print(f"  methods       : {', '.join(cfg.method_label(m) for m in cfg.methods)}")
    print("=" * 72)

    corr = None
    ef = None
    pb_summary = None
    pb_checks = None
    tests = None
    per_target_tests = None

    if cfg.run_correlations:
        print("[1/5] Correlations (Pearson, Spearman, Kendall + BH-FDR)")
        corr = run_correlations(cfg)
        corr = add_bootstrap_permutation(cfg, corr)
        print(f"      -> {cfg.tables_dir / 'correlations' / 'summary_all_proteins.csv'}")
        baselines = run_property_baselines(cfg)
        print(
            f"      -> {cfg.tables_dir / 'correlations' / 'property_null_baselines.csv'}"
            f" ({len(baselines)} rows)"
        )
        range_audit = run_pki_range_audit(cfg, corr)
        print(
            f"      -> {cfg.tables_dir / 'correlations' / 'pki_range_by_target.csv'}"
            f" ({len(range_audit)} targets)"
        )

    if cfg.run_enrichment:
        print("[2/5] Enrichment (EF1/5/10, nEF active + inactive)")
        ef = run_enrichment(cfg)
        print(f"      -> {cfg.tables_dir / 'enrichment' / 'ef_summary_all_proteins.csv'}")

    if cfg.run_posebusters:
        print("[3/5] PoseBusters pass rates (100%, 90%, 50%)")
        pb_summary, pb_checks = run_posebusters(cfg)
        print(f"      -> {cfg.tables_dir / 'posebusters' / 'pass_rates_by_target_method.csv'}")

    if cfg.run_inferential:
        if corr is None:
            p = cfg.tables_dir / "correlations" / "summary_all_proteins.csv"
            if p.exists():
                corr = pd.read_csv(p)
        if ef is None:
            p = cfg.tables_dir / "enrichment" / "ef_summary_all_proteins.csv"
            if p.exists():
                ef = pd.read_csv(p)
    if cfg.run_inferential and corr is not None and ef is not None:
        print("[4/5] Pairwise method tests (Wilcoxon + Holm/BH)")
        tests = run_pairwise_method_tests(
            corr, ef, cfg.methods, METHOD_LABELS, PRIMARY_METRICS
        )
        out = cfg.tables_dir / "inferential" / "pairwise_method_tests.csv"
        tests.to_csv(out, index=False)
        print(f"      -> {out}")
        per_target_tests = run_per_target_pairwise_tests(
            cfg, cfg.methods, METHOD_LABELS, PRIMARY_METRICS
        )
        out_pt = cfg.tables_dir / "inferential" / "per_target_pairwise_tests.csv"
        per_target_tests.to_csv(out_pt, index=False)
        print(f"      -> {out_pt}")
        omnibus = run_friedman_tests(corr, ef, cfg.methods, PRIMARY_METRICS)
        out_omni = cfg.tables_dir / "inferential" / "friedman_omnibus_tests.csv"
        omnibus.to_csv(out_omni, index=False)
        print(f"      -> {out_omni}")
        if pb_summary is None:
            ps = cfg.tables_dir / "posebusters" / "pass_rates_by_target_method.csv"
            if ps.exists():
                pb_summary = pd.read_csv(ps)
        if pb_summary is not None and not pb_summary.empty:
            pb_tests = run_posebusters_pairwise_tests(pb_summary)
            out_pb = cfg.tables_dir / "inferential" / "posebusters_pairwise_tests.csv"
            pb_tests.to_csv(out_pb, index=False)
            print(f"      -> {out_pb}")

    if cfg.run_figures:
        print("[5/5] Figures (PNG + SVG, dpi={})".format(cfg.figure_dpi))
        if corr is None:
            p = cfg.tables_dir / "correlations" / "summary_all_proteins.csv"
            if p.exists():
                corr = pd.read_csv(p)
        if ef is None:
            p = cfg.tables_dir / "enrichment" / "ef_summary_all_proteins.csv"
            if p.exists():
                ef = pd.read_csv(p)
        if pb_summary is None or pb_checks is None:
            ps = cfg.tables_dir / "posebusters" / "pass_rates_by_target_method.csv"
            pc = cfg.tables_dir / "posebusters" / "pass_rates_by_check_target_method.csv"
            if ps.exists():
                pb_summary = pd.read_csv(ps)
            if pc.exists():
                pb_checks = pd.read_csv(pc)
        if tests is None:
            pt = cfg.tables_dir / "inferential" / "pairwise_method_tests.csv"
            if pt.exists():
                tests = pd.read_csv(pt)
        if per_target_tests is None:
            pt2 = cfg.tables_dir / "inferential" / "per_target_pairwise_tests.csv"
            if pt2.exists():
                per_target_tests = pd.read_csv(pt2)

        if corr is not None:
            plot_correlation_heatmaps(corr, cfg)
            plot_all_scatters(cfg)
        if ef is not None:
            plot_nef_violins(ef, cfg)
        if pb_summary is not None and pb_checks is not None:
            plot_posebusters_summary(pb_summary, cfg)
            plot_posebusters_combined(pb_summary, pb_checks, cfg)
            plot_posebusters_per_check(pb_checks, cfg)
            plot_posebusters_pass_count_curves(cfg)
            plot_posebusters_at_least_all_tests_combined(pb_summary, pb_checks, cfg)
            plot_posebusters_at_least_all_tests_combined_var2(pb_summary, pb_checks, cfg)
        if tests is not None:
            plot_inferential_heatmaps(tests, cfg)
        if per_target_tests is not None:
            plot_per_target_inferential_heatmaps(per_target_tests, cfg)
        print(f"      -> {cfg.figures_dir}")

    print("Done.")


if __name__ == "__main__":
    main()
