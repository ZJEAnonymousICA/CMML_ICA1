#!/usr/bin/env python3
"""Export assessor-facing table data and a paper asset manifest."""

from __future__ import annotations

from pathlib import Path
import json

import pandas as pd


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"


def build_table_s1() -> pd.DataFrame:
    rows = [
        {
            "area": "Branch-choice semantics",
            "workshop_scaffold": "Divergent-junction router via branch_rule IDs",
            "expanded_model_where": "Convergent-junction decision with explicit BR1/BR3/BR5 rules (migration_rules.py)",
            "why_it_matters": "The biological decision problem is upstream choice at a convergent bifurcation",
        },
        {
            "area": "Endpoint definition",
            "workshop_scaffold": "No convergent-terminal loss metric or joint-success readout",
            "expanded_model_where": "Loss scored on convergent terminal segments; joint success reported (simulation_support.py)",
            "why_it_matters": "Topology-only stability can miss hierarchy failure",
        },
        {
            "area": "True occlusion",
            "workshop_scaffold": "Minimum diameter floor prevents complete pruning",
            "expanded_model_where": "Conductance can reach zero when segments empty (haemodynamics.py)",
            "why_it_matters": "Branch loss must be representable rather than clipped",
        },
        {
            "area": "Cell conservation",
            "workshop_scaffold": "Cells can effectively leave the system",
            "expanded_model_where": "Outlet-to-inlet recycling preserved (simulation_updated_model.py)",
            "why_it_matters": "Recycling preserves depletion dynamics and competition strength",
        },
        {
            "area": "Turn-aware transport",
            "workshop_scaffold": "Polarity transported in global coordinates",
            "expanded_model_where": "Segment-aware polarity rotation through corners/junctions (network_geometry.py)",
            "why_it_matters": "Keeps migration direction consistent across network turns",
        },
        {
            "area": "Intercalation smoothing",
            "workshop_scaffold": "Not implemented",
            "expanded_model_where": "Redistribution away from bifurcations (migration_rules.py)",
            "why_it_matters": "Suppresses non-biological pile-ups that distort calibre updates",
        },
        {
            "area": "Polarity realignment",
            "workshop_scaffold": "Simple heuristic update",
            "expanded_model_where": "Workshop-style mixed update used in main text (polarity_mixed.py)",
            "why_it_matters": "Tests whether the BR5 mechanism persists under richer polarity dynamics",
        },
    ]
    return pd.DataFrame(rows)


def build_asset_manifest() -> pd.DataFrame:
    rows = [
        {
            "paper_asset": "Main Figure 1",
            "scripts_to_run": "02_run_sim_full_analysis.py",
            "primary_outputs": "data/*.csv (SIM panel datasets) ; figures/Figure1.pdf (final assembled figure)",
        },
        {
            "paper_asset": "Main Figure 2",
            "scripts_to_run": "02_run_sim_full_analysis.py",
            "primary_outputs": "data/*.csv (SIM mechanism datasets) ; figures/Figure2.pdf (final assembled figure)",
        },
        {
            "paper_asset": "Workshop benchmark (Fig. S1)",
            "scripts_to_run": "01_run_workshop_benchmark.py",
            "primary_outputs": "data/workshop_scaffold_benchmark_timecourse.csv",
        },
        {
            "paper_asset": "Supplementary Figure 2",
            "scripts_to_run": "(precomputed only)",
            "primary_outputs": "(see submitted Supporting Materials PDF)",
        },
        {
            "paper_asset": "Supplementary Figure 3",
            "scripts_to_run": "02_run_sim_full_analysis.py",
            "primary_outputs": "data/br4_compare_multi_diam.csv ; data/rule_compare_runs.csv ; data/rule_compare_summary.csv",
        },
        {
            "paper_asset": "Supplementary Figure 4",
            "scripts_to_run": "02_run_sim_full_analysis.py",
            "primary_outputs": "data/alpha_plateau_summary.csv ; data/sweep_summary.csv ; data/sweep_branch_collapse_pct.csv ; data/sweep_regression.csv ; data/review_br5_runs.csv",
        },
        {
            "paper_asset": "Table S1 (CSV)",
            "scripts_to_run": "03_export_review_tables.py",
            "primary_outputs": "data/table_s1_scaffold_expansions.csv",
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    DATA_DIR.mkdir(exist_ok=True)

    table_s1 = build_table_s1()
    manifest = build_asset_manifest()

    table_s1.to_csv(DATA_DIR / "table_s1_scaffold_expansions.csv", index=False)
    manifest.to_csv(DATA_DIR / "paper_asset_manifest.csv", index=False)

    summary = {
        "generated_files": [
            "data/table_s1_scaffold_expansions.csv",
            "data/paper_asset_manifest.csv",
        ]
    }
    (DATA_DIR / "paper_asset_manifest.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
