#!/usr/bin/env python3
"""Generate all panel data for the current updated mixed-polarity ICA model."""

from __future__ import annotations

import math
import os
from multiprocessing import Pool, get_all_start_methods, get_context
import importlib as _il

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

# Force the script to run relative to its own directory so all CSV outputs land
# in the expected `data/` folder no matter where the user invokes it from.
BASE_DIR = os.path.dirname(__file__)
os.chdir(BASE_DIR)

import numpy as np
import pandas as pd

_net = _il.import_module("network_geometry")
_flow = _il.import_module("haemodynamics")
_sim = _il.import_module("simulation_updated_model")
_ana = _il.import_module("analysis_updated_model")

# Re-export the key model helpers under local names to keep the long analysis
# script readable without repeated module-qualified lookups.
NSEG = _net.NSEG
NNODE = _net.NNODE
run_simulation = _sim.run_simulation
DEFAULT_PARAMS = _sim.DEFAULT_PARAMS
compute_mean_diameter = _sim.compute_mean_diameter
alpha_sweep = _ana.alpha_sweep
generate_seeds = _ana.generate_seeds
run_multiple_seeds = _ana.run_multiple_seeds
wilson_ci = _ana.wilson_ci
compute_conductance = _flow.compute_conductance

DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR = os.path.join(BASE_DIR, "figures")

BR5_ALPHA = 0.45
COMPARE_MASTER_SEED = 20260413
SWEEP_MASTER_SEED = 42
TRACE_MASTER_SEED = 42
UNSTABLE_MASTER_SEED = 99

N_COMPARE = 100
SWEEP_SEEDS = 50
TIME_SEEDS = 60
N_WORKERS = min(8, os.cpu_count() or 1)

FINE_ALPHAS = np.round(np.arange(0.0, 1.0001, 0.01), 2)
BIAS_ALPHAS = np.round(np.arange(0.0, 1.0001, 0.05), 2)
SELECTED_ALPHAS = (0.10, 0.45, 1.00)

TOLERANCE_PCT_POINTS = 15.0
HIERARCHY_THRESHOLD_UM = 4.0
RUNAWAY_PROB_THRESHOLD = 0.80
SUSTAIN_STEPS = 2
FINAL_STEP = 36
DT_DAYS = DEFAULT_PARAMS["dt_hours"] / 24.0


def ensure_dirs() -> None:
    # All downstream writers assume these directories already exist.
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)


def _pool_factory(n_workers: int):
    # Match the multiprocessing policy used elsewhere in the repo.
    if "fork" in get_all_start_methods():
        return get_context("fork").Pool(processes=n_workers)
    return Pool(processes=n_workers)


def export_segment_coords() -> None:
    # Export a plotting-friendly geometry table so figure assembly can happen in
    # R or other tools without reimplementing the network layout logic.
    lengths = np.ones(NSEG) * DEFAULT_PARAMS["l_seg"]
    raw = _net.make_segment_coords(lengths)
    lower = raw["lower"]
    upper = raw["upper"]

    rows = []
    for seg in range(20):
        rows.append({
            "seg": seg,
            "x1": lower[seg, 0] * 1e6,
            "y1": lower[seg, 1] * 1e6,
            "x2": lower[seg + 1, 0] * 1e6,
            "y2": lower[seg + 1, 1] * 1e6,
        })
    for idx, seg in enumerate(range(20, 40)):
        rows.append({
            "seg": seg,
            "x1": upper[idx, 0] * 1e6,
            "y1": upper[idx, 1] * 1e6,
            "x2": upper[idx + 1, 0] * 1e6,
            "y2": upper[idx + 1, 1] * 1e6,
        })

    vessel_map = {}
    for s in range(5):
        vessel_map[s] = "feeding"
    for s in range(5, 15):
        vessel_map[s] = "proximal"
    for s in range(15, 20):
        vessel_map[s] = "draining"
    for s in range(20, 40):
        vessel_map[s] = "distal"
    for row in rows:
        row["vessel"] = vessel_map[row["seg"]]

    pd.DataFrame(rows).to_csv(os.path.join(DATA_DIR, "segment_coords.csv"), index=False)


def export_result_timeseries(result: dict, prefix: str) -> None:
    # This helper writes both branch-mean summaries and selected segment-level
    # snapshots for a single representative simulation run.
    n_steps = result["Ncell"].shape[0]
    time_days = result["time_days"]

    mean_df = pd.DataFrame({
        "step": range(n_steps),
        "time_days": time_days,
        "mean_diam_prox": result["mean_diam_prox"],
        "mean_diam_dist": result["mean_diam_dist"],
        "bifurcation_lost": result["bifurcation_lost"],
        "branch_collapse": result.get("branch_collapse", np.zeros(n_steps, dtype=int)),
    })
    mean_df.to_csv(os.path.join(DATA_DIR, f"{prefix}_mean_diam.csv"), index=False)

    rows = []
    for t_idx in [0, n_steps - 1]:
        # Only the initial and final segment snapshots are needed for the paper figures.
        for seg in range(NSEG):
            rows.append({
                "step": t_idx,
                "time_days": time_days[t_idx],
                "seg": seg,
                "Ncell": result["Ncell"][t_idx, seg],
                "D_um": result["D"][t_idx, seg],
                "Q_uLday": abs(result["Q"][t_idx, seg]) * 1e9 * 86400,
                "tau_Pa": result["tau"][t_idx, seg],
            })
    pd.DataFrame(rows).to_csv(os.path.join(DATA_DIR, f"{prefix}_segments.csv"), index=False)


def export_multiseed_diameters(results_list: list[dict], prefix: str) -> None:
    # Flatten many simulation histories into one long table for line plots with
    # seed-wise trajectories or ribbons.
    rows = []
    for seed_idx, result in enumerate(results_list):
        for step, time_days in enumerate(result["time_days"]):
            rows.append({
                "seed_idx": seed_idx,
                "step": step,
                "time_days": time_days,
                "mean_diam_prox": result["mean_diam_prox"][step],
                "mean_diam_dist": result["mean_diam_dist"][step],
                "loss_time": result["loss_time"],
                "stable": int(result["loss_time"] == 0),
                "hierarchy_success": int(result["hierarchy_success"]),
            })
    pd.DataFrame(rows).to_csv(os.path.join(DATA_DIR, f"{prefix}_multi_diam.csv"), index=False)


def export_prob_history(result: dict, prefix: str) -> None:
    if "prob_history" not in result or result["prob_history"] is None:
        return
    # BR5 probability histories are saved step-by-step so the emergence of
    # one-sided routing can be plotted directly.
    rows = []
    for step, probs in enumerate(result["prob_history"]):
        rows.append({
            "step": step,
            "P_tau1": probs["P_tau1"],
            "P_tau2": probs["P_tau2"],
            "P_n1": probs["P_n1"],
            "P_n2": probs["P_n2"],
            "P1": probs["P1"],
            "P2": probs["P2"],
        })
    pd.DataFrame(rows).to_csv(os.path.join(DATA_DIR, f"{prefix}_probs.csv"), index=False)


def export_pressure_cellcount(result: dict, prefix: str) -> None:
    # This export tracks the haemodynamic pressure drop against branch occupancy
    # to connect flow cues with cell redistribution over time.
    rows = []
    for step, time_days in enumerate(result["time_days"]):
        rows.append({
            "step": step,
            "time_days": time_days,
            "dp_bifurc": result["P"][step, 5] - result["P"][step, 15],
            "n_prox_mean": np.mean(result["Ncell"][step, 5:15]),
            "n_dist_mean": np.mean(result["Ncell"][step, 20:40]),
            "total_cells": np.sum(result["Ncell"][step]),
        })
    pd.DataFrame(rows).to_csv(
        os.path.join(DATA_DIR, f"{prefix}_pressure_cells.csv"), index=False
    )


def export_conductance_data(result: dict, prefix: str) -> None:
    # Recompute conductance from stored cell counts at a few representative time
    # points rather than storing yet another full history in every run.
    rows = []
    n_steps = result["Ncell"].shape[0] - 1
    for t_idx in [0, n_steps // 3, 2 * n_steps // 3, n_steps]:
        ncell = result["Ncell"][t_idx]
        diam, conductance, _ = compute_conductance(
            ncell,
            DEFAULT_PARAMS["cell_width"],
            DEFAULT_PARAMS["mu"],
            DEFAULT_PARAMS["l_seg"],
        )
        for seg in range(NSEG):
            if ncell[seg] > 0:
                rows.append({
                    "step": t_idx,
                    "seg": seg,
                    "Ncell": ncell[seg],
                    "G": conductance[seg],
                    "D_um": diam[seg] * 1e6,
                })
    pd.DataFrame(rows).to_csv(os.path.join(DATA_DIR, f"{prefix}_conductance.csv"), index=False)


def fisher_exact_twosided(a: int, b: int, c: int, d: int):
    # The counts correspond to a 2x2 contingency table:
    # [[success_a, failure_a], [success_b, failure_b]].
    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2

    lo = max(0, col1 - row2)
    hi = min(row1, col1)

    def prob(x):
        # Enumerate the hypergeometric probability of each feasible table with
        # the same margins as the observed one.
        return (
            math.comb(col1, x)
            * math.comb(total - col1, row1 - x)
            / math.comb(total, row1)
        )

    observed = prob(a)
    p_value = 0.0
    for x in range(lo, hi + 1):
        px = prob(x)
        if px <= observed + 1e-12:
            p_value += px

    if b * c == 0:
        odds_ratio = math.inf if a * d > 0 else 0.0
    else:
        odds_ratio = (a * d) / (b * c)

    return odds_ratio, min(p_value, 1.0)


def results_to_seed_df(results: list[dict], rule: str, alpha: float | None = None) -> pd.DataFrame:
    # Convert a list of per-seed result dictionaries into a tidy one-row-per-seed table.
    rows = []
    for seed_idx, result in enumerate(results):
        rows.append({
            "rule": rule,
            "alpha": alpha if alpha is not None else float("nan"),
            "seed_idx": seed_idx,
            "stable": int(result["loss_time"] == 0),
            "hierarchy_success": int(result["hierarchy_success"]),
            "joint_success": int((result["loss_time"] == 0) and result["hierarchy_success"]),
            "loss_time": int(result["loss_time"]),
            "final_prox_diam": float(result["mean_diam_prox"][-1]),
            "final_dist_diam": float(result["mean_diam_dist"][-1]),
            "final_delta": float(result["mean_diam_prox"][-1] - result["mean_diam_dist"][-1]),
        })
    return pd.DataFrame(rows)


def summarise_rule(seed_df: pd.DataFrame, rule: str) -> list[dict]:
    # Report three outcomes separately because the paper distinguishes topology
    # preservation, hierarchy conditional on survival, and their joint success.
    sub = seed_df.loc[seed_df["rule"] == rule].copy()
    n_total = len(sub)
    n_stable = int(sub["stable"].sum())
    n_joint = int(sub["joint_success"].sum())
    stable_pct, stable_lo, stable_hi = wilson_ci(n_stable, n_total)
    joint_pct, joint_lo, joint_hi = wilson_ci(n_joint, n_total)

    stable_subset = sub.loc[sub["stable"] == 1]
    if len(stable_subset) > 0:
        n_hier = int(stable_subset["hierarchy_success"].sum())
        hier_pct, hier_lo, hier_hi = wilson_ci(n_hier, len(stable_subset))
    else:
        n_hier = 0
        hier_pct = hier_lo = hier_hi = float("nan")

    return [
        {
            "rule": rule,
            "metric": "Topology preserved",
            "n_success": n_stable,
            "n_total": n_total,
            "pct": round(stable_pct, 2),
            "ci_lo": round(stable_lo, 2),
            "ci_hi": round(stable_hi, 2),
        },
        {
            "rule": rule,
            "metric": "Hierarchy among survivors",
            "n_success": n_hier,
            "n_total": len(stable_subset),
            "pct": round(hier_pct, 2) if pd.notna(hier_pct) else float("nan"),
            "ci_lo": round(hier_lo, 2) if pd.notna(hier_lo) else float("nan"),
            "ci_hi": round(hier_hi, 2) if pd.notna(hier_hi) else float("nan"),
        },
        {
            "rule": rule,
            "metric": "Joint success",
            "n_success": n_joint,
            "n_total": n_total,
            "pct": round(joint_pct, 2),
            "ci_lo": round(joint_lo, 2),
            "ci_hi": round(joint_hi, 2),
        },
    ]


def sweep_to_runs_df(sweep: dict, rule_label: str) -> pd.DataFrame:
    # Expand the alpha-sweep tensors into a long-form run table for grouped summaries.
    rows = []
    alpha_values = sweep["alpha_values"]
    seeds = sweep["seeds"]
    for i, alpha in enumerate(alpha_values):
        for j, seed in enumerate(seeds):
            rows.append({
                "rule": rule_label,
                "alpha": round(float(alpha), 2),
                "seed": int(seed),
                "loss_time": int(sweep["loss_time"][i, j]),
                "loss_branch": int(sweep["loss_branch"][i, j]),
                "branch_collapse_time": int(sweep["branch_collapse_time"][i, j]),
                "branch_collapse_branch": int(sweep["branch_collapse_branch"][i, j]),
                "final_prox_diam": float(sweep["final_prox_diam"][i, j]),
                "final_dist_diam": float(sweep["final_dist_diam"][i, j]),
                "hierarchy_success": int(sweep["hierarchy_success"][i, j]),
                "joint_success": int(sweep["joint_success"][i, j]),
                "stable": int(sweep["loss_time"][i, j] == 0),
            })
    return pd.DataFrame(rows)


def summarise_runs(df: pd.DataFrame) -> pd.DataFrame:
    # Summaries are computed alpha-by-alpha because each alpha defines a distinct
    # mechanistic balance between shear and occupancy cues.
    rows = []
    for alpha, grp in df.groupby("alpha", sort=True):
        n = len(grp)
        n_lost = int((grp["loss_time"] > 0).sum())
        n_collapsed = int((grp["branch_collapse_time"] > 0).sum())
        n_hierarchy = int(grp["hierarchy_success"].sum())
        n_joint = int(grp["joint_success"].sum())
        n_stable = int(grp["stable"].sum())
        stable_hierarchy = (
            100.0 * grp.loc[grp["stable"] == 1, "hierarchy_success"].mean()
            if n_stable > 0 else np.nan
        )
        # Wilson intervals are used consistently across the analysis outputs.
        loss_pct, loss_lo, loss_hi = wilson_ci(n_lost, n)
        collapse_pct, collapse_lo, collapse_hi = wilson_ci(n_collapsed, n)
        hier_pct, hier_lo, hier_hi = wilson_ci(n_hierarchy, n)
        joint_pct, joint_lo, joint_hi = wilson_ci(n_joint, n)
        rows.append({
            "alpha": round(float(alpha), 2),
            "n_total": n,
            "n_stable": n_stable,
            "loss_pct": round(loss_pct, 2),
            "loss_ci_lo": round(loss_lo, 2),
            "loss_ci_hi": round(loss_hi, 2),
            "branch_collapse_pct": round(collapse_pct, 2),
            "branch_collapse_ci_lo": round(collapse_lo, 2),
            "branch_collapse_ci_hi": round(collapse_hi, 2),
            "hierarchy_pct": round(hier_pct, 2),
            "hierarchy_ci_lo": round(hier_lo, 2),
            "hierarchy_ci_hi": round(hier_hi, 2),
            "joint_success_pct": round(joint_pct, 2),
            "joint_success_ci_lo": round(joint_lo, 2),
            "joint_success_ci_hi": round(joint_hi, 2),
            "stable_hierarchy_pct": round(stable_hierarchy, 2) if not np.isnan(stable_hierarchy) else np.nan,
            "mean_final_prox_diam": round(grp["final_prox_diam"].mean(), 3),
            "mean_final_dist_diam": round(grp["final_dist_diam"].mean(), 3),
            "mean_final_delta": round((grp["final_prox_diam"] - grp["final_dist_diam"]).mean(), 3),
        })
    return pd.DataFrame(rows).sort_values("alpha").reset_index(drop=True)


def export_time_surface(sweep: dict, filename: str, field: str) -> None:
    # Convert alpha x time matrices into long-form tables that contour/heatmap
    # plotting code can consume directly.
    rows = []
    for i, alpha in enumerate(sweep["alpha_values"]):
        for step in range(sweep["Nt"] + 1):
            rows.append({
                "alpha": round(float(alpha), 2),
                "step": step,
                "time_days": step * DT_DAYS,
                field: float(sweep[field][i, step]),
            })
    pd.DataFrame(rows).to_csv(os.path.join(DATA_DIR, filename), index=False)


def export_regression_breakdown(sweep: dict, filename: str) -> None:
    # This helper isolates which branch fails at each alpha, collapsing the full
    # seed-level tensor into one row per alpha.
    rows = []
    for i, alpha in enumerate(sweep["alpha_values"]):
        loss_time = sweep["loss_time"][i]
        loss_branch = sweep["loss_branch"][i]
        n = len(loss_time)
        rows.append({
            "alpha": round(float(alpha), 2),
            "total_pct": round(100.0 * np.mean(loss_time > 0), 2),
            "proximal_pct": round(100.0 * np.mean(loss_branch == 1), 2),
            "distal_pct": round(100.0 * np.mean(loss_branch == 2), 2),
            "n_total": n,
        })
    pd.DataFrame(rows).to_csv(os.path.join(DATA_DIR, filename), index=False)


def contiguous_peak_plateau(summary_df: pd.DataFrame):
    # Identify the broad "good-enough" alpha window around the optimum rather
    # than treating a single grid-point maximum as uniquely meaningful.
    summary_df = summary_df.sort_values("alpha").reset_index(drop=True)
    peak_idx = int(summary_df["joint_success_pct"].idxmax())
    peak_alpha = float(summary_df.loc[peak_idx, "alpha"])
    peak_value = float(summary_df.loc[peak_idx, "joint_success_pct"])
    threshold = peak_value - TOLERANCE_PCT_POINTS

    plateau = summary_df.loc[summary_df["joint_success_pct"] >= threshold].copy()
    plateau = plateau.sort_values("alpha").reset_index(drop=True)

    step = float(np.median(np.diff(summary_df["alpha"].to_numpy(dtype=float))))
    peak_pos = int(plateau.index[np.isclose(plateau["alpha"], peak_alpha)][0])
    left = peak_pos
    right = peak_pos

    while left > 0 and np.isclose(plateau.loc[left, "alpha"] - plateau.loc[left - 1, "alpha"], step):
        left -= 1
    while right < len(plateau) - 1 and np.isclose(plateau.loc[right + 1, "alpha"] - plateau.loc[right, "alpha"], step):
        right += 1

    return plateau.iloc[left:right + 1].copy(), peak_alpha, peak_value, threshold, step


def first_sustained_true(mask: np.ndarray, sustain_steps: int) -> int:
    # This compresses a noisy Boolean trace into the first time an effect stays
    # present for `sustain_steps` consecutive steps.
    run = 0
    for idx, flag in enumerate(mask):
        run = run + 1 if flag else 0
        if run >= sustain_steps:
            return idx - sustain_steps + 1
    return 0


def _run_bias_single(args):
    alpha, seed = args
    # This worker records only the first loss event, which is enough for the
    # topology-bias analysis and keeps the result payload small.
    result = run_simulation(branch_rule=5, alpha=float(alpha), seed=int(seed))
    loss_time = int(result["loss_time"])
    first_loss_branch = 0
    if loss_time > 0:
        first_loss_branch = int(result["bifurcation_lost"][loss_time])
    return {
        "alpha": round(float(alpha), 2),
        "seed": int(seed),
        "loss_time": loss_time,
        "first_loss_branch": first_loss_branch,
        "stable": int(loss_time == 0),
    }


def _run_imbalance_single(args):
    alpha, seed = args
    result = run_simulation(branch_rule=5, alpha=float(alpha), seed=int(seed))
    # Hierarchy onset is defined as a sustained diameter gap above a fixed threshold.
    delta = result["mean_diam_prox"] - result["mean_diam_dist"]
    mask = delta >= HIERARCHY_THRESHOLD_UM
    step = first_sustained_true(mask, SUSTAIN_STEPS)
    return {
        "alpha": round(float(alpha), 2),
        "seed": int(seed),
        "event_step": int(step),
        "event_day": round(step * DT_DAYS, 6),
        "event_observed": int(step > 0),
    }


def _run_runaway_single(args):
    alpha, seed = args
    result = run_simulation(
        branch_rule=5,
        alpha=float(alpha),
        seed=int(seed),
        record_probabilities=True,
    )
    probs = pd.DataFrame(result["prob_history"])
    # One-sided routing corresponds to either branch probability dominating strongly.
    dominance = np.maximum(probs["P1"], probs["P2"])
    mask = dominance >= RUNAWAY_PROB_THRESHOLD
    step = first_sustained_true(mask.to_numpy(), SUSTAIN_STEPS)
    if step > 0:
        dominant_branch = (
            "proximal-dominant" if probs["P1"].iloc[step] >= probs["P2"].iloc[step]
            else "distal-dominant"
        )
    else:
        dominant_branch = "none"
    return {
        "alpha": round(float(alpha), 2),
        "seed": int(seed),
        "event_step": int(step),
        "event_day": round(step * DT_DAYS, 6),
        "event_observed": int(step > 0),
        "dominant_branch": dominant_branch,
        "max_probability": float(dominance.max()),
    }


def cumulative_curve(runs: pd.DataFrame, event_label: str) -> pd.DataFrame:
    # Build a simple cumulative-incidence style curve from step-level event times.
    rows = []
    for alpha, grp in runs.groupby("alpha", sort=True):
        n = len(grp)
        for step in range(FINAL_STEP + 1):
            observed = ((grp["event_observed"] == 1) & (grp["event_step"] <= step)).sum()
            rows.append({
                "alpha": alpha,
                "event": event_label,
                "step": step,
                "time_days": step * DT_DAYS,
                "cumulative_pct": 100.0 * observed / n,
            })
    return pd.DataFrame(rows)


def kaplan_meier_curve(times, events, alpha):
    # This implementation is intentionally explicit so the output table carries
    # `n_at_risk` and `n_events` at each event time for plotting/inspection.
    times = pd.Series(times, dtype=int)
    events = pd.Series(events, dtype=int)
    unique_steps = sorted(times.loc[events == 1].unique())
    at_risk = len(times)
    surv = 1.0
    rows = [{
        "alpha": alpha,
        "step": 0,
        "time_days": 0.0,
        "survival_pct": 100.0,
        "n_at_risk": at_risk,
        "n_events": 0,
    }]
    for step in unique_steps:
        d = int(((times == step) & (events == 1)).sum())
        c = int(((times == step) & (events == 0)).sum())
        if at_risk > 0:
            # Standard Kaplan-Meier multiplicative survival update.
            surv *= (1.0 - d / at_risk)
        rows.append({
            "alpha": alpha,
            "step": int(step),
            "time_days": step * DT_DAYS,
            "survival_pct": 100.0 * surv,
            "n_at_risk": at_risk,
            "n_events": d,
        })
        at_risk -= d + c
    if rows[-1]["step"] != FINAL_STEP:
        rows.append({
            "alpha": alpha,
            "step": FINAL_STEP,
            "time_days": FINAL_STEP * DT_DAYS,
            "survival_pct": 100.0 * surv,
            "n_at_risk": at_risk,
            "n_events": 0,
        })
    return pd.DataFrame(rows)


def log_rank_test(df: pd.DataFrame, alpha_a: float, alpha_b: float):
    # Compare two survival curves via the standard observed-vs-expected event count.
    sub = df.loc[df["alpha"].isin([alpha_a, alpha_b])].copy()
    event_steps = sorted(sub.loc[sub["event"] == 1, "step"].unique())
    observed_a = expected_a = variance_a = 0.0
    for step in event_steps:
        n_a = int((sub["alpha"].eq(alpha_a) & (sub["step"] >= step)).sum())
        n_b = int((sub["alpha"].eq(alpha_b) & (sub["step"] >= step)).sum())
        d_a = int((sub["alpha"].eq(alpha_a) & sub["step"].eq(step) & sub["event"].eq(1)).sum())
        d_b = int((sub["alpha"].eq(alpha_b) & sub["step"].eq(step) & sub["event"].eq(1)).sum())
        n = n_a + n_b
        d = d_a + d_b
        if n <= 1 or d == 0:
            continue
        exp_a = d * (n_a / n)
        var_a = (n_a * n_b * d * (n - d)) / (n * n * (n - 1))
        observed_a += d_a
        expected_a += exp_a
        variance_a += var_a
    chi2 = ((observed_a - expected_a) ** 2 / variance_a) if variance_a > 0 else 0.0
    p_value = math.erfc(math.sqrt(chi2 / 2.0))
    return chi2, p_value


def run_matched_rule_comparison() -> pd.DataFrame:
    print(f"Running matched updated-model comparison ({N_COMPARE} seeds per rule)...", flush=True)
    # The same master seed is used across rules so differences arise from the
    # branch-choice rule, not from unmatched random seed sets.
    br3 = run_multiple_seeds(branch_rule=3, n_seeds=N_COMPARE, master_seed=COMPARE_MASTER_SEED)
    br4 = run_multiple_seeds(branch_rule=4, n_seeds=N_COMPARE, master_seed=COMPARE_MASTER_SEED)
    br5 = run_multiple_seeds(
        branch_rule=5,
        n_seeds=N_COMPARE,
        alpha=BR5_ALPHA,
        master_seed=COMPARE_MASTER_SEED,
    )

    export_multiseed_diameters(br3, "br3_compare")
    export_multiseed_diameters(br4, "br4_compare")
    export_multiseed_diameters(br5, "br5_compare")

    seed_df = pd.concat(
        [
            results_to_seed_df(br3, "BR3"),
            results_to_seed_df(br4, "BR4"),
            results_to_seed_df(br5, "BR5", alpha=BR5_ALPHA),
        ],
        ignore_index=True,
    )
    seed_df.to_csv(os.path.join(DATA_DIR, "rule_compare_runs.csv"), index=False)

    summary_rows = []
    for rule in ["BR3", "BR4", "BR5"]:
        summary_rows.extend(summarise_rule(seed_df, rule))
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(DATA_DIR, "rule_compare_summary.csv"), index=False)

    br3_joint = int(seed_df.loc[seed_df["rule"] == "BR3", "joint_success"].sum())
    br5_joint = int(seed_df.loc[seed_df["rule"] == "BR5", "joint_success"].sum())
    # Fisher's exact test is reserved for the headline BR3-vs-BR5 comparison.
    odds_ratio, p_value = fisher_exact_twosided(
        br5_joint, N_COMPARE - br5_joint, br3_joint, N_COMPARE - br3_joint
    )
    test_df = pd.DataFrame([{
        "metric": "Joint success",
        "rule_a": "BR5",
        "rule_b": "BR3",
        "n_a": N_COMPARE,
        "n_b": N_COMPARE,
        "success_a": br5_joint,
        "success_b": br3_joint,
        "risk_diff_pct": round(100.0 * ((br5_joint / N_COMPARE) - (br3_joint / N_COMPARE)), 2),
        "odds_ratio": round(odds_ratio, 3) if math.isfinite(odds_ratio) else odds_ratio,
        "p_value": p_value,
    }])
    test_df.to_csv(os.path.join(DATA_DIR, "br3_br5_joint_test.csv"), index=False)

    panel_df = summary_df.loc[summary_df["rule"].isin(["BR3", "BR5"])].copy()
    panel_df["metric_order"] = panel_df["metric"].map({
        "Topology preserved": 1,
        "Hierarchy among survivors": 2,
        "Joint success": 3,
    })
    panel_df["label"] = panel_df.apply(lambda r: f"{int(r['n_success'])}/{int(r['n_total'])}", axis=1)
    panel_df = panel_df.sort_values(["metric_order", "rule"]).reset_index(drop=True)
    panel_df.to_csv(os.path.join(DATA_DIR, "figure1d_endpoint_panel.csv"), index=False)

    annot_df = pd.DataFrame([
        {
            "metric": "Joint success",
            "x": 3.0,
            "y": float(panel_df.loc[panel_df["metric"] == "Joint success", "ci_hi"].max()) + 8.0,
            "label": f"Fisher exact p = {p_value:.2e}",
        },
        {
            "metric": "Hierarchy among survivors",
            "x": 2.0,
            "y": float(
                panel_df.loc[panel_df["metric"] == "Hierarchy among survivors", "ci_hi"].max()
            ) + 8.0,
            "label": "among stable seeds only",
        },
    ])
    annot_df.to_csv(os.path.join(DATA_DIR, "figure1d_endpoint_annotations.csv"), index=False)

    rule_outcome_rows = []
    for rule in ["BR3", "BR4", "BR5"]:
        sub = seed_df.loc[seed_df["rule"] == rule]
        rule_outcome_rows.append({
            "rule": rule,
            "n_total": len(sub),
            "topology_preserved_n": int(sub["stable"].sum()),
            "topology_preserved_pct": round(100.0 * sub["stable"].mean(), 2),
            "joint_success_n": int(sub["joint_success"].sum()),
            "joint_success_pct": round(100.0 * sub["joint_success"].mean(), 2),
            "mean_final_prox_um": round(sub["final_prox_diam"].mean(), 3),
            "mean_final_dist_um": round(sub["final_dist_diam"].mean(), 3),
        })
    pd.DataFrame(rule_outcome_rows).to_csv(
        os.path.join(DATA_DIR, "rule_outcome_summary.csv"), index=False
    )
    return seed_df


def run_representative_traces() -> None:
    print("Running representative BR1 and BR5 traces for updated model...", flush=True)
    # BR1 provides a deliberately simple comparator trace.
    br1 = run_simulation(branch_rule=1, seed=42)
    export_result_timeseries(br1, "br1")

    stable_result = None
    for seed in generate_seeds(150, master_seed=TRACE_MASTER_SEED):
        result = run_simulation(
            branch_rule=5,
            alpha=BR5_ALPHA,
            seed=int(seed),
            record_probabilities=True,
        )
        if result["loss_time"] == 0 and result["hierarchy_success"]:
            stable_result = result
            break
    if stable_result is None:
        # Fall back to a deterministic seed so the script always emits outputs.
        stable_result = run_simulation(
            branch_rule=5,
            alpha=BR5_ALPHA,
            seed=42,
            record_probabilities=True,
        )

    export_prob_history(stable_result, "br5_stable")
    export_pressure_cellcount(stable_result, "br5_stable")
    export_conductance_data(stable_result, "br5_stable")

    unstable_result = None
    for seed in generate_seeds(150, master_seed=UNSTABLE_MASTER_SEED):
        result = run_simulation(
            branch_rule=5,
            alpha=1.0,
            seed=int(seed),
            record_probabilities=True,
        )
        if result["loss_time"] > 0:
            unstable_result = result
            break
    if unstable_result is not None:
        export_result_timeseries(unstable_result, "br5_unstable")
        export_prob_history(unstable_result, "br5_unstable")


def run_fine_br5_sweep() -> pd.DataFrame:
    print(
        f"Running updated-model BR5 fine sweep ({len(FINE_ALPHAS)} alpha x {SWEEP_SEEDS} seeds)...",
        flush=True,
    )
    seeds = generate_seeds(SWEEP_SEEDS, master_seed=SWEEP_MASTER_SEED)
    n_alpha = len(FINE_ALPHAS)
    n_seeds = len(seeds)

    loss_time = np.zeros((n_alpha, n_seeds), dtype=int)
    loss_branch = np.zeros((n_alpha, n_seeds), dtype=int)
    branch_collapse_time = np.zeros((n_alpha, n_seeds), dtype=int)
    branch_collapse_branch = np.zeros((n_alpha, n_seeds), dtype=int)
    final_prox_diam = np.zeros((n_alpha, n_seeds))
    final_dist_diam = np.zeros((n_alpha, n_seeds))
    hierarchy_success = np.zeros((n_alpha, n_seeds), dtype=int)
    joint_success = np.zeros((n_alpha, n_seeds), dtype=int)

    for i, alpha in enumerate(FINE_ALPHAS):
        if i % 10 == 0:
            print(f"  alpha {i + 1}/{n_alpha} ({alpha:.2f})", flush=True)
        for j, seed in enumerate(seeds):
            # This sweep is intentionally serial: it is long-running but simpler
            # to reason about and reproduce than a more aggressively parallel driver.
            result = run_simulation(branch_rule=5, alpha=float(alpha), seed=int(seed))
            loss_time[i, j] = int(result["loss_time"])
            loss_branch[i, j] = int(result["bifurcation_lost"][-1])
            branch_collapse_time[i, j] = int(result["branch_collapse_time"])
            branch_collapse_branch[i, j] = int(result["branch_collapse"][-1])
            final_prox_diam[i, j] = float(result["mean_diam_prox"][-1])
            final_dist_diam[i, j] = float(result["mean_diam_dist"][-1])
            hierarchy_success[i, j] = int(result["hierarchy_success"])
            joint_success[i, j] = int((result["loss_time"] == 0) and result["hierarchy_success"])

    sweep = {
        "alpha_values": FINE_ALPHAS,
        "seeds": seeds,
        "loss_time": loss_time,
        "loss_branch": loss_branch,
        "branch_collapse_time": branch_collapse_time,
        "branch_collapse_branch": branch_collapse_branch,
        "final_prox_diam": final_prox_diam,
        "final_dist_diam": final_dist_diam,
        "hierarchy_success": hierarchy_success,
        "joint_success": joint_success,
        "Nt": FINAL_STEP,
    }

    sweep["loss_pct"] = np.zeros((n_alpha, FINAL_STEP + 1))
    sweep["branch_collapse_pct"] = np.zeros((n_alpha, FINAL_STEP + 1))
    for i in range(n_alpha):
        for t in range(1, FINAL_STEP + 1):
            # Convert exact event times into cumulative proportions by time.
            lost_by_t = np.sum((loss_time[i] > 0) & (loss_time[i] <= t))
            sweep["loss_pct"][i, t] = 100.0 * lost_by_t / n_seeds
            collapsed_by_t = np.sum(
                (branch_collapse_time[i] > 0) & (branch_collapse_time[i] <= t)
            )
            sweep["branch_collapse_pct"][i, t] = 100.0 * collapsed_by_t / n_seeds

    runs = sweep_to_runs_df(sweep, "BR5")
    summary = summarise_runs(runs)
    runs.to_csv(os.path.join(DATA_DIR, "review_br5_runs.csv"), index=False)
    summary.to_csv(os.path.join(DATA_DIR, "sweep_summary.csv"), index=False)
    export_time_surface(sweep, "sweep_loss_pct.csv", "loss_pct")
    export_time_surface(sweep, "sweep_branch_collapse_pct.csv", "branch_collapse_pct")
    export_regression_breakdown(sweep, "sweep_regression.csv")

    plateau, peak_alpha, peak_value, threshold, step = contiguous_peak_plateau(summary)
    # Export the plateau as an explicit artifact so the review text can cite it directly.
    pd.DataFrame([{
        "peak_alpha": round(peak_alpha, 2),
        "peak_joint_success_pct": round(peak_value, 2),
        "threshold_pct": round(threshold, 2),
        "plateau_alpha_lo": round(float(plateau["alpha"].min()), 2),
        "plateau_alpha_hi": round(float(plateau["alpha"].max()), 2),
        "grid_step": round(step, 4),
        "definition": "contiguous alphas around the maximum within 15 percentage points of peak joint success",
    }]).to_csv(os.path.join(DATA_DIR, "alpha_plateau_summary.csv"), index=False)

    out_df = summary.loc[:, [
        "alpha",
        "loss_pct",
        "stable_hierarchy_pct",
        "joint_success_pct",
    ]].copy()
    out_df["topology_preserved_pct"] = 100.0 - out_df["loss_pct"]
    out_df = out_df.rename(columns={"stable_hierarchy_pct": "hierarchy_among_survivors_pct"})
    out_df = out_df.loc[:, [
        "alpha",
        "topology_preserved_pct",
        "hierarchy_among_survivors_pct",
        "joint_success_pct",
        "loss_pct",
    ]]
    out_df.to_csv(os.path.join(DATA_DIR, "alpha_failure_decomposition.csv"), index=False)
    return summary


def run_loss_topology_bias() -> pd.DataFrame:
    print("Running updated-model first-loss topology bias sweep...", flush=True)
    seeds = generate_seeds(TIME_SEEDS, master_seed=SWEEP_MASTER_SEED)
    rows = []
    total = len(BIAS_ALPHAS) * len(seeds)
    counter = 0
    for alpha in BIAS_ALPHAS:
        for seed in seeds:
            counter += 1
            if counter % 200 == 0:
                print(f"  bias runs {counter}/{total}", flush=True)
            rows.append(_run_bias_single((alpha, seed)))

    runs = pd.DataFrame(rows)
    runs.to_csv(os.path.join(DATA_DIR, "loss_topology_bias_runs.csv"), index=False)

    summary_rows = []
    for alpha, grp in runs.groupby("alpha", sort=True):
        n = len(grp)
        n_lost = int((grp["loss_time"] > 0).sum())
        n_prox = int((grp["first_loss_branch"] == 1).sum())
        n_dist = int((grp["first_loss_branch"] == 2).sum())

        total_pct, total_lo, total_hi = wilson_ci(n_lost, n)
        prox_pct, prox_lo, prox_hi = wilson_ci(n_prox, n)
        dist_pct, dist_lo, dist_hi = wilson_ci(n_dist, n)

        prox_share = 100.0 * n_prox / n_lost if n_lost > 0 else np.nan
        dist_share = 100.0 * n_dist / n_lost if n_lost > 0 else np.nan

        # Report both absolute percentages over all runs and conditional shares
        # among only the runs that actually lost a branch.
        summary_rows.append({
            "alpha": round(float(alpha), 2),
            "n_total": n,
            "n_lost": n_lost,
            "n_proximal_first_loss": n_prox,
            "n_distal_first_loss": n_dist,
            "total_loss_pct": round(total_pct, 2),
            "total_loss_ci_lo": round(total_lo, 2),
            "total_loss_ci_hi": round(total_hi, 2),
            "proximal_first_loss_pct": round(prox_pct, 2),
            "proximal_first_loss_ci_lo": round(prox_lo, 2),
            "proximal_first_loss_ci_hi": round(prox_hi, 2),
            "distal_first_loss_pct": round(dist_pct, 2),
            "distal_first_loss_ci_lo": round(dist_lo, 2),
            "distal_first_loss_ci_hi": round(dist_hi, 2),
            "proximal_share_of_losses_pct": round(prox_share, 2),
            "distal_share_of_losses_pct": round(dist_share, 2),
        })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(DATA_DIR, "loss_topology_bias_summary.csv"), index=False)
    return runs


def run_temporal_diagnostics() -> pd.DataFrame:
    print("Running updated-model temporal onset diagnostics...", flush=True)
    seeds = generate_seeds(TIME_SEEDS, master_seed=SWEEP_MASTER_SEED)
    args = [(alpha, seed) for alpha in SELECTED_ALPHAS for seed in seeds]

    imbalance_rows = []
    for idx, arg in enumerate(args, start = 1):
        if idx % 60 == 0:
            print(f"  hierarchy-onset runs {idx}/{len(args)}", flush=True)
        imbalance_rows.append(_run_imbalance_single(arg))
    imbalance_runs = pd.DataFrame(imbalance_rows)
    imbalance_runs.to_csv(os.path.join(DATA_DIR, "time_to_imbalance_runs.csv"), index=False)
    imbalance_curve = cumulative_curve(imbalance_runs, "Hierarchy onset")
    imbalance_curve.to_csv(os.path.join(DATA_DIR, "time_to_imbalance_cumulative.csv"), index=False)

    runaway_rows = []
    for idx, arg in enumerate(args, start = 1):
        if idx % 60 == 0:
            print(f"  one-sided-routing runs {idx}/{len(args)}", flush=True)
        runaway_rows.append(_run_runaway_single(arg))
    runaway_runs = pd.DataFrame(runaway_rows)
    runaway_runs.to_csv(os.path.join(DATA_DIR, "time_to_runaway_runs.csv"), index=False)
    runaway_curve = cumulative_curve(runaway_runs, "One-sided routing")
    runaway_curve.to_csv(os.path.join(DATA_DIR, "time_to_runaway_cumulative.csv"), index=False)

    return pd.concat([imbalance_curve, runaway_curve], ignore_index=True)


def run_survival_from_bias(bias_runs: pd.DataFrame) -> None:
    # Reinterpret stable runs as right-censored at the final simulated step.
    runs = bias_runs.loc[bias_runs["alpha"].isin(SELECTED_ALPHAS)].copy()
    runs["step"] = runs["loss_time"].where(runs["loss_time"] > 0, FINAL_STEP).astype(int)
    runs["event"] = (runs["loss_time"] > 0).astype(int)

    curve = pd.concat([
        kaplan_meier_curve(
            runs.loc[runs["alpha"] == alpha, "step"],
            runs.loc[runs["alpha"] == alpha, "event"],
            alpha,
        )
        for alpha in SELECTED_ALPHAS
    ], ignore_index=True)
    curve.to_csv(os.path.join(DATA_DIR, "loss_survival_curves.csv"), index=False)

    test_rows = []
    for alpha_a, alpha_b in [(0.10, 0.45), (0.45, 1.00), (0.10, 1.00)]:
        chi2, p_value = log_rank_test(runs, alpha_a, alpha_b)
        test_rows.append({
            "alpha_a": alpha_a,
            "alpha_b": alpha_b,
            "chi2": chi2,
            "p_value": p_value,
        })
    pd.DataFrame(test_rows).to_csv(
        os.path.join(DATA_DIR, "loss_survival_tests.csv"), index=False
    )


def main() -> None:
    # The top-level workflow is ordered so foundational exports exist before the
    # more expensive simulation batches start writing their derived summaries.
    ensure_dirs()
    export_segment_coords()

    run_representative_traces()
    seed_df = run_matched_rule_comparison()
    sweep_summary = run_fine_br5_sweep()
    bias_runs = run_loss_topology_bias()
    run_temporal_diagnostics()
    run_survival_from_bias(bias_runs)

    peak_row = sweep_summary.loc[sweep_summary["joint_success_pct"].idxmax()]
    br3_joint = int(seed_df.loc[seed_df["rule"] == "BR3", "joint_success"].sum())
    br4_joint = int(seed_df.loc[seed_df["rule"] == "BR4", "joint_success"].sum())
    br5_joint = int(seed_df.loc[seed_df["rule"] == "BR5", "joint_success"].sum())

    # Final console output gives a quick human-readable digest after the long run.
    print("\nUpdated-model analysis complete.", flush=True)
    print(
        f"  BR3 joint success: {br3_joint}/{N_COMPARE} ({100.0 * br3_joint / N_COMPARE:.1f}%)",
        flush=True,
    )
    print(
        f"  BR4 joint success: {br4_joint}/{N_COMPARE} ({100.0 * br4_joint / N_COMPARE:.1f}%)",
        flush=True,
    )
    print(
        f"  BR5 joint success: {br5_joint}/{N_COMPARE} ({100.0 * br5_joint / N_COMPARE:.1f}%)",
        flush=True,
    )
    print(
        f"  BR5 peak alpha: {peak_row['alpha']:.2f} with joint success {peak_row['joint_success_pct']:.1f}%",
        flush=True,
    )


if __name__ == "__main__":
    main()
