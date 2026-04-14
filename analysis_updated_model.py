"""Analysis helpers for the current updated mixed-polarity model."""

import numpy as np
from multiprocessing import Pool, get_all_start_methods, get_context
import importlib as _il

_sim = _il.import_module("simulation_updated_model")
run_simulation = _sim.run_simulation
DEFAULT_PARAMS = _sim.DEFAULT_PARAMS


def generate_seeds(n_seeds=1000, master_seed=42):
    """Generate a list of random seeds for reproducible stochastic runs."""
    rng = np.random.RandomState(master_seed)
    return rng.randint(1, 10**9, size=n_seeds)


def wilson_ci(k, n, z=1.96):
    """Wilson score interval for a binomial proportion, returned in percent."""
    if n == 0:
        return 0.0, 0.0, 0.0
    phat = k / n
    denom = 1 + z**2 / n
    centre = (phat + z**2 / (2 * n)) / denom
    half_width = z * np.sqrt(phat * (1 - phat) / n + z**2 / (4 * n**2)) / denom
    return (
        phat * 100.0,
        max(0.0, centre - half_width) * 100.0,
        min(1.0, centre + half_width) * 100.0,
    )


def _pool_factory(n_workers):
    """Create a Pool, preferring fork on macOS for numpy-heavy workloads."""
    if "fork" in get_all_start_methods():
        return get_context("fork").Pool(processes=n_workers)
    return Pool(processes=n_workers)


def _run_single(args):
    """Worker function for parallel alpha sweep."""
    alpha, seed, params, Nt, branch_rule = args
    p = params.copy()
    p["Nt"] = Nt
    result = run_simulation(branch_rule=branch_rule, alpha=alpha, seed=seed,
                            params=p, record_probabilities=False)
    return {
        "alpha": alpha,
        "seed": seed,
        "loss_time": result["loss_time"],
        "loss_branch": result["bifurcation_lost"][-1],
        "branch_collapse_time": result["branch_collapse_time"],
        "branch_collapse_branch": result["branch_collapse"][-1],
        "final_prox_diam": result["mean_diam_prox"][-1],
        "final_dist_diam": result["mean_diam_dist"][-1],
        "hierarchy_success": int(result["hierarchy_success"]),
        "joint_success": int((result["loss_time"] == 0) and result["hierarchy_success"]),
    }


def alpha_sweep(alpha_values=None, seeds=None, params=None, Nt=36,
                n_workers=None, verbose=True, branch_rule=5):
    """Run exploratory-variant simulations across alpha values and seeds."""
    if alpha_values is None:
        alpha_values = np.arange(0, 1.01, 0.01)
    if seeds is None:
        seeds = generate_seeds(1000)
    if params is None:
        params = DEFAULT_PARAMS.copy()

    n_alpha = len(alpha_values)
    n_seeds = len(seeds)

    args_list = []
    for alpha in alpha_values:
        for seed in seeds:
            args_list.append((alpha, int(seed), params, Nt, branch_rule))

    total = len(args_list)
    if verbose:
        print(
            f"Running updated-model alpha sweep: "
            f"{n_alpha} alpha values x {n_seeds} seeds = {total} simulations"
        )

    if n_workers == 1:
        raw_results = []
        for idx, args in enumerate(args_list):
            if verbose and idx % 1000 == 0:
                print(f"  Progress: {idx}/{total}")
            raw_results.append(_run_single(args))
    else:
        with _pool_factory(n_workers) as pool:
            raw_results = pool.map(_run_single, args_list)

    loss_time = np.zeros((n_alpha, n_seeds), dtype=int)
    loss_branch = np.zeros((n_alpha, n_seeds), dtype=int)
    branch_collapse_time = np.zeros((n_alpha, n_seeds), dtype=int)
    branch_collapse_branch = np.zeros((n_alpha, n_seeds), dtype=int)
    final_prox_diam = np.zeros((n_alpha, n_seeds))
    final_dist_diam = np.zeros((n_alpha, n_seeds))
    hierarchy_success = np.zeros((n_alpha, n_seeds), dtype=int)
    joint_success = np.zeros((n_alpha, n_seeds), dtype=int)

    for idx, res in enumerate(raw_results):
        i = idx // n_seeds
        j = idx % n_seeds
        loss_time[i, j] = res["loss_time"]
        loss_branch[i, j] = res["loss_branch"]
        branch_collapse_time[i, j] = res["branch_collapse_time"]
        branch_collapse_branch[i, j] = res["branch_collapse_branch"]
        final_prox_diam[i, j] = res["final_prox_diam"]
        final_dist_diam[i, j] = res["final_dist_diam"]
        hierarchy_success[i, j] = res["hierarchy_success"]
        joint_success[i, j] = res["joint_success"]

    loss_pct = np.zeros((n_alpha, Nt + 1))
    branch_collapse_pct = np.zeros((n_alpha, Nt + 1))
    for i in range(n_alpha):
        for t in range(1, Nt + 1):
            lost_by_t = np.sum((loss_time[i] > 0) & (loss_time[i] <= t))
            loss_pct[i, t] = 100.0 * lost_by_t / n_seeds
            collapsed_by_t = np.sum(
                (branch_collapse_time[i] > 0) & (branch_collapse_time[i] <= t)
            )
            branch_collapse_pct[i, t] = 100.0 * collapsed_by_t / n_seeds

    return {
        "alpha_values": alpha_values,
        "seeds": seeds,
        "loss_time": loss_time,
        "loss_branch": loss_branch,
        "branch_collapse_time": branch_collapse_time,
        "branch_collapse_branch": branch_collapse_branch,
        "loss_pct": loss_pct,
        "branch_collapse_pct": branch_collapse_pct,
        "final_prox_diam": final_prox_diam,
        "final_dist_diam": final_dist_diam,
        "hierarchy_success": hierarchy_success,
        "joint_success": joint_success,
        "Nt": Nt,
    }


def compute_regression_stats(loss_time, loss_branch):
    """Compute regression statistics for a set of simulations at one alpha value."""
    n = len(loss_time)
    total_lost = np.sum(loss_time > 0)
    prox_lost = np.sum(loss_branch == 1)
    dist_lost = np.sum(loss_branch == 2)
    return {
        "total_pct": 100.0 * total_lost / n,
        "proximal_pct": 100.0 * prox_lost / n,
        "distal_pct": 100.0 * dist_lost / n,
    }


def run_multiple_seeds(branch_rule, n_seeds=10, alpha=0.45, params=None,
                       Nt=None, master_seed=42):
    """Run a bifurcation rule with multiple seeds under workshop-style polarity."""
    if params is None:
        params = DEFAULT_PARAMS.copy()
    if Nt is not None:
        params["Nt"] = Nt

    seeds = generate_seeds(n_seeds, master_seed)
    results = []
    for seed in seeds:
        result = run_simulation(
            branch_rule=branch_rule, alpha=alpha, seed=int(seed), params=params,
            record_probabilities=(branch_rule == 5)
        )
        results.append(result)
    return results
