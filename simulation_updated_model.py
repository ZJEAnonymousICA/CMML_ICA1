"""Updated mixed-polarity ABM used in the current ICA paper."""

import numpy as np
import importlib as _il

_flow = _il.import_module("haemodynamics")
solve_for_flow = _flow.solve_for_flow
compute_conductance = _flow.compute_conductance

_pol = _il.import_module("polarity_mixed")
realign_polarity = _pol.realign_polarity

_mig = _il.import_module("migration_rules")
cell_migration = _mig.cell_migration
compute_br5_probabilities = _mig.compute_br5_probabilities

_base = _il.import_module("simulation_support")
initialize_cells = _base.initialize_cells
deep_copy_cells = _base.deep_copy_cells
compute_mean_diameter = _base.compute_mean_diameter
check_bifurcation_loss = _base.check_bifurcation_loss
check_branch_collapse = _base.check_branch_collapse
DEFAULT_PARAMS = _base.DEFAULT_PARAMS.copy()

# Keep all model defaults except the mixed-polarity weights.
DEFAULT_PARAMS.update({
    "w1": 0.30,  # persistence
    "w2": 0.30,  # flow alignment (with flow)
    "w3": 0.10,  # same-segment neighbour alignment
    "w4": 0.30,  # random walk
})

_net = _il.import_module("network_geometry")
NNODE = _net.NNODE


def run_simulation(branch_rule, alpha=0.45, seed=None, params=None,
                   record_probabilities=False, snapshot_steps=None,
                   verbose=False):
    """Run the updated BR1-BR5 model with mixed polarity."""
    if params is None:
        params = DEFAULT_PARAMS.copy()

    if seed is not None:
        np.random.seed(seed)

    Nseg = params["Nseg"]
    n0 = params["n0"]
    cell_width = params["cell_width"]
    l_seg = params["l_seg"]
    mu = params["mu"]
    Pin = params["Pin"]
    Pout = params["Pout"]
    Nt = params["Nt"]
    dt_hours = params["dt_hours"]
    w1, w2, w3, w4 = params["w1"], params["w2"], params["w3"], params["w4"]
    epsilon = params.get("epsilon", 0.1)
    beta = params.get("beta", 0.3)
    gradient_span = params.get("gradient_span", 2)

    seg_cells = initialize_cells(Nseg, n0)

    Ncell = np.array([seg_cells[seg]["num"] for seg in range(Nseg)], dtype=float)
    D, G, H = compute_conductance(Ncell, cell_width, mu, l_seg)
    P, Q, tau = solve_for_flow(G, Pin, Pout, H)

    tau_ema = tau.copy() if branch_rule == 8 else None

    new_seg_cells = deep_copy_cells(seg_cells)
    for seg in range(Nseg):
        realign_polarity(seg, Q, seg_cells, new_seg_cells, w1, w2, w3, w4)

    time_days = np.zeros(Nt + 1)
    Ncell_history = np.zeros((Nt + 1, Nseg))
    D_history = np.zeros((Nt + 1, Nseg))
    Q_history = np.zeros((Nt + 1, Nseg))
    tau_history = np.zeros((Nt + 1, Nseg))
    P_history = np.zeros((Nt + 1, NNODE))
    mean_diam_prox = np.zeros(Nt + 1)
    mean_diam_dist = np.zeros(Nt + 1)
    bifurcation_lost = np.zeros(Nt + 1, dtype=int)
    branch_collapse = np.zeros(Nt + 1, dtype=int)

    Ncell_history[0] = Ncell
    D_history[0] = D * 1e6
    Q_history[0] = Q
    tau_history[0] = tau
    P_history[0] = P
    mean_diam_prox[0] = compute_mean_diameter(D * 1e6, "proximal")
    mean_diam_dist[0] = compute_mean_diameter(D * 1e6, "distal")

    prob_history = [] if record_probabilities else None
    snapshot_steps = set(snapshot_steps or [])
    cell_snapshots = {} if snapshot_steps else None
    if record_probabilities:
        prob_history.append(compute_br5_probabilities(seg_cells, tau, alpha))
    if cell_snapshots is not None and 0 in snapshot_steps:
        cell_snapshots[0] = deep_copy_cells(seg_cells)

    loss_time = 0
    branch_collapse_time = 0

    for t in range(1, Nt + 1):
        if verbose:
            print(f"Time step {t}/{Nt}")

        time_days[t] = t * dt_hours / 24.0

        migrate = np.zeros(Nseg)
        new_seg_cells = deep_copy_cells(seg_cells)

        for seg in range(Nseg):
            realign_polarity(seg, Q, seg_cells, new_seg_cells, w1, w2, w3, w4)
            cell_migration(
                seg, seg_cells, new_seg_cells, migrate, Q, tau,
                branch_rule, alpha, epsilon=epsilon,
                tau_ema=tau_ema, gradient_span=gradient_span
            )

            new_polarity = []
            for p in new_seg_cells[seg]["polarity"]:
                pv = np.array(p)
                if np.linalg.norm(pv) > 0:
                    new_polarity.append(p)
            new_seg_cells[seg]["polarity"] = new_polarity

        seg_cells = new_seg_cells
        for seg in range(Nseg):
            Ncell[seg] = seg_cells[seg]["num"]
            seg_cells[seg]["migration"] = [0] * int(Ncell[seg])

        D, G, H = compute_conductance(Ncell, cell_width, mu, l_seg)
        P, Q, tau = solve_for_flow(G, Pin, Pout, H)

        if branch_rule == 8 and tau_ema is not None:
            tau_ema = beta * tau + (1 - beta) * tau_ema

        Ncell_history[t] = Ncell.copy()
        D_history[t] = D * 1e6
        Q_history[t] = Q
        tau_history[t] = tau
        P_history[t] = P
        mean_diam_prox[t] = compute_mean_diameter(D * 1e6, "proximal")
        mean_diam_dist[t] = compute_mean_diameter(D * 1e6, "distal")

        loss = check_bifurcation_loss(seg_cells)
        bifurcation_lost[t] = loss
        if loss > 0 and loss_time == 0:
            loss_time = t

        collapse = check_branch_collapse(seg_cells)
        branch_collapse[t] = collapse
        if collapse > 0 and branch_collapse_time == 0:
            branch_collapse_time = t

        if record_probabilities:
            prob_history.append(compute_br5_probabilities(seg_cells, tau, alpha))
        if cell_snapshots is not None and t in snapshot_steps:
            cell_snapshots[t] = deep_copy_cells(seg_cells)

    results = {
        "time_days": time_days,
        "Ncell": Ncell_history,
        "D": D_history,
        "Q": Q_history,
        "tau": tau_history,
        "P": P_history,
        "mean_diam_prox": mean_diam_prox,
        "mean_diam_dist": mean_diam_dist,
        "bifurcation_lost": bifurcation_lost,
        "branch_collapse": branch_collapse,
        "loss_time": loss_time,
        "branch_collapse_time": branch_collapse_time,
        "branch_rule": branch_rule,
        "alpha": alpha,
        "seed": seed,
        "params": params,
        "hierarchy_success": bool(mean_diam_prox[-1] > mean_diam_dist[-1]),
    }

    if record_probabilities:
        results["prob_history"] = prob_history
    if cell_snapshots is not None:
        results["cell_snapshots"] = cell_snapshots

    return results
