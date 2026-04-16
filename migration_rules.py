"""
Cell migration logic for the updated A-branch ABM.

This module carries the BR1-BR5 logic used by the current paper workflow.
Legacy exploratory branches 6-8 remain implemented but are not used by the
current reproduction scripts.

Migration process:
1. Determine migration direction from polarity (threshold: component >= 0.5)
2. Apply intercalation/smoothing (diffusion scheme)
3. Move cells to target segments with polarity rotation at direction changes
4. Apply bifurcation rules at convergent (node 15) and divergent (node 5) junctions
"""

import numpy as np
import importlib as _il

# Centralize all geometry-aware behaviour through `network_geometry` so the
# migration code can reason in terms of named vessel regions and shared helpers.
_net = _il.import_module("network_geometry")
get_vessel_axis = _net.get_vessel_axis
rotate_polarity = _net.rotate_polarity
FEEDING = _net.FEEDING
PROXIMAL = _net.PROXIMAL
DRAINING = _net.DRAINING
DISTAL_UP = _net.DISTAL_UP
DISTAL_HOR = _net.DISTAL_HOR
DISTAL_DOWN = _net.DISTAL_DOWN


def _migration_direction(seg, polar_vect):
    """Determine migration direction for a cell based on polarity alignment.

    Returns:
        +1: downstream (with flow)
        -1: upstream (against flow)
         0: no migration (polarity not sufficiently aligned)
    """
    axis = get_vessel_axis(seg)

    # migrate_vect component along vessel axis
    # In MATLAB: migrate_vect = cell_size * polar_vect, threshold = cell_size/2
    # Equivalent to: polar_component >= 0.5
    component = np.dot(polar_vect, axis)

    # The explanatory block below preserves the reasoning that reconciles the
    # MATLAB coordinate conventions with the current segment-axis convention.
    # Direction meaning depends on vessel orientation relative to flow:
    # - Feeding (up) & Distal_up (up): positive component = downstream (with flow)
    # - Proximal (right) & Distal_hor (right): positive component = downstream
    # - Draining (down) & Distal_down (down): positive component = UPSTREAM (against flow)
    #   because the axis points down but flow goes down too, so going "with axis" = with flow = downstream
    #   Wait — in draining, flow goes DOWN (from node 15 to node 20). The axis is [0,-1].
    #   A cell with polarity aligned with [0,-1] (positive dot product) is pointing down = with flow.
    #   The MATLAB code for draining segments (16-20):
    #     if migrate_vect(2) >= cell_size/2 → migrate = -1 (upstream, meaning UP = against flow)
    #     if migrate_vect(2) <= -cell_size/2 → migrate = +1 (downstream, meaning DOWN = with flow)
    #   So for vertical-down segments, y >= 0.5 → upstream, y <= -0.5 → downstream
    #   In our convention: axis = [0, -1], component = -y. So component >= 0.5 means -y >= 0.5 means y <= -0.5 → downstream.
    #   And component <= -0.5 means y >= 0.5 → upstream.
    #   Actually let me just match the MATLAB exactly.

    if seg in FEEDING or seg in DISTAL_UP:
        # Vertical up: y >= 0.5 → downstream (+1), y <= -0.5 → upstream (-1)
        y = polar_vect[1]
        if y >= 0.5:
            return 1
        elif y <= -0.5:
            return -1
    elif seg in PROXIMAL or seg in DISTAL_HOR:
        # Horizontal right: x >= 0.5 → downstream (+1), x <= -0.5 → upstream (-1)
        x = polar_vect[0]
        if x >= 0.5:
            return 1
        elif x <= -0.5:
            return -1
    elif seg in DRAINING or seg in DISTAL_DOWN:
        # Vertical down: y >= 0.5 → upstream (-1), y <= -0.5 → downstream (+1)
        # (MATLAB: migrate_vect(2) >= cell_size/2 → seg_cells{seg,3} = -1)
        y = polar_vect[1]
        if y >= 0.5:
            return -1
        elif y <= -0.5:
            return 1
    return 0


def _apply_bifurcation_rule_downstream(seg, cell_vect, Q, tau, seg_cells,
                                        branch_rule, alpha, epsilon=0.1,
                                        tau_ema=None, gradient_span=2):
    """Apply bifurcation rule for a cell moving DOWNSTREAM at seg 4 (flow-divergent).

    Cell at seg 4 can go to:
      seg 5 (proximal, MATLAB seg 6) or seg 20 (distal, MATLAB seg 21)
    """
    if seg != 4:
        # This helper is only meaningful at the divergent branch point.
        return None

    if branch_rule == 1:
        # BR1: higher flow magnitude
        if abs(Q[5]) > abs(Q[20]):
            return 5
        else:
            return 20

    elif branch_rule == 2:
        # BR2: least direction change (polarity direction)
        # Compare dot product with each branch's direction
        # Proximal direction: [1, 0] (right), Distal direction: [0, 1] (up)
        if np.dot(cell_vect, [0, 1]) > np.dot(cell_vect, [1, 0]):
            return 20  # distal (continues up)
        else:
            return 5   # proximal (turns right)

    elif branch_rule == 3:
        # BR3: random 50/50
        if np.random.rand() < 0.5:
            return 20
        else:
            return 5

    elif branch_rule == 4:
        # BR4: biased random (P_proximal=0.7, P_distal=0.3)
        if np.random.rand() < 0.3:
            return 20
        else:
            return 5

    elif branch_rule == 5:
        # BR5: mechanistic (tau-based, MATLAB case 7)
        # Sum tau and cell counts over 3 segments near bifurcation
        tau1 = tau[5] + tau[6] + tau[7]     # proximal
        tau2 = tau[20] + tau[21] + tau[22]   # distal
        tau_sum = tau1 + tau2
        if tau_sum > 0:
            tau_ratio = tau2 / tau_sum
        else:
            tau_ratio = 0.5

        n1 = seg_cells[5]['num'] + seg_cells[6]['num'] + seg_cells[7]['num']
        n2 = seg_cells[20]['num'] + seg_cells[21]['num'] + seg_cells[22]['num']
        n_sum = n1 + n2
        if n_sum > 0:
            n_ratio = n2 / n_sum
        else:
            n_ratio = 0.5

        # `P2` is the probability of choosing the distal branch; the proximal
        # probability is its complement.
        P2 = alpha * tau_ratio + (1 - alpha) * n_ratio

        if np.random.rand() < P2:
            return 20  # distal
        else:
            return 5   # proximal

    elif branch_rule == 6:
        # BR6: Gradient-Sensing Rule (NOVEL)
        # Compute spatial WSS gradient near the bifurcation.
        # gradient_span controls how many segments apart the measurement points are.
        # Positive gradient = WSS increasing deeper into the branch = attractive.
        grad1 = tau[5 + gradient_span] - tau[5]    # proximal
        grad2 = tau[20 + gradient_span] - tau[20]  # distal
        g1 = max(grad1, 0)  # only positive gradients attract
        g2 = max(grad2, 0)
        g_sum = g1 + g2
        grad_ratio = g2 / g_sum if g_sum > 0 else 0.5

        n1 = seg_cells[5]['num'] + seg_cells[6]['num'] + seg_cells[7]['num']
        n2 = seg_cells[20]['num'] + seg_cells[21]['num'] + seg_cells[22]['num']
        n_sum = n1 + n2
        n_ratio = n2 / n_sum if n_sum > 0 else 0.5

        P2 = alpha * grad_ratio + (1 - alpha) * n_ratio
        if np.random.rand() < P2:
            return 20
        else:
            return 5

    elif branch_rule == 7:
        # BR7: Threshold Discrimination Rule (NOVEL)
        # Compute BR5 probabilities, then apply sensory deadband: if the
        # probability difference |P1 - P2| is below epsilon, cells cannot
        # distinguish branches and default to random choice (P = 0.5).
        tau1 = tau[5] + tau[6] + tau[7]
        tau2 = tau[20] + tau[21] + tau[22]
        tau_sum = tau1 + tau2
        tau_ratio = tau2 / tau_sum if tau_sum > 0 else 0.5

        n1 = seg_cells[5]['num'] + seg_cells[6]['num'] + seg_cells[7]['num']
        n2 = seg_cells[20]['num'] + seg_cells[21]['num'] + seg_cells[22]['num']
        n_sum = n1 + n2
        n_ratio = n2 / n_sum if n_sum > 0 else 0.5

        P2_raw = alpha * tau_ratio + (1 - alpha) * n_ratio
        P1_raw = 1.0 - P2_raw
        # Apply sensory threshold
        P2 = 0.5 if abs(P1_raw - P2_raw) < epsilon else P2_raw

        if np.random.rand() < P2:
            return 20
        else:
            return 5

    elif branch_rule == 8:
        # BR8: Temporal Integration Rule (NOVEL)
        # Use time-averaged (EMA) WSS instead of instantaneous WSS.
        # tau_ema is updated in the simulation loop: tau_ema = beta*tau + (1-beta)*tau_ema
        tau_src = tau_ema if tau_ema is not None else tau
        tau1 = tau_src[5] + tau_src[6] + tau_src[7]
        tau2 = tau_src[20] + tau_src[21] + tau_src[22]
        tau_sum = tau1 + tau2
        tau_ratio = tau2 / tau_sum if tau_sum > 0 else 0.5

        n1 = seg_cells[5]['num'] + seg_cells[6]['num'] + seg_cells[7]['num']
        n2 = seg_cells[20]['num'] + seg_cells[21]['num'] + seg_cells[22]['num']
        n_sum = n1 + n2
        n_ratio = n2 / n_sum if n_sum > 0 else 0.5

        P2 = alpha * tau_ratio + (1 - alpha) * n_ratio
        if np.random.rand() < P2:
            return 20
        else:
            return 5

    return 5  # default: proximal


def _apply_bifurcation_rule_upstream(seg, cell_vect, Q, tau, seg_cells,
                                      branch_rule, alpha, epsilon=0.1,
                                      tau_ema=None, gradient_span=2):
    """Apply bifurcation rule for a cell moving UPSTREAM at seg 15 (flow-convergent).

    Cell at seg 15 can go to:
      seg 14 (proximal end, MATLAB seg 15) or seg 39 (distal end, MATLAB seg 40)

    This is the PRIMARY decision point discussed in the paper.
    """
    if seg != 15:
        # This upstream decision exists only at the convergent junction.
        return None

    if branch_rule == 1:
        # BR1: higher flow magnitude
        if abs(Q[14]) > abs(Q[39]):
            return 14
        else:
            return 39

    elif branch_rule == 2:
        # BR2: least direction change (polarity direction)
        # From draining (vertical down), migrating up:
        # Proximal end (seg 14) is horizontal → cell needs to turn left ([-1, 0])
        # Distal end (seg 39) is vertical down → cell continues up ([0, 1])
        if np.dot(cell_vect, [0, 1]) > np.dot(cell_vect, [-1, 0]):
            return 39  # distal (continues up, least direction change)
        else:
            return 14  # proximal (turns left)

    elif branch_rule == 3:
        # BR3: random 50/50
        if np.random.rand() < 0.5:
            return 39
        else:
            return 14

    elif branch_rule == 4:
        # BR4: biased (P_proximal=0.7, P_distal=0.3)
        if np.random.rand() < 0.3:
            return 39
        else:
            return 14

    elif branch_rule == 5:
        # BR5: mechanistic (tau-based, MATLAB case 7)
        # Sum over 3 segments nearest convergence in each branch
        # MATLAB (1-indexed): tau(15)+tau(14)+tau(13) and tau(40)+tau(39)+tau(38)
        # Python (0-indexed): tau[14]+tau[13]+tau[12] and tau[39]+tau[38]+tau[37]
        tau1 = tau[14] + tau[13] + tau[12]  # proximal
        tau2 = tau[39] + tau[38] + tau[37]  # distal
        tau_sum = tau1 + tau2
        if tau_sum > 0:
            tau_ratio = tau2 / tau_sum
        else:
            tau_ratio = 0.5

        # MATLAB (1-indexed): seg_cells{15,1}+{14,1}+{13,1} and {40,1}+{39,1}+{38,1}
        # Python (0-indexed): seg 14,13,12 and 39,38,37
        n1 = seg_cells[14]['num'] + seg_cells[13]['num'] + seg_cells[12]['num']
        n2 = seg_cells[39]['num'] + seg_cells[38]['num'] + seg_cells[37]['num']
        n_sum = n1 + n2
        if n_sum > 0:
            n_ratio = n2 / n_sum
        else:
            n_ratio = 0.5

        # Again, `P2` represents the distal-choice probability at this junction.
        P2 = alpha * tau_ratio + (1 - alpha) * n_ratio

        if np.random.rand() < P2:
            return 39  # distal
        else:
            return 14  # proximal

    elif branch_rule == 6:
        # BR6: Gradient-Sensing Rule (NOVEL) — upstream version
        grad1 = tau[14 - gradient_span] - tau[14]  # proximal
        grad2 = tau[39 - gradient_span] - tau[39]  # distal
        g1 = max(grad1, 0)
        g2 = max(grad2, 0)
        g_sum = g1 + g2
        grad_ratio = g2 / g_sum if g_sum > 0 else 0.5

        n1 = seg_cells[14]['num'] + seg_cells[13]['num'] + seg_cells[12]['num']
        n2 = seg_cells[39]['num'] + seg_cells[38]['num'] + seg_cells[37]['num']
        n_sum = n1 + n2
        n_ratio = n2 / n_sum if n_sum > 0 else 0.5

        P2 = alpha * grad_ratio + (1 - alpha) * n_ratio
        if np.random.rand() < P2:
            return 39
        else:
            return 14

    elif branch_rule == 7:
        # BR7: Threshold Discrimination Rule (NOVEL) — upstream version
        tau1 = tau[14] + tau[13] + tau[12]
        tau2 = tau[39] + tau[38] + tau[37]
        tau_sum = tau1 + tau2
        tau_ratio = tau2 / tau_sum if tau_sum > 0 else 0.5

        n1 = seg_cells[14]['num'] + seg_cells[13]['num'] + seg_cells[12]['num']
        n2 = seg_cells[39]['num'] + seg_cells[38]['num'] + seg_cells[37]['num']
        n_sum = n1 + n2
        n_ratio = n2 / n_sum if n_sum > 0 else 0.5

        P2_raw = alpha * tau_ratio + (1 - alpha) * n_ratio
        P1_raw = 1.0 - P2_raw
        P2 = 0.5 if abs(P1_raw - P2_raw) < epsilon else P2_raw

        if np.random.rand() < P2:
            return 39
        else:
            return 14

    elif branch_rule == 8:
        # BR8: Temporal Integration Rule (NOVEL) — upstream version
        tau_src = tau_ema if tau_ema is not None else tau
        tau1 = tau_src[14] + tau_src[13] + tau_src[12]
        tau2 = tau_src[39] + tau_src[38] + tau_src[37]
        tau_sum = tau1 + tau2
        tau_ratio = tau2 / tau_sum if tau_sum > 0 else 0.5

        n1 = seg_cells[14]['num'] + seg_cells[13]['num'] + seg_cells[12]['num']
        n2 = seg_cells[39]['num'] + seg_cells[38]['num'] + seg_cells[37]['num']
        n_sum = n1 + n2
        n_ratio = n2 / n_sum if n_sum > 0 else 0.5

        P2 = alpha * tau_ratio + (1 - alpha) * n_ratio
        if np.random.rand() < P2:
            return 39
        else:
            return 14

    return 14  # default


def cell_migration(seg, seg_cells, new_seg_cells, migrate, Q, tau,
                   branch_rule, alpha, epsilon=0.1, tau_ema=None,
                   gradient_span=2):
    """Perform cell migration for one segment.

    Modifies new_seg_cells in place. Follows MATLAB cell_migration_Ub1.m exactly.

    Args:
        seg: current segment index (0-indexed)
        seg_cells: current cell state (list of dicts, read-only)
        new_seg_cells: target state being built (modified in place)
        migrate: migration counter array (NSEG,)
        Q: flow rates (NSEG,)
        tau: wall shear stress (NSEG,)
        branch_rule: 1–5
        alpha: BR5 weight parameter

    Returns:
        probabilities: dict with BR5 probability components (or None)
    """
    ncells = seg_cells[seg]['num']
    br5_probs = None

    if ncells == 0:
        return br5_probs

    # Step 1: Determine migration flags for each cell
    for cell_idx in range(ncells):
        polar_vect = np.array(seg_cells[seg]['polarity'][cell_idx], dtype=float)
        direction = _migration_direction(seg, polar_vect)
        # The current-step migration decision is stored alongside the cell so
        # later logic can inspect or suppress it before transport happens.
        seg_cells[seg]['migration'][cell_idx] = direction
        if direction != 0:
            migrate[seg] += 1

    # Step 2: Intercalation / diffusion smoothing
    # Find the "downstream" neighbour for smoothing
    down_seg = None
    if seg not in (19, 39):
        down_seg = seg + 1
    elif seg == 19:
        down_seg = 0  # periodic wrap
    elif seg == 39:
        # Choose the branch with fewer cells
        if seg_cells[14]['num'] < seg_cells[15]['num'] and seg_cells[14]['num'] != 0:
            down_seg = 14
        else:
            down_seg = 15

    # Special case for bifurcation segment
    if seg == 4:
        if seg_cells[5]['num'] < seg_cells[20]['num'] and seg_cells[5]['num'] != 0:
            down_seg = 5
        else:
            down_seg = 20

    if down_seg is not None and ncells > 0:
        if seg_cells[down_seg]['num'] < seg_cells[seg]['num'] and seg_cells[down_seg]['num'] != 0:
            # Randomly cancel one cell's migration
            random_cell = np.random.randint(0, ncells)
            seg_cells[seg]['migration'][random_cell] = 0
            migrate[seg] -= 1

    # Step 3: Move cells
    if migrate[seg] <= 0:
        # If every candidate move was suppressed, the segment contributes no transport this step.
        return br5_probs

    for cell_idx in range(ncells):
        direction = seg_cells[seg]['migration'][cell_idx]

        if direction == 1:
            # --- DOWNSTREAM movement ---
            # Remove the cell from the source segment immediately in the target structure.
            new_seg_cells[seg]['num'] -= 1
            cell_vect = np.array(seg_cells[seg]['polarity'][cell_idx], dtype=float)
            new_seg_cells[seg]['polarity'][cell_idx] = [0.0, 0.0]  # placeholder

            # Determine target
            if seg <= 18:
                target = seg + 1
            elif seg == 19:
                target = 0  # periodic BC
            elif seg <= 38:
                target = seg + 1
            elif seg == 39:
                target = 15

            # Apply bifurcation rule at flow-divergent bifurcation (seg 4)
            if seg == 4:
                result = _apply_bifurcation_rule_downstream(
                    seg, cell_vect, Q, tau, seg_cells, branch_rule, alpha,
                    epsilon=epsilon, tau_ema=tau_ema,
                    gradient_span=gradient_span)
                if result is not None:
                    target = result

            # Polarity rotation at direction changes
            # Rotations are discrete because the schematic network turns only in
            # right angles or in the periodic wrap-around.
            if seg == 4 and target == 5:
                cell_vect = rotate_polarity(cell_vect, -90)
            elif seg == 4 and target == 20:
                pass  # no rotation (same vertical direction)
            elif seg == 24 and target == 25:
                cell_vect = rotate_polarity(cell_vect, -90)
            elif seg == 14 and (target == 15 or target == 39):
                cell_vect = rotate_polarity(cell_vect, -90)
            elif seg == 34 and target == 35:
                cell_vect = rotate_polarity(cell_vect, -90)
            elif seg == 19 and target == 0:
                cell_vect = rotate_polarity(cell_vect, -180)

            # Move cell to target
            new_seg_cells[target]['num'] += 1
            new_seg_cells[target]['polarity'].append(cell_vect.tolist())

        elif direction == -1:
            # --- UPSTREAM movement ---
            # Upstream movement mirrors the same bookkeeping but uses different
            # routing logic at the convergent junction.
            new_seg_cells[seg]['num'] -= 1
            cell_vect = np.array(seg_cells[seg]['polarity'][cell_idx], dtype=float)
            new_seg_cells[seg]['polarity'][cell_idx] = [0.0, 0.0]  # placeholder

            # Determine target
            if seg >= 1 and seg <= 19:
                target = seg - 1
            elif seg == 0:
                target = 19  # periodic BC
            elif seg >= 21 and seg <= 39:
                target = seg - 1
            elif seg == 20:
                target = 4  # back to feeding

            # Apply bifurcation rule at flow-convergent bifurcation (seg 15)
            if seg == 15:
                result = _apply_bifurcation_rule_upstream(
                    seg, cell_vect, Q, tau, seg_cells, branch_rule, alpha,
                    epsilon=epsilon, tau_ema=tau_ema,
                    gradient_span=gradient_span)
                if result is not None:
                    target = result

            # Also handle upstream bifurcation at seg 5 (flow-divergent, cells merging)
            # MATLAB: seg==6 going upstream → target=5 (or 21 via branch rule)
            if seg == 5:
                # Cells from proximal going upstream enter feeding vessel
                target = 4
            if seg == 20:
                # Cells from distal going upstream enter feeding vessel
                target = 4

            # Polarity rotation at direction changes
            if (seg == 15 or seg == 39) and target == 14:
                cell_vect = rotate_polarity(cell_vect, 90)
            elif seg == 35 and target == 34:
                cell_vect = rotate_polarity(cell_vect, 90)
            elif seg == 25 and target == 24:
                cell_vect = rotate_polarity(cell_vect, 90)
            elif seg == 5 and (target == 4 or target == 20):
                cell_vect = rotate_polarity(cell_vect, 90)
            elif seg == 0 and target == 19:
                cell_vect = rotate_polarity(cell_vect, 180)
            elif seg == 20 and target == 4:
                cell_vect = rotate_polarity(cell_vect, 90)

            # Move cell to target
            new_seg_cells[target]['num'] += 1
            new_seg_cells[target]['polarity'].append(cell_vect.tolist())

    return br5_probs


def compute_br5_probabilities(seg_cells, tau, alpha):
    """Compute BR5 probability components for tracking.

    Returns dict with:
        P_tau1, P_tau2: shear stress probability components
        P_n1, P_n2: cell number probability components
        P1, P2: total branch probabilities
    For both the flow-convergent (upstream) and flow-divergent (downstream) bifurcations.
    """
    # Flow-convergent bifurcation (the main one in the paper)
    # Proximal (near convergence): segs 14, 13, 12
    # Distal (near convergence): segs 39, 38, 37
    tau1_conv = tau[14] + tau[13] + tau[12]
    tau2_conv = tau[39] + tau[38] + tau[37]
    tau_sum = tau1_conv + tau2_conv

    n1_conv = seg_cells[14]['num'] + seg_cells[13]['num'] + seg_cells[12]['num']
    n2_conv = seg_cells[39]['num'] + seg_cells[38]['num'] + seg_cells[37]['num']
    n_sum = n1_conv + n2_conv

    if tau_sum > 0:
        P_tau1 = tau1_conv / tau_sum
        P_tau2 = tau2_conv / tau_sum
    else:
        P_tau1 = P_tau2 = 0.5

    if n_sum > 0:
        P_n1 = n1_conv / n_sum
        P_n2 = n2_conv / n_sum
    else:
        P_n1 = P_n2 = 0.5

    # The final branch probabilities are convex combinations of haemodynamic
    # evidence and local occupancy, weighted by `alpha`.
    P1 = alpha * P_tau1 + (1 - alpha) * P_n1
    P2 = alpha * P_tau2 + (1 - alpha) * P_n2

    return {
        'P_tau1': P_tau1, 'P_tau2': P_tau2,
        'P_n1': P_n1, 'P_n2': P_n2,
        'P1': P1, 'P2': P2,
        'tau1': tau1_conv, 'tau2': tau2_conv,
        'n1': n1_conv, 'n2': n2_conv,
    }
