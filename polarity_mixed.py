"""
Workshop-style mixed polarity for the current updated ABM.

This module changes only the polarity update law. The rest of the model keeps
the paper-style geometry, flow solve, migration rules, recycling, and endpoint
definitions.
"""

import numpy as np
import importlib as _il

_net = _il.import_module("network_geometry")
get_vessel_axis = _net.get_vessel_axis


def realign_polarity(seg, Q, seg_cells, new_seg_cells, w1, w2, w3, w4):
    """Realign cell polarity vectors using the workshop weighted-sum heuristic.

    Args:
        seg: segment index
        Q: flow array (NSEG,)
        seg_cells: current cell state (list of dicts)
        new_seg_cells: updated cell state (modified in place)
        w1: persistence weight
        w2: flow-alignment weight
        w3: same-segment neighbour-alignment weight
        w4: random-walk weight
    """
    ncells = seg_cells[seg]["num"]
    if ncells == 0:
        return

    axis = get_vessel_axis(seg)
    flow_dir = np.sign(Q[seg]) * axis

    neigh_dir = np.array([0.0, 0.0])
    if w3 > 0:
        all_pols = np.array(seg_cells[seg]["polarity"], dtype=float)
        sum_pol = np.sum(all_pols, axis=0)
        norm_sum = np.linalg.norm(sum_pol)
        if norm_sum > 0:
            neigh_dir = sum_pol / norm_sum

    for cell_idx in range(ncells):
        old_pol = np.array(seg_cells[seg]["polarity"][cell_idx], dtype=float)

        rand_dir = np.array([0.0, 0.0])
        if w4 > 0:
            rand_dir = np.random.randn(2)
            norm_rand = np.linalg.norm(rand_dir)
            if norm_rand > 0:
                rand_dir /= norm_rand

        new_pol = w1 * old_pol + w2 * flow_dir + w3 * neigh_dir + w4 * rand_dir

        norm_new = np.linalg.norm(new_pol)
        if norm_new > 0:
            new_pol /= norm_new
        else:
            new_pol = old_pol

        seg_cells[seg]["polarity"][cell_idx] = new_pol.tolist()
        new_seg_cells[seg]["polarity"][cell_idx] = new_pol.tolist()
