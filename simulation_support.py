"""Shared helpers and defaults for the updated mixed-polarity ABM."""

from __future__ import annotations

import numpy as np


DEFAULT_PARAMS = {
    "Nseg": 40,
    "n0": 8,
    "cell_width": 5e-6,
    "l_seg": 10e-6,
    "mu": 0.0035,
    "Pin": 100.0,
    "Pout": 0.0,
    "w1": 0.0,
    "w2": 1.0,
    "w3": 0.0,
    "w4": 0.0,
    "Nt": 36,
    "dt_hours": 10 / 3,
    "epsilon": 0.1,
    "beta": 0.3,
    "gradient_span": 2,
}


def initialize_cells(nseg: int, n0: int) -> list[dict]:
    """Create the initial per-segment cell state with random unit polarity."""
    seg_cells = []
    for _seg in range(nseg):
        polarities = []
        for _ in range(n0):
            vec = np.random.randn(2)
            vec /= np.linalg.norm(vec)
            polarities.append(vec.tolist())
        seg_cells.append(
            {
                "num": n0,
                "polarity": polarities,
                "migration": [0] * n0,
            }
        )
    return seg_cells


def deep_copy_cells(seg_cells: list[dict]) -> list[dict]:
    """Clone the mutable cell-state structure."""
    copied = []
    for seg_data in seg_cells:
        copied.append(
            {
                "num": seg_data["num"],
                "polarity": [p[:] for p in seg_data["polarity"]],
                "migration": list(seg_data["migration"]),
            }
        )
    return copied


def compute_mean_diameter(diam_um: np.ndarray, branch: str = "proximal") -> float:
    """Compute mean diameter over the branch subsets used in the paper."""
    if branch == "proximal":
        return float(np.mean(diam_um[5:15]))
    if branch == "distal":
        return float(np.mean(diam_um[20:30]))
    if branch == "feeding":
        return float(np.mean(diam_um[0:5]))
    if branch == "draining":
        return float(np.mean(diam_um[15:20]))
    return float(np.mean(diam_um))


def check_bifurcation_loss(seg_cells: list[dict]) -> int:
    """Score Edgar et al.'s terminal-segment branch-loss endpoint."""
    if seg_cells[14]["num"] == 0:
        return 1
    if seg_cells[39]["num"] == 0:
        return 2
    return 0


def check_branch_collapse(seg_cells: list[dict]) -> int:
    """Score the stricter whole-branch collapse endpoint."""
    for seg in range(5, 15):
        if seg_cells[seg]["num"] == 0:
            return 1
    for seg in range(20, 40):
        if seg_cells[seg]["num"] == 0:
            return 2
    return 0
