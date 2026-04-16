#!/usr/bin/env python3
"""
Run the provided workshop scaffold once (default settings) and record a small set
of comparable diagnostics for the Supporting Materials.

This script does NOT modify the workshop code. It imports the workshop modules and
replays the same 40-step loop while logging summary state at each step.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class WorkshopParams:
    # These defaults mirror the scaffold script so that the benchmark reflects
    # the behaviour reviewers would see if they ran the provided workshop code.
    nt: int = 40
    pin: float = 4 * 98.0
    pout: float = 1 * 98.0
    mu: float = 3.5e-3
    nseg: int = 40
    num_cell: int = 10
    cell_size: float = 5e-6

    branch_rule: int = 3
    branch_alpha: float = 1.0

    # Polarity re-alignment weights (as in abm_ec_simulation_v2.py)
    w2: float = 0.30  # flow component
    w3: float = 0.10  # neighbour alignment
    w4: float = 0.30  # random re-alignment

    @property
    def w1(self) -> float:
        # The persistence weight is whatever mass remains after allocating the
        # explicit flow, neighbour, and random components.
        return 1.0 - self.w2 - self.w3 - self.w4


def _import_workshop_modules(workshop_dir: Path) -> None:
    # The scaffold modules use plain local imports rather than a package
    # structure, so we prepend their directory to `sys.path` temporarily.
    sys.path.insert(0, str(workshop_dir))

    # Import side effects are fine; we only need functions.
    global solve_for_flow, cell_migration, realign_polarity, make_segments
    from solve_for_flow import solve_for_flow  # type: ignore
    from cell_migration import cell_migration  # type: ignore
    from realign_polarity import realign_polarity  # type: ignore
    from make_segments import make_segments  # type: ignore


def initialize_segments(nseg: int, num_cell: int) -> list[dict]:
    # Recreate the scaffold's per-segment mutable state:
    # - `num` tracks occupancy
    # - `polarity` stores one 2D unit vector per cell
    # - `migration` stores a per-cell move flag for the current step
    seg_cells: list[dict] = [{} for _ in range(nseg)]
    for seg in range(nseg):
        seg_cells[seg]["num"] = int(num_cell)
        seg_cells[seg]["polarity"] = [np.random.randn(2) for _ in range(num_cell)]
        for v in seg_cells[seg]["polarity"]:
            # Normalizing here makes each random polarity a direction only,
            # independent of the initial Gaussian draw magnitude.
            v /= np.linalg.norm(v)
        seg_cells[seg]["migration"] = np.zeros(num_cell)
    return seg_cells


def compute_conductance(
    nseg: int,
    ncell: np.ndarray,
    cell_size: float,
    mu: float,
    length_m: np.ndarray,
    min_d: float = 1e-7,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Mirror workshop compute_conductance() including the minimum diameter floor.
    """
    # `d` is segment diameter, `g` is Poiseuille conductance, and `h` is the
    # factor used later to convert flow into wall shear stress.
    d = np.zeros(nseg)
    g = np.zeros(nseg)
    h = np.zeros(nseg)
    for seg in range(nseg):
        # In both scaffold and updated models, diameter is derived from the
        # number of endothelial cells lining the segment circumference.
        calculated_d = ncell[seg] * cell_size / np.pi
        # The workshop scaffold enforces a hard minimum calibre so segments
        # never truly occlude numerically.
        d[seg] = max(calculated_d, min_d)
        g[seg] = (np.pi * d[seg] ** 4) / (128 * mu * length_m[seg])
        h[seg] = (32 * mu) / (np.pi * d[seg] ** 3)
    return d, g, h


def main() -> None:
    root = Path(__file__).resolve().parent
    workshop_dir = root / "provided_workshop_scaffold"
    out_csv = root / "data" / "workshop_scaffold_benchmark_timecourse.csv"

    # Load the original workshop helper functions without modifying them.
    _import_workshop_modules(workshop_dir)

    params = WorkshopParams()
    # A fixed seed makes the benchmark deterministic and therefore comparable
    # across reruns and against the submitted supporting material.
    np.random.seed(123456789)

    # The scaffold uses uniform segment lengths and a uniform initial cell count.
    length_m = np.ones(params.nseg) * 10e-6
    ncell = np.ones(params.nseg) * params.num_cell

    # Geometry is created in the workshop scaffold; included here for parity.
    _segments = make_segments(length_m)  # noqa: F841

    seg_cells = initialize_segments(params.nseg, params.num_cell)

    # Solve the initial haemodynamic state before any remodelling happens.
    d, g, h = compute_conductance(params.nseg, ncell, params.cell_size, params.mu, length_m)
    p, q, tau = solve_for_flow(g, params.pin, params.pout, h)

    records: list[dict] = []

    def record_state(step: int) -> None:
        # These slices reproduce the branch groupings used elsewhere in the
        # project so the scaffold benchmark can be compared directly.
        # Branch labels match the manuscript convention used elsewhere.
        prox_idx = slice(5, 15)   # segments 5-14
        dist_idx = slice(20, 40)  # segments 20-39

        # The scaffold's branch-routing behaviour at the divergent junction can
        # be summarized by the relative shear magnitude in the two daughter arms.
        val_lower = float(abs(tau[5]))
        val_upper = float(abs(tau[20]))
        p_lower = val_lower / (val_lower + val_upper) if (val_lower + val_upper) > 0 else np.nan

        records.append(
            {
                "step": step,
                "total_cells": float(np.sum(ncell)),
                "prox_cells": float(np.sum(ncell[prox_idx])),
                "dist_cells": float(np.sum(ncell[dist_idx])),
                "prox_mean_diam_um": float(np.mean(d[prox_idx]) * 1e6),
                "dist_mean_diam_um": float(np.mean(d[dist_idx]) * 1e6),
                "terminal14_cells": float(ncell[14]),
                "terminal39_cells": float(ncell[39]),
                "p_lower_divergent_tau": p_lower,
            }
        )

    record_state(step=0)

    for t in range(params.nt):
        # `migrate` is a per-segment counter that the scaffold migration routine
        # mutates as it marks cells for movement.
        migrate = np.zeros(params.nseg)

        # The scaffold performs its polarity update into a separate structure so
        # all cells observe the same pre-update neighbourhood state.
        realigned_cells = [dct.copy() for dct in seg_cells]
        for i in range(len(seg_cells)):
            realigned_cells[i]["polarity"] = list(seg_cells[i]["polarity"])
            realigned_cells[i]["migration"] = list(seg_cells[i]["migration"])

        for seg in range(params.nseg):
            # Realignment uses the current flow field and same-step neighbour
            # state, but it does not yet move any cells between segments.
            realign_polarity(seg, q, seg_cells, realigned_cells,
                             params.w1, params.w2, params.w3, params.w4)

        # Migration consumes the newly realigned polarity vectors and returns an
        # updated segment-state list after routing cells through the network.
        seg_cells = cell_migration(realigned_cells, migrate, q,
                                   params.branch_rule, params.branch_alpha, tau)

        for seg in range(params.nseg):
            # The haemodynamic solve only depends on segment occupancy, so we
            # collapse the richer cell state back to a simple cell-count vector.
            ncell[seg] = seg_cells[seg]["num"]

        # Recompute diameters, flows, and shear stress after remodelling.
        d, g, h = compute_conductance(params.nseg, ncell, params.cell_size, params.mu, length_m)
        p, q, tau = solve_for_flow(g, params.pin, params.pout, h)

        record_state(step=t + 1)

    # Export as a plain table so downstream review figures can reuse the same data.
    df = pd.DataFrame.from_records(records)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv} ({len(df)} rows).")


if __name__ == "__main__":
    main()
