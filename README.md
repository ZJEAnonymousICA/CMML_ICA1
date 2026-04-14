<p align="center">
  <img src="./assets/ZJE_logo.png" alt="ZJU-UoE Institute logo" width="780">
</p>

<h1 align="center">CMML3 ICA1: Improved Workshop's Model and BR Implementation </h1>

<p align="center">
  Reviewer-facing reproduction package (WBM vs SIM).
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.x-3776AB?logo=python&logoColor=white" alt="Python badge">
  <img src="https://img.shields.io/badge/Scope-WBM%20%2B%20SIM-1f883d" alt="Scope badge">
  <img src="https://img.shields.io/badge/Outputs-Data%20%7C%20Figures-6f42c1" alt="Outputs badge">
</p>

---

## Overview

This folder is a compact, reviewer-facing reproduction package. It contains:

- the original workshop base model (WBM) under [`provided_workshop_scaffold/`](./provided_workshop_scaffold/)
- the self-implemented model (SIM) used to generate the precomputed outputs in [`data/`](./data/) and [`figures/`](./figures/)
- minimal scripts to regenerate the WBM benchmark and SIM analysis outputs

No absolute paths are required. All scripts resolve paths relative to this folder.

## Repository Layout

```text
Improved_Model_BRs/
  README.md
  CHANGES.md

  01_run_workshop_benchmark.py        # WBM benchmark run (imports provided_workshop_scaffold/)
  02_run_sim_full_analysis.py         # SIM full analysis (writes data/*.csv)
  03_export_review_tables.py          # exports data/table_s1_scaffold_expansions.csv + manifest

  provided_workshop_scaffold/         # original teaching code (unchanged)

  network_geometry.py                 # A-branch topology + turn-aware polarity rotation
  haemodynamics.py                    # conductance + Kirchhoff flow solve
  polarity_mixed.py                   # workshop-style mixed polarity update
  migration_rules.py                  # BR routing + intercalation (BR1/BR3/BR5 used)
  simulation_support.py               # defaults + endpoints (bifurcation loss, etc.)
  simulation_updated_model.py         # simulation engine (SIM)
  analysis_updated_model.py           # alpha sweep + multi-seed utilities

  data/                               # precomputed and regenerated CSV/JSON outputs
  figures/                            # precomputed Figure1.pdf/Figure2.pdf + panels/ PNGs
  visualization/                      # optional panel-by-panel R scripts (PNG outputs)
  assets/                             # logo
```

## Core Code

| File | Purpose |
| --- | --- |
| [`network_geometry.py`](./network_geometry.py) | A-branch geometry and turn logic |
| [`haemodynamics.py`](./haemodynamics.py) | Kirchhoff flow solve and conductance remodelling |
| [`polarity_mixed.py`](./polarity_mixed.py) | Workshop-style persistence/flow/neighbour/random polarity update |
| [`migration_rules.py`](./migration_rules.py) | BR1-BR5 bifurcation-routing mechanics used by the improved model |
| [`simulation_support.py`](./simulation_support.py) | Shared defaults and endpoint helpers |
| [`simulation_updated_model.py`](./simulation_updated_model.py) | Updated simulation engine |
| [`analysis_updated_model.py`](./analysis_updated_model.py) | Multi-seed and alpha-sweep utilities |

## Output Locations

| Location | Contents |
| --- | --- |
| [`data/`](./data/) | Paper-facing CSV and JSON outputs |
| [`figures/`](./figures/) | Rendered PNG and PDF figure assets |

## What Changed (WBM → SIM)


| Area | Workshop base model (WBM) | SIM implementation (where) | Why it matters |
| --- | --- | --- | --- |
| Branch-choice semantics | Divergent-junction router via `branch_rule` IDs | Convergent-junction decision with explicit BR1/BR3/BR5 rules ([migration_rules.py](./migration_rules.py)) | The decision problem of interest is upstream choice at a *convergent* bifurcation |
| Endpoint definition | No convergent-terminal loss metric or joint-success readout | Loss scored on convergent terminal segments; joint success reported ([simulation_support.py](./simulation_support.py)) | Topology-only stability can miss hierarchy failure |
| True occlusion | Minimum diameter floor prevents complete pruning | Conductance can reach zero when segments empty ([haemodynamics.py](./haemodynamics.py)) | Branch loss must be representable rather than clipped |
| Cell conservation | Cells can effectively leave the system | Outlet-to-inlet recycling preserved ([simulation_updated_model.py](./simulation_updated_model.py)) | Recycling preserves depletion dynamics and competition strength |
| Turn-aware transport | Polarity transported in global coordinates | Segment-aware polarity rotation through corners/junctions ([network_geometry.py](./network_geometry.py)) | Keeps migration direction consistent across network turns |
| Intercalation smoothing | Not implemented | Redistribution away from bifurcations ([migration_rules.py](./migration_rules.py)) | Suppresses non-biological pile-ups that distort calibre updates |
| Polarity realignment | Simple heuristic update | Workshop-style mixed update used in SIM ([polarity_mixed.py](./polarity_mixed.py)) | Tests whether BR feedback persists under the workshop’s richer polarity design |

The same information is also available as a CSV for reviewers who prefer a copy/paste table:
`data/table_s1_scaffold_expansions.csv`.

## File Relationships (Call Graph)

### WBM benchmark

```text
01_run_workshop_benchmark.py
  -> provided_workshop_scaffold/solve_for_flow.py
  -> provided_workshop_scaffold/cell_migration.py
  -> provided_workshop_scaffold/realign_polarity.py
  -> provided_workshop_scaffold/make_segments.py
  -> writes data/workshop_scaffold_benchmark_timecourse.csv
```

### SIM analysis

```text
02_run_sim_full_analysis.py
  -> simulation_updated_model.py
     -> haemodynamics.py
     -> polarity_mixed.py
     -> migration_rules.py
     -> simulation_support.py
     -> network_geometry.py
  -> analysis_updated_model.py
     -> simulation_updated_model.py
  -> writes data/*.csv (SIM panel + mechanism datasets)
```

### Review tables

```text
03_export_review_tables.py
  -> writes data/table_s1_scaffold_expansions.csv
  -> writes data/paper_asset_manifest.csv + data/paper_asset_manifest.json
```

## Quick Start

To reproduce the analysis outputs, run exactly this from `Workshop/Improved_Model_BRs`:

```bash
cd /Users/coellearth/Desktop/CMML3/ICA1/Workshop/Improved_Model_BRs

# 1. WBM benchmark (runs the provided workshop code once and logs a timecourse CSV)
/Users/coellearth/.virtualenvs/r-reticulate/bin/python 01_run_workshop_benchmark.py

# 2. SIM full analysis (writes all SIM datasets into data/*.csv)
/Users/coellearth/.virtualenvs/r-reticulate/bin/python 02_run_sim_full_analysis.py

# 3. Export the WBM→SIM change table and a small asset manifest
/Users/coellearth/.virtualenvs/r-reticulate/bin/python 03_export_review_tables.py
```

## Key Parameters (As Used In The Precomputed Outputs)

### WBM benchmark defaults (01_run_workshop_benchmark.py)

- `nt=40` steps; `pin=4*98`, `pout=1*98`; `mu=3.5e-3`
- `nseg=40`; `num_cell=10`; `cell_size=5e-6`
- `branch_rule=3`; `branch_alpha=1.0`
- polarity weights (workshop heuristic): `w2=0.30` (flow), `w3=0.10` (neighbour), `w4=0.30` (random), `w1=1-w2-w3-w4` (persistence)

### SIM run and analysis defaults (02_run_sim_full_analysis.py)

Core simulation defaults come from `simulation_support.DEFAULT_PARAMS` and the mixed-polarity override in `simulation_updated_model.DEFAULT_PARAMS`:

- geometry: `Nseg=40` segments; `NNODE=40` nodes (see `network_geometry.py`)
- initialisation: `n0=8` cells/segment, random unit polarity
- flow/remodelling: `cell_width=5e-6`, `l_seg=10e-6`, `mu=0.0035`, `Pin=100`, `Pout=0`
- duration: `Nt=36` steps, `dt_hours=10/3` (5 days total)
- mixed polarity weights (SIM): `w1=0.30`, `w2=0.30`, `w3=0.10`, `w4=0.30`
- BR parameters: `alpha` (swept for BR5), `epsilon=0.1`, `beta=0.3`, `gradient_span=2`

The analysis script uses the following fixed seeds and sample sizes (for reproducibility):

- matched comparison: `N_COMPARE=100`, `COMPARE_MASTER_SEED=20260413`
- BR5 fine sweep: `FINE_ALPHAS=0.00..1.00 (step 0.01)`, `SWEEP_SEEDS=50`, `SWEEP_MASTER_SEED=42`
- topology-bias sweep: `BIAS_ALPHAS=0.00..1.00 (step 0.05)`, `TIME_SEEDS=60`
- multiprocessing: `N_WORKERS=min(8, cpu_count)`

## ICA Figures

- Main figures (final, as submitted): [`figures/Figure1.pdf`](./figures/Figure1.pdf), [`figures/Figure2.pdf`](./figures/Figure2.pdf)
- Panel PNGs used to build figures: `figures/panels/*.png`
