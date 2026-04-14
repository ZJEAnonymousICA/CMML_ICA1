# What Changed Relative To The Workshop Base Model (WBM)

The original workshop teaching code is preserved in [`provided_workshop_scaffold/`](./provided_workshop_scaffold/).

This repository implements a self-implemented model (SIM) on top of that base, so the convergent-bifurcation branching-rule comparison is well-defined and reproducible.

The same content is also exported as a CSV: `data/table_s1_scaffold_expansions.csv`.

| Area | Workshop base model (WBM) | SIM implementation (where) | Why it matters |
| --- | --- | --- | --- |
| Branch-choice semantics | Divergent-junction router via `branch_rule` IDs | Convergent-junction decision with explicit BR rules ([migration_rules.py](./migration_rules.py)) | The decision problem of interest is upstream choice at a *convergent* bifurcation |
| Endpoint definition | No convergent-terminal loss metric or joint-success readout | Loss scored on convergent terminal segments; joint success reported ([simulation_support.py](./simulation_support.py)) | Topology-only stability can miss hierarchy failure |
| True occlusion | Minimum diameter floor prevents complete pruning | Conductance can reach zero when segments empty ([haemodynamics.py](./haemodynamics.py)) | Branch loss must be representable rather than clipped |
| Cell conservation | Cells can effectively leave the system | Outlet-to-inlet recycling preserved ([simulation_updated_model.py](./simulation_updated_model.py)) | Recycling preserves depletion dynamics and competition strength |
| Turn-aware transport | Polarity transported in global coordinates | Segment-aware polarity rotation through corners/junctions ([network_geometry.py](./network_geometry.py)) | Keeps migration direction consistent across network turns |
| Intercalation smoothing | Not implemented | Redistribution away from bifurcations ([migration_rules.py](./migration_rules.py)) | Suppresses non-biological pile-ups that distort calibre updates |
| Polarity realignment | Simple heuristic update | Workshop-style mixed update used in SIM ([polarity_mixed.py](./polarity_mixed.py)) | Tests whether BR feedback persists under the workshop’s richer polarity design |

