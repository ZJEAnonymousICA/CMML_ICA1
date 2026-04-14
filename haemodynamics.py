"""
Pressure, flow, and wall shear stress solver for the A-branch vessel network.

Translates MATLAB solve_for_flow_Ub1.m to Python.
Uses Kirchhoff's current law at each node to assemble a linear system,
then solves for nodal pressures and computes segment flows and WSS.

Node indexing (0-based, MATLAB is 1-based):
  Lower path: nodes 0(inlet)–20(outlet), segments 0–19
  Upper path: nodes 21–39, segments 20–39
  Bifurcation: node 5 (segs 4, 5, 20 meet)
  Convergence: node 15 (segs 14, 15, 39 meet)
"""

import numpy as np

NNODE = 40
NSEG = 40


def solve_for_flow(G, Pin, Pout, H):
    """Solve for pressures, flows, and wall shear stress in the network.

    Args:
        G: conductance array (NSEG,) in m^3/(Pa·s)
        Pin: inlet pressure (Pa)
        Pout: outlet pressure (Pa)
        H: shear stress coefficient array (NSEG,), H = 32μ/(πD³)

    Returns:
        P: nodal pressures (NNODE,) in Pa
        Q: segment flow rates (NSEG,) in m³/s
        tau: wall shear stress (NSEG,) in Pa
    """
    G = G.copy()
    G[G == 0] = 1e-25  # prevent singular matrix

    C = np.zeros((NNODE, NNODE))
    B = np.zeros(NNODE)

    # --- Node 0 (inlet boundary): P[0] = Pin ---
    C[0, 0] = G[0]
    B[0] = G[0] * Pin

    # --- Node 1: adjacent to inlet ---
    # G[0]*(P[0]-P[1]) = G[1]*(P[1]-P[2])
    # Since P[0]=Pin is known, move to RHS:
    # (G[0]+G[1])*P[1] - G[1]*P[2] = G[0]*Pin
    C[1, 1] = G[0] + G[1]
    C[1, 2] = -G[1]
    B[1] = G[0] * Pin

    # --- Nodes 2–4: standard interior nodes on feeding vessel ---
    for node in range(2, 5):
        seg_in = node - 1   # segment entering this node
        seg_out = node       # segment leaving this node
        C[node, node-1] = -G[seg_in]
        C[node, node] = G[seg_in] + G[seg_out]
        C[node, node+1] = -G[seg_out]

    # --- Node 5: bifurcation ---
    # Seg 4 (in from node 4), Seg 5 (out to node 6), Seg 20 (out to node 21)
    C[5, 4] = -G[4]
    C[5, 5] = G[4] + G[5] + G[20]
    C[5, 6] = -G[5]
    C[5, 21] = -G[20]

    # --- Nodes 6–14: interior nodes on proximal branch ---
    for node in range(6, 15):
        seg_in = node - 1
        seg_out = node
        C[node, node-1] = -G[seg_in]
        C[node, node] = G[seg_in] + G[seg_out]
        C[node, node+1] = -G[seg_out]

    # --- Node 15: convergence ---
    # Seg 14 (in from node 14), Seg 39 (in from node 39), Seg 15 (out to node 16)
    C[15, 14] = -G[14]
    C[15, 15] = G[14] + G[15] + G[39]
    C[15, 16] = -G[15]
    C[15, 39] = -G[39]

    # --- Nodes 16–18: interior nodes on draining vessel ---
    for node in range(16, 19):
        seg_in = node - 1
        seg_out = node
        C[node, node-1] = -G[seg_in]
        C[node, node] = G[seg_in] + G[seg_out]
        C[node, node+1] = -G[seg_out]

    # --- Node 19: adjacent to outlet ---
    # G[18]*(P[18]-P[19]) = G[19]*(P[19]-P[20])
    # P[20]=Pout is known, move to RHS:
    C[19, 18] = -G[18]
    C[19, 19] = G[18] + G[19]
    B[19] = G[19] * Pout

    # --- Node 20: outlet boundary: P[20] = Pout ---
    C[20, 20] = G[19]
    B[20] = G[19] * Pout

    # --- Node 21: first upper-path node ---
    # Seg 20 comes from node 5, goes to node 21
    # Seg 21 goes from node 21 to node 22
    C[21, 5] = -G[20]
    C[21, 21] = G[20] + G[21]
    C[21, 22] = -G[21]

    # --- Nodes 22–38: interior nodes on upper path ---
    for node in range(22, 39):
        seg_in = node - 1   # segments 21–37
        seg_out = node       # segments 22–38
        C[node, node-1] = -G[seg_in]
        C[node, node] = G[seg_in] + G[seg_out]
        C[node, node+1] = -G[seg_out]

    # --- Node 39: last upper-path node ---
    # Seg 38 (in from node 38), Seg 39 (out to node 15)
    C[39, 15] = -G[39]
    C[39, 38] = -G[38]
    C[39, 39] = G[38] + G[39]

    # Solve the linear system
    P = np.linalg.solve(C, B)

    # --- Calculate flow rates ---
    Q = np.zeros(NSEG)

    # Lower path segments 0–19: seg i connects node i → node i+1
    for seg in range(20):
        n_up, n_dn = seg, seg + 1
        Q[seg] = -G[seg] * (P[n_dn] - P[n_up])

    # Seg 20: node 5 → node 21
    Q[20] = -G[20] * (P[21] - P[5])

    # Segs 21–38: node seg → node seg+1
    for seg in range(21, 39):
        Q[seg] = -G[seg] * (P[seg+1] - P[seg])

    # Seg 39: node 39 → node 15
    Q[39] = -G[39] * (P[15] - P[39])

    # --- Calculate wall shear stress ---
    tau = H * Q

    return P, Q, tau


def compute_conductance(Ncell, cell_width, mu, l_seg):
    """Compute diameter, conductance, and shear stress coefficient from cell counts.

    Args:
        Ncell: cell counts per segment (NSEG,)
        cell_width: lateral width of an EC (m), w in paper
        mu: dynamic viscosity (Pa·s)
        l_seg: segment length (m)

    Returns:
        D: lumen diameters (NSEG,) in m
        G: conductances (NSEG,) in m³/(Pa·s)
        H: shear stress coefficients (NSEG,)
    """
    D = np.zeros(NSEG)
    G = np.zeros(NSEG)
    H = np.zeros(NSEG)

    for seg in range(NSEG):
        if Ncell[seg] >= 1:
            D[seg] = Ncell[seg] * cell_width / np.pi
        else:
            D[seg] = 0.0

        if D[seg] > 0:
            G[seg] = np.pi * D[seg]**4 / (128 * mu * l_seg)
            H[seg] = 32 * mu / (np.pi * D[seg]**3)
        else:
            G[seg] = 0.0
            H[seg] = 0.0

    return D, G, H
