'''
This file has been changed a lot, and the relied functions from other files are also modified. 
Main changes in this file include:
    1. Add a minimum diameter to prevent total occlusion, if diameter equal to 0, use minimum diameter.
    2. Modify the main time stepping loop to separate realign and migration steps.
--By Jinliang
'''

import numpy as np
import matplotlib.pyplot as plt
from solve_for_flow import solve_for_flow
from cell_migration import cell_migration
from realign_polarity import realign_polarity
from plot_network import plot_network
from make_segments import make_segments

# Set random seed for reproducibility
np.random.seed(123456789)

# Input parameters
Nt = 40  # Number of time steps
Pin = 4 * 98  # Inlet pressure (Pa)
Pout = 1 * 98  # Outlet pressure (Pa)

mu = 3.5e-3  # Dynamic viscosity of blood (Pa-s)
Nn = 40  # Number of nodes
Nseg = 40  # Number of segments
num_cell = 10  # Initial number of cells per segment
cell_size = 5e-6  # Size of each cell (m)

branch_rule = 3  # Branching rule - new
branch_alpha = 1.0  # Branching parameter

# Polarization re-alignment weights
w2 = 0.30  # Flow component weight
w3 = 0.10  # Neighbor re-alignment weight
w4 = 0.30  # Random re-alignment weight
w1 = 1 - w2 - w3 - w4  # Persistence component

# Initialize segment properties
L = np.ones(Nseg) * 10e-6  # Segment lengths (m)
Ncell = np.ones(Nseg) * num_cell  # Segment cell number array
D = np.zeros(Nseg)  # Segment diameters (m)
G = np.zeros(Nseg)  # Segment conductance array (m^4/Pa-s-m)
H = np.zeros(Nseg)  # Shear stress calculation factor

tau = np.zeros(Nseg)  # Shear stress array
segments = make_segments(L)  # Generate segment structure

# Initialize segment cell structures
def initialize_segments(Nseg, num_cell):
    seg_cells = [{} for _ in range(Nseg)]
    for seg in range(Nseg):
        seg_cells[seg]['num'] = int(num_cell)  # Number of cells
        seg_cells[seg]['polarity'] = [np.random.randn(2) for _ in range(num_cell)]  # Random polarity vectors
        for v in seg_cells[seg]['polarity']:
            v /= np.linalg.norm(v)  # Normalize to unit vectors
        seg_cells[seg]['migration'] = np.zeros(num_cell)  # Migration indicator
    return seg_cells

seg_cells = initialize_segments(Nseg, num_cell)

print("Property of the segment 0:", seg_cells[0]) # To be honest, you shoud use a class to store these properties.

# Compute initial segment conductance and shear stress
def compute_conductance(Nseg, Ncell, cell_size, mu, L):
    D = np.zeros(Nseg)
    G = np.zeros(Nseg)
    H = np.zeros(Nseg)
    min_D = 1e-7 # Minimum diameter (0.1 micron) to prevent total occlusion
    
    for seg in range(Nseg):
        # Calculate diameter based on cell count
        # Modified to allow non-integer effective Ncell or just using the min_D
        calculated_D = Ncell[seg] * cell_size / np.pi
        
        if calculated_D < min_D:
            D[seg] = min_D
        else:
            D[seg] = calculated_D
            
        G[seg] = (np.pi * D[seg]**4) / (128 * mu * L[seg])
        H[seg] = (32 * mu) / (np.pi * D[seg]**3)
    return D, G, H

D, G, H = compute_conductance(Nseg, Ncell, cell_size, mu, L)

# Solve for initial flow
P, Q, tau = solve_for_flow(G, Pin, Pout, H)

plot_network(segments, D, P, Q, seg_cells, tau)

# Time stepping for migration process
for t in range(Nt):
    print(f'Time step {t+1}/{Nt}')
    
    migrate = np.zeros(Nseg)
    new_seg_cells = [entry.copy() for entry in seg_cells] # Deep copy might be safer here
    
    # 1. Realign Polarity (rotate in origin)
    realigned_cells = [d.copy() for d in seg_cells] # Shallow copies of dicts
    # Deep copy the lists inside
    for i in range(len(seg_cells)):
        realigned_cells[i]['polarity'] = list(seg_cells[i]['polarity']) 
        realigned_cells[i]['migration'] = list(seg_cells[i]['migration'])

    for seg in range(Nseg):
        realign_polarity(seg, Q, seg_cells, realigned_cells, w1, w2, w3, w4)
    
    # 2. Cell Migration
    final_cells = cell_migration(realigned_cells, migrate, Q, branch_rule, branch_alpha, tau)
    
    # 3. Update State
    seg_cells = final_cells
    
    # Update Ncell
    for seg in range(Nseg):
        Ncell[seg] = seg_cells[seg]['num']

    # Update conductance and shear stress with NEW Ncell
    D, G, H = compute_conductance(Nseg, Ncell, cell_size, mu, L)
    
    # Solve for flow with updated parameters
    P, Q, tau = solve_for_flow(G, Pin, Pout, H)
    
    # Plot only every 20 time steps
    if (t + 1) % 20 == 0:
        plot_network(segments, D, P, Q, seg_cells, tau)