'''
Totally reconstructed.
    1. Write the moving cells into the corresponding segments instead of simply marking them.
    2. Provided the cells with a navigation, enabling them to know the corresponding downstream routes.
    3. Cells move according to the flow or shear force.
--By Jinliang
'''

import numpy as np

def cell_migration(seg_cells, migrate, Q, branch_rule, branch_alpha, tau, cell_size=5e-6):
    """
    Handle cellular migration with full topological routing and physical transport.
    
    Args:
        seg_cells (list): Current state of cells in each segment.
        migrate (array): Migration flags (mostly for statistics).
        Q (array): Flow rates.
        branch_rule (int): Rule ID for migration logic.
        branch_alpha (float): Threshold parameter.
        tau (array): Shear stress.
        cell_size (float): Size of cell (m).
        
    Returns:
        new_seg_cells (list): Updated cellular state after migration.
    """
    
    Nseg = len(seg_cells)
    
    # Initialize empty structure for the NEXT time step
    # We cannot modify in-place because a cell moving from 4 to 5 shouldn't be processed again in 5.
    new_seg_cells = [{'num': 0, 'polarity': [], 'migration': []} for _ in range(Nseg)]
    
    # Probability of migration event (can be tuned)
    # In this scaffold rewrite, most cells do not attempt to move on a given step.
    mchance = 0.01
    
    # === 1. Define Topology / Routing Table ===
    # Helper to find neighbors
    def get_downstream_candidates(current_seg):
        # Downstream follows the imposed vessel orientation used throughout the scaffold.
        # Specific Bifurcation Logic
        if current_seg == 4: return [5, 20]     # Bifurcation: Lower(5), Upper(20)
        if current_seg == 14: return [15]       # Lower Reunion
        if current_seg == 39: return [15]       # Upper Reunion
        if current_seg == 19: return [-1]       # Exit (Outlet)
        
        # General Linear Logic
        # Lower Path
        if 0 <= current_seg < 19: return [current_seg + 1]
        # Upper Path
        if 20 <= current_seg < 39: return [current_seg + 1]
        
        # Fallback (should not happen)
        return []

    def get_upstream_candidates(current_seg):
        # Upstream reverses those same topological connections.
        # Specific Logic
        if current_seg == 5: return [4]         # From Lower back to Main
        if current_seg == 20: return [4]        # From Upper back to Main
        if current_seg == 15: return [14, 39]   # Reunion: Back to Lower(14) or Upper(39)
        if current_seg == 0: return [-1]        # Exit (Inlet)

        # General Linear Logic
        if 0 < current_seg <= 19: return [current_seg - 1]
        if 20 < current_seg <= 39: return [current_seg - 1]
        
        return []

    # === 2. Process Every Segment ===
    for seg in range(Nseg):
        current_data = seg_cells[seg]
        num_cells = current_data['num']
        
        if num_cells == 0:
            continue
            
        for i in range(num_cells):
            # Extract cell properties
            pol = current_data['polarity'][i] # Vector [x, y]
            
            # --- Decision Logic (Simplified from v1/v2 concepts) ---
            # Try to migrate based on polarity orientation relative to "Flow Direction"
            # But here we simplify: If vector Y > threshold -> move Up/Down? 
            # Better: Project polarity onto vessel axis.
            
            # Determine Vessel Axis Direction (Geometric)
            # This is hardcoded based on make_segments logic
            axis_vec = np.array([0, 0])
            if (0 <= seg <= 4) or (20 <= seg <= 24): axis_vec = np.array([0, 1]) # Up
            elif (5 <= seg <= 14) or (25 <= seg <= 34): axis_vec = np.array([1, 0]) # Right
            elif (15 <= seg <= 19) or (35 <= seg <= 39): axis_vec = np.array([0, -1]) # Down
            
            # Dot product to check alignment
            # alignment > 0: Downstream
            # alignment < 0: Upstream
            alignment = np.dot(pol, axis_vec)
            
            move_decision = 0 # 0: Stay, 1: Downstream, -1: Upstream
            
            # Randomized migration check
            if np.random.rand() <= mchance:
                # Threshold for movement (e.g. projection > 0.5 cell size)
                # Normalized vector dot product is just cos(theta). 
                # Let's say if aligned enough, we move.
                if alignment > 0.5: 
                    move_decision = 1
                elif alignment < -0.5:
                    move_decision = -1
            
            # --- Routing Logic ---
            target_seg = seg # Default: Stay
            
            if move_decision == 1: # Downstream
                candidates = get_downstream_candidates(seg)
                if candidates == [-1]: # Exit system
                    target_seg = -1
                elif len(candidates) == 1:
                    target_seg = candidates[0]
                elif len(candidates) > 1:
                    # Bifurcation Decision (Node 5 -> 5 or 20)
                    # The scaffold exposes several simple routing heuristics
                    # that choose between the lower and upper branches.
                    prob_lower = 0.5 # Default Rule 1: Random

                    # Identify indices for lower vs upper branch start
                    # Based on structure: Lower starts at 5, Upper starts at 20
                    idx_lower = 5
                    idx_upper = 20

                    if branch_rule == 2 or branch_rule == 6: # Flow based (Q)
                        val_lower = abs(Q[idx_lower])
                        val_upper = abs(Q[idx_upper]) if Nseg > idx_upper else 0
                        total = val_lower + val_upper
                        if total > 0:
                            prob_lower = val_lower / total
                            
                    elif branch_rule == 3: # Shear Stress based (tau)
                        val_lower = abs(tau[idx_lower])
                        val_upper = abs(tau[idx_upper]) if Nseg > idx_upper else 0
                        total = val_lower + val_upper
                        if total > 0:
                            prob_lower = val_lower / total
                    
                    # Apply decision
                    if np.random.rand() < prob_lower:
                        target_seg = idx_lower # Lower path start
                    else:
                        target_seg = idx_upper # Upper path start
            
            elif move_decision == -1: # Upstream
                candidates = get_upstream_candidates(seg)
                if candidates == [-1]: # Exit system
                    target_seg = -1
                elif len(candidates) == 1:
                    target_seg = candidates[0]
                elif len(candidates) > 1:
                    # Reunion Backwards (Node 15 -> 14 or 39)
                    # Random choice for now
                    target_seg = np.random.choice(candidates)
            
            # --- Execute Transport ---
            if target_seg != -1:
                # Append to the NEW list of the target segment
                # Note: We need to handle Polarity Rotation if vessel changes direction
                # For simplicity now: Keep polarity global (world coordinates)
                # But in reality, cells align to local vessel. 
                # Let's just copy it for now.
                
                new_seg_cells[target_seg]['num'] += 1
                new_seg_cells[target_seg]['polarity'].append(pol)
                new_seg_cells[target_seg]['migration'].append(0) # Reset flag
            else:
                # Cell exited the system (Pruned/Lost)
                pass

    return new_seg_cells
