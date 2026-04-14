"""
Realign the polarity vectors of all cells in a segment based on weight factors.
Correct the direction logic by updating segment ranges.
Uses vector addition for robust alignment.
"""

import numpy as np

def realign_polarity(seg, Q, seg_cells, realigned_cells, w1, w2, w3, w4):

    if seg_cells[seg]['num'] != 0:
        
        # Calculate average polarity of the segment (Neighbor Alignment)
        neigh_dir = np.array([0., 0.])
        if w3 > 0:
            # Sum all polarity vectors in this segment
            # seg_cells[seg]['polarity'] is a list of lists/arrays
            all_pols = np.array(seg_cells[seg]['polarity'])
            sum_pol = np.sum(all_pols, axis=0)
            norm_sum = np.linalg.norm(sum_pol)
            if norm_sum > 0:
                neigh_dir = sum_pol / norm_sum

        for i in range(seg_cells[seg]['num']):
            
            # 1. Current Polarity (Persistence)
            old_pol = np.array(seg_cells[seg]['polarity'][i])
            
            # 2. Flow Direction
            # Map segment ID to geometric direction (Consistent with make_segments.py)
            flow_dir = np.array([0., 0.])
            
            if (0 <= seg <= 4) or (20 <= seg <= 24):
                flow_dir = np.array([0., 1.]) # Up
            elif (5 <= seg <= 14) or (25 <= seg <= 34):
                flow_dir = np.array([1., 0.]) # Right
            elif (15 <= seg <= 19) or (35 <= seg <= 39):
                flow_dir = np.array([0., -1.]) # Down
            
            # Adjust for flow direction
            if Q[seg] < 0:
                flow_dir = -flow_dir
                
            # 3. Random Walk
            rand_dir = np.random.randn(2)
            norm_rand = np.linalg.norm(rand_dir)
            if norm_rand > 0:
                rand_dir /= norm_rand
                
            # 4. Weighted Vector Sum (Added neighbor component)
            new_pol = w1 * old_pol + w2 * flow_dir + w3 * neigh_dir + w4 * rand_dir
            
            # Normalize
            norm_new = np.linalg.norm(new_pol)
            if norm_new > 0:
                new_pol /= norm_new
            else:
                new_pol = old_pol
            
            # Update realigned_cells structure
            realigned_cells[seg]['polarity'][i] = list(new_pol)
            
    return realigned_cells
