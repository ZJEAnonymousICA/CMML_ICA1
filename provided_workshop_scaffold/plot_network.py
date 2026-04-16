'''
Just add some labels to the network plot for better visualization.
--By Jinliang
'''

import numpy as np
import matplotlib.pyplot as plt

def plot_network(segments, D, P, Q, seg_cells, tau=None):
    """
    Plot the bifurcating vessel network with numerical labels.
    """
    
    plt.figure(figsize=(14, 6))
    
    # === Subplot 1: Network Structure ===
    plt.subplot(1, 2, 1)
    plt.title('Network: Q [blue, pL/s] / Tau [red, Pa]')
    
    def plot_seg(seg_idx, pt1, pt2):
        # Color encodes the sign of the solved flow rate for this segment.
        if Q[seg_idx] > 0:
            color = "red"  # Flow is positive
        else:
            color = "blue"
        
        # Draw the vessel shape
        lw = D[seg_idx] * 1e6 / 2
        plt.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], 
                 color=color, linewidth=lw, alpha=0.6)

        # Draw Text Labels, only label the middle segment of each straight section to avoid clutter
        
        should_label = False
        # Vertical Up (0-4), center is 2
        if seg_idx == 2: should_label = True
        # Horizontal Right (5-14), centers near 9
        if seg_idx == 9: should_label = True
        # Vertical Down (15-19), center is 17
        if seg_idx == 17: should_label = True
        
        # Upper Loop
        if seg_idx == 22: should_label = True
        if seg_idx == 29: should_label = True
        if seg_idx == 37: should_label = True
        
        if should_label:
            mid_x = (pt1[0] + pt2[0]) / 2 + 2 # Offset slightly
            mid_y = (pt1[1] + pt2[1]) / 2
            
            # Format numbers
            # Flow is converted from SI units into pL/s for readability.
            q_val = Q[seg_idx] * 1e12 # Convert to pL/s
            tau_val = tau[seg_idx] if tau is not None else 0
            
            label_text = f"Q:{q_val:.2f}\nT:{tau_val:.2f}"
            
            plt.text(mid_x, mid_y, label_text, 
                     fontsize=8, color='black', 
                     bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # Draw Lower Loop (Seg 0-19)
    for seg in range(20):
        pt1 = segments[seg]
        pt2 = segments[seg + 1]
        plot_seg(seg, pt1, pt2)

    # Draw Upper Loop (Seg 20-39)
    offset = 21 
    for seg in range(20, 40):
        idx = offset + (seg - 20)
        pt1 = segments[idx]
        pt2 = segments[idx + 1]
        plot_seg(seg, pt1, pt2)
    
    plt.grid(True, alpha=0.3)
    plt.axis('equal') 
    plt.xlabel('x (microns)')
    plt.ylabel('y (microns)')

    # === Subplot 2: Polarity ===
    plt.subplot(1, 2, 2)
    plt.title('Distribution of Cell Polarity')
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.grid(True)
    
    theta = np.linspace(0, 2*np.pi, 100)
    plt.plot(np.cos(theta), np.sin(theta), 'k:', alpha=0.3)
    
    for seg in range(len(seg_cells)):
        if seg_cells[seg]['num'] > 0:
            for cell in range(seg_cells[seg]['num']):
                # Each line shows one cell's polarity as a vector from the origin.
                polarity = seg_cells[seg]['polarity'][cell]
                plt.plot([0, polarity[0]], [0, polarity[1]], 'b-', alpha=0.1)
    
    # Save current frame to verify without display
    # plt.savefig('current_frame.png') 
    plt.show()
