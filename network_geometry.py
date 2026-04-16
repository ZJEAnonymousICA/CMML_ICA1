"""
Network topology and geometry for the A-branch vessel model.

The network consists of 40 segments and 40 nodes arranged as:
- Lower path: Feeding (seg 0-4, up) → Proximal (seg 5-14, right) → Draining (seg 15-19, down)
- Upper path: Distal (seg 20-24, up; seg 25-34, right; seg 35-39, down)

Bifurcation at node 5 (flow-divergent), Convergence at node 15 (flow-convergent).
Inlet at node 0, Outlet at node 20.

Reference: Edgar et al. (2021), PLoS Comput Biol.
"""

import numpy as np

# --- Constants ---
# Segment and node counts are shared across the entire updated model, so this
# module acts as the single source of truth for the network layout.
NSEG = 40
NNODE = 40  # 21 lower-path nodes + 19 upper-path nodes (nodes 5 and 15 shared)

# Segment group ranges (0-indexed)
FEEDING = range(0, 5)       # 5 segments, vertical up
PROXIMAL = range(5, 15)     # 10 segments, horizontal right
DRAINING = range(15, 20)    # 5 segments, vertical down
DISTAL_UP = range(20, 25)   # 5 segments, vertical up
DISTAL_HOR = range(25, 35)  # 10 segments, horizontal right
DISTAL_DOWN = range(35, 40) # 5 segments, vertical down
DISTAL = range(20, 40)      # all 20 distal segments

# Key nodes
NODE_INLET = 0
NODE_OUTLET = 20
NODE_BIFURCATION = 5   # flow-divergent (left side)
NODE_CONVERGENCE = 15  # flow-convergent (right side)


def get_vessel_axis(seg):
    """Return the geometric direction vector for a segment."""
    # Each contiguous straight run in the A-branch network has a fixed axis, so
    # segment identity alone is enough to determine its local coordinate frame.
    if seg in FEEDING or seg in DISTAL_UP:
        return np.array([0.0, 1.0])   # up
    elif seg in PROXIMAL or seg in DISTAL_HOR:
        return np.array([1.0, 0.0])   # right
    elif seg in DRAINING or seg in DISTAL_DOWN:
        return np.array([0.0, -1.0])  # down
    else:
        raise ValueError(f"Invalid segment index: {seg}")


def get_seg_nodes(seg):
    """Return (upstream_node, downstream_node) for a segment following flow direction.

    Flow direction: inlet → feeding → bifurcation → branches → convergence → draining → outlet

    Convention: Q > 0 means flow from upstream_node to downstream_node.
    """
    # The lower path uses consecutive node numbering, while the upper path is
    # spliced onto shared bifurcation/convergence nodes.
    if 0 <= seg <= 19:
        # Lower path: seg i connects node i to node i+1
        return (seg, seg + 1)
    elif seg == 20:
        # First distal segment: node 5 (bifurcation) to node 21
        return (5, 21)
    elif 21 <= seg <= 38:
        # Upper path interior: node seg to node seg+1
        return (seg, seg + 1)
    elif seg == 39:
        # Last distal segment: node 39 to node 15 (convergence)
        return (39, 15)
    else:
        raise ValueError(f"Invalid segment index: {seg}")


def get_downstream_target(seg):
    """Return list of candidate target segments when a cell moves DOWNSTREAM (with flow).

    Downstream = increasing segment number along the flow direction.
    """
    # Special cases
    if seg == 4:
        # Bifurcation: cell can go to proximal (seg 5) or distal (seg 20)
        return [5, 20]
    if seg == 19:
        # Outlet: cell exits (or wraps to inlet via periodic BC)
        return [0]  # periodic BC: wrap to segment 0
    if seg == 39:
        # End of distal branch: enters convergence → draining
        return [15]
    # Default: next segment
    # Outside junctions, downstream movement is just "advance one segment".
    if 0 <= seg <= 18:
        return [seg + 1]
    if 20 <= seg <= 38:
        return [seg + 1]
    return []


def get_upstream_target(seg):
    """Return list of candidate target segments when a cell moves UPSTREAM (against flow).

    This is the primary migration direction for ECs during remodelling.
    """
    # Special cases
    if seg == 0:
        # Inlet: cell wraps to outlet via periodic BC
        return [19]  # periodic BC: wrap to segment 19
    if seg == 5:
        # Start of proximal: goes back to feeding vessel
        return [4]
    if seg == 20:
        # Start of distal: goes back to feeding vessel (bifurcation)
        return [4]
    if seg == 15:
        # Start of draining: cell at convergence can go to proximal (seg 14) or distal (seg 39)
        return [14, 39]
    # Default: previous segment
    # Outside junctions, upstream movement is just "step to the previous segment".
    if 1 <= seg <= 19:
        return [seg - 1]
    if 21 <= seg <= 39:
        return [seg - 1]
    return []


def get_rotation_angle_downstream(seg, target):
    """Return polarity rotation angle (degrees) for a cell moving downstream from seg to target."""
    # These hard-coded turns encode the actual corners in the network so a
    # cell's polarity can stay consistent with its new segment orientation.
    # Horizontal → Vertical down (right-angle turn)
    if seg == 14 and target == 15:
        return -90.0
    if seg == 34 and target == 35:
        return -90.0
    # Vertical up → Horizontal right (right-angle turn at bifurcation)
    if seg == 4 and target == 5:
        return -90.0
    if seg == 24 and target == 25:
        return -90.0
    # Periodic BC wrap (outlet → inlet): 180° turn
    if seg == 19 and target == 0:
        return -180.0
    # Vertical up → Vertical up (no turn, same direction at bifurcation)
    if seg == 4 and target == 20:
        return 0.0
    # Vertical down → Horizontal (convergence to draining continuation)
    if seg == 39 and target == 15:
        return -90.0
    return 0.0


def get_rotation_angle_upstream(seg, target):
    """Return polarity rotation angle (degrees) for a cell moving upstream from seg to target."""
    # The upstream mapping is not simply the negative of the downstream mapping,
    # because the same junction can imply different geometric continuations.
    # Horizontal left ← Vertical down (right-angle turn)
    if seg == 15 and target == 14:
        return 90.0
    if seg == 35 and target == 34:
        return 90.0
    # Vertical down ← Horizontal right (right-angle turn at bifurcation/convergence)
    if seg == 5 and (target == 4 or target == 20):
        return 90.0
    if seg == 25 and target == 24:
        return 90.0
    # Periodic BC wrap (inlet → outlet): 180° turn
    if seg == 0 and target == 19:
        return 180.0
    # Convergence upstream: draining → proximal or distal
    if seg == 15 and target == 39:
        return 0.0  # straight continuation (vertical down → vertical down)
    # Vertical down going up
    if seg == 20 and target == 4:
        return 90.0
    return 0.0


def rotate_polarity(polar_vect, angle_deg):
    """Apply a 2D rotation to a polarity vector."""
    if angle_deg == 0.0:
        return polar_vect
    # A standard 2D rotation matrix is enough because polarity lives in the
    # plane of the schematic vessel network.
    theta = np.radians(angle_deg)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])
    rotated = R @ polar_vect
    norm = np.linalg.norm(rotated)
    if norm > 0:
        # Renormalization removes any accumulated floating-point drift.
        rotated /= norm
    return rotated


def make_segment_coords(L):
    """Generate segment endpoint coordinates for plotting.

    Args:
        L: array of segment lengths (m), length NSEG

    Returns:
        segments: dict with 'lower' and 'upper' coordinate arrays
                  lower: (21, 2) array for nodes 0-20
                  upper: (20, 2) array for nodes 5, 21-39
    """
    # Geometry is stored in micrometres because the coordinates are used
    # primarily for plotting rather than for physical calculations.
    L_um = L * 1e6  # convert to micrometres for plotting

    v_up = np.array([0, 1])
    v_right = np.array([1, 0])
    v_down = np.array([0, -1])

    # Lower path: nodes 0-20
    lower = np.zeros((21, 2))
    # Feeding vessel (seg 0-4): 5 segments going up
    for i in range(5):
        # Each node is constructed cumulatively from the previous node plus one
        # segment-length step in the appropriate direction.
        lower[i+1] = lower[i] + v_up * L_um[i]
    # Proximal branch (seg 5-14): 10 segments going right
    for i in range(10):
        lower[5+i+1] = lower[5+i] + v_right * L_um[5+i]
    # Draining vessel (seg 15-19): 5 segments going down
    for i in range(5):
        lower[15+i+1] = lower[15+i] + v_down * L_um[15+i]

    # Upper path: start from node 5 (bifurcation)
    upper = np.zeros((21, 2))
    upper[0] = lower[5]  # node 5 = bifurcation
    # Distal up (seg 20-24): 5 segments going up
    for i in range(5):
        upper[i+1] = upper[i] + v_up * L_um[20+i]
    # Distal horizontal (seg 25-34): 10 segments going right
    for i in range(10):
        upper[5+i+1] = upper[5+i] + v_right * L_um[25+i]
    # Distal down (seg 35-39): 5 segments going down
    for i in range(5):
        upper[15+i+1] = upper[15+i] + v_down * L_um[35+i]

    return {'lower': lower, 'upper': upper}
