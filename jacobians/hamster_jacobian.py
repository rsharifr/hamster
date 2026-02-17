"""
Auto-generated from jacobians/symbolic_hamster.py. Do not edit by hand.
Run build_hamster_jacobian.py to regenerate.
"""

import numpy as np


def hamster_jacobian(D: float, r: float, theta: float) -> np.ndarray:
    """
    Forward Jacobian: body velocity [v_x, v_y, Omega] -> wheel rates [omega1, omega2, omega3].
    Parameters: D = distance center to wheel, r = wheel radius, theta = robot orientation (rad).
    Returns (3, 3) array.
    """
    J_flat = np.array([
        -2/3*r*np.cos(theta),
        -2/3*r*np.sin(theta),
        (1/3)*r/D,
        (2/3)*r*np.sin(theta + (1/6)*np.pi),
        -2/3*r*np.cos(theta + (1/6)*np.pi),
        (1/3)*r/D,
        (1/3)*r*(-np.sqrt(3)*np.sin(theta) + np.cos(theta)),
        (2/3)*r*np.sin(theta + (1/3)*np.pi),
        (1/3)*r/D,
    ], dtype=float)
    return J_flat.reshape((3, 3), order="F")
