"""
Auto-generated from jacobians/symbolic_fsquirrle.py. Do not edit by hand.
Run build_fsquirrle_jacobian.py to regenerate.
"""

import numpy as np
from numpy.typing import ArrayLike


def fsquirrle_jacobian_inverse(
    R: ArrayLike,
    A1: ArrayLike,
    A2: ArrayLike,
    A3: ArrayLike,
    C1_u: ArrayLike,
    C2_u: ArrayLike,
    C3_u: ArrayLike,
    r_p: float,
    p_s: float,
) -> np.ndarray:
    """
    Inverse Jacobian: dR -> dq (actuator rates). Returns (4, 3) array.
    """
    Rx, Ry, Rz = np.asarray(R, dtype=float)
    A1x, A1y, _ = np.asarray(A1, dtype=float)
    A2x, A2y, _ = np.asarray(A2, dtype=float)
    A3x, A3y, _ = np.asarray(A3, dtype=float)
    C1_ux, C1_uy, C1_uz = np.asarray(C1_u, dtype=float)
    C2_ux, C2_uy, C2_uz = np.asarray(C2_u, dtype=float)
    C3_ux, C3_uy, C3_uz = np.asarray(C3_u, dtype=float)

    J_inv_flat = np.array([
        -((-A1x + C1_ux + Rx)/np.sqrt((C1_uz + Rz)**2 + (-A1x + C1_ux + Rx)**2 + (-A1y + C1_uy + Ry)**2) + (-A1x + C1_ux + Rx)/np.sqrt((-A1x + C1_ux + Rx)**2 + (-A1y + C1_uy + Ry)**2))/r_p,
        -((-A2x + C2_ux + Rx)/np.sqrt((C2_uz + Rz)**2 + (-A2x + C2_ux + Rx)**2 + (-A2y + C2_uy + Ry)**2) + (-A2x + C2_ux + Rx)/np.sqrt((-A2x + C2_ux + Rx)**2 + (-A2y + C2_uy + Ry)**2))/r_p,
        -((-A3x + C3_ux + Rx)/np.sqrt((C3_uz + Rz)**2 + (-A3x + C3_ux + Rx)**2 + (-A3y + C3_uy + Ry)**2) + (-A3x + C3_ux + Rx)/np.sqrt((-A3x + C3_ux + Rx)**2 + (-A3y + C3_uy + Ry)**2))/r_p,
        0,
        -((-A1y + C1_uy + Ry)/np.sqrt((C1_uz + Rz)**2 + (-A1x + C1_ux + Rx)**2 + (-A1y + C1_uy + Ry)**2) + (-A1y + C1_uy + Ry)/np.sqrt((-A1x + C1_ux + Rx)**2 + (-A1y + C1_uy + Ry)**2))/r_p,
        -((-A2y + C2_uy + Ry)/np.sqrt((C2_uz + Rz)**2 + (-A2x + C2_ux + Rx)**2 + (-A2y + C2_uy + Ry)**2) + (-A2y + C2_uy + Ry)/np.sqrt((-A2x + C2_ux + Rx)**2 + (-A2y + C2_uy + Ry)**2))/r_p,
        -((-A3y + C3_uy + Ry)/np.sqrt((C3_uz + Rz)**2 + (-A3x + C3_ux + Rx)**2 + (-A3y + C3_uy + Ry)**2) + (-A3y + C3_uy + Ry)/np.sqrt((-A3x + C3_ux + Rx)**2 + (-A3y + C3_uy + Ry)**2))/r_p,
        0,
        (-C1_uz - Rz)/(r_p*np.sqrt((C1_uz + Rz)**2 + (-A1x + C1_ux + Rx)**2 + (-A1y + C1_uy + Ry)**2)),
        (-C2_uz - Rz)/(r_p*np.sqrt((C2_uz + Rz)**2 + (-A2x + C2_ux + Rx)**2 + (-A2y + C2_uy + Ry)**2)),
        (-C3_uz - Rz)/(r_p*np.sqrt((C3_uz + Rz)**2 + (-A3x + C3_ux + Rx)**2 + (-A3y + C3_uy + Ry)**2)),
        2*np.pi/p_s,
    ], dtype=float)
    return J_inv_flat.reshape((4, 3), order="F")


def fsquirrle_jacobian(
    R: ArrayLike,
    A1: ArrayLike,
    A2: ArrayLike,
    A3: ArrayLike,
    C1_u: ArrayLike,
    C2_u: ArrayLike,
    C3_u: ArrayLike,
    r_p: float,
    p_s: float,
) -> np.ndarray:
    """Forward Jacobian (pseudoinverse of J_inv). Returns (3, 4) array."""
    J_inv = fsquirrle_jacobian_inverse(R, A1, A2, A3, C1_u, C2_u, C3_u, r_p, p_s)
    return np.linalg.pinv(J_inv)


def fsquirrle_bases(
    R: ArrayLike,
    A1: ArrayLike,
    A2: ArrayLike,
    A3: ArrayLike,
    C1_u: ArrayLike,
    C2_u: ArrayLike,
    C3_u: ArrayLike,
):
    """
    Basis matrices for cable forces. Returns (B, B_u, B_l), each (3, 4).
    """
    Rx, Ry, Rz = np.asarray(R, dtype=float)
    A1x, A1y, _ = np.asarray(A1, dtype=float)
    A2x, A2y, _ = np.asarray(A2, dtype=float)
    A3x, A3y, _ = np.asarray(A3, dtype=float)
    C1_ux, C1_uy, C1_uz = np.asarray(C1_u, dtype=float)
    C2_ux, C2_uy, C2_uz = np.asarray(C2_u, dtype=float)
    C3_ux, C3_uy, C3_uz = np.asarray(C3_u, dtype=float)

    B_flat = np.array([
        (A1x - C1_ux - Rx)/np.sqrt((-C1_uz - Rz)**2 + (A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2) + (A1x - C1_ux - Rx)/np.sqrt((A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2),
        (A1y - C1_uy - Ry)/np.sqrt((-C1_uz - Rz)**2 + (A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2) + (A1y - C1_uy - Ry)/np.sqrt((A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2),
        (-C1_uz - Rz)/np.sqrt((-C1_uz - Rz)**2 + (A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2),
        (A2x - C2_ux - Rx)/np.sqrt((-C2_uz - Rz)**2 + (A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2) + (A2x - C2_ux - Rx)/np.sqrt((A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2),
        (A2y - C2_uy - Ry)/np.sqrt((-C2_uz - Rz)**2 + (A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2) + (A2y - C2_uy - Ry)/np.sqrt((A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2),
        (-C2_uz - Rz)/np.sqrt((-C2_uz - Rz)**2 + (A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2),
        (A3x - C3_ux - Rx)/np.sqrt((-C3_uz - Rz)**2 + (A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2) + (A3x - C3_ux - Rx)/np.sqrt((A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2),
        (A3y - C3_uy - Ry)/np.sqrt((-C3_uz - Rz)**2 + (A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2) + (A3y - C3_uy - Ry)/np.sqrt((A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2),
        (-C3_uz - Rz)/np.sqrt((-C3_uz - Rz)**2 + (A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2),
        0,
        0,
        1,
    ], dtype=float)
    B = B_flat.reshape((3, 4), order="F")

    B_u_flat = np.array([
        (A1x - C1_ux - Rx)/np.sqrt((-C1_uz - Rz)**2 + (A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2),
        (A1y - C1_uy - Ry)/np.sqrt((-C1_uz - Rz)**2 + (A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2),
        (-C1_uz - Rz)/np.sqrt((-C1_uz - Rz)**2 + (A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2),
        (A2x - C2_ux - Rx)/np.sqrt((-C2_uz - Rz)**2 + (A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2),
        (A2y - C2_uy - Ry)/np.sqrt((-C2_uz - Rz)**2 + (A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2),
        (-C2_uz - Rz)/np.sqrt((-C2_uz - Rz)**2 + (A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2),
        (A3x - C3_ux - Rx)/np.sqrt((-C3_uz - Rz)**2 + (A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2),
        (A3y - C3_uy - Ry)/np.sqrt((-C3_uz - Rz)**2 + (A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2),
        (-C3_uz - Rz)/np.sqrt((-C3_uz - Rz)**2 + (A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2),
        0,
        0,
        0,
    ], dtype=float)
    B_u = B_u_flat.reshape((3, 4), order="F")

    B_l_flat = np.array([
        (A1x - C1_ux - Rx)/np.sqrt((A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2),
        (A1y - C1_uy - Ry)/np.sqrt((A1x - C1_ux - Rx)**2 + (A1y - C1_uy - Ry)**2),
        0,
        (A2x - C2_ux - Rx)/np.sqrt((A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2),
        (A2y - C2_uy - Ry)/np.sqrt((A2x - C2_ux - Rx)**2 + (A2y - C2_uy - Ry)**2),
        0,
        (A3x - C3_ux - Rx)/np.sqrt((A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2),
        (A3y - C3_uy - Ry)/np.sqrt((A3x - C3_ux - Rx)**2 + (A3y - C3_uy - Ry)**2),
        0,
        0,
        0,
        1,
    ], dtype=float)
    B_l = B_l_flat.reshape((3, 4), order="F")

    return B, B_u, B_l
