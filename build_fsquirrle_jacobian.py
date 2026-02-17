"""
Build script: read symbolic Jacobians and bases from jacobians/symbolic_fsquirrle.py,
then generate efficient numerical functions in jacobians/fsquirrle_jacobian.py.

Run this after editing symbolic_fsquirrle.py (e.g. when robot physics change).
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from jacobians.symbolic_fsquirrle import (
    J_inv_flat,
    B_flat,
    B_u_flat,
    B_l_flat,
)
from jacobians.codegen import exprs_to_numpy_code


def _code_block(exprs):
    return "\n".join(f"        {c}," for c in exprs_to_numpy_code(exprs))

J_inv_lines = _code_block(J_inv_flat)
B_lines = _code_block(B_flat)
B_u_lines = _code_block(B_u_flat)
B_l_lines = _code_block(B_l_flat)

content = f'''"""
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
{J_inv_lines}
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
{B_lines}
    ], dtype=float)
    B = B_flat.reshape((3, 4), order="F")

    B_u_flat = np.array([
{B_u_lines}
    ], dtype=float)
    B_u = B_u_flat.reshape((3, 4), order="F")

    B_l_flat = np.array([
{B_l_lines}
    ], dtype=float)
    B_l = B_l_flat.reshape((3, 4), order="F")

    return B, B_u, B_l
'''

content = content.replace("from np.typing", "from numpy.typing")

out_dir = os.path.join(os.path.dirname(__file__), "jacobians")
out_path = os.path.join(out_dir, "fsquirrle_jacobian.py")
with open(out_path, "w", encoding="utf-8") as f:
    f.write(content)

print(f"Wrote {out_path}")
