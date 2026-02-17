"""
Symbolic kinematics and Jacobians for the FSquirrle cable-driven robot.

Edit this file when the robot physics change (anchor layout, cable attachment,
pulley radius, screw pitch, number of cables). After editing, run
build_fsquirrle_jacobian.py to regenerate the efficient numerical functions
in fsquirrle_jacobian.py.
"""

import sympy as sp

# =============================================================================
# State and parameters (edit to match your robot)
# =============================================================================
# R     : robot position (3,) in world frame
# A1,A2,A3 : anchor positions (3,) in world; z-component is set to 0 below
# C1_u, C2_u, C3_u : upper cable attachment points relative to R
# r_p   : pulley radius (cable length rate -> actuator rate)
# p_s   : screw pitch (vertical motion -> 4th actuator rate)

Rx, Ry, Rz = sp.symbols("Rx Ry Rz", real=True)
dRx, dRy, dRz = sp.symbols("dRx dRy dRz", real=True)
R = sp.Matrix([Rx, Ry, Rz])
dR = sp.Matrix([dRx, dRy, dRz])

A1x, A1y = sp.symbols("A1x A1y", real=True)
A2x, A2y = sp.symbols("A2x A2y", real=True)
A3x, A3y = sp.symbols("A3x A3y", real=True)
A1 = sp.Matrix([A1x, A1y, sp.Integer(0)])
A2 = sp.Matrix([A2x, A2y, sp.Integer(0)])
A3 = sp.Matrix([A3x, A3y, sp.Integer(0)])

C1_ux, C1_uy, C1_uz = sp.symbols("C1_ux C1_uy C1_uz", real=True)
C2_ux, C2_uy, C2_uz = sp.symbols("C2_ux C2_uy C2_uz", real=True)
C3_ux, C3_uy, C3_uz = sp.symbols("C3_ux C3_uy C3_uz", real=True)
C1_u = sp.Matrix([C1_ux, C1_uy, C1_uz])
C2_u = sp.Matrix([C2_ux, C2_uy, C2_uz])
C3_u = sp.Matrix([C3_ux, C3_uy, C3_uz])

r_p, p_s = sp.symbols("r_p p_s", real=True, positive=True)

# =============================================================================
# Geometry: global positions and cable vectors
# =============================================================================
# Upper attachment in world: C_u_global = C_u + R
# Lower section at anchor height (z=0): C_l_global = [C_u_global.x, C_u_global.y, 0]

C1_u_global = C1_u + R
C2_u_global = C2_u + R
C3_u_global = C3_u + R
C1_l_global = sp.Matrix([C1_u_global[0], C1_u_global[1], 0])
C2_l_global = sp.Matrix([C2_u_global[0], C2_u_global[1], 0])
C3_l_global = sp.Matrix([C3_u_global[0], C3_u_global[1], 0])

# Cable vectors: from attachment (global) to anchor
l1_u = A1 - C1_u_global
l2_u = A2 - C2_u_global
l3_u = A3 - C3_u_global
l1_l = A1 - C1_l_global
l2_l = A2 - C2_l_global
l3_l = A3 - C3_l_global


def _norm_sym(v):
    return sp.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


# Unit directions
l1_u_dir = l1_u / _norm_sym(l1_u)
l2_u_dir = l2_u / _norm_sym(l2_u)
l3_u_dir = l3_u / _norm_sym(l3_u)
l1_l_dir = l1_l / _norm_sym(l1_l)
l2_l_dir = l2_l / _norm_sym(l2_l)
l3_l_dir = l3_l / _norm_sym(l3_l)

# =============================================================================
# Actuator rates: cable length rate and screw rate
# =============================================================================
# Cable length rate = d(l)/dt = direction · dR; actuator rate = (length rate) / r_p
# Screw: dq4 = dRz / (p_s / (2*pi))

dl1_u = l1_u_dir.dot(dR)
dl2_u = l2_u_dir.dot(dR)
dl3_u = l3_u_dir.dot(dR)
dl1_l = l1_l_dir.dot(dR)
dl2_l = l2_l_dir.dot(dR)
dl3_l = l3_l_dir.dot(dR)
dl1 = dl1_u + dl1_l
dl2 = dl2_u + dl2_l
dl3 = dl3_u + dl3_l

dq1 = dl1 / r_p
dq2 = dl2 / r_p
dq3 = dl3 / r_p
dq4 = dRz / (p_s / (2 * sp.pi))

# =============================================================================
# Inverse Jacobian: dq = J_inv @ dR  (4x3)
# =============================================================================
J_inv_sym = sp.Matrix([
    [dq1.diff(dRx), dq1.diff(dRy), dq1.diff(dRz)],
    [dq2.diff(dRx), dq2.diff(dRy), dq2.diff(dRz)],
    [dq3.diff(dRx), dq3.diff(dRy), dq3.diff(dRz)],
    [dq4.diff(dRx), dq4.diff(dRy), dq4.diff(dRz)],
])
J_inv_sym = sp.simplify(J_inv_sym)

# =============================================================================
# Basis matrices for force allocation (columns = cable/screw directions)
# =============================================================================
B_u_sym = sp.Matrix.hstack(l1_u_dir, l2_u_dir, l3_u_dir, sp.Matrix([0, 0, 0]))
B_l_sym = sp.Matrix.hstack(l1_l_dir, l2_l_dir, l3_l_dir, sp.Matrix([0, 0, 1]))
B_sym = B_l_sym + B_u_sym

# For codegen: flatten column-major
J_inv_flat = [J_inv_sym[i, j] for j in range(3) for i in range(4)]
B_flat = [B_sym[i, j] for j in range(4) for i in range(3)]
B_u_flat = [B_u_sym[i, j] for j in range(4) for i in range(3)]
B_l_flat = [B_l_sym[i, j] for j in range(4) for i in range(3)]

# Input symbols for numerical functions (R, A1, A2, A3, C1_u, C2_u, C3_u, r_p, p_s)
# Scalar names for codegen
_input_scalars = (
    Rx, Ry, Rz,
    A1x, A1y, A2x, A2y, A3x, A3y,
    C1_ux, C1_uy, C1_uz, C2_ux, C2_uy, C2_uz, C3_ux, C3_uy, C3_uz,
    r_p, p_s,
)
