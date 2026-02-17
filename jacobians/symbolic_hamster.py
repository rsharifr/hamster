"""
Symbolic kinematics and Jacobian for the HAMSTER omniwheel robot.

Edit this file when the robot physics change (wheel layout, number of wheels,
wheel radius, geometry). After editing, run build_hamster_jacobian.py to
regenerate the efficient numerical Jacobian in hamster_jacobian.py.
"""

import sympy as sp

# =============================================================================
# Parameters (edit to match your robot)
# =============================================================================
# theta : robot orientation in the plane (rad)
# r     : wheel radius
# D     : distance from robot center to each wheel
# Body velocity: v_x, v_y (m/s), Omega (rad/s)

theta = sp.Symbol("theta", real=True)
r = sp.Symbol("r", real=True, positive=True)
D = sp.Symbol("D", real=True, positive=True)
v_x, v_y, Omega = sp.symbols("v_x v_y Omega", real=True)

# =============================================================================
# Wheel geometry (edit if you change the number or arrangement of wheels)
# =============================================================================
# Three wheels at 120°: angles relative to global X are
#   phi1 = pi + theta,   phi2 = 5*pi/3 + theta,   phi3 = pi/3 + theta

phi1_expr = sp.pi + theta
phi2_expr = sp.Integer(5) * sp.pi / 3 + theta
phi3_expr = sp.pi / 3 + theta

# =============================================================================
# Kinematics: wheel angular rate vs body velocity
# =============================================================================
# Each wheel: omega_i = (1/r) * (v_x*cos(phi_i) + v_y*sin(phi_i) + D*Omega)

omega1 = (1 / r) * (v_x * sp.cos(phi1_expr) + v_y * sp.sin(phi1_expr) + D * Omega)
omega2 = (1 / r) * (v_x * sp.cos(phi2_expr) + v_y * sp.sin(phi2_expr) + D * Omega)
omega3 = (1 / r) * (v_x * sp.cos(phi3_expr) + v_y * sp.sin(phi3_expr) + D * Omega)

# =============================================================================
# Inverse Jacobian: maps [v_x, v_y, Omega] -> [omega1, omega2, omega3]
# =============================================================================
J_inv_sym = sp.Matrix([
    [omega1.diff(v_x), omega1.diff(v_y), omega1.diff(Omega)],
    [omega2.diff(v_x), omega2.diff(v_y), omega2.diff(Omega)],
    [omega3.diff(v_x), omega3.diff(v_y), omega3.diff(Omega)],
])
J_inv_sym = sp.simplify(J_inv_sym)

# =============================================================================
# Forward Jacobian: maps [omega1, omega2, omega3] -> [v_x, v_y, Omega]
# =============================================================================
J_sym = sp.simplify(J_inv_sym.inv())

# Inputs for the numerical function (state/parameters)
# Order must match the generated function signature: hamster_jacobian(D, r, theta)
_input_symbols = (D, r, theta)
_output_matrix = J_sym
_output_shape = (3, 3)
