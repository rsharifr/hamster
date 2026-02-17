"""
FSquirrle cable-driven robot motor selection and workspace analysis.

Uses generated Jacobians and bases from jacobians (see build_fsquirrle_jacobian.py).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import lsq_linear

from jacobians import fsquirrle_jacobian_inverse, fsquirrle_bases


def fibonacci_sphere(n: int) -> np.ndarray:
    """
    Python port of `fibonacci_sphere.m`.

    Returns
    -------
    points : (n, 3) ndarray
        Approximately uniform unit vectors on the sphere.
    """
    points = np.zeros((n, 3), dtype=float)
    phi = np.pi * (np.sqrt(5.0) - 1.0)  # golden angle

    if n == 1:
        points[0] = np.array([0.0, 0.0, 1.0])
        return points

    for i in range(n):
        z = 1.0 - (i / (n - 1)) * 2.0
        radius = np.sqrt(max(0.0, 1.0 - z * z))
        theta = phi * i
        x = np.cos(theta) * radius
        y = np.sin(theta) * radius
        points[i] = np.array([x, y, z])

    return points


def run_fsquirrle_motor_selection(
    v: float = 1.0,
    F: float = 20.0,
    AnchorCircleRadius: float = 0.6,
    RobotCircleRadius: float = 0.1,
    r_pulley: float = 0.005,
    p_screw: float = 0.010,
    eff_screw: float = 0.5,
    gear_ratio: float = 1.0,
    number_of_vectors: int = 1000,
    x_range=None,
    y_range=None,
    z_range=None,
    make_plots: bool = True,
):
    """
    Python port (numerically equivalent) of `run__FSquirrle_motorSelection.m`.

    Notes
    -----
    - Acceleration-related terms are commented out in the MATLAB script; we
      similarly ignore them here.
    """
    if x_range is None:
        x_range = np.arange(-0.3, 0.3 + 1e-12, 0.05)
    if y_range is None:
        y_range = np.arange(-0.3, 0.1 + 1e-12, 0.05)
    if z_range is None:
        z_range = np.arange(0.1, 0.5 + 1e-12, 0.1)

    # Anchor points
    A1 = AnchorCircleRadius * np.array(
        [np.cos(np.deg2rad(90.0)), np.sin(np.deg2rad(90.0)), 0.0]
    )
    A2 = AnchorCircleRadius * np.array(
        [np.cos(np.deg2rad(210.0)), np.sin(np.deg2rad(210.0)), 0.0]
    )
    A3 = AnchorCircleRadius * np.array(
        [np.cos(np.deg2rad(-30.0)), np.sin(np.deg2rad(-30.0)), 0.0]
    )

    # Robot anchor points relative to handle
    C1_u = np.array(
        [
            RobotCircleRadius * np.cos(np.deg2rad(90.0)),
            RobotCircleRadius * np.sin(np.deg2rad(90.0)),
            0.1,
        ]
    )
    C2_u = np.array(
        [
            RobotCircleRadius * np.cos(np.deg2rad(210.0)),
            RobotCircleRadius * np.sin(np.deg2rad(210.0)),
            0.1,
        ]
    )
    C3_u = np.array(
        [
            RobotCircleRadius * np.cos(np.deg2rad(-30.0)),
            RobotCircleRadius * np.sin(np.deg2rad(-30.0)),
            0.1,
        ]
    )

    unit_vectors = fibonacci_sphere(number_of_vectors)

    tau_F_only_line = np.full(
        (len(x_range), len(y_range), len(z_range)), np.nan, dtype=float
    )
    tau_F_only_screw = np.full_like(tau_F_only_line, np.nan)
    omega_pulley_line = np.full_like(tau_F_only_line, np.nan)
    omega_pulley_screw = np.full_like(tau_F_only_line, np.nan)

    for i, x in enumerate(x_range):
        for j, y in enumerate(y_range):
            for k, z in enumerate(z_range):
                R = np.array([x, y, z], dtype=float)

                J_inv = fsquirrle_jacobian_inverse(
                    R, A1, A2, A3, C1_u, C2_u, C3_u, r_pulley, p_screw
                )
                B, _, B_l = fsquirrle_bases(R, A1, A2, A3, C1_u, C2_u, C3_u)

                # Check B_l angle condition as in MATLAB
                dots = np.array(
                    [
                        np.dot(B_l[:, 0], B_l[:, 1]),
                        np.dot(B_l[:, 0], B_l[:, 2]),
                        np.dot(B_l[:, 1], B_l[:, 2]),
                    ]
                )
                if np.min(dots) < np.cos(np.deg2rad(165.0)):
                    continue

                omega_pulley_list = []
                tau_F_list = []

                for u in unit_vectors:
                    hand_force = -u * F

                    # lsqnonneg(B, handForce)
                    res = lsq_linear(B, hand_force, bounds=(0.0, np.inf))
                    line_forces = res.x

                    scale = np.array(
                        [
                            r_pulley / gear_ratio,
                            r_pulley / gear_ratio,
                            r_pulley / gear_ratio,
                            p_screw / (2.0 * np.pi * eff_screw),
                        ]
                    )
                    tau_F_list.append(line_forces * scale)

                    Vel = u * v
                    omega = J_inv @ Vel
                    omega_scaled = omega * np.array(
                        [gear_ratio, gear_ratio, gear_ratio, 1.0]
                    )
                    omega_pulley_list.append(omega_scaled)

                tau_F_arr = np.stack(tau_F_list, axis=1)
                omega_arr = np.stack(omega_pulley_list, axis=1)

                omega_pulley_line[i, j, k] = np.max(omega_arr[0:3, :])
                omega_pulley_screw[i, j, k] = np.max(omega_arr[3, :])
                tau_F_only_line[i, j, k] = np.max(tau_F_arr[0:3, :])
                tau_F_only_screw[i, j, k] = np.max(tau_F_arr[3, :])

    if make_plots:
        from mpl_toolkits.mplot3d import Axes3D

        xx, yy = np.meshgrid(x_range, y_range, indexing="ij")

        fig1 = plt.figure(1)
        fig1.clf()
        for idx in range(omega_pulley_line.shape[2]):
            ax1 = fig1.add_subplot(1, 2, 1, projection="3d")
            ax1.set_title("max torque - line motors")
            ax1.plot_surface(
                xx,
                yy,
                tau_F_only_line[:, :, idx],
                cmap="viridis",
                edgecolor="none",
            )

            ax2 = fig1.add_subplot(1, 2, 2, projection="3d")
            ax2.set_title("max torque - screw motor")
            ax2.plot_surface(
                xx,
                yy,
                tau_F_only_screw[:, :, idx] / gear_ratio,
                cmap="viridis",
                edgecolor="none",
            )

        fig2 = plt.figure(2)
        fig2.clf()
        for idx in range(omega_pulley_line.shape[2]):
            ax1 = fig2.add_subplot(1, 2, 1, projection="3d")
            ax1.set_title("max speed - line motors")
            ax1.plot_surface(
                xx,
                yy,
                omega_pulley_line[:, :, idx] * gear_ratio * 60.0 / (2.0 * np.pi),
                cmap="viridis",
                edgecolor="none",
            )

            ax2 = fig2.add_subplot(1, 2, 2, projection="3d")
            ax2.set_title("max speed - screw motor")
            ax2.plot_surface(
                xx,
                yy,
                omega_pulley_screw[:, :, idx] * 60.0 / (2.0 * np.pi),
                cmap="viridis",
                edgecolor="none",
            )

    return {
        "x_range": np.asarray(x_range),
        "y_range": np.asarray(y_range),
        "z_range": np.asarray(z_range),
        "tau_F_only_line": tau_F_only_line,
        "tau_F_only_screw": tau_F_only_screw,
        "omega_pulley_line": omega_pulley_line,
        "omega_pulley_screw": omega_pulley_screw,
    }


if __name__ == "__main__":
    print("Running FSquirrle motor selection analysis...")
    print("This may take a while due to the workspace sweep...")
    results = run_fsquirrle_motor_selection()
    print("Analysis complete. Close plots to exit.")
    plt.show()
