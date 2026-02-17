"""
HAMSTER omniwheel robot motor selection analysis.

Uses the generated Jacobian from jacobians (see build_hamster_jacobian.py).
"""

import numpy as np
import matplotlib.pyplot as plt

from jacobians import hamster_jacobian


def run_hamster_motor_selection(
    v: float = 1.0,
    a: float = 2.0,
    F: float = 20.0,
    T: float = 0.0,
    Omega: float = 0.0,
    phaseshift: float = np.pi,
    D: float = 0.15,
    r: float = 0.04,
    theta: float = 0.0,
    m: float = 5.0,
    I: float = 0.4,
    gear_ratio: float = 1.0,
    phi_step: float = 0.1,
    make_plots: bool = True,
):
    """
    Python port of `run__HAMSTER_motorSelection.m`.

    Returns
    -------
    results : dict
        Dictionary containing:
        - phi_range
        - tau_F_only
        - tau_a_only
        - tau_combined
        - omega_wheel
    """
    J = hamster_jacobian(D, r, theta)
    M = np.diag([m, m, I])

    phi_range = np.arange(0.0, 2.0 * np.pi + 1e-12, phi_step)

    tau_F_only = []
    tau_a_only = []
    omega_wheel = []

    for phi in phi_range:
        Fx = F * np.cos(phi)
        Fy = F * np.sin(phi)

        hand_force = -np.array([Fx, Fy, T])
        tau_F_only.append(J.T @ hand_force)

        ax = a * np.cos(phi + phaseshift)
        ay = a * np.sin(phi + phaseshift)
        robot_acc = np.array([ax, ay, 0.0])
        tau_a_only.append(J.T @ (M @ robot_acc))

        vx = v * np.cos(phi)
        vy = v * np.sin(phi)
        V = np.array([vx, vy, Omega])

        omega_wheel.append(np.linalg.solve(J, V))

    tau_F_only = np.stack(tau_F_only, axis=1)
    tau_a_only = np.stack(tau_a_only, axis=1)
    omega_wheel = np.stack(omega_wheel, axis=1)

    tau_combined = tau_F_only + tau_a_only

    if make_plots:
        deg = np.rad2deg(phi_range)
        plt.figure(1)
        plt.clf()

        plt.subplot(3, 1, 1)
        plt.plot(deg, tau_F_only.T / gear_ratio)
        plt.title("motor torque for hand force only")

        plt.subplot(3, 1, 2)
        plt.plot(deg, tau_a_only.T / gear_ratio)
        plt.title("motor torque for linear acceleration only")

        plt.subplot(3, 1, 3)
        plt.plot(deg, tau_combined.T / gear_ratio)
        plt.title("combined effect of acceleration and force")

        plt.figure(2)
        plt.clf()
        plt.plot(deg, (omega_wheel * gear_ratio).T)
        plt.title("motor speed")

        plt.tight_layout()

    return {
        "phi_range": phi_range,
        "tau_F_only": tau_F_only,
        "tau_a_only": tau_a_only,
        "tau_combined": tau_combined,
        "omega_wheel": omega_wheel,
    }


if __name__ == "__main__":
    print("Running HAMSTER motor selection analysis...")
    results = run_hamster_motor_selection()
    print("Analysis complete. Close plots to exit.")
    plt.show()
