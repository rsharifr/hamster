import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are


def lqr(A, B, Q, R):
    """
    Simple continuous-time LQR, equivalent to MATLAB lqr(A,B,Q,R).
    """
    # Solve continuous-time algebraic Riccati equation
    P = solve_continuous_are(A, B, Q, R)
    # Compute LQR gain
    K = np.linalg.inv(R) @ (B.T @ P)
    return K


def run_1d_lqr_test(
    M: float = 5.0,
    M_ad: float = 0.5,
    Kd_ad: float = 0.0,
    F_des: float = 0.0,
    dt: float = 0.001,
    Kp_hand: float = 50.0,
    Kd_hand: float = 50.0,
    hand_frequency: float = 0.5 * np.pi,
    hand_amp: float = 0.3,
    Q_scale: float = 1e3,
    R_val: float = 1.0,
    gain_p: float = 10.0,
    gain_d: float = 10.0,
    t_end: float = 10.0,
    robot_noise: float = 0.001,
    force_noise: float = 0.05,
    delay_sec: float = 0.01,
    make_plots: bool = True,
):
    """
    Python port of `run__1D_lqr_test.m`.

    Returns
    -------
    results : dict with time vectors and logs.
    """
    A = np.array([[0.0, 1.0], [0.0, 0.0]])
    B = np.array([[0.0], [1.0 / M]])

    A_ad = np.array([[0.0, 1.0], [0.0, -Kd_ad]])
    B_ad = np.array([[0.0], [1.0 / M_ad]])

    Q = Q_scale * np.diag([1.0, 1.0])
    R = np.array([[R_val]])

    K_lqr = lqr(A, B, Q, R)

    N = int(round(t_end / dt))
    delay = int(round(delay_sec / dt))

    IC = np.array([0.0, hand_frequency * hand_amp])
    X = np.zeros((2, N + 1))
    X[:, 0] = IC
    X_ad = np.copy(X)

    F_measured_log = np.zeros(N)
    hand_log = np.zeros((N, 2))
    robot_log = np.zeros((N, 2))
    model_log = np.zeros((N, 2))

    for i in range(N):
        t = i * dt

        hand_p = hand_amp * np.sin(hand_frequency * t)
        hand_v = hand_amp * hand_frequency * np.cos(hand_frequency * t)

        robot_p = X[0, i] + robot_noise * np.random.randn()
        robot_v = X[1, i] + robot_noise * np.random.randn()

        model_p = X_ad[0, i]
        model_v = X_ad[1, i]

        F_int = (hand_p - robot_p) * Kp_hand + (hand_v - robot_v) * Kd_hand
        F_measured = F_int + force_noise * np.random.randn()

        F_measured_log[i] = F_measured
        hand_log[i] = np.array([hand_p, hand_v])
        robot_log[i] = np.array([robot_p, robot_v])
        model_log[i] = np.array([model_p, model_v])

        if i <= delay:
            F_feedback = 0.0
            robot_feedback = np.array([0.0, 0.0])
        else:
            F_feedback = F_measured_log[i - delay]
            robot_feedback = robot_log[i - delay]

        Xdot_ad = (A_ad @ X_ad[:, i]) + (B_ad.flatten() * (F_des + F_feedback))
        X_ad[:, i + 1] = X_ad[:, i] + dt * Xdot_ad

        u_pid = gain_p * (model_p - robot_feedback[0]) + gain_d * (
            model_v - robot_feedback[1]
        )

        # Note: MATLAB code computes u_lqr and u_k but does not use them in dynamics.
        u_lqr = -K_lqr @ (robot_feedback - X_ad[:, i])
        u_k = -12.0 * (F_des - F_feedback)
        _ = (u_lqr, u_k)  # keep for parity, though unused

        Xdot = (A @ X[:, i]) + B.flatten() * (u_pid + F_measured)
        X[:, i + 1] = X[:, i] + dt * Xdot

    t_states = np.linspace(0.0, t_end + dt, N)
    t_forces = np.linspace(0.0, t_end, N)

    if make_plots:
        plt.figure(1)
        plt.clf()

        plt.subplot(1, 2, 1)
        plt.plot(
            t_states,
            robot_log[:, 0],
            t_states,
            model_log[:, 0],
            t_states,
            hand_log[:, 0],
            linewidth=2,
        )
        plt.title("states")
        plt.legend(["robot x", "model x", "hand x"])

        plt.subplot(1, 2, 2)
        plt.plot(t_forces, F_measured_log)
        plt.ylabel("f measured")

        plt.tight_layout()

    return {
        "t_states": t_states,
        "t_forces": t_forces,
        "F_measured_log": F_measured_log,
        "hand_log": hand_log,
        "robot_log": robot_log,
        "model_log": model_log,
    }


if __name__ == "__main__":
    print("Running 1D LQR admittance control test...")
    results = run_1d_lqr_test()
    print("Simulation complete. Close plots to exit.")
    plt.show()

