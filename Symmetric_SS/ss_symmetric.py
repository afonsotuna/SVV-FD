import math
import numpy as np
import control


def ss_sym(rho=1.225, theta_0=0, m=4157.174, v=80, C_x_q = -0.2817, C_z_q = -5.6629, C_m_alpha=-0.7249, C_m_delta_e = -1.4968, C_m_q = -8.7941):
    # dimensions
    W = m * 9.81
    S = 30  # m^2
    c = 2.0569  # m

    # aerodynamic coefficients
    u_c = m / (rho * S * c)

    # inertia
    K_yy = 1.3925 ** 0.5

    # long. force deriv.
    C_x_u = -0.095
    C_x_alpha = -0.4797
    C_x_delta_e = -0.0373
    C_x_q = C_x_q
    C_x_0 = 2 * W * math.sin(theta_0) / (rho * v ** 2 * S)

    # normal force deriv.
    C_z_u = -0.3762
    C_z_alpha = -5.7434
    C_z_alphad = -0.0035
    C_z_q = C_z_q
    C_z_delta_e = -0.6961
    C_z_0 = -W * 2 * math.cos(theta_0) / (rho * v ** 2 * S)

    # pitch moment deriv.
    C_m_u = 0.0699
    C_m_alpha = C_m_alpha  # UPDATE, previous -0.5626
    C_m_alphad = 0.1780
    C_m_q = C_m_q
    C_m_delta_e = C_m_delta_e  # UPDATED, previous was -1.1642

    # MATRICES
    P = c / v * np.array(
        [[-2 * u_c, 0, 0, 0],
         [0, C_z_alphad - 2 * u_c, 0, 0],
         [0, 0, -1, 0],
         [0, C_m_alphad, 0, -2 * u_c * K_yy ** 2]])

    Q = ([[-C_x_u, -C_x_alpha, -C_z_0, C_x_q],
          [-C_z_u, -C_z_alpha, C_x_0, -(C_z_q + 2 * u_c)],
          [0, 0, 0, -1],
          [-C_m_u, -C_m_alpha, 0, -C_m_q]])

    R = ([[-C_x_delta_e], [-C_z_delta_e], [0], [-C_m_delta_e]])

    A = np.dot(np.linalg.inv(P), Q)
    B = np.dot(np.linalg.inv(P), R)
    C = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    D = np.array([[0], [0], [0], [0]])

    sys = control.StateSpace(A, B, C, D)

    return sys
