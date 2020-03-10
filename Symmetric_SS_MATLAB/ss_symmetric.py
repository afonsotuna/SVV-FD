import math
import numpy as np
import control


def ss_sym(rho=1.225, theta_0=0, m=4157.174, V=80):
    # dimensions
    W = m * 9.81
    S = 30  # m^2
    c = 2.0569  # m
    b = 15.911  # m

    # aerodynamic coefficients
    C_D_0 = 0.04
    C_L_alpha = 5.084
    e = 0.8
    u_c = m / (rho * S * c)

    # inertia
    K_xx = (0.019) ** 0.5
    K_yy = (1.3925) ** 0.5
    K_zz = (0.042) ** 0.5
    K_xz = 0.002

    # long. force deriv.
    C_x_u = -0.0279
    C_x_alpha = -0.4797
    C_x_alphad = 0.0833
    C_x_q = -0.2817
    C_x_delta_e = -0.0373
    C_x_0 = 2 * W * math.sin(theta_0) / (rho * V ** 2 * S)

    # normal force deriv.
    C_z_u = -0.3762
    C_z_alpha = -5.7434
    C_z_alphad = -0.0035
    C_z_q = -5.6629
    C_z_delta_e = -0.6961
    C_z_0 = -W * 2 * math.cos(theta_0) / (rho * V ** 2 * S)

    # pitch moment deriv.
    C_m_0 = 0.0297
    C_m_u = 0.0699
    C_m_alpha = -0.5626
    C_m_alphad = 0.1780
    C_m_q = -8.7941
    C_m_delta_e = -1.1642
    C_m_T_c = -0.0064

    # MATRIX INPUTS

    x_u = V / c * C_x_u / (2 * u_c)
    x_alpha = V / c * C_x_alpha / (2 * u_c)
    x_theta = V / c * C_z_0 / (2 * u_c)

    z_u = V / c * C_z_u / (2 * u_c - C_z_alphad)
    z_alpha = V / c * C_z_alpha / (2 * u_c - C_z_alphad)
    z_theta = V / c * C_x_0 / (2 * u_c - C_z_alphad)
    z_q = V / c * (2 * u_c + C_z_q) / (2 * u_c - C_z_alphad)

    m_u = V / c * (C_m_u + C_z_u * C_m_alphad / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)
    m_alpha = V / c * (C_m_alpha + C_z_alpha * C_m_alphad / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)
    m_theta = V / c * (C_x_0 * C_m_alphad / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)
    m_q = V / c * (C_m_q + C_m_alphad * (2 * u_c + C_z_q) / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)

    x_delta_e = V / c * C_x_delta_e / (2 * u_c)
    z_delta_e = V / c * C_z_delta_e / (2 * u_c - C_z_alphad)
    m_delta_e = V / c * (C_m_delta_e + C_z_delta_e * C_m_alphad / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)

    # MATRICES
    P = c / V * np.array(
        [[-2 * u_c, 0, 0, 0], [0, C_z_alphad - 2 * u_c, 0, 0], [0, 0, -1, 0], [0, C_m_alphad, 0, -2 * u_c * K_yy ** 2]])

    Q = ([[-C_x_u, -C_x_alpha, -C_z_0, 0], [-C_z_u, -C_z_alpha, C_x_0, -(C_z_q + 2 * u_c)], [0, 0, 0, -1],
          [-C_m_u, -C_m_alpha, 0, -C_m_q]])

    R = ([[-C_x_delta_e], [-C_z_delta_e], [0], [-C_m_delta_e]])

    A = np.dot(np.linalg.inv(P), Q)
    B = np.dot(np.linalg.inv(P), R)
    C = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    D = np.array([[0], [0], [0], [0]])

    sys = control.StateSpace(A, B, C, D)

    return sys
