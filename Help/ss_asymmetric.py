import scipy as sp
import numpy as np
import control


def ss_asym(rho=1.225, theta_0=0, m=4157.174, v=80, CY_b=-0.7500, Cn_r=-0.2061, Cn_p=-0.0602, Cl_r=0.2376, Cl_p=-0.7108):
    # dimensions
    W = m * 9.81
    S = 30  # m^2
    c = 2.0569  # m
    b = 15.911  # m

    # aerodynamic coefficients
    u_b = m / (rho * S * b)

    # inertia
    K_xx = (0.019) ** 0.5
    K_zz = (0.042) ** 0.5
    K_xz = 0.002

    # LATERAL FORCE DERIVATIVES
    CY_b = CY_b
    CY_b_d = 0  # bd stands for derivative of beta
    CY_p = -0.0304
    CY_r = 0.8495
    CY_sa = -0.0400  # derivative to aileron deflection
    CY_sr = 0.2300  # derivative to rudder deflection

    # ROLL MOMENT DERIVATIVES
    Cl_b = -0.1026
    Cl_p = Cl_p  # bd stands for derivative of beta
    Cl_r = Cl_r
    Cl_sa = -0.2309  # derivative to aileron deflection
    Cl_sr = 0.0344  # derivative to rudder deflection

    # YAW MOMENT DERIVATIVES
    Cn_b = 0.1348
    Cn_b_d = 0  # bd stands for derivative of beta
    Cn_p = Cn_p
    Cn_r = Cn_r
    Cn_sa = -0.0120  # derivative to aileron deflection
    Cn_sr = -0.0939  # derivative to rudder deflection

    # other parameters only required in asym  matrix
    C_L = (W * np.cos(theta_0)) / (0.5 * rho * v ** 2 * S)  # component of weight along the Zs axis

    # BUILDING MATRICES

    P = np.zeros((4, 4))
    P[0, 0] = (CY_b_d - 2 * u_b) * b / v
    P[0, 1] = 0
    P[0, 2] = 0
    P[0, 3] = 0

    P[1, 0] = 0
    P[1, 1] = (-1 / 2) * (b / v)
    P[1, 2] = 0
    P[1, 3] = 0

    P[2, 0] = 0
    P[2, 1] = 0
    P[2, 2] = -4 * u_b * K_xx ** 2 * (b / v)
    P[2, 3] = 4 * u_b * K_xz * (b / v)
    P[2, 2] = -2 * u_b * K_xx ** 2 * (b / v) ** 2
    P[2, 3] = 2 * u_b * K_xz * (b / v) ** 2

    P[3, 0] = (Cn_b_d) * (b / v)
    P[3, 1] = 0
    P[3, 2] = 4 * u_b * K_xz * (b / v)
    P[3, 3] = -4 * u_b * K_zz ** 2 * (b / v)
    P[3, 2] = 2 * u_b * K_xz * (b / v) ** 2
    P[3, 3] = -2 * u_b * K_zz ** 2 * (b / v) ** 2

    Q = np.zeros((4, 4))
    Q[0, 0] = -CY_b
    Q[0, 1] = -C_L
    Q[0, 2] = -CY_p
    Q[0, 3] = -(CY_r - 4 * u_b)
    Q[0, 2] = -CY_p * b / (2 * v)
    Q[0, 3] = -(CY_r - 4 * u_b) * b / (2 * v)

    Q[1, 0] = 0
    Q[1, 1] = 0
    Q[1, 2] = -1
    Q[1, 3] = 0
    Q[1, 2] = -1 * b / (2 * v)
    Q[1, 3] = 0

    Q[2, 0] = -Cl_b
    Q[2, 1] = 0
    Q[2, 2] = -Cl_p
    Q[2, 3] = -Cl_r
    Q[2, 2] = -Cl_p * b / (2 * v)
    Q[2, 3] = -Cl_r * b / (2 * v)

    Q[3, 0] = -Cn_b
    Q[3, 1] = 0
    Q[3, 2] = -Cn_p
    Q[3, 3] = -Cn_r
    Q[3, 2] = -Cn_p * b / (2 * v)
    Q[3, 3] = -Cn_r * b / (2 * v)

    R = np.zeros((4, 2))
    R[0, 0] = -CY_sa
    R[0, 1] = -CY_sr

    R[1, 0] = 0
    R[1, 1] = 0

    R[2, 0] = -Cl_sa
    R[2, 1] = -Cl_sr

    R[3, 0] = -Cn_sa
    R[3, 1] = -Cn_sr

    A = np.dot(np.linalg.inv(P), Q)
    B = np.dot(np.linalg.inv(P), R)
    C = np.identity(4)
    D = np.zeros((4, 2))

    # print(P, Q, R, A, B, C, D)

    sys = control.StateSpace(A, B, C, D)  # the state vector = [beta, phi, pb/2V, rb/2V]
    # print(u_b, v)
    return sys
