import scipy as sp
import numpy as np
import control

# Paramaters

rho = 1.225
theta_0 = 0
m = 4157.174  # kg
W = m * 9.81
V = 300

# dimensions
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

# LATERAL FORCE DERIVATIVES
CY_b = -0.7500
CY_b_d = 0  # bd stands for derivative of beta
CY_p = -0.0304
CY_r = 0.8495
CY_sa = -0.0400  # derivative to aileron deflection
CY_sr = 0.2300  # derivative to rudder deflection

# ROLL MOMENT DERIVATIVES
Cl_b = -0.1026
Cl_p = -0.7108  # bd stands for derivative of beta
Cl_r = 0.2376
Cl_sa = -0.2309  # derivative to aileron deflection
Cl_sr = 0.0344  # derivative to rudder deflection

# YAW MOMENT DERIVATIVES
Cn_b = 0.1348
Cn_b_d = 0  # bd stands for derivative of beta
Cn_p = -0.0602
Cn_r = -0.2061
Cn_sa = -0.0120  # derivative to aileron deflection
Cn_sr = -0.0939  # derivative to rudder deflection

# other parameters only required in this matrix
C_L = (W * np.cos(theta_0)) / (0.5 * rho * V ** 2 * S)  # component of weight along the Zs axis

# BUILDING MATRIXES

P = np.zeros((4, 4))
P[0, 0] = (CY_b_d - 2 * u_c) * b / V
P[0, 1] = 0
P[0, 2] = 0
P[0, 3] = 0

P[1, 0] = 0
P[1, 1] = (-1 / 2) * (b / V)
P[1, 2] = 0
P[1, 3] = 0

P[2, 0] = 0
P[2, 1] = 0
P[2, 2] = -4 * u_c * K_xx ** 2 * (b / V)
P[2, 3] = 4 * u_c * K_xz * (b / V)

P[3, 0] = (Cn_b_d) * (b / V)
P[3, 1] = 0
P[3, 2] = 4 * u_c * K_xz * (b / V)
P[3, 3] = -4 * u_c * K_zz ** 2 * (b / V)

Q = np.zeros((4, 4))
Q[0, 0] = -CY_b
Q[0, 1] = -C_L
Q[0, 2] = -CY_p
Q[0, 3] = -(CY_r - 4 * u_c)

Q[1, 0] = 0
Q[1, 1] = 0
Q[1, 2] = -1
Q[1, 3] = 0

Q[2, 0] = -Cl_b
Q[2, 1] = 0
Q[2, 2] = -Cl_p
Q[2, 3] = -Cl_r

Q[3, 0] = -Cn_b
Q[3, 1] = 0
Q[3, 2] = -Cn_p
Q[3, 3] = -Cn_r

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

# print(C)
# print(B)

# the state vector = [beta, phi, pb/2V, rb/2V]
# print(np.linalg.eigvals(A))
#print(P)
#print(R)
sys = control.StateSpace(A, B, C, D)
