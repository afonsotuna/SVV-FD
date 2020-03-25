C_D_0 = 0.04
C_L_alpha = 5.084
e = 0.8

b = 15.911  # m

K_xx = (0.019) ** 0.5
K_zz = (0.042) ** 0.5
K_xz = 0.002


C_x_alphad = 0.0833
C_x_q = -0.2817

C_m_0 = 0.0297
C_m_T_c = -0.0064

# MATRIX INPUTS

x_u = v / c * C_x_u / (2 * u_c)
x_alpha = v / c * C_x_alpha / (2 * u_c)
x_theta = v / c * C_z_0 / (2 * u_c)

z_u = v / c * C_z_u / (2 * u_c - C_z_alphad)
z_alpha = v / c * C_z_alpha / (2 * u_c - C_z_alphad)
z_theta = v / c * C_x_0 / (2 * u_c - C_z_alphad)
z_q = v / c * (2 * u_c + C_z_q) / (2 * u_c - C_z_alphad)

m_u = v / c * (C_m_u + C_z_u * C_m_alphad / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)
m_alpha = v / c * (C_m_alpha + C_z_alpha * C_m_alphad / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)
m_theta = v / c * (C_x_0 * C_m_alphad / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)
m_q = v / c * (C_m_q + C_m_alphad * (2 * u_c + C_z_q) / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)

x_delta_e = v / c * C_x_delta_e / (2 * u_c)
z_delta_e = v / c * C_z_delta_e / (2 * u_c - C_z_alphad)
m_delta_e = v / c * (C_m_delta_e + C_z_delta_e * C_m_alphad / (2 * u_c - C_z_alphad)) / (2 * u_c * K_yy ** 2)