import matplotlib.pyplot as plt
import control
import scipy.io
import numpy as np
from Symmetric_SS.ss_symmetric import ss_sym
import math as m


def num_model_sym_data(output=1, t_lookup=3717, t_limit=14, block_fuel=2700, passenger_weight=771, c=2.0569,
                       C_x_q=-0.2817, C_z_q=-5.6629, C_m_alpha=-0.7249, C_m_delta_e=-1.4968, C_m_q=-8.7941):
    # Outputs: 1 - alpha / 2 - theta / 3 - q

    t_interval = t_lookup + t_limit

    # Flight data imported
    mat = scipy.io.loadmat('flight_actual.mat')
    flight_data = mat['flightdata']
    flight_data = flight_data[0, 0]

    # Get data location
    index = int((t_lookup - flight_data['time'][0][0][0][0][0]) / 0.1)
    n_points = int(t_limit / 0.1) + 1

    # Obtain correct weight (manoeuvre start) and velocity, get system
    used_fuel = flight_data['lh_engine_FU'][0][0][0][index] + flight_data['rh_engine_FU'][0][0][0][index]
    mass_event = (block_fuel - used_fuel + 9165) * 0.453592 + passenger_weight
    tas_event = flight_data['Dadc1_tas'][0][0][0][index] * 0.514444

    # obtain correct rho
    h_p = flight_data['Dadc1_alt'][0][0][0][index] * 0.3048
    p = 101325 * (1 + (-0.0065 * h_p / 288.15)) ** (-9.81 / (-0.0065 * 287.05))
    T = flight_data['Dadc1_sat'][0][0][0][index] + 273.15
    rho = p / (287.05 * T)

    # Obtain correspondent flight data
    data_event = np.zeros((n_points, 2))
    if output == 0:
        for i in range(n_points):
            data_event[i, 0] = flight_data['time'][0][0][0][0][index + i] - flight_data['time'][0][0][0][0][index]
            data_event[i, 1] = flight_data['Dadc1_tas'][0][0][0][index + i] * 0.514444
    if output == 1:
        for i in range(n_points):
            data_event[i, 0] = flight_data['time'][0][0][0][0][index + i] - flight_data['time'][0][0][0][0][index]
            data_event[i, 1] = (flight_data['vane_AOA'][0][0][0][index + i] - flight_data['vane_AOA'][0][0][0][
                index] * 0) * m.pi / 180  # Output alpha
    elif output == 2:
        for i in range(n_points):
            data_event[i, 0] = flight_data['time'][0][0][0][0][index + i] - flight_data['time'][0][0][0][0][index]
            data_event[i, 1] = (flight_data['Ahrs1_Pitch'][0][0][0][index + i] - 0 *
                                flight_data['Ahrs1_Pitch'][0][0][0][
                                    index]) * m.pi / 180  # Output theta
    elif output == 3:
        for i in range(n_points):
            data_event[i, 0] = flight_data['time'][0][0][0][0][index + i] - flight_data['time'][0][0][0][0][index]
            data_event[i, 1] = ((flight_data['Ahrs1_bPitchRate'][0][0][0][index + i] * c) / \
                                flight_data['Dadc1_tas'][0][0][0][index] * 0.514444) * m.pi / 180  # Output q
            data_event[i, 1] = flight_data['Ahrs1_bPitchRate'][0][0][0][index + i] * m.pi / 180  # Output q
    t1 = data_event[:, 0]
    y1 = data_event[:, 1]

    # Define initial conditions
    u_hat_0 = 0
    alpha_0 = flight_data['vane_AOA'][0][0][0][index] * m.pi / 180  # flight_data[index, 0]
    theta_0 = flight_data['Ahrs1_Pitch'][0][0][0][index] * m.pi / 180
    qc_v_0 = (flight_data['Ahrs1_bPitchRate'][0][0][0][index] * c) / flight_data['Dadc1_tas'][0][0][0][index] / 0.514444
    initial_cond = np.array([[u_hat_0], [0], [0], [qc_v_0[0]]])

    # Obtain impulse response
    input_delta_e = (flight_data['delta_e'][0][0][0][index:index + n_points] - flight_data['delta_e'][0][0][0][
        index]) * m.pi / 180
    input_tot = np.array([input_delta_e[:, 0]], dtype='float')

    sys = ss_sym(rho=float(rho), m=float(mass_event), theta_0=float(theta_0), v=float(tas_event), C_x_q=C_x_q,
                 C_z_q=C_z_q, C_m_alpha=C_m_alpha, C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
    t2, out, p2 = control.forced_response(sys, T=t1, U=input_tot)

    out[0, :] *= flight_data['Dadc1_tas'][0][0][0][index] * 0.514444
    out[0, :] += flight_data['Dadc1_tas'][0][0][0][index] * 0.514444
    out[1, :] += alpha_0
    out[2, :] += theta_0

    y2 = out[output, :]

    if output == 3:
        # y1 *= flight_data['Dadc1_tas'][0][0][0][index]*0.514444 / c
        y2 *= flight_data['Dadc1_tas'][0][0][0][index] * 0.514444 / c

    return y1, y2, t1, t2, input_delta_e, t_lookup, t_interval


def make_plot_sym_data(output=1, t_lookup=3717, t_limit=14, block_fuel=2700, passenger_weight=771, c=2.0569,
                       C_x_q=-0.2817, C_z_q=-5.6629, C_m_alpha=-0.7249, C_m_delta_e=-1.4968, C_m_q=-8.7941):
    y1, y2, t1, t2, input_delta_e, t_lookup, t_interval = num_model_sym_data(output=output, t_lookup=t_lookup,
                                                                             t_limit=t_limit,
                                                                             block_fuel=block_fuel,
                                                                             passenger_weight=passenger_weight,
                                                                             c=c, C_x_q=C_x_q, C_z_q=C_z_q,
                                                                             C_m_alpha=C_m_alpha,
                                                                             C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)

    if output == 0:
        plt.plot(t1, y1, label=r'Reference data - $u$')
        plt.plot(t2, y2, label=r'System response - $u$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Velocity [m/s]')
        plt.title(
            'Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

    if output == 1:
        plt.plot(t1, y1, label=r'Reference data - $\alpha$')
        plt.plot(t2, y2, label=r'System response - $\alpha$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Angle of attack [rad]')
        plt.title(
            'Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

    elif output == 2:
        plt.plot(t1, y1, label=r'Reference data - $\theta$')
        plt.plot(t2, y2, label=r'System response - $\theta$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Pitch angle [rad]')
        plt.title(
            'Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

    elif output == 3:
        plt.plot(t1, y1, label=r'Reference data - $q$')
        plt.plot(t2, y2, label=r'System response - $q$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Pitch rate [rad/s]')
        plt.title(
            'Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

        # elif output == 4:
        plt.plot(t2, input_delta_e, label='Elevator input')
        plt.legend()
        plt.xlabel('Elevator input between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.ylabel('Deflection [rad]')
        plt.show()

    return


t_rn = 3124.6
t_lim = 3.5
C_x_q, C_z_q, C_m_alpha, C_m_delta_e, C_m_q = -0.2817, -5.6629, -0.7249, -1.4968, -8.7941
# C_x_q, C_z_q, C_m_alpha, C_m_delta_e,C_m_q = -13.193219977520021, -7.130847541981485, -0.21272996346475043, -0.6919484170523861,-11.722501242413275

make_plot_sym_data(output=0, t_lookup=t_rn, t_limit=t_lim, C_x_q=C_x_q, C_z_q=C_z_q, C_m_alpha=C_m_alpha,
                   C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
make_plot_sym_data(output=1, t_lookup=t_rn, t_limit=t_lim, C_x_q=C_x_q, C_z_q=C_z_q, C_m_alpha=C_m_alpha,
                   C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
make_plot_sym_data(output=2, t_lookup=t_rn, t_limit=t_lim, C_x_q=C_x_q, C_z_q=C_z_q, C_m_alpha=C_m_alpha,
                   C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
# make_plot_sym_data(output=3, t_lookup=t_rn, t_limit=t_lim, C_x_q=C_x_q, C_z_q=C_z_q, C_m_alpha=C_m_alpha, C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)

t_rn = 3219
t_lim = 170
make_plot_sym_data(output=0, t_lookup=t_rn, t_limit=t_lim, C_x_q=C_x_q, C_z_q=C_z_q, C_m_alpha=C_m_alpha,
                   C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
make_plot_sym_data(output=1, t_lookup=t_rn, t_limit=t_lim, C_x_q=C_x_q, C_z_q=C_z_q, C_m_alpha=C_m_alpha,
                   C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
make_plot_sym_data(output=2, t_lookup=t_rn, t_limit=t_lim, C_x_q=C_x_q, C_z_q=C_z_q, C_m_alpha=C_m_alpha,
                   C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
# make_plot_sym_data(output=3, t_lookup=t_rn, t_limit=t_lim, C_x_q=C_x_q, C_z_q=C_z_q, C_m_alpha=C_m_alpha, C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
