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
    mat = scipy.io.loadmat('clean_flight_data.mat')
    flight_data = mat['clean_data']

    # Get data location
    index = int((t_lookup - flight_data[0, 48]) / 0.1)
    n_points = int(t_limit / 0.1) + 1

    # Obtain correct weight (manoeuvre start) and velocity, get system
    used_fuel = flight_data[index, 14] + flight_data[index, 15]
    mass_event = (block_fuel - used_fuel + 9165) * 0.453592 + passenger_weight
    tas_event = flight_data[index, 42] * 0.514444

    # obtain correct rho
    h_p = flight_data[index, 37] * 0.3048
    p = 101325 * (1 + (-0.0065 * h_p / 288.15)) ** (-9.81 / (-0.0065 * 287.05))
    T = flight_data[index, 35] + 273.15
    rho = p / (287.05 * T)

    # Obtain correspondent flight data
    data_event = np.zeros((n_points, 2))
    if output == 1:
        for i in range(n_points):
            data_event[i, 0] = flight_data[index + i, 48] - flight_data[index, 48]
            data_event[i, 1] = flight_data[index + i, 0] - flight_data[index, 0]  # Output alpha
    elif output == 2:
        for i in range(n_points):
            data_event[i, 0] = flight_data[index + i, 48] - flight_data[index, 48]
            data_event[i, 1] = (flight_data[index + i, 22]) - (flight_data[index, 22])  # Output theta
    elif output == 3:
        for i in range(n_points):
            data_event[i, 0] = flight_data[index + i, 48] - flight_data[index, 48]
            data_event[i, 1] = ((flight_data[index + i, 27] * c) / flight_data[index, 42])  # Output q
    t1 = data_event[:, 0]
    y1 = data_event[:, 1] * m.pi / 180

    # Define initial conditions
    u_hat_0 = 0
    alpha_0 = 0  # flight_data[index, 0]
    theta_0 = flight_data[index, 22] * m.pi / 180
    qc_v_0 = (flight_data[index, 27] * c) / flight_data[index, 42]
    initial_cond = np.array([[u_hat_0], [alpha_0], [theta_0], [qc_v_0]])

    # Obtain impulse response
    input_delta_e = -flight_data[index:index + n_points, 16] * m.pi / 180
    sys = ss_sym(rho=rho, m=mass_event, theta_0=theta_0, v=tas_event, C_x_q=C_x_q, C_z_q=C_z_q, C_m_alpha=C_m_alpha,
                 C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
    t2, out, p2 = control.forced_response(sys, T=t1, U=input_delta_e)

    y2 = out[output, :]

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

    if output == 1:
        plt.plot(t1, y1, label=r'Reference data - $\alpha$')
        plt.plot(t2, y2, label=r'System response - $\alpha$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Angle of attack [deg]')
        plt.title(
            'Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

    elif output == 2:
        plt.plot(t1, y1, label=r'Reference data - $\theta$')
        plt.plot(t2, y2, label=r'System response - $\theta$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Pitch angle [deg]')
        plt.title(
            'Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

    elif output == 3:
        plt.plot(t1, y1, label=r'Reference data - $q$')
        plt.plot(t2, y2, label=r'System response - $q$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Pitch rate [deg/s]')
        plt.title(
            'Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

    elif output == 4:
        plt.plot(t2, input_delta_e, label='Elevator input')
        plt.legend()
        plt.xlabel('Aileron and rudder input between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.ylabel('Deflection [rad]')
        plt.show()

    return

# t_rn = 3125
# t_lim = 30
#
# make_plot_sym_data(output=1,t_lookup=t_rn,t_limit=t_lim)
# make_plot_sym_data(output=2,t_lookup=t_rn,t_limit=t_lim)
# make_plot_sym_data(output=3,t_lookup=t_rn,t_limit=t_lim)
# make_plot_sym(output=3, t_lookup=3229, t_limit=200, C_x_q=-13.193219977520021, C_z_q=-7.130847541981485, C_m_alpha=-0.21272996346475043, C_m_delta_e=-0.6919484170523861,C_m_q=-11.722501242413275)
