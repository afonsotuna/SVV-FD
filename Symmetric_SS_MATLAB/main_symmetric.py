import matplotlib.pyplot as plt
import control
import scipy.io
import numpy as np
from Symmetric_SS_MATLAB.ss_symmetric import ss_sym
import math as m


def plot_compare():  # OUTPUTS: 0 - u (not working) / 1 - alpha / 2 - theta / 3 - q
    reference = int(input("Flight data (0) or reference data (1)?"))
    motion_str = input("What's the motion? (no capitals)")
    t_0 = int(input("Starting time? (s)"))
    delta_t = int(input("Event duration? (s)"))
    output = int(input("Output? 1 - alpha / 2 - theta / 3 - q"))
    if reference:
        # Flight data imported
        mat = scipy.io.loadmat('reference_clean.mat')
        flight_data = mat['clean_data']

        # Define event (from flight test sheet)
        t_lookup = t_0 - flight_data[0, 47]  # time (in seconds) at which event happens
        t_limit = delta_t  # length of interval (in seconds, multiples of 0.1)
        t_interval = t_lookup + t_limit
        block_fuel = 4050  # block fuel (lbs)
        passenger_weight = 695  # sum of passenger weight (kgs)
        c = 2.0569
        motion = '(' + motion_str + ')'

        # Get data location
        index = int((t_lookup - flight_data[0, 47]) / 0.1)
        n_points = int(t_limit / 0.1) + 1

        # Obtain correct weight (manoeuvre start) and velocity, get system
        used_fuel = flight_data[index, 13] + flight_data[index, 14]
        mass_event = (block_fuel - used_fuel + 9165) * 0.453592 + passenger_weight
        # print(mass_event)
        tas_event = flight_data[index, 41] * 0.514444

        # Obtain correspondent flight data
        data_event = np.zeros((n_points, 2))
        for i in range(n_points):
            data_event[i, 0] = flight_data[index + i, 47] - flight_data[index, 47]
            if output == 1:
                column_index_data = 0
            elif output == 2:
                column_index_data = 21
            elif output == 3:
                column_index_data = 26
            data_event[i, 1] = flight_data[
                index + i, column_index_data]  # Outputs: ?? - u / 0 - alpha / 21 - theta / 26 - q
        t1 = data_event[:, 0]
        y1 = data_event[:, 1] * m.pi / 180

        # Define initial conditions
        u_hat_0 = 0
        alpha_0 = 0  # flight_data[index, 0]
        theta_0 = flight_data[index, 21] * m.pi / 180
        qc_v_0 = (flight_data[index, 26] * c) / flight_data[index, 41]
        initial_cond = np.array([[u_hat_0], [alpha_0], [theta_0], [qc_v_0]])

        # Obtain impulse response
        input_delta_e = flight_data[index:index + n_points, 16] * m.pi / 180
        sys = ss_sym(rho=1.028, m=mass_event, theta_0=theta_0, v=tas_event)
        t2, out, p2 = control.forced_response(sys, T=t1, U=input_delta_e)
        if output == 1:
            y2 = out[1, :] + initial_cond[1]  # Outputs: 0 - u / 1 - alpha / 2 - theta / 3 - q
        elif output == 2:
            y2 = out[2, :] + initial_cond[2]
        elif output == 3:
            y2 = out[2, :] * (flight_data[index, 41] / c) + initial_cond[2]

        # IN DEBUGGING - DON'T TOUCH (currently looking at phugoid for reference data)
        plt.plot(t2, y2, label='System response')
        plt.plot(t1, y1, label='Reference data')
        # plt.plot(t2, input_delta_e, label='Elevator input')
        plt.legend()
        plt.xlabel('Reference data vs system response between ' + str(int(t_lookup)) + ' [s] and ' + str(
            int(t_interval)) + ' [s] ' + motion + '.')
        if output == 1:
            plt.ylabel('Angle of attack \u03B1 [radians]')
        elif output == 2:
            plt.ylabel('Pitch angle \u03B8 [radians]')
        elif output == 3:
            plt.ylabel('Pitch rate q [radians/s]')
        plt.show()
    else:
        file_name = input("What is the .mat file name? (no extension, must be in this directory)")
        plot_compare_flight(t_0, delta_t, motion_str, output, file_name)
    return


def plot_compare_flight(t_0, delta_t, motion_str, output, file_name):
    # Flight data imported
    import_file = file_name + '.mat'
    mat = scipy.io.loadmat(import_file)
    flight_data = mat['clean_data']

    # Define event (from flight test sheet)
    t_lookup = t_0 - flight_data[0, 48]  # time (in seconds) at which event happens
    t_limit = delta_t  # length of interval (in seconds, multiples of 0.1)
    t_interval = t_lookup + t_limit
    block_fuel = 2700  # block fuel (lbs)
    passenger_weight = 771  # sum of passenger weight (kgs)
    c = 2.0569
    motion = '(' + motion_str + ')'

    # Get data location
    index = int((t_lookup - flight_data[0, 48]) / 0.1)
    n_points = int(t_limit / 0.1) + 1

    # Obtain correct weight (manoeuvre start) and velocity, get system
    used_fuel = flight_data[index, 14] + flight_data[index, 15]
    mass_event = (block_fuel - used_fuel + 9165) * 0.453592 + passenger_weight
    # print(mass_event)
    tas_event = flight_data[index, 42] * 0.514444

    # Obtain correspondent flight data
    data_event = np.zeros((n_points, 2))
    for i in range(n_points):
        data_event[i, 0] = flight_data[index + i, 48] - flight_data[index, 48]
        if output == 1:
            column_index_data = 0
        elif output == 2:
            column_index_data = 22
        elif output == 3:
            column_index_data = 27
        data_event[i, 1] = flight_data[
            index + i, column_index_data]  # Outputs: ?? - u / 0 - alpha / 22 - theta / 27 - q
    t1 = data_event[:, 0]
    y1 = data_event[:, 1] * m.pi / 180

    # Define initial conditions
    u_hat_0 = 0
    alpha_0 = 0
    theta_0 = flight_data[index, 22] * m.pi / 180
    qc_v_0 = (flight_data[index, 27] * c) / flight_data[index, 42]
    initial_cond = np.array([[u_hat_0], [alpha_0], [theta_0], [qc_v_0]])

    # Obtain impulse response
    input_delta_e = flight_data[index:index + n_points, 17] * m.pi / 180
    sys = ss_sym(rho=1.028, m=mass_event, theta_0=theta_0, v=tas_event)
    t2, out, p2 = control.forced_response(sys, T=t1, U=input_delta_e)
    if output == 1:
        y2 = out[1, :] + initial_cond[1]
    elif output == 2:
        y2 = out[2, :] + initial_cond[2]
    elif output == 3:
        y2 = out[2, :] * (flight_data[index, 42] / c) + initial_cond[2]

    # Plot making
    plt.plot(t2, y2, label='System response')
    plt.plot(t1, y1, label='Reference data')
    # plt.plot(t2, input_delta_e, label='Elevator input')
    plt.legend()
    plt.xlabel('Reference data vs system response between ' + str(int(t_lookup)) + ' [s] and ' + str(
        int(t_interval)) + ' [s] ' + motion + '.')
    if output == 1:
        plt.ylabel('Angle of attack \u03B1 [radians]')
    elif output == 2:
        plt.ylabel('Pitch angle \u03B8 [radians]')
    elif output == 3:
        plt.ylabel('Pitch rate q [radians/s]')
    plt.show()
    return


run = True
while run:
    plot_compare()
    repeat_status = int(input("Repeat? 1/0"))
    if repeat_status:
        run = True
    else:
        run = False
