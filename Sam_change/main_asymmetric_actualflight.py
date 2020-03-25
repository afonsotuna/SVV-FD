import matplotlib.pyplot as plt
import control
import scipy.io
import numpy as np
from Sam_change.ss_asymmetric import ss_asym
import math as m


def num_model_asym_reference(output=1, t_lookup=3717, t_limit=14, eigenmotion = "dutch roll", block_fuel=4050, passenger_weight=695, b=15.911,
                             CY_b=-0.7500, Cn_r=-0.2061, Cn_p=-0.0602, Cl_r=0.2376, Cl_p=-0.7108):
    # Outputs: 1 - phi / 2 - pb/2V / 3 - rb/2V

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
    T = flight_data['Dadc1_sat'][0][0][0][index]+ 273.15
    rho = p / (287.05 * T)

    # Obtain correspondent flight data
    data_event = np.zeros((n_points, 2))
    if output == 1:
        for i in range(n_points):
            data_event[i, 0] = flight_data['time'][0][0][0][0][index + i] - flight_data['time'][0][0][0][0][index]
            data_event[i, 1] = flight_data['Ahrs1_Roll'][0][0][0][index + i] # Output phi
    elif output == 2:
        for i in range(n_points):
            data_event[i, 0] = flight_data['time'][0][0][0][0][index + i] - flight_data['time'][0][0][0][0][index]
            data_event[i, 1] = flight_data['Ahrs1_bRollRate'][0][0][0][index + i] # Output pb/2V
    elif output == 3:
        for i in range(n_points):
            data_event[i, 0] = flight_data['time'][0][0][0][0][index + i] - flight_data['time'][0][0][0][0][index]
            data_event[i, 1] = flight_data['Ahrs1_bYawRate'][0][0][0][index + i] # Output rb/2V
    t1 = data_event[:, 0]
    y1 = data_event[:, 1] * m.pi / 180

    # Obtain impulse response
    # input_delta_a = flight_data[index:index + n_points, 15] * m.pi / 180 + 0.006983947075106035
    # input_delta_r = flight_data[index:index + n_points, 17] * m.pi / 180 + 0 * 0.022041530548663223
    # input_tot = np.array([input_delta_a, input_delta_r])

    if eigenmotion == "dutch roll":
        input_delta_a = (flight_data['delta_a'][0][0][0][index:index + n_points] * m.pi / 180 - flight_data['delta_a'][0][0][0][index+n_points] * m.pi/180)
        input_delta_r = -flight_data['delta_r'][0][0][0][index:index + n_points] * m.pi / 180

    if eigenmotion == "aperiodic":
        input_delta_a = -(flight_data['delta_a'][0][0][0][index:index + n_points] * m.pi / 180 - flight_data['delta_a'][0][0][0][index+n_points] * m.pi/180)
        input_delta_r = -flight_data['delta_r'][0][0][0][index:index + n_points] * m.pi / 180

    if eigenmotion == "spiral":
        input_delta_a = -(flight_data['delta_a'][0][0][0][index:index + n_points] * m.pi / 180 - flight_data['delta_a'][0][0][0][index+n_points] * m.pi/180)
        input_delta_r = -(flight_data['delta_r'][0][0][0][index:index + n_points] * m.pi / 180 - flight_data['delta_r'][0][0][0][index+n_points] * m.pi/180)

    input_tot = np.array([input_delta_a[:,0], input_delta_r[:,0]])

    sys = ss_asym(rho=rho, m=mass_event, theta_0=flight_data['Ahrs1_Pitch'][0][0][0][index] * m.pi / 180, v=tas_event, CY_b=CY_b, Cn_r=Cn_r,
                  Cn_p=Cn_p, Cl_r=Cl_r, Cl_p=Cl_p)
    # t2, out, p2 = control.forced_response(sys, T=t1, U=input_tot,X0=[0., -flight_data['Ahrs1_Roll'][0][0][0][index][0]* m.pi / 180,  -flight_data['Ahrs1_bRollRate'][0][0][0][index][0]* m.pi / 180, -flight_data['Ahrs1_bYawRate'][0][0][0][index][0]* m.pi / 180])
    t2, out, p2 = control.forced_response(sys, T=t1, U=input_tot,X0=[0., flight_data['Ahrs1_Roll'][0][0][0][index][0]* m.pi / 180,  flight_data['Ahrs1_bRollRate'][0][0][0][index][0]* m.pi / 180, flight_data['Ahrs1_bYawRate'][0][0][0][index][0]* m.pi / 180])

    if eigenmotion == "dutch roll":
        out = out

    y2 = out[output, :]  # Outputs: 1 - phi / 2 - pb/2V / 3 - rb/2V

    # flight data, model response, time vectors and input vectors
    return y1, y2, t1, t2, input_delta_a, input_delta_r, t_lookup, t_interval

def make_plot_asym(output=1, eigenmotion = "dutch roll", t_lookup=3717, t_limit=14, block_fuel=4050, passenger_weight=695, b=15.911, CY_b=-0.7500,
                   Cn_r=-0.2061, Cn_p=-0.0602, Cl_r=0.2376, Cl_p=-0.7108):
    y1, y2, t1, t2, input_delta_a, input_delta_r, t_lookup, t_interval = num_model_asym_reference(output=output,
                                                                                                  t_lookup=t_lookup,
                                                                                                  t_limit=t_limit,
                                                                                                  eigenmotion=eigenmotion,
                                                                                                  block_fuel=block_fuel,
                                                                                                  passenger_weight=passenger_weight,
                                                                                                  b=b, CY_b=CY_b,
                                                                                                  Cn_r=Cn_r, Cn_p=Cn_p,
                                                                                                  Cl_r=Cl_r, Cl_p=Cl_p)

    if output == 1:
        plt.plot(t1, y1, label=r'Reference data - $\phi$')
        plt.plot(t2, y2, label=r'System response - $\phi$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Roll angle [rad]')
        plt.title('Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

    elif output == 2:
        plt.plot(t1, y1, label=r'Reference data - $p$')
        plt.plot(t2, y2, label=r'System response - $p$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Roll rate [rad/s]')
        plt.title('Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

    elif output == 3:
        plt.plot(t1, y1, label=r'Reference data - $r$')
        plt.plot(t2, y2, label=r'System response - $r$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Yaw rate [rad/s]')
        plt.title('Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.show()

    elif output == 4:
        plt.plot(t2, input_delta_a, label='Aileron input')
        plt.plot(t2, input_delta_r, label='Rudder input')
        plt.legend()
        plt.xlabel('Aileron and rudder input between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
        plt.ylabel('Deflection [rad]')
        plt.show()

    return


#make_plot_asym(output=3, CY_b=-0.7419644366092554,Cn_r=-0.1552558675361675, Cn_p=0.0, Cl_r=0, Cl_p=-0.6898931437580158)
make_plot_asym(output=1, eigenmotion = "aperiodic", t_lookup=3050, t_limit=15)
make_plot_asym(output=2, eigenmotion = "aperiodic", t_lookup=3050, t_limit=15)
make_plot_asym(output=3, eigenmotion = "aperiodic", t_lookup=3050, t_limit=15)
# make_plot_asym(output=1, eigenmotion = "dutch roll", t_lookup=3455, t_limit=15, CY_b = -2.560698700667486, Cn_r = 0.0, Cn_p = 0.0, Cl_r = 0.09170271668625395, Cl_p = -0.7930823749878441)
# make_plot_asym(output=2, eigenmotion = "dutch roll", t_lookup=3455, t_limit=15, CY_b = -2.560698700667486, Cn_r = 0.0, Cn_p = 0.0, Cl_r = 0.09170271668625395, Cl_p = -0.7930823749878441)
# make_plot_asym(output=3, eigenmotion = "dutch roll", t_lookup=3455, t_limit=15, CY_b = -2.560698700667486, Cn_r = 0.0, Cn_p = 0.0, Cl_r = 0.09170271668625395, Cl_p = -0.7930823749878441)
# make_plot_asym(output=1, eigenmotion = "aperiodic", t_lookup=3050, t_limit=20, CY_b = -2.560698700667486, Cn_r = 0.0, Cn_p = 0.0, Cl_r = 0.09170271668625395, Cl_p = -0.7930823749878441)
# make_plot_asym(output=2, eigenmotion = "aperiodic", t_lookup=3050, t_limit=20, CY_b = -2.560698700667486, Cn_r = 0.0, Cn_p = 0.0, Cl_r = 0.09170271668625395, Cl_p = -0.7930823749878441)
# make_plot_asym(output=3, eigenmotion = "aperiodic", t_lookup=3050, t_limit=20, CY_b = -2.560698700667486, Cn_r = 0.0, Cn_p = 0.0, Cl_r = 0.09170271668625395, Cl_p = -0.7930823749878441)
