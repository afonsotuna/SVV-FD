import matplotlib.pyplot as plt
import control
import scipy.io
import numpy as np
from Asymmetric_SS.ss_asymmetric import ss_asym
import math as m
# Define event (from flight test sheet)
t_lookup = 3717                 # time (in seconds) at which event happens
t_limit = 14                    # length of interval (in seconds, multiples of 0.1)
t_interval = t_lookup + t_limit
block_fuel = 4050               # block fuel (lbs)
passenger_weight = 695          # sum of passenger weight (kgs)
c = 2.0569
b = 15.911

# Flight data imported
mat = scipy.io.loadmat('reference_clean.mat')
flight_data = mat['clean_data']

# Get data location
index = int((t_lookup - flight_data[0, 47]) / 0.1)
n_points = int(t_limit / 0.1) + 1

# Obtain correct weight (manoeuvre start) and velocity, get system
used_fuel = flight_data[index, 13] + flight_data[index, 14]
mass_event = (block_fuel - used_fuel + 9165) * 0.453592 + passenger_weight
tas_event = flight_data[index, 41] * 0.514444

# Obtain correspondent flight data
data_event = np.zeros((n_points, 2))
for i in range(n_points):
    data_event[i, 0] = flight_data[index + i, 47] - flight_data[index, 47]
    data_event[i, 1] = flight_data[index + i, 20] - flight_data[
        index, 20]  # Outputs: ?? - beta / 21 - phi / (26*b)/(2*V) - pb/2V / (28*b)/(2*V) - rb/2V
t1 = data_event[:, 0]
y1 = data_event[:, 1]*m.pi/180

# Define initial conditions
beta_0 = 0
phi_0 = flight_data[index, 20] * m.pi/180
pb_2v_0 = (flight_data[index, 26] * b) / (2 * tas_event)
rb_2v_0 = (flight_data[index, 28] * b) / (2 * tas_event)
pb_2v_0 = flight_data[index, 26]
rb_2v_0 = flight_data[index, 28]
initial_cond = np.array([[beta_0], [phi_0], [pb_2v_0], [rb_2v_0]])

# Obtain impulse response
input_delta_a = flight_data[index:index + n_points, 15] * m.pi/180 + 0.006983947075106035
input_delta_r = flight_data[index:index + n_points, 17] * m.pi/180 + 0*0.022041530548663223
input_tot = np.array([input_delta_a, input_delta_r])
sys = ss_asym(rho=1.028,m=mass_event, theta_0=flight_data[index, 21] * m.pi/180, v=tas_event)
t2, out, p2 = control.forced_response(sys, T=t1, U=input_tot)
# t2, out = control.impulse_response(sys, T=t1)

# If it's Dutch roll:
out = -out

# t2, out = control.impulse_response(sys, T=t1)
y2 = out[1, :]  # Outputs: 0 - beta / 1 - phi / 2 - pb/2V / 3 - rb/2V
#hehehe
# IN DEBUGGING - DON'T TOUCH (currently looking at dutch roll 1 for reference data)
# plt.plot(t1, y1, label='Reference data - phi')
# plt.plot(t2, y2, label='System response - phi')
plt.plot(t2, input_delta_a, label='Aileron input')
plt.plot(t2, input_delta_r, label='Rudder input')
plt.legend()
plt.xlabel('Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
plt.ylabel('Radians')
plt.show()
