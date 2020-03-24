import matplotlib.pyplot as plt
import control
import scipy.io
import numpy as np
from Symmetric_SS_MATLAB.ss_symmetric import ss_sym
import math as m
# Define event (from flight test sheet)
t_lookup = 3630  # time (in seconds) at which event happens
t_limit = 40  # length of interval (in seconds, multiples of 0.1)
t_interval = t_lookup + t_limit
block_fuel = 4050  # block fuel (lbs)
passenger_weight = 695  # sum of passenger weight (kgs)
c = 2.0569

# Flight data imported
mat = scipy.io.loadmat('reference_clean.mat')
flight_data = mat['clean_data']

# Get data location
index = int((t_lookup - flight_data[0, 47]) / 0.1)
n_points = int(t_limit / 0.1) + 1

# Obtain correct weight (manoeuvre start) and velocity, get system
used_fuel = flight_data[index, 13] + flight_data[index, 14]
mass_event = (block_fuel - used_fuel + 9165) * 0.453592 + passenger_weight
print(mass_event)
tas_event = flight_data[index, 41] * 0.514444

# Obtain correspondent flight data
data_event = np.zeros((n_points, 2))
for i in range(n_points):
    data_event[i, 0] = flight_data[index + i, 47] - flight_data[index, 47]
    data_event[i, 1] = (flight_data[index + i, 21])#*c)/flight_data[index, 41]  #- flight_data[index, 0]# Outputs: ?? - u / 0 - alpha / 21 - theta / 26 - q
t1 = data_event[:, 0]
y1 = data_event[:, 1]* m.pi/180

# Define initial conditions
u_hat_0 = 0
alpha_0 = 0 #flight_data[index, 0]
theta_0 = flight_data[index, 21] * m.pi/180
qc_v_0 = (flight_data[index, 26] * c) / flight_data[index, 41]
initial_cond = np.array([[u_hat_0], [alpha_0], [theta_0], [qc_v_0]])

# Obtain impulse response
input_delta_e = flight_data[index:index + n_points, 16]* m.pi/180
sys = ss_sym(rho=1.028,m=mass_event, theta_0=theta_0, v=tas_event)
t2, out, p2 = control.forced_response(sys, T=t1, U=input_delta_e)
y2 = out[2, :] + initial_cond[2] # Outputs: 0 - u / 1 - alpha / 2 - theta / 3 - q

# IN DEBUGGING - DON'T TOUCH (currently looking at phugoid for reference data)
plt.plot(t2, y2, label='System response - AoA')
plt.plot(t1, y1, label='Reference data - pitch angle')
# plt.plot(t2, input_delta_e, label='Elevator input')
plt.legend()
plt.xlabel('Reference data vs system response between ' + str(t_lookup) + ' [s] and ' + str(t_interval) + ' [s].')
plt.ylabel('Degrees')
plt.show()
