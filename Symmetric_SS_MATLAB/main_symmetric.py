import matplotlib.pyplot as plt
import control
import scipy.io
import numpy as np
from Symmetric_SS_MATLAB.ss_symmetric import ss_sym

# Define event (from flight test sheet)
t_lookup = 3237  # time (in seconds) at which event happens
t_limit = 120  # length of interval (in seconds, multiples of 0.1)
block_fuel = 4050  # block fuel (lbs)
passenger_weight = 695  # sum of passenger weight (kgs)

# Flight data imported and SS system computed
mat = scipy.io.loadmat('reference_clean.mat')
flight_data = mat['clean_data']

# Get data location
index = int((t_lookup - flight_data[0, 47]) / 0.1)
n_points = int(t_limit / 0.1) + 1

# Obtain correct weight (manoeuvre start) and velocity, get system
used_fuel = flight_data[index, 13] + flight_data[index, 14]
mass_event = (block_fuel - used_fuel + 9165) * 0.453592 + passenger_weight
tas_event = flight_data[index, 41] * 0.514444
sys = ss_sym(m=mass_event, v=tas_event)

# Obtain correspondent flight data
data_event = np.zeros((n_points, 2))
for i in range(n_points):
    data_event[i, 0] = flight_data[index + i, 47] - flight_data[index, 47]
    data_event[i, 1] = flight_data[index + i, 0] - flight_data[
        index, 0]  # Outputs: ?? - u / 0 - alpha / ?? - theta / ?? - q
t1 = data_event[:, 0]
y1 = data_event[:, 1]

# Obtain impulse response
input_delta_e = flight_data[index:index + n_points, 16]
# t2, y2 = control.impulse_response(sys, T=t1, output=2)
t2, out, p2 = control.forced_response(sys, T=t1, U=input_delta_e)
y2 = out[1, :]  # Outputs: 0 - u / 1 - alpha / 2 - theta / 3 - q

# IN DEBUGGING - DON'T TOUCH (currently looking at phugoid for reference data)
plt.plot(t1, y1)
plt.plot(t2, y2)
plt.plot(t2, input_delta_e)
plt.show()