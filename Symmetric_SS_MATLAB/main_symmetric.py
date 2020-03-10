import matplotlib.pyplot as plt
import control
import scipy.io
import numpy as np
from Symmetric_SS_MATLAB.ss_symmetric import ss_sym

# Flight data imported and SS system computed
mat = scipy.io.loadmat('flight1_clean.mat')
flight_data = mat['clean_data']
sys = ss_sym(m=6324, V=80)  # for blank arguments => rho=1.225, theta_0=0, m=4157.174, V=80

# Define event (from flight test sheet)
t_lookup = 3020  # time (in seconds) at which event happens
t_limit = 140  # length of interval (in seconds, multiples of 0.1)

# Obtain correspondent flight data
index = int((t_lookup - flight_data[0, 48]) / 0.1)
n_points = int(t_limit / 0.1) + 1
data_event = np.zeros((n_points, 2))
for i in range(n_points):
    data_event[i, 0] = flight_data[index + i, 48] - flight_data[index, 48]
    data_event[i, 1] = flight_data[index + i, 0] - flight_data[index, 0]  # Outputs: ?? - u / 0 - alpha / ?? - theta / ?? - q
t1 = data_event[:, 0]
y1 = data_event[:, 1]

# Obtain impulse response
t2, y2 = control.impulse_response(sys, T=t1, output=2)  # Outputs: 1 - u / 2 - alpha / 3 - theta / 4 - q


# IN DEBUGGING - DON'T TOUCH (currently looking at phugoid for flight 1)
plt.plot(t1, y1)
plt.plot(t2, y2)
plt.show()
