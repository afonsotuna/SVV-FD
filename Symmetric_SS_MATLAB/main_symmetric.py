import matplotlib.pyplot as plt
import control
import scipy.io
from Symmetric_SS_MATLAB.ss_symmetric import ss_sym

# Flight data imported and SS system computed
mat = scipy.io.loadmat('flight1_clean.mat')
flight_data = mat['clean_data']
sys = ss_sym()  # for blank arguments => rho=1.225, theta_0=0, m=4157.174, V=300

# Obtain impulse response
t, y = control.impulse_response(sys, output=2)  # Outputs: 1 - u / 2 - alpha / 3 - theta / 4 - q


plt.plot(t, y)
plt.show()
