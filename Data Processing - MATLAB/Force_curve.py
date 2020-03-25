import pandas as Pd
import numpy as np
import matplotlib.pyplot as plt
import xlrd
from scipy import stats, optimize
from Reduction_function import reduce

rho0 = 1.2250
g = 9.81
S = 30.00
C_m_T_c = -0.0064 # From appendix C tables
BEM = 9165    # Basic empty mass, taken from mass report for 2020 [lbs]
c = 2.0569    # MAC, taken from cit_par.py [m]
d = 0.686     # diameter of engine [m], from online source
W_s = 60500   # standard weight [N]

ref = 0
plot = 1
# reading data unto pandas dataframe:
if ref == 1:
    exc_data = 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx'
else:
    exc_data = 'C:/Users/nano2598\Documents/actual college ugh/Third year/SVV-FD/Data_repo/20200311_V4.xlsx'

data = Pd.read_excel(exc_data, header=55, usecols='B, D:M', skiprows=[57], nrows=8)

Fe = data['Fe'].iloc[1:].to_numpy(dtype=float)  # Putting pressure alt (hp) values from the dataframe column to np array

hp = data['hp'].iloc[1:].to_numpy(dtype=float)  # Putting pressure alt (hp) values from the dataframe column to np array
V_ias = data['IAS'].iloc[1:].to_numpy(dtype=float)  # Putting V_IAS values from the dataframe column to np array
Tm = data['TAT'].iloc[1:].to_numpy(dtype=float)  # Putting TAT values from the dataframe column to np array

V_e, T, M, Re_range = reduce(hp, V_ias, Tm)
F_used = data['F. used'].iloc[1:].to_numpy(dtype=float)  # Same as what's done above but for fuel used

pax_masses = np.array(Pd.read_excel(exc_data, header=None, usecols='H', skiprows=7, nrows=9))
init_fuel = xlrd.open_workbook(exc_data).sheet_by_index(0).cell_value(17, 3)
ramp_mass = (BEM + init_fuel) * 0.453592 + np.sum(
    pax_masses)  # Calculating total ramp mass, converting some lbs terms to kg
mass = ramp_mass - F_used * 0.453592  # Mass at each data point, also converting to kg

V_e_red = V_e * np.sqrt(W_s / (mass * g))
Fe_red = Fe * (W_s/(mass * 9.81))

coeffs = np.polyfit(V_e_red, Fe_red, 2)

if plot == 1:
    x1 = np.linspace(0, 120, 120)
    y1 = coeffs[0] * x1 ** 2 + coeffs[1] * x1 + coeffs[2]
    xlab = ('Reduced velocity (\u1E7C$_e$) [m/s]')
    ylab = ('Reduced stick force (F$_e^*$) [N]')
    plt.axes(xlim=(60, 105), ylim=(100, -60), xlabel=xlab, ylabel=ylab)
    plt.scatter(V_e_red, Fe_red, s=70, c='r', marker='x', label='Reduced stick force data')
    plt.plot(x1, y1, 'darkturquoise', label='Quadratic least-squares fit')
    plt.grid()
    plt.legend(loc='best')
    plt.savefig('figs/Reduced_Force_Curve')
    plt.show()
