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
Cm0 = 0.0297

ref = 1
# reading data unto pandas dataframe:
if ref == 1:
    exc_data = 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx'
else:
    exc_data = 'C:/Users/nano2598\Documents/actual college ugh/Third year/SVV-FD/Data_repo/20200311_V4.xlsx'

data = Pd.read_excel(exc_data, header=72, usecols='B, D:M', nrows=3)

if ref == 1:
    cg1 = 282.1501
    cg2 = 280.0019
else:
    cg1 = 280.20982
    cg2 = 277.60802

hp = data['hp'].iloc[1:].to_numpy(dtype=float)  # Putting pressure alt (hp) values from the dataframe column to np array
V_ias = data['IAS'].iloc[1:].to_numpy(dtype=float)    # Putting V_IAS values from the dataframe column to np array
Tm = data['TAT'].iloc[1:].to_numpy(dtype=float)     # Putting TAT values from the dataframe column to np array
delta_e = data['de'].iloc[1:].to_numpy(dtype=float)

V_e, T, M, Re = reduce(hp, V_ias, Tm)
F_used = data['F. used'].iloc[1:].to_numpy(dtype=float)  # Same as what's done above but for fuel used

pax_masses = np.array(Pd.read_excel(exc_data, header=None, usecols='H', skiprows=7, nrows=9))
init_fuel = xlrd.open_workbook(exc_data).sheet_by_index(0).cell_value(17, 3)
ramp_mass = (BEM + init_fuel)*0.453592 + np.sum(pax_masses)  # Calculating total ramp mass, converting some lbs terms to kg
mass = ramp_mass - F_used * 0.453592                         # Mass at each data point, also converting to kg

C_L = (mass * g) / (0.5 * S * rho0 * V_e ** 2 )
avgC_L = np.average(C_L) # The two are actually almost the same

C_m_deltae = -(avgC_L * (cg2 - cg1) * 0.0254) / ((delta_e[1] - delta_e[0]) * c)
C_m_deltae_rad = C_m_deltae * (180/np.pi)

#print(C_m_deltae)
print('C_m_delta_e in radians is', C_m_deltae_rad)

### Now the second block of data (line 56 in excel) has to be read to find the gradient of the trim curve to find C_m_alpha

data2 = Pd.read_excel(exc_data, header=55, usecols='B, D:M', skiprows=[57], nrows=8)


hp2 = data2['hp'].iloc[1:].to_numpy(dtype=float)  # Putting pressure alt (hp) values from the dataframe column to np array
V_ias2 = data2['IAS'].iloc[1:].to_numpy(dtype=float)    # Putting V_IAS values from the dataframe column to np array
Tm2 = data2['TAT'].iloc[1:].to_numpy(dtype=float)     # Putting TAT values from the dataframe column to np array

V_e2, T2, M2, Re2 = reduce(hp2, V_ias2, Tm2)
F_used = data2['F. used'].iloc[1:].to_numpy(dtype=float)  # Same as what's done above but for fuel used

pax_masses = np.array(Pd.read_excel(exc_data, header=None, usecols='H', skiprows=7, nrows=9))
init_fuel = xlrd.open_workbook(exc_data).sheet_by_index(0).cell_value(17, 3)
ramp_mass = (BEM + init_fuel)*0.453592 + np.sum(pax_masses)  # Calculating total ramp mass, converting some lbs terms to kg
mass = ramp_mass - F_used * 0.453592                         # Mass at each data point, also converting to kg

V_e_red2 = V_e2 * np.sqrt(W_s/(mass * g))

aoa = data2['a'].iloc[1:].to_numpy(dtype=float)
defl = data2['de'].iloc[1:].to_numpy(dtype=float)

gradient, intercept, r_value, p_value, std_err = stats.linregress(aoa, defl)
#mn = np.min(aoa)
#mx = np.max(aoa)
#x1 = np.linspace(mn,mx,500)
#y1 = gradient*x1+intercept
#plt.scatter(aoa, defl, marker='x')
#plt.plot(x1, y1,'-g')
#plt.show()

C_m_alpha = -C_m_deltae_rad * gradient
print('C_m_alpha in radians is', C_m_alpha)


f = open('thrust_dats/thr_trim_ours.dat', 'r')
thrust = f.read().split()
thrust = np.array(thrust, dtype=float)
f.close()

f = open('thrust_dats/thr_trim_ours_std.dat', 'r')
thrust_std = f.read().split()
thrust_std = np.array(thrust_std, dtype=float)
f.close()

thr_coeff = np.zeros(int(0.5 * len(thrust)))
thr_coeff_std = np.zeros(int(0.5 * len(thrust)))
for i in range(int(len(thrust) * 0.5)):
    thr_coeff[i] = (thrust[2*i] + thrust[2*i+1]) / (0.5 * rho0 * (V_e_red2[i] ** 2) * (d ** 2))
    thr_coeff_std[i] = (thrust_std[2*i] + thrust_std[2*i+1]) / (0.5 * rho0 * (V_e_red2[i] ** 2) * (d ** 2))

delta_e_red = defl - (C_m_T_c/C_m_deltae) * (thr_coeff_std - thr_coeff)
print(delta_e_red)
print(V_e_red2)

def inv_f (x, a, b):
    return (a + (b/x ** 2))

popt, pcov = optimize.curve_fit(inv_f, V_e_red2, delta_e_red)
xlab = ('Reduced velocity (\u1E7C$_e$) [m/s]')
ylab = ('Reduced deflection ($\delta_{e_{eq}}^*$) [deg]')
plt.axes(xlim=(50, 200), ylim=(4, -3), xlabel=xlab, ylabel=ylab)
x1 = np.linspace(1, 250, 250)
plt.plot(x1, inv_f(x1, *popt), 'darkturquoise', label='Inverse square least-squares fit')
plt.scatter(V_e_red2, delta_e_red, s=70, c='red', marker='x', label='Reduced elevator deflection data')
plt.plot(x1, (-Cm0/C_m_deltae) * np.ones(250) * 180/np.pi, '-g')
plt.plot(x1, popt[0] * np.ones(250))
plt.grid()
plt.legend(loc='best')
plt.show()
