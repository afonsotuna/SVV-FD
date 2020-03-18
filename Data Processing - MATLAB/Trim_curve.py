import pandas as Pd
import numpy as np
import matplotlib.pyplot as plt
import xlrd
from scipy import stats
from Reduction_function import reduce

rho0 = 1.2250
g = 9.81
S = 30.00
BEM = 9165    # Basic empty mass, taken from mass report for 2020 [lbs]
c = 2.0569    # MAC, taken from cit_par.py

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

V_e, T, M = reduce(hp, V_ias, Tm)
F_used = data['F. used'].iloc[1:].to_numpy(dtype=float)  # Same as what's done above but for fuel used

pax_masses = np.array(Pd.read_excel(exc_data, header=None, usecols='H', skiprows=7, nrows=9))
init_fuel = xlrd.open_workbook(exc_data).sheet_by_index(0).cell_value(17, 3)
ramp_mass = (BEM + init_fuel)*0.453592 + np.sum(pax_masses)  # Calculating total ramp mass, converting some lbs terms to kg
mass = ramp_mass - F_used * 0.453592                         # Mass at each data point, also converting to kg

defl = (mass * g) / (0.5 * S * rho0 * V_e ** 2 )
avgdefl = np.average(defl) # The two are actually almost the same

C_m_deltae = -(avgdefl * (cg2 - cg1) * 0.0254) / ((delta_e[1] - delta_e[0]) * c)
C_m_deltae_rad = C_m_deltae * (180/np.pi)

#print(C_m_deltae)
print('C_m_delta_e in radians is', C_m_deltae_rad)

### Now the second block of data (line 56 in excel) has to be read to find the gradient of the trim curve to find C_m_alpha

data2 = Pd.read_excel(exc_data, header=55, usecols='B, D:M', skiprows=[57], nrows=8)

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