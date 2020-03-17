import pandas as Pd
import numpy as np
import matplotlib.pyplot as plt
import xlrd
from scipy import stats
from Reduction_function import reduce
from Weight_and_Balance.weight_cg import cg

rho0 = 1.2250
g = 9.81
S = 30.00
BEM = 9165    # Basic empty mass, taken from mass report for 2020 [lbs]
c = 2.0569    # MAC, taken from cit_par.py

ref_data = 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx'
data = Pd.read_excel(ref_data, header=72, usecols='B, D:M', nrows=3)
cg1 = xlrd.open_workbook(ref_data).sheet_by_index(0).cell_value(70, 2)
cg2 = xlrd.open_workbook(ref_data).sheet_by_index(0).cell_value(70, 7)

hp = data['hp'].iloc[1:].to_numpy(dtype=float)  # Putting pressure alt (hp) values from the dataframe column to np array
V_ias = data['IAS'].iloc[1:].to_numpy(dtype=float)    # Putting V_IAS values from the dataframe column to np array
Tm = data['TAT'].iloc[1:].to_numpy(dtype=float)     # Putting TAT values from the dataframe column to np array
delta_e = data['de'].iloc[1:].to_numpy(dtype=float)

V_e, T, M = reduce(hp, V_ias, Tm)
F_used = data['F. used'].iloc[1:].to_numpy(dtype=float)  # Same as what's done above but for fuel used

pax_masses = np.array(Pd.read_excel(ref_data, header=None, usecols='H', skiprows=7, nrows=9))
init_fuel = xlrd.open_workbook(ref_data).sheet_by_index(0).cell_value(17, 3)
ramp_mass = (BEM + init_fuel)*0.453592 + np.sum(pax_masses)  # Calculating total ramp mass, converting some lbs terms to kg
mass = ramp_mass - F_used * 0.453592                         # Mass at each data point, also converting to kg

C_L = (mass * g) / (0.5 * S * rho0 * V_e ** 2 )
avgC_L = np.average(C_L) # The two are actually almost the same

C_m_deltae = (avgC_L * (cg2 - cg1) * 0.0254) / ((delta_e[1] - delta_e[0]) * c)
C_m_deltae_rad = C_m_deltae * (180/np.pi)

print(cg(46.9*60))
print(cg(47.5*60))

