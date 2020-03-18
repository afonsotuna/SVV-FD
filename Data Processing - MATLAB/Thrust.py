import pandas as Pd
import numpy as np
from Reduction_function import reduce
#This file creates the matlab.dat needed for the thrust.exe file

rho0   = 1.2250          # air density at sea level [kg/m^3]
lamb   = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
p0 = 101325              # Pressure at sea level in ISA [Pa]
gamm = 1.4               # Ratio of specific heats [-]
m_f_std = 0.048          # Standard mass fuel flow [kg/s]

ref = 0      # change to one to run for reference data
std_ff = 1   # change to one to run with standard fuel flow

# reading data unto pandas dataframe:
if ref == 1:
    exc_data = 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx'
else:
    exc_data = 'C:/Users/nano2598\Documents/actual college ugh/Third year/SVV-FD/Data_repo/20200311_V4.xlsx'

data = Pd.read_excel(exc_data, header=24, usecols='B, D:J', skiprows=[26], nrows=7)
print(data)

hp = data['hp'].iloc[1:].to_numpy(dtype=float) * 0.3048

V_ias = data['IAS'].iloc[1:].to_numpy(dtype=float)    # Putting V_IAS values from the dataframe column to np array

FFl = data['FFl'].iloc[1:].to_numpy(dtype=float) * (0.453592/3600)
FFr = data['FFr'].iloc[1:].to_numpy(dtype=float) * (0.453592/3600)

Tm = data['TAT'].iloc[1:].to_numpy(dtype=float)     # Putting TAT values from the dataframe column to np array

V_e, T, M = reduce(hp, V_ias, Tm)
delta_T = T - (Temp0 + lamb * hp)

if std_ff == 1:
    for i in range(np.shape(FFr)[0]):
        FFl[i] = m_f_std
        FFr[i] = m_f_std

f = open("matlab.dat", "w+")
for i in range (np.shape(M)[0]):
    f.write(str(hp[i]) + ' ' + str(M[i]) + ' ' + str(delta_T[i]) + ' ' + str(FFl[i]) + ' ' + str(FFr[i]) + '\n')
