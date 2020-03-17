import pandas as Pd
import numpy as np

#This file creates the matlab.dat needed for the thrust.exe file

rho0   = 1.2250          # air density at sea level [kg/m^3]
lamb   = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
p0 = 101325              # Pressure at sea level in ISA [Pa]
gamm = 1.4               # Ratio of specific heats [-]

ref_data = 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx'
data = Pd.read_excel(ref_data, header=55, usecols='B, D:M', skiprows=[57], nrows=8)
print(data)

hp = data['hp'].iloc[1:].to_numpy(dtype=float)
hp = hp * 0.3048                                      # Converting to meters
p = p0 * (1 + (lamb * hp / Temp0)) ** (-g / (lamb * R))  # Calculating air pressure values from pressure altitude hp


V_ias = data['IAS'].iloc[1:].to_numpy(dtype=float)    # Putting V_IAS values from the dataframe column to np array
V_cas = V_ias - 2                                     # Conversion from IAS to CAS according to appendix A
V_cas_SI = V_cas * 0.51444444                         # Converting to m/s

# Huge formula for calculating the mach number, also from the assignment
M = np.sqrt( (2/(gamm - 1) * ( ( 1 + (p0/p) * ((1 + ((gamm - 1) / (2*gamm)) * rho0/p0 * V_cas_SI ** 2) ** (gamm/(gamm - 1)) - 1 )) ** ((gamm-1)/gamm) - 1) ) )
FFl = data['FFl'].iloc[1:].to_numpy(dtype=float) * (0.453592/3600)
FFr = data['FFr'].iloc[1:].to_numpy(dtype=float) * (0.453592/3600)

Tm = data['TAT'].iloc[1:].to_numpy(dtype=float)     # Putting TAT values from the dataframe column to np array
T = (Tm + 273.15) / (1 + (M ** 2) * (gamm-1)/2)     # Correcting the measured TAT for ram rise to find SAT, important
                                                    # to convert to Kelvin here
delta_T = T - (Temp0 + lamb * hp)


f = open("matlab.dat", "w+")
for i in range (np.shape(M)[0]):
    f.write(str(hp[i]) + ' ' + str(M[i]) + ' ' + str(delta_T[i]) + ' ' + str(FFl[i]) + ' ' + str(FFr[i]) + '\n')
