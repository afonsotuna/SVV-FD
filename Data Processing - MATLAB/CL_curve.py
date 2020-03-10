import pandas as Pd
import numpy as np
import matplotlib.pyplot as plt
import xlrd

# Some needed constants, copy-pasted from the provided cit_par.py file

rho0   = 1.2250          # air density at sea level [kg/m^3]
lamb   = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
p0 = 101325              # Pressure at sea level in ISA [Pa]
gamm = 1.4               # Ratio of specific heats [-]
S      = 30.00	         # wing area [m^2]

BEM = 9165               # Basic empty mass, taken from mass report for 2020 [lbs]

# reading data unto pandas dataframe:
ref_data = 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx'
data = Pd.read_excel(ref_data, header=24, usecols='B, D:J', skiprows=[26], nrows=7)
#print(data)

hp = data['hp'].iloc[1:].to_numpy(dtype=float)  # Putting pressure alt (hp) values from the dataframe column to np array
hp = hp * 0.3048                                      # Converting to meters
p = p0 * (1 + (lamb * hp / Temp0)) ** (-g / (lamb * R))  # Calculating air pressure values from pressure altitude hp

V_ias = data['IAS'].iloc[1:].to_numpy(dtype=float)    # Putting V_IAS values from the dataframe column to np array
V_cas = V_ias - 2                                     # Conversion from IAS to CAS according to appendix A
V_cas_SI = V_cas * 0.51444444                         # Converting to m/s

# Huge formula for calculating the mach number, also from the assignment
M = np.sqrt( (2/(gamm - 1) * ( ( 1 + (p0/p) * ((1 + ((gamm - 1) / (2*gamm)) * rho0/p0 * V_cas_SI ** 2) ** (gamm/(gamm - 1)) - 1 )) ** ((gamm-1)/gamm) - 1) ) )

Tm = data['TAT'].iloc[1:].to_numpy(dtype=float)     # Putting TAT values from the dataframe column to np array
T = (Tm + 273.15) / (1 + (M ** 2) * (gamm-1)/2)     # Correcting the measured TAT for ram rise to find SAT, important
                                                    # to convert to Kelvin here

a = np.sqrt(gamm * R * T)                           # Calculating the speed of sound at each point
V_tru = M*a                                         # Calculating true airspeed
V_e = V_tru * np.sqrt( (p/(R*T))/rho0 )             # Calculating equivalent airspeed

F_used = data['F. used'].iloc[1:].to_numpy(dtype=float)  # Same as what's done above but for fuel used

pax_masses = np.array(Pd.read_excel(ref_data, header=None, usecols='H', skiprows=7, nrows=9))
init_fuel = xlrd.open_workbook(ref_data).sheet_by_index(0).cell_value(17, 3)
ramp_mass = (BEM + init_fuel)*0.453592 + np.sum(pax_masses)  # Calculating total ramp mass, converting some lbs terms to kg
mass = ramp_mass - F_used * 0.453592                         # Mass at each data point, also converting

C_L = (mass * g) / (0.5 * S * rho0 * V_e ** 2 )
aoa = data['a'].iloc[1:].to_numpy(dtype=float)

plt.scatter(aoa, C_L)
plt.show()