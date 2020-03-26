import pandas as Pd
import numpy as np
import matplotlib.pyplot as plt
import xlrd
from scipy import stats
# Some needed constants, copy-pasted from the provided cit_par.py file

rho0   = 1.2250          # air density at sea level [kg/m^3]
lamb   = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
p0 = 101325              # Pressure at sea level in ISA [Pa]
gamm = 1.4               # Ratio of specific heats [-]
S      = 30.00	         # wing area [m^2]
c = 2.0569    # MAC, taken from cit_par.py [m]
mu_0 = 1.7894 * 10 ** -5 # viscosity at ISA sea-level, [kg/(m*s)]
W_s = 60500

# Function that takes in values of pressure altitude, indicated airspeed, and TAT,
# and outputs the corresponding equivalent airspeed (not reduced), mach, and static temp.
def reduce (hp, V_ias, Tm):

    hp = hp * 0.3048                                      # Converting to meters
    p = p0 * (1 + (lamb * hp / Temp0)) ** (-g / (lamb * R))  # Calculating air pressure values from pressure altitude hp

    V_cas = V_ias - 2                                     # Conversion from IAS to CAS according to appendix A
    V_cas_SI = V_cas * 0.51444444                         # Converting to m/s

    # Huge formula for calculating the mach number, also from the assignment
    M = np.sqrt( (2/(gamm - 1) * ( ( 1 + (p0/p) * ((1 + ((gamm - 1) / (2*gamm)) * rho0/p0 * V_cas_SI ** 2) ** (gamm/(gamm - 1)) - 1 )) ** ((gamm-1)/gamm) - 1) ) )

    T = (Tm + 273.15) / (1 + (M ** 2) * (gamm-1)/2)     # Correcting the measured TAT for ram rise to find SAT, important
                                                        # to convert to Kelvin here

    a = np.sqrt(gamm * R * T)                           # Calculating the speed of sound at each point
    V_tru = M*a                                         # Calculating true airspeed
    V_e = V_tru * np.sqrt( (p/(R*T))/rho0 )             # Calculating equivalent airspeed

    mu = mu_0 * ((T / Temp0) ** (3 / 2)) * ((Temp0 + 110) / (T + 110))
    Re = ((p/(R * T)) * V_tru * c) / mu
    Re_range = [np.format_float_positional(min(Re), 2, fractional=False), np.format_float_positional(max(Re), 3, fractional=False)]

    return (V_e, T, M, Re_range)