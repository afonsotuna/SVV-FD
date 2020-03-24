import pandas as Pd
import numpy as np
import matplotlib.pyplot as plt
import xlrd
from scipy import stats, optimize
from Reduction_function import reduce

# Some needed constants, copy-pasted from the provided cit_par.py file

rho0   = 1.2250          # air density at sea level [kg/m^3]
lamb   = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
p0 = 101325              # Pressure at sea level in ISA [Pa]
gamm = 1.4               # Ratio of specific heats [-]
S      = 30.00	         # wing area [m^2]
W_s = 60500              # standard weight [N]

BEM = 9165               # Basic empty mass, taken from mass report for 2020 [lbs]

def CD_curve (ref, data=None, datathr=None, plot=0):

    if ref == 1 and data is None:
        exc_data = 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx'
        f = open('thr_CLCD_ref.dat')
    elif data is None:
        exc_data = 'C:/Users/nano2598\Documents/actual college ugh/Third year/SVV-FD/Data_repo/20200311_V4.xlsx'
        f = open('thr_CLCD_ours.dat', 'r')
    else:
        exc_data = data
        f = open(datathr, 'r')
    data = Pd.read_excel(exc_data, header=24, usecols='B, D:J', skiprows=[26], nrows=7)

    hp = data['hp'].iloc[1:].to_numpy(dtype=float)  # Putting pressure alt (hp) values from the dataframe column to np array
    V_ias = data['IAS'].iloc[1:].to_numpy(dtype=float)  # Putting V_IAS values from the dataframe column to np array
    Tm = data['TAT'].iloc[1:].to_numpy(dtype=float)  # Putting TAT values from the dataframe column to np array

    V_e, T, M = reduce(hp, V_ias, Tm)

    thrust = f.read().split()
    thrust = np.array(thrust, dtype=float)
    print(thrust)
    f.close()

    CD = np.zeros(int(0.5 * len(thrust)))
    for i in range(int(0.5 * len(thrust))):
        CD[i] = (thrust[2 * i] + thrust[2 * i + 1]) / (0.5 * rho0 * (V_e[i] ** 2) * S)

    aoa = data['a'].iloc[1:].to_numpy(dtype=float)
    # plt.scatter(aoa, CD)
    # plt.show()

    return aoa, CD
