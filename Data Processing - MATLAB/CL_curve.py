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
mu_0 = 1.7894 * 10 ** -5
BEM = 9165               # Basic empty mass, taken from mass report for 2020 [lbs]

def CL_curve (ref, data=None, plot=0):  ## Use ref = 1 for ref data, any other number for our data, you can also override with
                                ## your own set of data

    # reading data unto pandas dataframe:
    if ref == 1 and data is None:
        exc_data = 'Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx'
    elif data is None:
        exc_data = 'C:/Users/nano2598\Documents/actual college ugh/Third year/SVV-FD/Data_repo/20200311_V4.xlsx'
    else:
        exc_data = data
    data = Pd.read_excel(exc_data, header=24, usecols='B, D:J', skiprows=[26], nrows=7)

    #print(data)
    hp = data['hp'].iloc[1:].to_numpy(dtype=float)  # Putting pressure alt (hp) values from the dataframe column to np array
    V_ias = data['IAS'].iloc[1:].to_numpy(dtype=float)    # Putting V_IAS values from the dataframe column to np array
    Tm = data['TAT'].iloc[1:].to_numpy(dtype=float)     # Putting TAT values from the dataframe column to np array

    V_e, T, M, Re_range = reduce(hp, V_ias, Tm)
    F_used = data['F. used'].iloc[1:].to_numpy(dtype=float)  # Same as what's done above but for fuel used

    pax_masses = np.array(Pd.read_excel(exc_data, header=None, usecols='H', skiprows=7, nrows=9))
    init_fuel = xlrd.open_workbook(exc_data).sheet_by_index(0).cell_value(17, 3)
    ramp_mass = (BEM + init_fuel)*0.453592 + np.sum(pax_masses)  # Calculating total ramp mass, converting some lbs terms to kg
    mass = ramp_mass - F_used * 0.453592                         # Mass at each data point, also converting to kg

    V_e_red = V_e * np.sqrt(W_s/(mass * g))
    C_L = (mass * g) / (0.5 * S * rho0 * V_e ** 2 )
    aoa = data['a'].iloc[1:].to_numpy(dtype=float)

    gradient, intercept, r_value, p_value, std_err = stats.linregress(aoa, C_L)

    if plot == 1:
        x1 = np.linspace(0, 14, 500)
        y1 = gradient * x1 + intercept
        fig = plt.figure()
        ax = fig.add_subplot(111)
        xlab = ('Angle of attack (AoA) [deg]')
        ylab = ('Lift coefficient (C$_L$) [-]')
        plt.axes(xlim=(0,12) , ylim=(0,1.2), xlabel=xlab, ylabel=ylab)
        plt.xticks(np.arange(0,13,step=1))
        ylist = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
        plt.yticks(np.arange(0,1.3,0.1), labels=ylist)
        plt.scatter(aoa, C_L, s=70, c='r', marker='x', label='Calculated C$_L$ values')
        plt.plot(x1, y1, 'darkturquoise', label='Line of best fit')
        plt.grid()
        plt.legend(loc='best')
        if ref == 1:
            plt.savefig('figs/CL_alpha_ref')
        else:
            plt.savefig('figs/CL_alpha_our_data')
        plt.show()

    #print('C_L_alpha in 1/deg is', gradient)
    print('C_L_alpha in 1/rad is', gradient * (180 / np.pi))
    #print('r_value of the regression is', r_value)
    print("Re range is from " + np.format_float_scientific(float(Re_range[0])) + ' to ' + np.format_float_scientific(float(Re_range[1])))
    print('Mach range is from ' + str(round(min(M), 2)) + ' to ' + str(round(max(M), 2)))


    return(gradient*(180/np.pi), aoa, C_L)

CL_curve(1, plot=1)
