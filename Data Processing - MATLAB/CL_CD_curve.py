import matplotlib.pyplot as plt
from scipy import stats, optimize
import numpy as np
from CL_curve import CL_curve
from CD_curve import CD_curve

c = 2.0569    # MAC, taken from cit_par.py [m]
b      = 15.911	          # wing span [m]
S = 30.0          # wing area [m^2]

def CLCD (ref, plot=0):
    CL = CL_curve(ref)[2]
    CD = CD_curve(ref)[1]
    CLsq = CL ** 2
    plt.scatter(CL, CD, marker='x')
    plt.scatter(CLsq, CD, marker='o')

    gradient, intercept, r_value, p_value, std_err = stats.linregress(CLsq, CD)
    plt.plot(np.linspace(0,1,10), intercept + gradient * np.linspace(0,1,10))
    #plt.show()

    AR = b ** 2 / S
    e1 = 1 / (np.pi * AR * gradient)
    CD01 = intercept

    if plot == 1:
        x1 = np.linspace(0, 1, 100)
        y1 = intercept + gradient * x1 ** 2
        xlab = ('Lift coefficient (C$_L$) [-]')
        ylab = ('Drag coefficient (C$_D$) [-]')
        plt.axes(xlim=(0,1) , ylim=(0,0.08), xlabel=xlab, ylabel=ylab)
        xlist = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        plt.xticks(np.arange(0,1.1,step=0.1), labels=xlist)
        ylist = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08]
        plt.yticks(np.arange(0,0.09,0.01), labels=ylist)
        plt.scatter(CL, CD, s=70, c='r', marker='x', label='Calculated C$_L$ values')
        plt.plot(x1, y1, 'darkturquoise', label='Least-squares fit')
        plt.grid()
        plt.legend(loc='best')
        if ref == 1:
            plt.savefig('figs/CL_CD_ref')
        else:
            plt.savefig('figs/CL_CD_our_data')
        plt.show()

    print('CD_0 is ' ,CD01, ' and oswald is ', e1)
    return (CL)

CLCD(0, 1)