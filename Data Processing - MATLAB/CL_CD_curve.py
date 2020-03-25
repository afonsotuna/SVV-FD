import matplotlib.pyplot as plt
from scipy import stats, optimize
import numpy as np
from CL_curve import CL_curve
from CD_curve import CD_curve

c = 2.0569    # MAC, taken from cit_par.py [m]
b      = 15.911	          # wing span [m]

def CLCD (ref):
    CL = CL_curve(ref)[2]
    CD = CD_curve(ref)[1]
    CLsq = CL ** 2
    plt.scatter(CL, CD, marker='x')
    plt.scatter(CLsq, CD, marker='o')

    gradient, intercept, r_value, p_value, std_err = stats.linregress(CLsq, CD)
    plt.plot(np.linspace(0,1,10), intercept + gradient * np.linspace(0,1,10))
    plt.show()
    print(gradient, intercept, r_value)
    return (CL)
print('AR' , 30 / c ** 2, b ** 2/30)
CLCD(0)