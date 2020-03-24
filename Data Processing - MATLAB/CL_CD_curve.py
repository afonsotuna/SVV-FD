import matplotlib.pyplot as plt
from scipy import stats, optimize
import numpy as np
from CL_curve import CL_curve
from CD_curve import CD_curve

def CLCD (ref):
    CL = CL_curve(ref)[2]
    aoa, CD = CD_curve(ref)
    CLsq = CL ** 2
    plt.scatter(CL, CD, marker='x')
    plt.scatter(CLsq, CD, marker='o')

    gradient, intercept, r_value, p_value, std_err = stats.linregress(CLsq, CD)
    plt.plot(np.linspace(0,1,10), intercept + gradient * np.linspace(0,1,10))
    plt.show()
    return (CL)

CLCD(0)