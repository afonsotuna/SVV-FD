import matplotlib.pyplot as plt
from CL_curve import CL_curve
from CD_curve import CD_curve

def CLCD (ref):
    CL = CL_curve(ref)[2]
    aoa, CD = CD_curve(ref)
    CLsq = CL ** 2
    plt.scatter(CL, CD, marker='x')
    plt.scatter(CLsq, CD, marker='o')
    plt.show()
    return (CL)

CLCD(0)