import numpy as np
import matplotlib.pyplot as plt

#======================================= WEIGHT AND CENTRE OF GRAVITY AS FUNCTIONS OF TIME =======================================

# UNITS:
# Masses in [lbs]
# Moments in [lbs-in]

# CONVERSION FACTORS
in_to_m = 2.54E-2
lbs_to_kg = 0.453592

# BASIC EMPTY MASS SPECIFICATION
masses = {'BEM': 9165.0}
moments = {'BEM': 9165.0 * 291.65 / 100}

# INITIAL PAYLOAD SPECIFICATION
with open('Weight and Balance/payload.txt') as f:
    # INPUT: file with position, arm, mass [kg]
    # OUTPUT: moments and masses dictionaries in lbs-in and lbs respectively
    for line in f:
        line = line.strip()
        line = line.split(',')
        label = line[0]
        moments[label] = int(line[1]) * float(line[2]) / lbs_to_kg / 100
        masses[label] = float(line[2]) / lbs_to_kg

# PRE-PROCESSING FOR FUEL FLOW [lbs/s] AS A FUNCTION OF TIME [s]
time = list(np.arange(0, 9, 0.1))
right_engine = list(np.zeros(90))
left_engine = list(np.zeros(90))

with open('Weight and Balance/fuelflow.csv') as f:
    for line in f:
        line = line.strip()
        line = line.split(',')
        time.append(float(line[0]))
        right_engine.append(float(line[1]) / 3600)
        left_engine.append(float(line[2]) / 3600)

rff = lambda t: np.interp(t, time, right_engine)
lff = lambda t: np.interp(t, time, right_engine)

# INTERPOLATOR FOR FUEL MOMENT AS A FUNCTION OF FUEL WEIGHT
x = []
y = []

with open('Weight and Balance/fuel.txt') as f:
    for line in f:
        line = line.strip()
        line = line.split(',')
        x.append(int(line[0]))
        y.append(float(line[1]))

fuel_moment = lambda m: np.interp(m, x, y)

# INITIAL FUEL SPECIFICATION
masses['fuel'] = 750.0
moments['fuel'] = fuel_moment(750.0)

# RAMP MASS AND CG LOCATION
RM = sum(masses.values())
RCG = sum(moments.values()) * 100 / sum(masses.values())

# DEFINE INTEGRATOR FOR TOTAL FUEL CONSUMED
def integrate(f, a, b):
    x = np.arange(a, b + 0.1, 0.1)
    y = f(x)
    return np.trapz(y, x)

# DEFINE MASS [lbs] AS A FUNCTION OF TIME [s]
def mass(t):
    masses['fuel'] = masses['fuel'] - integrate(rff, 0, t) - integrate(lff, 0 , t)
    assert fuel >= 0, 'Aircraft fuel depleted at given time.'
    return RM - integrate(rff, 0, t) - integrate(lff, 0 , t)

def cg(t):
    masses['fuel'] = masses['fuel'] - integrate(rff, 0, t) - integrate(lff, 0 , t)
    assert fuel >= 0, 'Aircraft fuel depleted at given time.'
    moments['fuel'] = fuel_moment(fuel)
    return sum(moments.values()) * 100 / sum(masses.values())













