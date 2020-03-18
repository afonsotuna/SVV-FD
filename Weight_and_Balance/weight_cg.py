import numpy as np
import matplotlib.pyplot as plt
from scipy import io

# ======================================= WEIGHT AND CENTRE OF GRAVITY AS FUNCTIONS OF TIME ============================

# UNITS:
# Masses in [lbs]
# Moments in [lbs-in/100]

# CONVERSION FACTORS
in_to_m = 2.54E-2
lbs_to_kg = 0.453592

# BASIC EMPTY MASS SPECIFICATION
masses = {'BEM': 9165.0}
moments = {'BEM': 9165.0 * 291.65 / 100}

# INITIAL PAYLOAD SPECIFICATION
with open('payload.txt') as f:
    # INPUT: file with position, arm, mass [kg]
    # OUTPUT: moments and masses dictionaries in lbs-in/100 and lbs respectively
    for line in f:
        line = line.strip()
        line = line.split(',')
        label = line[0]
        moments[label] = int(line[1]) * float(line[2]) / lbs_to_kg / 100
        masses[label] = float(line[2]) / lbs_to_kg

# INTERPOLATOR FOR FUEL MOMENT AS A FUNCTION OF FUEL WEIGHT
x = []
y = []

with open('fuel.txt') as f:
    for line in f:
        line = line.strip()
        line = line.split(',')
        x.append(int(line[0]))
        y.append(float(line[1]))

fuel_moment = lambda m: np.interp(m, x, y)

# INITIAL FUEL SPECIFICATION
fuel_init = 2700.0
masses['fuel'] = fuel_init
moments['fuel'] = fuel_moment(fuel_init)

# IMPORTING FUEL USED DATA
mat = io.loadmat('clean_flight_data.mat')
flight_data = mat['clean_data']
left_FU = flight_data[:, 14]
right_FU = flight_data[:, 15]
time = flight_data[:, 48]

fuel_used = lambda t: np.interp(t, time, left_FU) + np.interp(t, time, right_FU)

# DEFINITION OF XLEMAC AND XTEMAC to calculate MAC [in]
xlemac = 261.45
xtemac = 342.43
mac = xtemac - xlemac

# DEFINE MASS [kg] AS A FUNCTION OF TIME [s]
def mass(t):
    masses['fuel'] = fuel_init
    moments['fuel'] = fuel_moment(fuel_init)
    if t >= 9:
        fuel = fuel_init - fuel_used(t)
        assert fuel >= 0, 'Block fuel depleted at given time.'
        masses['fuel'] = fuel
        moments['fuel'] = fuel_moment(fuel)
        return f'{round(sum(masses.values()) * lbs_to_kg, 3)} kg'
    else:
        return f'{round(sum(masses.values()) * lbs_to_kg, 3)} kg'


# DEFINE CENTRE OF GRAVITY LOCATION [% MAC] AS A FUNCTION OF TIME [s]
# DATUM: LEMAC
def cg(t):
    masses['fuel'] = fuel_init
    moments['fuel'] = fuel_moment(fuel_init)
    moments['seat8'] = 86 * 288 / lbs_to_kg / 100
    if t >= 9:
        fuel = fuel_init - fuel_used(t)
        assert fuel >= 0, 'Block fuel depleted at given time.'
        masses['fuel'] = fuel
        moments['fuel'] = fuel_moment(fuel)
    if 47 * 60 <= t < 49 * 60:
        moments['seat8'] = 86 * 131 / lbs_to_kg / 100
    cg = sum(moments.values()) * 100 / sum(masses.values())
    return f'{round((cg - xlemac) / mac * 100, 3)} % MAC'