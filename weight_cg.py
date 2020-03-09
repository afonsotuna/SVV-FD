from scipy import interpolate
import matplotlib.pyplot as plt

#======================================= WEIGHT AND CENTRE OF GRAVITY AS FUNCTIONS OF TIME =======================================

# CONVERSION FACTORS
in_to_m = 2.54E-2
lbs_to_kg = 0.453592

# BASIC EMPTY MASS SPECIFICATION
masses = {'BEM': 9165.0}
moments = {'BEM': 9165.0 * 291.65 / 100}

# PAYLOAD SPECIFICATION
with open('Weight and Balance/payload.txt') as f:
    # Input file: position, arm, mass [kg]
    # Output: moments and masses dictionaries in lbs-in and lbs respectively
    for line in f:
        line = line.strip()
        line = line.split(',')
        label = line[0]
        moments[label] = int(line[1]) * float(line[2]) / lbs_to_kg / 100
        masses[label] = float(line[2]) / lbs_to_kg

# INITIAL FUEL SPECIFICATION
masses['fuel'] = 750.0
moments['fuel'] = 750.0 * 287.58 / 100

# INTERPOLATOR FOR FUEL MOMENT AS A FUNCTION OF FUEL WEIGHT
x = []
y = []

with open('Weight and Balance/fuel.txt') as f:
    for line in f:
        line = line.strip()
        line = line.split(',')
        x.append(int(line[0]))
        y.append(float(line[1]))

fuel = interpolate.interp1d(x, y)












