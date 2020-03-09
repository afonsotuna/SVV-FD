#======================================= WEIGHT AND CENTRE OF GRAVITY AS FUNCTIONS OF TIME =======================================

# BASIC EMPTY MASS SPECIFICATION
BEM_w = 9165.0
BEM_arm = 291.65
BEM_m = BEM_w * BEM_arm

# INITIAL FUEL SPECIFICATION
fuel_w = 750.0
fuel_arm = 287.58
fuel_m = fuel_w * fuel_arm

# CONVERSION FACTORS
in_to_m = 2.54E-2
lbs_to_kg = 0.453592

# PAYLOAD SPECIFICATION
with open('payload.txt') as f:
    for line in f:
        line = line.strip()
        line = line.split()
        print(line)





