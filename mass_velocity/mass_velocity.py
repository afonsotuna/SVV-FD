import scipy.io


def mass_velocity_reference(t):
    # Flight data imported
    mat = scipy.io.loadmat('reference_clean.mat')
    flight_data = mat['clean_data']
    block_fuel = 4050  # block fuel (lbs)
    passenger_weight = 695  # sum of passenger weight (kgs)

    # Get data location
    t_lookup = t - flight_data[0, 47]
    index = int((t_lookup - flight_data[0, 47]) / 0.1)

    # Obtain correct weight (manoeuvre start) and velocity, get system
    used_fuel = flight_data[index, 13] + flight_data[index, 14]
    mass_event = (block_fuel - used_fuel + 9165) * 0.453592 + passenger_weight
    tas_event = flight_data[index, 41] * 0.514444

    return mass_event, tas_event


def mass_velocity_flight(t):
    # Flight data imported
    mat = scipy.io.loadmat('clean_flight_data.mat')
    flight_data = mat['clean_data']
    block_fuel = 2700  # block fuel (lbs)
    passenger_weight = 771  # sum of passenger weight (kgs)

    # Get data location
    t_lookup = t - flight_data[0, 48]
    index = int((t_lookup - flight_data[0, 48]) / 0.1)

    # Obtain correct weight (manoeuvre start) and velocity, get system
    used_fuel = flight_data[index, 14] + flight_data[index, 15]
    mass_event = (block_fuel - used_fuel + 9165) * 0.453592 + passenger_weight
    tas_event = flight_data[index, 42] * 0.514444

    return mass_event, tas_event

short_period_mass, short_period_vel = mass_velocity_reference(3635)
print("short period =", short_period_mass,short_period_vel)

phugoid_mass, phugoid_vel = mass_velocity_reference(3237)
print("phugoid =", phugoid_mass, phugoid_vel)

dutch_roll_mass, dutch_roll_vel = mass_velocity_reference(3717)
print("dutch roll =", dutch_roll_mass, dutch_roll_vel)

spiral_mass, spiral_vel= mass_velocity_reference(3920)
print("spiral =", spiral_mass, spiral_vel)


aperiodic_mass, aperiodic_vel = mass_velocity_reference(3550)
print("aperiodic =", aperiodic_mass,aperiodic_vel)