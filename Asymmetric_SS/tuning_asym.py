import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
from Asymmetric_SS.main_asymmetric import num_model_asym_reference


def error_def(y1, y2):
    # error calculation
    error_raw = y2 - y1  # first diff

    data_range = np.amax(y1) - np.amin(y1)  # range of data magnitude
    error_norm = error_raw / data_range  # normalize for range and abs
    error_norm_sq = error_norm ** 2  # square normalized error
    error_sum = np.sum(error_norm_sq)  # sum normalized square errors
    error = error_sum ** 0.5  # square root normalized error

    error_std = error / len(y2)  # get standard error

    return error_std


def error_function(parameters, block_fuel=4050, passenger_weight=695, b=15.911):
    error_tot = 0
    CY_b = parameters[0]
    Cn_r = parameters[1]
    Cn_p = parameters[2]
    Cl_r = parameters[3]
    Cl_p = parameters[4]

    # Dutch Roll
    for i in range(3):
        output = i + 1
        y1_DR, y2_DR, _, _, _, _, _, _ = num_model_asym_reference(output=output, t_lookup=3717, t_limit=14,
                                                                  eigenmotion="spiral", block_fuel=block_fuel,
                                                                  passenger_weight=passenger_weight, b=b, CY_b=CY_b,
                                                                  Cn_r=Cn_r, Cn_p=Cn_p, Cl_r=Cl_r, Cl_p=Cl_p)
        error_tot += error_def(y1_DR, y2_DR)

    # Aperiodic Roll
    for i in range(3):
        output = i + 1
        y1_AR, y2_AR, _, _, _, _, _, _ = num_model_asym_reference(output=output, t_lookup=3550, t_limit=14,
                                                                  eigenmotion="aperiodic", block_fuel=block_fuel,
                                                                  passenger_weight=passenger_weight, b=b, CY_b=CY_b,
                                                                  Cn_r=Cn_r, Cn_p=Cn_p, Cl_r=Cl_r, Cl_p=Cl_p)
        error_tot += error_def(y1_AR, y2_AR)

    # Spiral
    for i in range(3):
        output = i + 1
        y1_SP, y2_SP, _, _, _, _, _, _ = num_model_asym_reference(output=output, t_lookup=3900, t_limit=120,
                                                                  eigenmotion="spiral", block_fuel=block_fuel,
                                                                  passenger_weight=passenger_weight, b=b, CY_b=CY_b,
                                                                  Cn_r=Cn_r, Cn_p=Cn_p, Cl_r=Cl_r, Cl_p=Cl_p)
        error_tot += error_def(y1_SP, y2_SP)

    return error_tot


def error_minimize(x_bounds):
    x0 = np.array(([-0.7500], [-0.2061], [-0.0602], [0.2376], [-0.7108]))
    result = minimize(error_function, x0, bounds=x_bounds, options={'disp': True})
    CY_b, Cn_r, Cn_p, Cl_r, Cl_p = result.x
    return CY_b, Cn_r, Cn_p, Cl_r, Cl_p


x_bounds = [[-5, 0], [-5, 0], [-5, 0], [0, 5], [-5, 0]]

print(error_minimize(x_bounds))
