import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
from Asymmetric_SS.main_asymmetric import num_model_asym_data


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


def error_function_asym(parameters, block_fuel=2700, passenger_weight=771):
    error_tot = 0
    CY_b = parameters[0]
    Cn_r = parameters[1]
    Cn_p = parameters[2]
    Cl_r = parameters[3]
    Cl_p = parameters[4]

    # Dutch Roll
    for i in range(3):
        output = i + 1
        y1_DR, y2_DR, _, _, _, _, _, _ = num_model_asym_data(output=output, t_lookup=3455, t_limit=15,
                                                             eigenmotion="dutch roll", block_fuel=block_fuel,
                                                             passenger_weight=passenger_weight, CY_b=CY_b,
                                                             Cn_r=Cn_r, Cn_p=Cn_p, Cl_r=Cl_r, Cl_p=Cl_p)
        errorDR = error_def(y1_DR, y2_DR)
        # print("Dutch roll error: ", errorDR)
        error_tot += errorDR

    # Aperiodic Roll
    for i in range(3):
        output = i + 1
        y1_AR, y2_AR, _, _, _, _, _, _ = num_model_asym_data(output=output, t_lookup=3050, t_limit=20,
                                                             eigenmotion="aperiodic", block_fuel=block_fuel,
                                                             passenger_weight=passenger_weight, CY_b=CY_b,
                                                             Cn_r=Cn_r, Cn_p=Cn_p, Cl_r=Cl_r, Cl_p=Cl_p)
        errorAR = error_def(y1_AR, y2_AR)
        # print("Aperiodic roll error: ", errorAR)
        error_tot += errorAR

    # Spiral
    for i in range(3):
        output = i + 1
        y1_SP, y2_SP, _, _, _, _, _, _ = num_model_asym_data(output=output, t_lookup=3590, t_limit=120,
                                                             eigenmotion="spiral", block_fuel=block_fuel,
                                                             passenger_weight=passenger_weight, CY_b=CY_b,
                                                             Cn_r=Cn_r, Cn_p=Cn_p, Cl_r=Cl_r, Cl_p=Cl_p)
        errorSP = error_def(y1_SP, y2_SP)
        # print("Spiral error: ", errorSP)
        error_tot += errorSP

    return error_tot


def error_minimize_asym(x_bounds):
    x0 = np.array(([-0.7500], [-0.2061], [-0.0602], [0.2376], [-0.7108]))
    result = minimize(error_function_asym, x0, bounds=x_bounds, options={'disp': True})
    CY_b, Cn_r, Cn_p, Cl_r, Cl_p = result.x
    return CY_b, Cn_r, Cn_p, Cl_r, Cl_p


x_bounds = [[-5, 0], [-5, 0], [-5, 0], [0, 5], [-5, 0]]

# print(error_minimize_asym(x_bounds))

# pars= np.array(([-2.662089857595346], [0.0], [0.0], [0.07580034708397379], [-0.7020261870931681]))
# pars = np.array(([-0.7500], [-0.2061], [-0.0602], [0.2376], [-0.7108]))
# error_function_asym(pars)
