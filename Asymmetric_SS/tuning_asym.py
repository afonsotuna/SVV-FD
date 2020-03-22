import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
from Asymmetric_SS.main_asymmetric import num_model_asym_reference


def error_def(y1, y2):
    # error calculation
    error_raw = y2 - y1  # first diff

    data_range = np.amax(y1) - np.amin(y1)  # range of data magnitude
    error_norm = abs(error_raw / data_range)  # normalize for range and abs

    error_sum = np.sum(error_norm)
    error_std = error_sum / len(y2)

    return error_std


def error_function(parameters, t_lookup=3717, t_limit=14, block_fuel=4050, passenger_weight=695, b=15.911, rho=1.208):
    Kzz = parameters[0]
    Cn_r = parameters[1]
    CY_b = parameters[2]

    # get data and response vectors - try for Dutch Roll first
    y1_DR, y2_DR, _, _, _, _, _, _ = num_model_asym_reference(t_lookup=t_lookup, t_limit=t_limit, block_fuel=block_fuel,
                                                              passenger_weight=passenger_weight, b=b, rho=rho, Kzz=Kzz,
                                                              Cn_r=Cn_r, CY_b=CY_b)

    error_DR = error_def(y1_DR, y2_DR)
    return error_DR


def error_minimize():
    x0 = np.array(([0.042 ** 0.5], [-0.2061], [-0.7500]))
    result = minimize(error_function, x0)
    Kzz, Cn_r, CY_b = result.x
    return Kzz, Cn_r, CY_b
