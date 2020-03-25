import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
from Symmetric_SS.main_symmetric import num_model_sym_data


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


def error_function_sym(parameters, block_fuel=2700, passenger_weight=771, c=2.0569):
    error_tot = 0
    C_x_q = parameters[0]
    C_z_q = parameters[1]
    C_m_alpha = parameters[2]
    C_m_delta_e = parameters[3]
    C_m_q = parameters[4]

    # Phugoid
    for output in range(4):
        y1_PG, y2_PG, _, _, _, _, _ = num_model_sym_data(output=output, t_lookup=3219, t_limit=170,
                                                         block_fuel=block_fuel,
                                                         passenger_weight=passenger_weight, c=c,
                                                         C_x_q=C_x_q, C_z_q=C_z_q,
                                                         C_m_alpha=C_m_alpha,
                                                         C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
        error_tot += error_def(y1_PG, y2_PG)

    # Short Period
    for output in range(4):
        y1_SP, y2_SP, _, _, _, _, _ = num_model_sym_data(output=output, t_lookup=3124.6, t_limit=3.5,
                                                         block_fuel=block_fuel,
                                                         passenger_weight=passenger_weight, c=c,
                                                         C_x_q=C_x_q, C_z_q=C_z_q,
                                                         C_m_alpha=C_m_alpha,
                                                         C_m_delta_e=C_m_delta_e, C_m_q=C_m_q)
        error_tot += error_def(y1_SP, y2_SP)

    return error_tot


def error_minimize_sym(x_bounds):
    x0 = np.array(([-0.2817], [-5.6629], [-0.7249], [-1.4968], [-8.7941]))
    result = minimize(error_function_sym, x0, bounds=x_bounds, options={'disp': True})
    C_x_q, C_z_q, C_m_alpha, C_m_delta_e, C_m_q = result.x
    return C_x_q, C_z_q, C_m_alpha, C_m_delta_e, C_m_q


x_bounds = [[-20, 0], [-20, 0], [-5, 0], [-5, 0], [-20, 0]]

print(error_minimize_sym(x_bounds))
