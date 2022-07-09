
import const as c

import numpy as np


def generate_Hf(mother_sample, fixed=None):
    return np.random.choice(mother_sample) if fixed is None else fixed

def Hf2Lu(gen_Hf):
    chur_hh_init = extrapolate(4.56, c.chur_hh, c.chur_lh)
    return (gen_Hf - chur_hh_init) / (np.exp(c.half_life * 4.56) - 1)

def lin_intrapolate(x, y1, y2, turn_points=c.turn_points):
    return ((y2 - y1)*x + (turn_points[1]*y1 - turn_points[0]*y2)) / np.ptp(turn_points)

def get_crust_Lu(oxy, y1_tuple=c.y1_tuple, y2_tuple=c.y2_tuple, ave_tuple=c.ave_tuple,
           Lu_crust_err=True):
    oxy_arr = np.ones(1) * oxy if type(oxy) == float else oxy
    low_branch = np.random.normal(y1_tuple[0], y1_tuple[1]*Lu_crust_err, len(oxy_arr))
    # Create an array whose length is indentical to oxy_arr.
    # Let them be y1+/-sigma_y1.
    # Use low_ind to select those low oxygen elements.
    # Assign them with the low_branch values, and mask the others to be 0.
    high_branch = np.random.normal(y2_tuple[0], y2_tuple[1] * Lu_crust_err, len(oxy_arr))
    low_ind, high_ind = oxy_arr < c.turn_points[0], oxy_arr > c.turn_points[1]
    # If nan, enter nan_branch
    nan_branch = np.random.normal(ave_tuple[0], ave_tuple[1] * Lu_crust_err, len(oxy_arr))
    nan_ind = np.isnan(oxy_arr)
    # For the middle branch, bootstrapping module has already been built.
    # Directly call lin_intrapolate function will automatically use it.
    mid_branch = lin_intrapolate(oxy_arr, np.random.normal(y1_tuple[0], y1_tuple[1] * Lu_crust_err),
                                 np.random.normal(y2_tuple[0], y2_tuple[1] * Lu_crust_err))
    res_arr = (np.where(low_ind, low_branch, 0.) + np.where(high_ind, high_branch, 0.) +
        np.where(nan_ind, nan_branch, 0.) + np.where(~low_ind & ~high_ind & ~nan_ind, mid_branch, 0.))
    return float(res_arr) if type(oxy) == float else res_arr

def extrapolate(delta_t, Hf, Lu):
    return Hf - Lu * (np.exp(c.half_life * (delta_t)) - 1)

def epsilon(t, Hf, Lu): 
    return (extrapolate(t, Hf, Lu) / extrapolate(t, c.chur_hh, c.chur_lh) - 1) * 1e4

def solve_T(t, Hf_arr, Lu_arr,
            Hf_am=c.am_hh, Lu_am=c.am_lh, Lu_crustal=0.021, method='direct'):
    '''
    Parameters
    ----------
    method: str. If method == 'direct', get Hf_am_T from 0 directly.
            Else, two steps. Two methods return the same result.
    '''
    # Method1: get Hf_am_T by T-> t
    if method != 'direct':
        Hf_sample_t = extrapolate(t, Hf_arr, Lu_arr)
        Hf_am_t = extrapolate(t, Hf_am, Lu_am)
        delta_Hf = Hf_sample_t - Hf_am_t
        delta_Lu = Lu_crustal - (Lu_am * np.exp(c.half_life * t))
        T = t + 1/c.half_life * np.log(1 + delta_Hf / delta_Lu)
        # For grains that lie above the AM curve, use U-Pb age to represent its model age
        return np.where(t > T, t, T)
    # Method2: get Hf_am_T by T-> 0
    else:
        Hf_sample_t = extrapolate(t, Hf_arr, Lu_arr)
        T = 1/c.half_life * np.log((np.exp(c.half_life * t) * 
            (Hf_am + Lu_am - Hf_sample_t - Lu_crustal)) / 
            (np.exp(c.half_life * t) * Lu_am - Lu_crustal) )
        return np.where(t > T, t, T)

def T_boot(t, Hf_arr, Lu_arr, oxy_arr,
          mother_sample=None, percentiles=(2.5, 97.5)):
    matrix = []
    for i in range(100):
        Hf_am = generate_Hf(mother_sample)
        Lu_am = Hf2Lu(Hf_am)
        matrix.append(solve_T(t, Hf_arr, Lu_arr, oxy_arr, Hf_am=Hf_am, Lu_am=Lu_am))

    print("Median age is ", np.percentile(matrix, 50, axis=0))
    print(percentiles[0], "% age is", np.percentile(matrix, percentiles[0], axis=0))
    print(percentiles[1], "% age is", np.percentile(matrix, percentiles[1], axis=0))


def calc_epsilon_error(data, times=100, meas_err=1, percentiles=(2.5, 97.5)):
    '''
    Calcuate error of epsilon_hf based on measurement errors
    '''
    t0, Hf0, Lu0 = data.u_pb_age / 1e3, data.hf_hf, data.lu_hf
    t_err = data.age_2se / 1e3 / 2
    Hf_err, Lu_err = data.hf_hf_2se / 2, data.lu_hf_2se / 2

    epsilon_list = []
    for i in range(times):
        t = np.random.normal(t0, t_err * meas_err)
        Hf = np.random.normal(Hf0, Hf_err * meas_err)
        Lu = np.random.normal(Lu0, Lu_err * meas_err)
        epsilon_hf = epsilon(t, Hf, Lu)
        epsilon_list.append(epsilon_hf)
        
    return (np.percentile(epsilon_list, 50, axis=0),
            np.percentile(epsilon_list, percentiles[0], axis=0),
            np.percentile(epsilon_list, percentiles[1], axis=0))

def calc_T(data, mother_sample, times=5000, percentiles=(2.5, 97.5), 
           meas_err=True, Hf_am_err=True, Lu_crust_err=True, fixed_Lu=False):
    t0, Hf0 = data.u_pb_age / 1e3, data.hf_hf
    Lu0, oxy0 = data.lu_hf, data.o
    t_err, Hf_err = data.age_2se / 1e3 / 2, data.hf_hf_2se / 2
    Lu_err, oxy_err = data.lu_hf_2se / 2, data.o_sed

    T_list = []
    for i in range(times):
        t = np.random.normal(t0, t_err * meas_err)
        Hf = np.random.normal(Hf0, Hf_err * meas_err)
        Lu = np.random.normal(Lu0, Lu_err * meas_err)
        oxy = np.random.normal(oxy0, oxy_err * meas_err)

        fixed = None if Hf_am_err else c.am_hh
        Hf_am = generate_Hf(mother_sample, fixed=fixed)
        Lu_am = Hf2Lu(Hf_am)

        if fixed_Lu:
            Lu_crustal = np.where(t <= 2.5, c.ave_tuple[0], 0.) + np.where(t > 2.5, c.y1_tuple[0], 0.)
        else:
            Lu_crustal = get_crust_Lu(oxy, Lu_crust_err=Lu_crust_err)

        T = solve_T(t, Hf, Lu, Hf_am=Hf_am, Lu_am=Lu_am, Lu_crustal=Lu_crustal)
        T_list.append(T)

    return (np.percentile(T_list, 50, axis=0),
            np.percentile(T_list, percentiles[0], axis=0),
            np.percentile(T_list, percentiles[1], axis=0))

