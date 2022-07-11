
import const as c

import numpy as np


def extrapolate(delta_t, Hf, Lu):
    return Hf - Lu * (np.exp(c.half_life * (delta_t)) - 1)

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


def epsilon(t, Hf, Lu): 
    return (extrapolate(t, Hf, Lu) / extrapolate(t, c.chur_hh, c.chur_lh) - 1) * 1e4

def solve_T(u_pb_age, hf_hf, lu_hf,
            Hf_am=c.am_hh, Lu_am=c.am_lh, Lu_crustal=c.ave_tuple[0]):
    # Method: get Hf_am_T by T-> 0
    Hf_sample_t = extrapolate(u_pb_age, hf_hf, lu_hf)
    T = 1/c.half_life * np.log((np.exp(c.half_life * u_pb_age) * 
        (Hf_am + Lu_am - Hf_sample_t - Lu_crustal)) / 
        (np.exp(c.half_life * u_pb_age) * Lu_am - Lu_crustal) )
    return np.where(u_pb_age > T, u_pb_age, T)
