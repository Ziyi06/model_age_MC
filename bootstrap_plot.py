
import const as c
from utils import epsilon, extrapolate, solve_T, get_crust_Lu, Hf2Lu

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def T_boot(u_pb_age, hf_hf, lu_hf, oxygen=0, 
           u_pb_age_err=0, hf_hf_err=0, lu_hf_err=0, oxygen_err=0,
           mantle_hh=c.am_hh, times=500, 
           meas_err=True, mantle_hh_err=True, Lu_crust_err=True, 
           plot=False):
    '''
    Calculate Hf model ages and errors using bootstrapping.

    Parameters
    ----------
    u_pb_age: float or array-like
        Measured zircon U-Pb age in a unit of Gyr.

    hf_hf: float or array-like
        Measured zircon 176Hf/177Hf.

    lu_hf: float or array-like
        Measured zircon 176Lu/177Hf.

    oxygen: float or array-like, optional (defaulted 0)
        Measured zircon δ18O in a unit of ‰.
        Note that if this value is 0, crustal Lu/Hf will be constant values; 
        otherwise, crustal Lu/Hf will depend on measured δ18O values. 

    u_pb_age_err, hf_hf_err, lu_hf_err, oxygen_err: float or array-like, optional 
        (defaulted 0)
        Measurement 1 standard errors. 

    mantle_hh: float, tuple, or array-like, optional (defaulted 0.283160 for 
        modern arc mantle)
        Modern mantle 176Hf/177Hf (either arc mantle or depleted mantle).
        Use float if there is only one preferred value.
        Use the tuple shape of (mean, std) if mean and standard deviation are known.
        Use array-like if a distribution is available (e.g., Iizuka et al. 2013).
    
    times: int, optional (defaulted 500)
        Repetition times for bootstrapping, defaulted to be 5000.

    meas_err: bool, optional (defaulted True)
        If True, the calculation of model age errors considers measurement error.

    mantle_hh_err: bool, optional (defaulted True)
        If True, the calculation of model age errors considers uncertainty in the
        Hf isotope of the mantle reservoir.

    Lu_crust_err: bool, optional (defaulted True)
        If True, the calculation of model age errors considers uncertainty in the 
        176Lu/177Hf for the crustal source region.

    plot: bool, optional (defaulted False)
        If True, plot εHf(t) versus U-Pb age for the multiple bootstrapping trials.

    Returns
    -------
    T_list: array
        Calculated Hf model ages, with a length of times (defaulted to 500).
    
    epsilon_list: array
        Calculated εHf(t), with a length of times (defaulted to 500). 

    '''
    if plot:
        plt.figure()
        ax = plt.subplot(111)
    
    T_arr = np.zeros(times)
    eps_arr = np.zeros(times)
    for i in range(times):
        t = np.random.normal(u_pb_age, u_pb_age_err * meas_err)
        Hf = np.random.normal(hf_hf, hf_hf_err * meas_err)
        Lu = np.random.normal(lu_hf, lu_hf_err * meas_err)
        oxy = np.random.normal(oxygen, oxygen_err * meas_err)

        if type(mantle_hh) == tuple:  
            assert len(mantle_hh) == 2
            Hf_am = np.random.normal(mantle_hh[0], mantle_hh[1] * mantle_hh_err)
        elif type(mantle_hh) == float:  
            Hf_am = mantle_hh
        else:  
            Hf_am = np.random.choice(mantle_hh) if mantle_hh_err else np.mean(mantle_hh)
        Lu_am = Hf2Lu(Hf_am)

        if np.any(oxy): 
            Lu_crustal = get_crust_Lu(oxy, Lu_crust_err=Lu_crust_err)
        else:
            Lu_crustal = np.where(t <= 2.5, c.ave_tuple[0], 0.) + np.where(t > 2.5, c.y1_tuple[0], 0.)

        T = solve_T(t, Hf, Lu, Hf_am=Hf_am, Lu_am=Lu_am, Lu_crustal=Lu_crustal)
        T_arr[i] = T
        eps_arr[i] = epsilon(t, Hf, Lu)

        if plot:
            x1 = np.arange(0, t, .001) # zircon 
            ax.plot(x1, epsilon(x1, Hf, Lu), color='#367db7', alpha=.01, lw=2) 

            x2 = np.arange(t, T, .001) # crustal 
            ax.plot(x2, (extrapolate(x2 - t, extrapolate(t, Hf, Lu), Lu_crustal) /
                        extrapolate(x2, c.chur_hh, c.chur_lh) - 1) * 1e4, 
                        color='g', alpha=.01, lw=2, zorder=6) 
            
            x3 = np.arange(0, 4., .001) # AM curve
            ax.plot(x3, epsilon(x3, Hf_am, Lu_am), color='#df1d27', alpha=.01, lw=2, zorder=2) 
            ax.vlines(x=T, ymin=-50., ymax=epsilon(T, Hf_am, Lu_am), 
                        color='dimgrey', alpha=.01, lw=2, linestyle="-", zorder=3) # Model age intercepts

    if plot:
        axins = inset_axes(ax, width="40%", height="40%", loc=4, borderpad=1.3)
        axins.hist(T_arr, bins=30, linewidth=.8, histtype='step', color='k')
        axins.axvline(x=np.percentile(T_arr, 50), ls='--', lw=1, color='k')
        axins.axvline(x=np.percentile(T_arr, 2.5), ls='--', lw=1, color='k')
        axins.axvline(x=np.percentile(T_arr, 97.5), ls='--', lw=1, color='k')
        axins.tick_params(axis='both', which='major', length=1, labelsize=6, pad=1)
    
        ax.set_xlabel('U-Pb age (Ga)')
        ax.set_ylabel('$\epsilon$Hf(t)')

        plt.show()
        plt.savefig('./figure/bootstrap.pdf')

    return T_arr, eps_arr
