
import const as c
from bootstrap import epsilon, generate_Hf, Hf2Lu
from bootstrap import solve_T, extrapolate, get_crust_Lu

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

matplotlib.rcParams['font.family'] = 'Arial'
fontsize = 11

def plot_example(data, mother_sample, ax, colors, ind=480, times=100,
                meas_err=1, Hf_am_err=1, Lu_crust_err=1):
    '''
    This is to plot an example with any index showing epsilon_hf versus U-Pb age
    '''
    t0, Hf0 = float(data[ind:ind+1].u_pb_age / 1e3), float(data[ind:ind+1].hf_hf)
    Lu0, oxy0 = float(data[ind:ind+1].lu_hf), float(data[ind:ind+1].o)
    t_err, Hf_err = float(data[ind:ind+1].age_2se / 1e3 / 2), float(data[ind:ind+1].hf_hf_2se / 2)
    Lu_err, oxy_err = float(data[ind:ind+1].lu_hf_2se / 2), float(data[ind:ind+1].o_sed)
    #! Note here O error is 1 sed 
   
    ax.set_xlim(0, 4.)
    ax.set_ylim(-30, 20)

    alpha = .01
    lw = 2

    T_list = []
    for i in range(times):
        # In each bootstrap, generate 1. random t, hf, lu, oxy based on measured errors
        t = np.random.normal(t0, t_err * meas_err)
        Hf = np.random.normal(Hf0, Hf_err * meas_err)
        Lu = np.random.normal(Lu0, Lu_err * meas_err)
        oxy = np.random.normal(oxy0, oxy_err * meas_err)

        # 2. a random AM 176Hf_177Hf from mother_sample 
        fixed = None if Hf_am_err else c.am_hh
        Hf_am = generate_Hf(mother_sample, fixed=fixed)  
        Lu_am = Hf2Lu(Hf_am)

        # 3. a random oxygen-based crustal 176Lu_177Hf from its error envelop
        Lu_crustal = get_crust_Lu(oxy, Lu_crust_err=Lu_crust_err)

        # Calculating model age T 
        T = solve_T(t, Hf, Lu, Hf_am=Hf_am, Lu_am=Lu_am, Lu_crustal=Lu_crustal)
        T_list.append(float(T))
        
        x1 = np.arange(0, t, .001) # zircon 
        ax.plot(x1, epsilon(x1, Hf, Lu), color=colors[0], alpha=alpha, lw=lw) 

        x2 = np.arange(t, T, .001) # crustal 
        ax.plot(x2, (extrapolate(x2 - t, extrapolate(t, Hf, Lu), Lu_crustal) /
                    extrapolate(x2, c.chur_hh, c.chur_lh) - 1) * 1e4, 
                    color=colors[1], alpha=alpha, lw=lw, zorder=6) 
        
        x3 = np.arange(0, 4., .001) # AM curve
        ax.plot(x3, epsilon(x3, Hf_am, Lu_am), color=colors[2], alpha=alpha, lw=lw, zorder=2) 
        ax.vlines(x=T, ymin=-50., ymax=epsilon(T, Hf_am, Lu_am), 
                    color=colors[3], alpha=alpha, lw=lw, linestyle="-", zorder=3) # Model age intercepts
    # print(np.percentile(T_list, 50), np.percentile(T_list, 2.5), np.percentile(T_list, 97.5))
    return T_list

def plot_steps(data, mother_sample, colors, 
                ind=242, times=1000, percentiles=(2.5, 97.5),
                save_fig=0, show=1):
    """
    A 2x2 plot to show the steps of calculating model age errors using bootstrap
    (1-3): individual effect; 
    (4): combined effect
    """

    _, ((a1, a2), (a3, a4)) = plt.subplots(2, 2, figsize=(6, 4.9), sharex=True, sharey=True)

    options = [[0, 1, 0], [0, 0, 1], [0, 1, 0], [1, 1, 1]]
    for ax, option, ilabel in zip([a1, a2, a3, a4], options, ['A', 'B', 'C', 'D']):
        T_list = plot_example(data, mother_sample, ax, colors, ind=ind, times=times, 
                meas_err=option[0], Hf_am_err=option[1], Lu_crust_err=option[2])

        ax.annotate(ilabel, xy=(0.92, 0.91), xycoords="axes fraction", 
                    fontsize=fontsize+2, fontweight='extra bold')
        ax.set_xlim(0, 4.)
        ax.set_ylim(-30, 20)
    
        # Insert a subplot within each panel 
        axins = inset_axes(ax, width="40%", height="40%", loc=4, borderpad=1.3)
        axins.hist(T_list, bins=30, linewidth=.8, histtype='step', color='k')
        axins.axvline(x=np.percentile(T_list, 50), ls='--', lw=1, color='k')
        axins.axvline(x=np.percentile(T_list, percentiles[0]), ls='--', lw=1, color='k')
        axins.axvline(x=np.percentile(T_list, percentiles[1]), ls='--', lw=1, color='k')
        axins.tick_params(axis='both', which='major', length=1, labelsize=fontsize-5, pad=1)
        axins.set_xlim(0.85, 1.45)

    a3.set_xlabel('U-Pb age (Ga)')
    a4.set_xlabel('U-Pb age (Ga)')
    a1.set_ylabel('$\epsilon$Hf(t)')
    a3.set_ylabel('$\epsilon$Hf(t)')
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1, wspace=0.05)

    if save_fig:
        plt.savefig('./figure/plot_'+ind+'_.pdf', bbox_inches='tight', dpi=600)
    if show:
        plt.show()
    
data = pd.read_excel('./data/all_europe_raw.xlsx', sheet_name="o_hf")
    
am_hf_hf_ori = pd.read_csv('./data/hf_dist_392.csv')['176Hf_177Hf']
mother_sample = am_hf_hf_ori[(am_hf_hf_ori < np.percentile(am_hf_hf_ori, 97.5))
                & (am_hf_hf_ori > np.percentile(am_hf_hf_ori, 2.5))]

blue, green, red, dimgrey = '#367db7', 'g', '#df1d27', 'dimgrey'
colors = [blue, green, red, dimgrey]

plot_steps(data, mother_sample, colors, 
            ind=242, times=100, percentiles=(2.5, 97.5),
            save_fig=0, show=1)

