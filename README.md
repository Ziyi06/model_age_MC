# model_age_MC
Estimates on zircon Hf model ages and their uncertainties using Monte Carlo bootstrapping, used in Zhu et al. in submission.

This module provides a model age and a εHf(t) along with their errors, by bootstrapping (i) uncertainty in the Hf isotopic composition of the mantle reservoir, (ii) uncertainty in the 176Lu/177Hf for the crustal source region, and (iii) uncertainty in measurements. One can decide which uncertainty(s) to be considered. The minimum requirement of data is a zircon U-Pb age, a 176Hf/177Hf, and a 176Lu/177Hf.

Run "python bootstrap_plot.py" as a demo. T_boot() is the core function. The detailed parameter description can be found in the comments. 

If you already have an array containing a 176Hf/177Hf distribution of arc/depleted mantle, load it and set it to the parameter "mantle_hh" in T_boot(). If not, you can use the recommended data compilation of modern arc mantle or leave it defaulted as a constant. 
Note that the parameter "plot" can only be True when the inputs are all float. If the inputs are array-like, remember to switch it to False.
